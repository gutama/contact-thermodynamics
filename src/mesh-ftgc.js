/**
 * FTGC (Fundamental Theorem of Geometric Calculus) on Triangle Meshes
 * 
 * Implements the discrete geometric derivative ∇ on triangle meshes using
 * the split differential operator: ∇ = Σᵢ eⁱ ⊗ ∂/∂xᵢ
 * 
 * Key insight: Cotan weights encode the mesh metric (reciprocal basis).
 * 
 * Staggered storage by grade:
 *   Grade 0 (scalars)   → Vertices
 *   Grade 1 (vectors)   → Edges
 *   Grade 2 (bivectors) → Faces
 * 
 * The geometric derivative unifies DEC operators:
 *   ∇ = ∇· + ∇∧  (inner + outer = divergence + curl)
 * 
 * FTGC: ∫_M ∇F dV = ∮_∂M F dS
 * 
 * @module mesh-ftgc
 * @license MIT
 */

(function (global) {
    'use strict';

    const EPSILON = 1e-10;

    // ============================================================================
    // SPARSE MATRIX (minimal COO/CSR-like for Laplacian etc.)
    // ============================================================================

    /**
     * Simple sparse matrix in CSR format for efficient matrix-vector products.
     */
    class SparseMatrix {
        /**
         * @param {number} rows 
         * @param {number} cols 
         */
        constructor(rows, cols) {
            this.rows = rows;
            this.cols = cols;
            // COO format for building
            this._entries = new Map(); // "row,col" -> value
        }

        /**
         * Set entry (accumulates if called multiple times).
         */
        set(row, col, value) {
            const key = `${row},${col}`;
            this._entries.set(key, (this._entries.get(key) || 0) + value);
        }

        /**
         * Get entry.
         */
        get(row, col) {
            return this._entries.get(`${row},${col}`) || 0;
        }

        /**
         * Build CSR format for efficient matvec.
         */
        buildCSR() {
            const entries = [];
            for (const [key, val] of this._entries) {
                const [r, c] = key.split(',').map(Number);
                entries.push({ r, c, v: val });
            }
            entries.sort((a, b) => a.r - b.r || a.c - b.c);

            this._csrRowPtr = new Uint32Array(this.rows + 1);
            this._csrColIdx = new Uint32Array(entries.length);
            this._csrValues = new Float64Array(entries.length);

            let ptr = 0;
            let currentRow = 0;
            for (let i = 0; i < entries.length; i++) {
                const { r, c, v } = entries[i];
                while (currentRow <= r) {
                    this._csrRowPtr[currentRow] = ptr;
                    currentRow++;
                }
                this._csrColIdx[ptr] = c;
                this._csrValues[ptr] = v;
                ptr++;
            }
            while (currentRow <= this.rows) {
                this._csrRowPtr[currentRow] = ptr;
                currentRow++;
            }
        }

        /**
         * Matrix-vector product: y = A * x
         * @param {Float64Array} x 
         * @returns {Float64Array}
         */
        matvec(x) {
            if (!this._csrRowPtr) this.buildCSR();

            const y = new Float64Array(this.rows);
            for (let row = 0; row < this.rows; row++) {
                let sum = 0;
                for (let ptr = this._csrRowPtr[row]; ptr < this._csrRowPtr[row + 1]; ptr++) {
                    sum += this._csrValues[ptr] * x[this._csrColIdx[ptr]];
                }
                y[row] = sum;
            }
            return y;
        }

        /**
         * Transpose matrix-vector product: y = Aᵀ * x
         * @param {Float64Array} x 
         * @returns {Float64Array}
         */
        matvecT(x) {
            if (!this._csrRowPtr) this.buildCSR();

            const y = new Float64Array(this.cols);
            for (let row = 0; row < this.rows; row++) {
                const xr = x[row];
                for (let ptr = this._csrRowPtr[row]; ptr < this._csrRowPtr[row + 1]; ptr++) {
                    y[this._csrColIdx[ptr]] += this._csrValues[ptr] * xr;
                }
            }
            return y;
        }

        /**
         * Get diagonal entries as Float64Array.
         */
        diagonal() {
            const diag = new Float64Array(Math.min(this.rows, this.cols));
            for (let i = 0; i < diag.length; i++) {
                diag[i] = this.get(i, i);
            }
            return diag;
        }
    }

    /**
     * Create a diagonal matrix from an array.
     * @param {Float64Array} diag 
     * @returns {SparseMatrix}
     */
    function diagMatrix(diag) {
        const n = diag.length;
        const M = new SparseMatrix(n, n);
        for (let i = 0; i < n; i++) {
            M.set(i, i, diag[i]);
        }
        M.buildCSR();
        return M;
    }

    // ============================================================================
    // MESH MULTIVECTOR FIELD (Staggered Storage)
    // ============================================================================

    /**
     * A multivector-valued field on a triangle mesh with staggered storage.
     * 
     * Different grades live on different mesh elements:
     *   Grade 0 (scalars)   → Vertices  (nVertices values)
     *   Grade 1 (vectors)   → Edges     (nEdges values)
     *   Grade 2 (bivectors) → Faces     (nFaces values)
     */
    class MeshMultivectorField {
        /**
         * @param {TriangleMesh} mesh 
         * @param {Float64Array|null} grade0 - Scalar field on vertices
         * @param {Float64Array|null} grade1 - Vector field on edges
         * @param {Float64Array|null} grade2 - Bivector field on faces
         */
        constructor(mesh, grade0 = null, grade1 = null, grade2 = null) {
            this.mesh = mesh;
            this.grade0 = grade0 || new Float64Array(mesh.nVertices);
            this.grade1 = grade1 || new Float64Array(mesh.nEdges);
            this.grade2 = grade2 || new Float64Array(mesh.nFaces);
        }

        /**
         * Add two multivector fields.
         */
        add(other) {
            const g0 = new Float64Array(this.grade0.length);
            const g1 = new Float64Array(this.grade1.length);
            const g2 = new Float64Array(this.grade2.length);
            for (let i = 0; i < g0.length; i++) g0[i] = this.grade0[i] + other.grade0[i];
            for (let i = 0; i < g1.length; i++) g1[i] = this.grade1[i] + other.grade1[i];
            for (let i = 0; i < g2.length; i++) g2[i] = this.grade2[i] + other.grade2[i];
            return new MeshMultivectorField(this.mesh, g0, g1, g2);
        }

        /**
         * Subtract two multivector fields.
         */
        sub(other) {
            const g0 = new Float64Array(this.grade0.length);
            const g1 = new Float64Array(this.grade1.length);
            const g2 = new Float64Array(this.grade2.length);
            for (let i = 0; i < g0.length; i++) g0[i] = this.grade0[i] - other.grade0[i];
            for (let i = 0; i < g1.length; i++) g1[i] = this.grade1[i] - other.grade1[i];
            for (let i = 0; i < g2.length; i++) g2[i] = this.grade2[i] - other.grade2[i];
            return new MeshMultivectorField(this.mesh, g0, g1, g2);
        }

        /**
         * Scalar multiplication.
         */
        scale(s) {
            const g0 = new Float64Array(this.grade0.length);
            const g1 = new Float64Array(this.grade1.length);
            const g2 = new Float64Array(this.grade2.length);
            for (let i = 0; i < g0.length; i++) g0[i] = this.grade0[i] * s;
            for (let i = 0; i < g1.length; i++) g1[i] = this.grade1[i] * s;
            for (let i = 0; i < g2.length; i++) g2[i] = this.grade2[i] * s;
            return new MeshMultivectorField(this.mesh, g0, g1, g2);
        }

        /**
         * Extract only grade-k component.
         */
        gradeSelect(k) {
            const result = new MeshMultivectorField(this.mesh);
            if (k === 0) result.grade0 = new Float64Array(this.grade0);
            else if (k === 1) result.grade1 = new Float64Array(this.grade1);
            else if (k === 2) result.grade2 = new Float64Array(this.grade2);
            return result;
        }

        /**
         * Clone the field.
         */
        clone() {
            return new MeshMultivectorField(
                this.mesh,
                new Float64Array(this.grade0),
                new Float64Array(this.grade1),
                new Float64Array(this.grade2)
            );
        }

        /**
         * Create a pure scalar (grade-0) field on vertices.
         */
        static scalarField(mesh, values) {
            return new MeshMultivectorField(mesh, new Float64Array(values), null, null);
        }

        /**
         * Create a pure vector (grade-1) field on edges.
         */
        static vectorField(mesh, values) {
            return new MeshMultivectorField(mesh, null, new Float64Array(values), null);
        }

        /**
         * Create a pure bivector (grade-2) field on faces.
         */
        static bivectorField(mesh, values) {
            return new MeshMultivectorField(mesh, null, null, new Float64Array(values));
        }
    }

    // ============================================================================
    // MESH GEOMETRIC DERIVATIVE ∇
    // ============================================================================

    /**
     * The discrete geometric derivative ∇ on triangle meshes.
     * 
     * Implements the split differential operator: ∇ = Σᵢ eⁱ ⊗ ∂/∂xᵢ
     * 
     * where:
     *   - eⁱ are reciprocal basis vectors (encoded in cotan weights)
     *   - ∂/∂xᵢ are partial derivatives (finite differences along edges)
     * 
     * The geometric derivative unifies what DEC calls separate operators:
     *   ∇ = ∇· + ∇∧  (inner + outer = divergence + curl)
     */
    class MeshGeometricDerivative {
        /**
         * @param {TriangleMesh} mesh 
         */
        constructor(mesh) {
            this.mesh = mesh;

            // Build split operator components (sparse matrices)
            this._wedge_0to1 = this._buildWedge0to1();  // scalar → vector
            this._wedge_1to2 = this._buildWedge1to2();  // vector → bivector

            // Metric tensors (diagonal, encode reciprocal basis via cotan weights)
            this._metricVertices = this._buildVertexMetric();  // dual cell areas
            this._metricEdges = this._buildEdgeMetric();       // COTAN WEIGHTS!
            this._metricFaces = this._buildFaceMetric();       // triangle areas

            // Inverse metrics (for raising indices)
            this._metricVerticesInv = this._invertDiag(this._metricVertices);
            this._metricEdgesInv = this._invertDiag(this._metricEdges);
            this._metricFacesInv = this._invertDiag(this._metricFaces);

            // Precomputed Laplacian
            this._laplacian = null;
            this._laplacianWeak = null;
        }

        _invertDiag(diag) {
            const inv = new Float64Array(diag.length);
            for (let i = 0; i < diag.length; i++) {
                inv[i] = Math.abs(diag[i]) > EPSILON ? 1 / diag[i] : 0;
            }
            return inv;
        }

        // ========================================================================
        // BUILD SPLIT OPERATOR COMPONENTS
        // ========================================================================

        /**
         * Build ∇∧ from grade 0 to grade 1 (outer derivative on scalars).
         * (∇∧f)[edge] = f[j] - f[i]
         * 
         * This is the discrete exterior derivative d₀ in DEC language.
         */
        _buildWedge0to1() {
            const nV = this.mesh.nVertices;
            const nE = this.mesh.nEdges;
            const edges = this.mesh.edges;

            const mat = new SparseMatrix(nE, nV);
            for (let e = 0; e < nE; e++) {
                const i = edges[e * 2];
                const j = edges[e * 2 + 1];
                mat.set(e, i, -1);
                mat.set(e, j, +1);
            }
            mat.buildCSR();
            return mat;
        }

        /**
         * Build ∇∧ from grade 1 to grade 2 (outer derivative on vectors).
         * (∇∧ω)[face] = circulation around face
         * 
         * This is the discrete exterior derivative d₁ in DEC language.
         * Note: ∇∧(∇∧f) = 0 identically (curl of gradient is zero)
         */
        _buildWedge1to2() {
            const nE = this.mesh.nEdges;
            const nF = this.mesh.nFaces;
            const faces = this.mesh.faces;

            const mat = new SparseMatrix(nF, nE);

            for (let t = 0; t < nF; t++) {
                const i = faces[t * 3];
                const j = faces[t * 3 + 1];
                const k = faces[t * 3 + 2];

                const triEdges = [[i, j], [j, k], [k, i]];
                for (const [a, b] of triEdges) {
                    const eIdx = this.mesh.getEdgeIndex(a, b);
                    // Sign depends on edge orientation vs triangle orientation
                    const sign = a < b ? 1 : -1;
                    mat.set(t, eIdx, sign);
                }
            }
            mat.buildCSR();
            return mat;
        }

        /**
         * Build metric at vertices (mixed Voronoi dual cell areas).
         */
        _buildVertexMetric() {
            return this.mesh.vertexDualAreas();
        }

        /**
         * Build metric at edges: THE COTAN WEIGHTS!
         * 
         * This is where the "Swiss Army knife" appears in geometric calculus:
         * The cotan weights encode the inner product of reciprocal basis vectors.
         */
        _buildEdgeMetric() {
            return this.mesh.cotanWeights();
        }

        /**
         * Build metric at faces (triangle areas).
         */
        _buildFaceMetric() {
            return this.mesh.faceAreas();
        }

        // ========================================================================
        // GEOMETRIC DERIVATIVE OPERATIONS
        // ========================================================================

        /**
         * Apply the outer derivative ∇∧ (grade-raising part of ∇).
         * ∇∧: grade k → grade k+1
         */
        wedge(F) {
            const result = new MeshMultivectorField(this.mesh);
            // ∇∧(scalar) → vector
            result.grade1 = this._wedge_0to1.matvec(F.grade0);
            // ∇∧(vector) → bivector
            result.grade2 = this._wedge_1to2.matvec(F.grade1);
            return result;
        }

        /**
         * Apply the inner derivative ∇· (grade-lowering part of ∇).
         * ∇·: grade k → grade k-1
         * 
         * Using the metric (cotan weights):
         *   ∇·V = M⁻¹_vertices (∇∧)ᵀ M_edges V
         */
        inner(F) {
            const result = new MeshMultivectorField(this.mesh);

            // ∇·(vector) → scalar (divergence)
            // δ₁ = M₀⁻¹ d₀ᵀ M₁
            const tmp1 = new Float64Array(this.mesh.nEdges);
            for (let i = 0; i < tmp1.length; i++) {
                tmp1[i] = this._metricEdges[i] * F.grade1[i];
            }
            const tmp2 = this._wedge_0to1.matvecT(tmp1);
            for (let i = 0; i < result.grade0.length; i++) {
                result.grade0[i] = this._metricVerticesInv[i] * tmp2[i];
            }

            // ∇·(bivector) → vector
            // δ₂ = M₁⁻¹ d₁ᵀ M₂
            const tmp3 = new Float64Array(this.mesh.nFaces);
            for (let i = 0; i < tmp3.length; i++) {
                tmp3[i] = this._metricFaces[i] * F.grade2[i];
            }
            const tmp4 = this._wedge_1to2.matvecT(tmp3);
            for (let i = 0; i < result.grade1.length; i++) {
                result.grade1[i] = this._metricEdgesInv[i] * tmp4[i];
            }

            return result;
        }

        /**
         * Apply the full geometric derivative ∇.
         * ∇F = ∇·F + ∇∧F
         */
        apply(F) {
            return this.wedge(F).add(this.inner(F));
        }

        /**
         * Gradient of scalar field (convenience method).
         * ∇f for scalar f gives a vector (grade 1).
         * @param {Float64Array} f - Scalar field on vertices
         * @returns {Float64Array} Vector field on edges
         */
        grad(f) {
            return this._wedge_0to1.matvec(f);
        }

        /**
         * Divergence of vector field (convenience method).
         * ∇·V for vector V gives a scalar (grade 0).
         * @param {Float64Array} V - Vector field on edges
         * @returns {Float64Array} Scalar field on vertices
         */
        div(V) {
            const tmp = new Float64Array(this.mesh.nEdges);
            for (let i = 0; i < tmp.length; i++) {
                tmp[i] = this._metricEdges[i] * V[i];
            }
            const tmp2 = this._wedge_0to1.matvecT(tmp);
            const result = new Float64Array(this.mesh.nVertices);
            for (let i = 0; i < result.length; i++) {
                result[i] = this._metricVerticesInv[i] * tmp2[i];
            }
            return result;
        }

        /**
         * Curl of vector field (convenience method).
         * ∇∧V for vector V gives a bivector (grade 2).
         * @param {Float64Array} V - Vector field on edges
         * @returns {Float64Array} Bivector field on faces
         */
        curl(V) {
            return this._wedge_1to2.matvec(V);
        }

        /**
         * Laplacian of scalar field.
         * ∇²f = ∇·(∇f) = div(grad(f))
         * 
         * This is the famous cotan Laplacian!
         * Convention: ∇² is negative semi-definite (negative at maxima).
         * 
         * @param {Float64Array} f - Scalar field on vertices
         * @param {Object} opts - Options
         * @param {Uint8Array} opts.dirichletMask - Mask for Dirichlet vertices (1 = fixed)
         * @param {Float64Array} opts.dirichletValues - Values at Dirichlet vertices
         * @returns {Float64Array} Laplacian on vertices
         */
        laplacian(f, opts = {}) {
            const { dirichletMask, dirichletValues } = opts;

            const gradF = this.grad(f);
            // Negate to get standard convention: ∇²f ≤ 0 at maxima
            const lapF = this.div(gradF);
            for (let i = 0; i < lapF.length; i++) {
                lapF[i] = -lapF[i];
            }

            // Apply Dirichlet constraints
            if (dirichletMask) {
                for (let i = 0; i < lapF.length; i++) {
                    if (dirichletMask[i]) {
                        lapF[i] = 0; // No contribution from fixed vertices
                    }
                }
            }

            return lapF;
        }

        /**
         * Build the Laplacian as a sparse matrix.
         * L = -M₀⁻¹ (∇∧)ᵀ M₁ (∇∧)
         * @returns {SparseMatrix}
         */
        laplacianMatrix() {
            if (this._laplacian) return this._laplacian;

            const nV = this.mesh.nVertices;
            const nE = this.mesh.nEdges;
            const edges = this.mesh.edges;

            const L = new SparseMatrix(nV, nV);

            // Cotan Laplacian: L_ij = -w_ij for neighbors, L_ii = Σ w_ij
            for (let e = 0; e < nE; e++) {
                const i = edges[e * 2];
                const j = edges[e * 2 + 1];
                const w = this._metricEdges[e];

                // Off-diagonal (negated for standard convention)
                L.set(i, j, w);
                L.set(j, i, w);

                // Diagonal
                L.set(i, i, -w);
                L.set(j, j, -w);
            }

            // Apply inverse vertex metric (mass matrix)
            // L_final = M⁻¹ L
            for (let i = 0; i < nV; i++) {
                const scale = this._metricVerticesInv[i];
                for (let j = 0; j < nV; j++) {
                    const val = L.get(i, j);
                    if (Math.abs(val) > EPSILON) {
                        L._entries.set(`${i},${j}`, val * scale);
                    }
                }
            }

            L.buildCSR();
            this._laplacian = L;
            return L;
        }

        /**
         * Build the "weak" Laplacian (without mass matrix inverse).
         * L_weak = -(∇∧)ᵀ M₁ (∇∧)
         * @returns {SparseMatrix}
         */
        laplacianWeak() {
            if (this._laplacianWeak) return this._laplacianWeak;

            const nV = this.mesh.nVertices;
            const nE = this.mesh.nEdges;
            const edges = this.mesh.edges;

            const L = new SparseMatrix(nV, nV);

            for (let e = 0; e < nE; e++) {
                const i = edges[e * 2];
                const j = edges[e * 2 + 1];
                const w = this._metricEdges[e];

                L.set(i, j, w);
                L.set(j, i, w);
                L.set(i, i, -w);
                L.set(j, j, -w);
            }

            L.buildCSR();
            this._laplacianWeak = L;
            return L;
        }

        // ========================================================================
        // VERIFICATION / SANITY CHECKS
        // ========================================================================

        /**
         * Verify ∇∧(∇∧f) = 0 (curl of gradient is zero).
         * @param {Float64Array} f - Test scalar field
         * @returns {number} Maximum absolute error
         */
        verifyCurlGradZero(f) {
            const gradF = this.grad(f);
            const curlGradF = this.curl(gradF);
            let maxErr = 0;
            for (let i = 0; i < curlGradF.length; i++) {
                maxErr = Math.max(maxErr, Math.abs(curlGradF[i]));
            }
            return maxErr;
        }

        /**
         * Verify ∇²(constant) ≈ 0.
         * @param {number} c - Constant value
         * @returns {number} Maximum absolute error
         */
        verifyLaplacianConstantZero(c = 1) {
            const f = new Float64Array(this.mesh.nVertices).fill(c);
            const lapF = this.laplacian(f);
            let maxErr = 0;
            for (let i = 0; i < lapF.length; i++) {
                maxErr = Math.max(maxErr, Math.abs(lapF[i]));
            }
            return maxErr;
        }
    }

    // ============================================================================
    // EXPORTS
    // ============================================================================

    const MeshFTGCModule = {
        SparseMatrix,
        diagMatrix,
        MeshMultivectorField,
        MeshGeometricDerivative
    };

    // CommonJS
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = MeshFTGCModule;
    }
    // AMD
    else if (typeof define === 'function' && define.amd) {
        define(function () { return MeshFTGCModule; });
    }
    // Browser global
    else {
        global.ContactThermo = global.ContactThermo || {};
        Object.assign(global.ContactThermo, MeshFTGCModule);
    }

})(typeof globalThis !== 'undefined' ? globalThis : typeof window !== 'undefined' ? window : this);
