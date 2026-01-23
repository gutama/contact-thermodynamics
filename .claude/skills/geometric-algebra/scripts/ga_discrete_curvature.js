/**
 * Discrete Differential Geometry via Geometric Algebra
 * 
 * Implements discrete curvature concepts from Keenan Crane's DDG course,
 * adapted for Geometric Algebra (PGA and CGA).
 * 
 * Key insight: The wedge product naturally computes determinants and signed areas/volumes,
 * which are fundamental to discrete curvature calculations.
 * 
 * This module provides both:
 * - INTEGRAL APPROACH: Angle defect, cotan-Laplace (local vertex quantities)
 * - VARIATIONAL APPROACH: Steiner's formula, curvature flows (global derivatives)
 * 
 * References:
 * - Keenan Crane, "Discrete Differential Geometry: An Applied Introduction"
 * - De Keninck, Roelfs, Dorst, Eelbode, "Clean up your Mesh!"
 * - Roelfs & De Keninck, "From Invariant Decomposition to Spinors"
 * 
 * @license MIT
 * @version 2.0.0
 */

(function(global) {
    'use strict';

    const EPSILON = 1e-10;
    const abs = Math.abs;
    const sqrt = Math.sqrt;
    const sin = Math.sin;
    const cos = Math.cos;
    const atan2 = Math.atan2;
    const acos = Math.acos;
    const PI = Math.PI;

    // ============================================================================
    // WEDGE PRODUCT AND DETERMINANTS
    // ============================================================================

    /**
     * The wedge product of n vectors in n-dimensional space equals the determinant
     * of the matrix formed by those vectors as columns, times the pseudoscalar.
     * 
     * For 2 vectors: a ∧ b = (a₁b₂ - a₂b₁) e₁₂
     * For 3 vectors: a ∧ b ∧ c = det([a,b,c]) e₁₂₃
     * 
     * This provides a coordinate-free way to compute signed areas and volumes.
     */
    const WedgeDeterminant = {
        /**
         * 2D determinant via wedge product
         * Returns signed area of parallelogram
         */
        det2D(a, b) {
            return a[0] * b[1] - a[1] * b[0];
        },

        /**
         * 3D determinant via wedge product (scalar triple product)
         * Returns signed volume of parallelepiped
         */
        det3D(a, b, c) {
            return a[0] * (b[1] * c[2] - b[2] * c[1])
                 - a[1] * (b[0] * c[2] - b[2] * c[0])
                 + a[2] * (b[0] * c[1] - b[1] * c[0]);
        },

        /**
         * Area of triangle via wedge: A = ½|v₀v₁ ∧ v₀v₂|
         */
        triangleArea(v0, v1, v2) {
            const e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
            const e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
            // Cross product magnitude = |e1 ∧ e2|
            const cx = e1[1] * e2[2] - e1[2] * e2[1];
            const cy = e1[2] * e2[0] - e1[0] * e2[2];
            const cz = e1[0] * e2[1] - e1[1] * e2[0];
            return 0.5 * sqrt(cx * cx + cy * cy + cz * cz);
        },

        /**
         * Signed volume of tetrahedron via wedge: V = ⅙|v₀v₁ ∧ v₀v₂ ∧ v₀v₃|
         * This is 1/6 the determinant of the edge matrix
         */
        tetraVolume(v0, v1, v2, v3) {
            const e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
            const e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
            const e3 = [v3[0] - v0[0], v3[1] - v0[1], v3[2] - v0[2]];
            return this.det3D(e1, e2, e3) / 6;
        },

        /**
         * Normal vector via wedge (un-normalized)
         * n = (v₁ - v₀) ∧ (v₂ - v₀) in ℝ³
         */
        triangleNormal(v0, v1, v2) {
            const e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
            const e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
            return [
                e1[1] * e2[2] - e1[2] * e2[1],
                e1[2] * e2[0] - e1[0] * e2[2],
                e1[0] * e2[1] - e1[1] * e2[0]
            ];
        }
    };

    // ============================================================================
    // DISCRETE CURVATURE (following Crane's DDG)
    // ============================================================================

    /**
     * Discrete Curvature Calculator
     * 
     * INTEGRAL APPROACH (local vertex quantities):
     * 1. Angle Defect (Gaussian Curvature): K_i = 2π - Σ θ_j
     * 2. Mean Curvature Normal (Laplace-Beltrami): H⃗_i = (1/2A_i) Σ (cot α + cot β)(v_j - v_i)
     * 3. Gauss-Bonnet: Σ K_i = 2π·χ
     * 
     * VARIATIONAL APPROACH (global derivatives):
     * 4. Steiner's Formula: V(t) = V + At + Ht² + (4π/3)Kt³
     * 5. Gradient Hierarchy: ∇Vol→Area, ∇Area→H, ∇H→K, ∇K→0
     * 6. Curvature Flows: Mean curvature flow, Willmore flow
     */
    class DiscreteCurvature {
        constructor() {
            this.vertices = [];
            this.faces = [];
            this.edges = [];
            this.vertexFaces = [];  // faces incident to each vertex
            this.vertexEdges = [];  // edges incident to each vertex
        }

        /**
         * Load mesh from vertex/face arrays
         * @param {Array} vertices - [[x,y,z], ...]
         * @param {Array} faces - [[i,j,k], ...] (triangles)
         */
        loadMesh(vertices, faces) {
            this.vertices = vertices.map(v => [...v]);
            this.faces = faces.map(f => [...f]);
            this._buildTopology();
        }

        _buildTopology() {
            const V = this.vertices.length;
            const F = this.faces.length;
            
            // Initialize per-vertex lists
            this.vertexFaces = Array.from({ length: V }, () => []);
            this.vertexEdges = Array.from({ length: V }, () => new Set());
            
            // Edge map: canonical key -> [v0, v1, f_left, f_right]
            const edgeMap = new Map();
            
            for (let fi = 0; fi < F; fi++) {
                const f = this.faces[fi];
                for (let j = 0; j < 3; j++) {
                    const v0 = f[j];
                    const v1 = f[(j + 1) % 3];
                    
                    this.vertexFaces[v0].push(fi);
                    
                    const key = v0 < v1 ? `${v0}-${v1}` : `${v1}-${v0}`;
                    if (!edgeMap.has(key)) {
                        edgeMap.set(key, { v0: Math.min(v0, v1), v1: Math.max(v0, v1), faces: [] });
                    }
                    edgeMap.get(key).faces.push(fi);
                    
                    this.vertexEdges[v0].add(key);
                    this.vertexEdges[v1].add(key);
                }
            }
            
            this.edges = Array.from(edgeMap.values());
            this.edgeMap = edgeMap;
        }

        // ========================================================================
        // BASIC GEOMETRY
        // ========================================================================

        /**
         * Edge length between two vertices
         */
        edgeLength(v0, v1) {
            const p0 = this.vertices[v0];
            const p1 = this.vertices[v1];
            const dx = p1[0] - p0[0];
            const dy = p1[1] - p0[1];
            const dz = p1[2] - p0[2];
            return sqrt(dx * dx + dy * dy + dz * dz);
        }

        /**
         * Tip angle at vertex v0 in triangle [v0, v1, v2]
         * Uses: cos(θ) = (e1·e2)/(|e1||e2|)
         */
        tipAngle(v0, v1, v2) {
            const p0 = this.vertices[v0];
            const p1 = this.vertices[v1];
            const p2 = this.vertices[v2];
            
            const e1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
            const e2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
            
            const dot = e1[0] * e2[0] + e1[1] * e2[1] + e1[2] * e2[2];
            const len1 = sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
            const len2 = sqrt(e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2]);
            
            if (len1 < EPSILON || len2 < EPSILON) return 0;
            
            const cosTheta = Math.max(-1, Math.min(1, dot / (len1 * len2)));
            return acos(cosTheta);
        }

        /**
         * Cotangent of angle opposite to edge in a triangle
         * cot(θ) = cos(θ)/sin(θ)
         */
        cotanOpposite(faceIdx, edgeV0, edgeV1) {
            const f = this.faces[faceIdx];
            // Find the vertex opposite to the edge
            let opposite = -1;
            for (let i = 0; i < 3; i++) {
                if (f[i] !== edgeV0 && f[i] !== edgeV1) {
                    opposite = f[i];
                    break;
                }
            }
            if (opposite < 0) return 0;
            
            const angle = this.tipAngle(opposite, edgeV0, edgeV1);
            const sinA = sin(angle);
            if (abs(sinA) < EPSILON) return 0;
            return cos(angle) / sinA;
        }

        /**
         * Vertex normal (area-weighted average of incident face normals)
         */
        vertexNormal(vertexIdx) {
            let nx = 0, ny = 0, nz = 0;
            
            for (const fi of this.vertexFaces[vertexIdx]) {
                const f = this.faces[fi];
                const n = WedgeDeterminant.triangleNormal(
                    this.vertices[f[0]], this.vertices[f[1]], this.vertices[f[2]]
                );
                nx += n[0];
                ny += n[1];
                nz += n[2];
            }
            
            const len = sqrt(nx * nx + ny * ny + nz * nz);
            if (len < EPSILON) return [0, 0, 1];
            return [nx / len, ny / len, nz / len];
        }

        // ========================================================================
        // INTEGRAL APPROACH: GAUSSIAN CURVATURE (ANGLE DEFECT)
        // ========================================================================

        /**
         * Gaussian Curvature at vertex (Angle Defect)
         * K = 2π - Σ θ_j for interior vertices
         * 
         * This is the integrated Gaussian curvature over the dual cell.
         * To get pointwise curvature, divide by the dual cell area.
         */
        gaussianCurvature(vertexIdx) {
            const faces = this.vertexFaces[vertexIdx];
            let angleSum = 0;
            
            for (const fi of faces) {
                const f = this.faces[fi];
                // Find position of vertex in face
                const idx = f.indexOf(vertexIdx);
                const v1 = f[(idx + 1) % 3];
                const v2 = f[(idx + 2) % 3];
                angleSum += this.tipAngle(vertexIdx, v1, v2);
            }
            
            // Check if boundary vertex
            const isBoundary = this._isBoundaryVertex(vertexIdx);
            return isBoundary ? (PI - angleSum) : (2 * PI - angleSum);
        }

        /**
         * Check if vertex is on boundary
         */
        _isBoundaryVertex(vertexIdx) {
            for (const edgeKey of this.vertexEdges[vertexIdx]) {
                const edge = this.edgeMap.get(edgeKey);
                if (edge.faces.length === 1) return true;
            }
            return false;
        }

        /**
         * Total Gaussian Curvature (should equal 2π·χ by Gauss-Bonnet)
         */
        totalGaussianCurvature() {
            let total = 0;
            for (let i = 0; i < this.vertices.length; i++) {
                total += this.gaussianCurvature(i);
            }
            return total;
        }

        // ========================================================================
        // INTEGRAL APPROACH: MEAN CURVATURE (COTAN-LAPLACE)
        // ========================================================================

        /**
         * Mean Curvature Normal at vertex (Discrete Laplace-Beltrami)
         * Hn = (1/2A) Σ (cot α + cot β)(v_j - v_i)
         * 
         * Returns [Hx, Hy, Hz] - the mean curvature times normal vector
         */
        meanCurvatureNormal(vertexIdx) {
            const vi = this.vertices[vertexIdx];
            let Hn = [0, 0, 0];
            
            // For each edge incident to this vertex
            for (const edgeKey of this.vertexEdges[vertexIdx]) {
                const edge = this.edgeMap.get(edgeKey);
                const vj_idx = edge.v0 === vertexIdx ? edge.v1 : edge.v0;
                const vj = this.vertices[vj_idx];
                
                // Sum cotangent weights from incident triangles
                let cotSum = 0;
                for (const fi of edge.faces) {
                    cotSum += this.cotanOpposite(fi, edge.v0, edge.v1);
                }
                
                // (v_j - v_i) weighted by cotangent sum
                const diff = [vj[0] - vi[0], vj[1] - vi[1], vj[2] - vi[2]];
                Hn[0] += cotSum * diff[0];
                Hn[1] += cotSum * diff[1];
                Hn[2] += cotSum * diff[2];
            }
            
            // Compute mixed area (Voronoi area for obtuse-safe)
            const area = this.mixedArea(vertexIdx);
            
            if (area > EPSILON) {
                Hn[0] /= (2 * area);
                Hn[1] /= (2 * area);
                Hn[2] /= (2 * area);
            }
            
            return Hn;
        }

        /**
         * Mean Curvature (scalar) = |Hn| / 2
         */
        meanCurvature(vertexIdx) {
            const Hn = this.meanCurvatureNormal(vertexIdx);
            return sqrt(Hn[0] * Hn[0] + Hn[1] * Hn[1] + Hn[2] * Hn[2]) / 2;
        }

        /**
         * Mixed Area (Voronoi region with obtuse handling)
         * From Meyer et al., "Discrete Differential-Geometry Operators"
         */
        mixedArea(vertexIdx) {
            let area = 0;
            
            for (const fi of this.vertexFaces[vertexIdx]) {
                const f = this.faces[fi];
                const idx = f.indexOf(vertexIdx);
                const v0 = f[idx];
                const v1 = f[(idx + 1) % 3];
                const v2 = f[(idx + 2) % 3];
                
                const p0 = this.vertices[v0];
                const p1 = this.vertices[v1];
                const p2 = this.vertices[v2];
                
                // Check for obtuse angles
                const angle0 = this.tipAngle(v0, v1, v2);
                const angle1 = this.tipAngle(v1, v2, v0);
                const angle2 = this.tipAngle(v2, v0, v1);
                
                const triArea = WedgeDeterminant.triangleArea(p0, p1, p2);
                
                if (angle0 > PI / 2) {
                    // Obtuse at our vertex: use half triangle area
                    area += triArea / 2;
                } else if (angle1 > PI / 2 || angle2 > PI / 2) {
                    // Obtuse at other vertex: use quarter triangle area
                    area += triArea / 4;
                } else {
                    // Non-obtuse: use Voronoi area
                    const cot1 = cos(angle1) / (sin(angle1) + EPSILON);
                    const cot2 = cos(angle2) / (sin(angle2) + EPSILON);
                    
                    const e1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
                    const e2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
                    const len1sq = e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2];
                    const len2sq = e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2];
                    
                    area += (cot2 * len1sq + cot1 * len2sq) / 8;
                }
            }
            
            return area;
        }

        /**
         * Principal Curvatures from Gaussian and Mean
         * κ₁, κ₂ = H ± √(H² - K)
         */
        principalCurvatures(vertexIdx) {
            const H = this.meanCurvature(vertexIdx);
            const K = this.gaussianCurvature(vertexIdx) / this.mixedArea(vertexIdx);
            
            const discriminant = H * H - K;
            if (discriminant < 0) {
                // Numerical error - return mean as both
                return [H, H];
            }
            
            const sqrtD = sqrt(discriminant);
            return [H + sqrtD, H - sqrtD];
        }

        // ========================================================================
        // TOPOLOGY
        // ========================================================================

        /**
         * Euler Characteristic χ = V - E + F
         */
        eulerCharacteristic() {
            return this.vertices.length - this.edges.length + this.faces.length;
        }

        /**
         * Verify Gauss-Bonnet: ∫K dA = 2πχ
         */
        verifyGaussBonnet() {
            const totalK = this.totalGaussianCurvature();
            const chi = this.eulerCharacteristic();
            const expected = 2 * PI * chi;
            return {
                totalGaussianCurvature: totalK,
                eulerCharacteristic: chi,
                expected: expected,
                error: abs(totalK - expected),
                verified: abs(totalK - expected) < 0.01
            };
        }

        // ========================================================================
        // VARIATIONAL APPROACH: DIHEDRAL ANGLES
        // ========================================================================

        /**
         * Dihedral angle at an edge (angle between face normals)
         * @param {Object} edge - Edge object with {v0, v1, faces}
         * @returns {number} Dihedral angle in radians [0, π]
         */
        dihedralAngleAtEdge(edge) {
            if (edge.faces.length !== 2) {
                return PI; // Boundary edge has "flat" dihedral
            }
            
            const f1 = this.faces[edge.faces[0]];
            const f2 = this.faces[edge.faces[1]];
            
            // Get the two normals
            const n1 = WedgeDeterminant.triangleNormal(
                this.vertices[f1[0]], this.vertices[f1[1]], this.vertices[f1[2]]
            );
            const n2 = WedgeDeterminant.triangleNormal(
                this.vertices[f2[0]], this.vertices[f2[1]], this.vertices[f2[2]]
            );
            
            // Normalize
            const len1 = sqrt(n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]);
            const len2 = sqrt(n2[0] * n2[0] + n2[1] * n2[1] + n2[2] * n2[2]);
            
            if (len1 < EPSILON || len2 < EPSILON) return 0;
            
            const dot = (n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]) / (len1 * len2);
            return acos(Math.max(-1, Math.min(1, dot)));
        }

        /**
         * Total Mean Curvature via Steiner's formula (edge-based)
         * H_total = (1/2) Σ_edges l_e × θ_e
         * 
         * where θ_e = π - (dihedral angle) is the exterior dihedral angle
         * 
         * NOTE: This is triangulation-dependent (unlike total Gaussian K)
         */
        totalMeanCurvatureSteiner() {
            let H = 0;
            
            for (const edge of this.edges) {
                const length = this.edgeLength(edge.v0, edge.v1);
                const dihedral = this.dihedralAngleAtEdge(edge);
                const exteriorAngle = PI - dihedral;
                
                H += 0.5 * length * exteriorAngle;
            }
            
            return H;
        }

        // ========================================================================
        // VARIATIONAL APPROACH: STEINER'S FORMULA
        // ========================================================================

        /**
         * Mesh volume using tetrahedra from origin
         */
        meshVolume() {
            let volume = 0;
            const origin = [0, 0, 0];
            
            for (const f of this.faces) {
                volume += WedgeDeterminant.tetraVolume(
                    origin,
                    this.vertices[f[0]],
                    this.vertices[f[1]],
                    this.vertices[f[2]]
                );
            }
            
            return abs(volume);
        }

        /**
         * Mesh surface area
         */
        meshArea() {
            let area = 0;
            for (const f of this.faces) {
                area += WedgeDeterminant.triangleArea(
                    this.vertices[f[0]],
                    this.vertices[f[1]],
                    this.vertices[f[2]]
                );
            }
            return area;
        }

        /**
         * Steiner's Formula coefficients
         * 
         * When a polyhedron P is "mollified" by Minkowski sum with ball B(t):
         *     V(P ⊕ B(t)) = V + A·t + H·t² + (4π/3)·t³
         * 
         * The coefficients are:
         *   - V: original volume
         *   - A: surface area
         *   - H: total mean curvature (edge-based)
         *   - K: total Gaussian curvature (= 4π for genus-0)
         * 
         * INSIGHT: Differentiating w.r.t. t recovers the hierarchy:
         *   dV/dt|₀ = A
         *   d²V/dt²|₀ = 2H
         *   d³V/dt³|₀ = 4πK (constant = 4π)
         * 
         * @returns {Object} {V, A, H, K, chi}
         */
        steinerCoefficients() {
            return {
                V: this.meshVolume(),
                A: this.meshArea(),
                H: this.totalMeanCurvatureSteiner(),
                K: this.totalGaussianCurvature(),
                chi: this.eulerCharacteristic()
            };
        }

        /**
         * Compute mollified volume at radius t using Steiner's formula
         * @param {number} t - Mollification radius
         * @returns {number} Volume of Minkowski sum with ball of radius t
         */
        mollifiedVolume(t) {
            const { V, A, H } = this.steinerCoefficients();
            return V + A * t + H * t * t + (4 * PI / 3) * t * t * t;
        }

        // ========================================================================
        // VARIATIONAL APPROACH: WILLMORE ENERGY
        // ========================================================================

        /**
         * Willmore Energy: W = Σ H² · A
         * 
         * Measures total bending of surface. Minimizing Willmore energy
         * produces "fair" surfaces with minimal bending variation.
         * 
         * For a sphere of radius r: W = 4π (independent of r!)
         * This is the minimum possible for any genus-0 surface.
         * 
         * @returns {number} Total Willmore energy
         */
        willmoreEnergy() {
            let W = 0;
            
            for (let i = 0; i < this.vertices.length; i++) {
                const H = this.meanCurvature(i);
                const A = this.mixedArea(i);
                W += H * H * A;
            }
            
            return W;
        }

        // ========================================================================
        // VARIATIONAL APPROACH: CURVATURE FLOWS
        // ========================================================================

        /**
         * Mean Curvature Flow: df/dt = -H⃗
         * 
         * Moves each vertex in the direction opposite its mean curvature normal.
         * This minimizes surface area (gradient descent on area functional).
         * 
         * Properties:
         * - Convex surfaces shrink to points
         * - Removes noise/high-frequency detail
         * - Can develop singularities (necks pinch off)
         * 
         * @param {number} dt - Time step (small, e.g., 0.001-0.01)
         * @returns {Array} New vertex positions
         */
        meanCurvatureFlowStep(dt = 0.01) {
            const newVertices = [];
            
            for (let i = 0; i < this.vertices.length; i++) {
                const Hn = this.meanCurvatureNormal(i);
                const v = this.vertices[i];
                
                // Move in negative gradient direction (towards minimum area)
                newVertices.push([
                    v[0] - dt * Hn[0],
                    v[1] - dt * Hn[1],
                    v[2] - dt * Hn[2]
                ]);
            }
            
            return newVertices;
        }

        /**
         * Apply mean curvature flow in-place
         * @param {number} dt - Time step
         * @param {number} steps - Number of iterations
         */
        applyMeanCurvatureFlow(dt = 0.01, steps = 1) {
            for (let s = 0; s < steps; s++) {
                this.vertices = this.meanCurvatureFlowStep(dt);
                this._buildTopology();
            }
            return this;
        }

        /**
         * Gaussian Curvature Flow: df/dt = -K·n
         * 
         * Moves vertices proportional to their Gaussian curvature.
         * 
         * @param {number} dt - Time step
         * @returns {Array} New vertex positions
         */
        gaussCurvatureFlowStep(dt = 0.01) {
            const newVertices = [];
            
            for (let i = 0; i < this.vertices.length; i++) {
                const K = this.gaussianCurvature(i);
                const n = this.vertexNormal(i);
                const v = this.vertices[i];
                
                newVertices.push([
                    v[0] - dt * K * n[0],
                    v[1] - dt * K * n[1],
                    v[2] - dt * K * n[2]
                ]);
            }
            
            return newVertices;
        }

        /**
         * Apply Gaussian curvature flow in-place
         */
        applyGaussCurvatureFlow(dt = 0.01, steps = 1) {
            for (let s = 0; s < steps; s++) {
                this.vertices = this.gaussCurvatureFlowStep(dt);
                this._buildTopology();
            }
            return this;
        }

        // ========================================================================
        // VARIATIONAL APPROACH: GRADIENT OF AREA
        // ========================================================================

        /**
         * Gradient of Triangle Area w.r.t vertex position
         * 
         * ∇_{vi} Area = (1/2) × (n × e_opposite)
         * 
         * where n is the face normal and e_opposite is the edge opposite to vi
         * 
         * @param {number} faceIdx - Face index
         * @param {number} vertexInFace - Which vertex (0, 1, or 2)
         * @returns {Array} [gx, gy, gz] gradient vector
         */
        gradTriangleArea(faceIdx, vertexInFace) {
            const f = this.faces[faceIdx];
            const v0 = this.vertices[f[0]];
            const v1 = this.vertices[f[1]];
            const v2 = this.vertices[f[2]];
            
            // Face normal (unnormalized, magnitude = 2 × area)
            const n = WedgeDeterminant.triangleNormal(v0, v1, v2);
            const len = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
            if (len < EPSILON) return [0, 0, 0];
            
            // Normalize
            n[0] /= len;
            n[1] /= len;
            n[2] /= len;
            
            // Edge opposite to vertex i
            let e;
            if (vertexInFace === 0) {
                e = [v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]];
            } else if (vertexInFace === 1) {
                e = [v0[0] - v2[0], v0[1] - v2[1], v0[2] - v2[2]];
            } else {
                e = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
            }
            
            // Gradient = (1/2) × (n × e)
            return [
                0.5 * (n[1] * e[2] - n[2] * e[1]),
                0.5 * (n[2] * e[0] - n[0] * e[2]),
                0.5 * (n[0] * e[1] - n[1] * e[0])
            ];
        }

        // ========================================================================
        // SCHLÄFLI FORMULA
        // ========================================================================

        /**
         * Schläfli's Formula: Σ l_e × dθ_e = 0
         * 
         * For any infinitesimal deformation of a polyhedron, the sum of
         * (edge length) × (change in dihedral angle) equals zero.
         * 
         * This is a fundamental constraint on polyhedral deformations.
         * 
         * @param {Array} dThetas - Array of dihedral angle changes for each edge
         * @returns {number} The Schläfli sum (should be zero for valid deformations)
         */
        schlafliSum(dThetas) {
            let sum = 0;
            
            for (let i = 0; i < this.edges.length; i++) {
                const edge = this.edges[i];
                const length = this.edgeLength(edge.v0, edge.v1);
                sum += length * (dThetas[i] || 0);
            }
            
            return sum;
        }

        /**
         * Get all current dihedral angles
         * @returns {Array} Dihedral angles for each edge
         */
        getAllDihedralAngles() {
            return this.edges.map(edge => this.dihedralAngleAtEdge(edge));
        }

        // ========================================================================
        // MESH UTILITIES
        // ========================================================================

        /**
         * Mesh center of mass (centroid)
         */
        meshCentroid() {
            let cx = 0, cy = 0, cz = 0;
            let totalVolume = 0;
            const origin = [0, 0, 0];
            
            for (const f of this.faces) {
                const v0 = this.vertices[f[0]];
                const v1 = this.vertices[f[1]];
                const v2 = this.vertices[f[2]];
                
                const vol = WedgeDeterminant.tetraVolume(origin, v0, v1, v2);
                // Centroid of tetrahedron is average of 4 vertices
                const tx = (v0[0] + v1[0] + v2[0]) / 4;
                const ty = (v0[1] + v1[1] + v2[1]) / 4;
                const tz = (v0[2] + v1[2] + v2[2]) / 4;
                
                cx += vol * tx;
                cy += vol * ty;
                cz += vol * tz;
                totalVolume += vol;
            }
            
            if (abs(totalVolume) < EPSILON) return [0, 0, 0];
            return [cx / totalVolume, cy / totalVolume, cz / totalVolume];
        }

        /**
         * Clone the mesh
         */
        clone() {
            const copy = new DiscreteCurvature();
            copy.loadMesh(this.vertices, this.faces);
            return copy;
        }
    }

    // ============================================================================
    // CGA DISCRETE CURVATURE
    // ============================================================================

    /**
     * CGA-based discrete curvature calculations
     * 
     * Key idea: In CGA, circles and spheres are first-class objects.
     * The osculating circle at a curve point encodes curvature directly.
     */
    class CGACurvature {
        constructor(cgaAlgebra) {
            this.cga = cgaAlgebra;
        }

        point(x, y, z) {
            return this.cga.point(x, y, z);
        }

        circleThrough3(p1, p2, p3) {
            const P1 = this.point(p1[0], p1[1], p1[2]);
            const P2 = this.point(p2[0], p2[1], p2[2]);
            const P3 = this.point(p3[0], p3[1], p3[2]);
            return P1.wedge(P2).wedge(P3);
        }

        sphereThrough4(p1, p2, p3, p4) {
            const P1 = this.point(p1[0], p1[1], p1[2]);
            const P2 = this.point(p2[0], p2[1], p2[2]);
            const P3 = this.point(p3[0], p3[1], p3[2]);
            const P4 = this.point(p4[0], p4[1], p4[2]);
            return P1.wedge(P2).wedge(P3).wedge(P4);
        }

        osculatingCircleCurvature(p0, p1, p2) {
            const d1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
            const d2 = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
            
            const cross = [
                d1[1] * d2[2] - d1[2] * d2[1],
                d1[2] * d2[0] - d1[0] * d2[2],
                d1[0] * d2[1] - d1[1] * d2[0]
            ];
            
            const crossMag = sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);
            const d1Mag = sqrt(d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]);
            
            if (d1Mag < EPSILON) return 0;
            return 2 * crossMag / (d1Mag * d1Mag * d1Mag);
        }
    }

    // ============================================================================
    // PGA DISCRETE CURVATURE
    // ============================================================================

    /**
     * PGA-based discrete curvature using join operations
     */
    class PGACurvature {
        constructor(pgaAlgebra) {
            this.pga = pgaAlgebra;
        }

        point(x, y, z) {
            return this.pga.point(x, y, z);
        }

        edgeLength(p0, p1) {
            const dx = p1[0] - p0[0];
            const dy = p1[1] - p0[1];
            const dz = p1[2] - p0[2];
            return sqrt(dx * dx + dy * dy + dz * dz);
        }

        triangleArea(p0, p1, p2) {
            return WedgeDeterminant.triangleArea(p0, p1, p2);
        }

        tetraVolume(p0, p1, p2, p3) {
            return WedgeDeterminant.tetraVolume(p0, p1, p2, p3);
        }

        triangleNormal(p0, p1, p2) {
            const n = WedgeDeterminant.triangleNormal(p0, p1, p2);
            const len = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
            if (len < EPSILON) return [0, 0, 1];
            return [n[0] / len, n[1] / len, n[2] / len];
        }

        meshVolume(faces, vertices) {
            let volume = 0;
            const origin = [0, 0, 0];
            
            for (const f of faces) {
                volume += this.tetraVolume(origin, vertices[f[0]], vertices[f[1]], vertices[f[2]]);
            }
            
            return abs(volume);
        }

        meshCentroid(faces, vertices) {
            let cx = 0, cy = 0, cz = 0;
            let totalVolume = 0;
            const origin = [0, 0, 0];
            
            for (const f of faces) {
                const v0 = vertices[f[0]];
                const v1 = vertices[f[1]];
                const v2 = vertices[f[2]];
                
                const vol = this.tetraVolume(origin, v0, v1, v2);
                const tx = (v0[0] + v1[0] + v2[0]) / 4;
                const ty = (v0[1] + v1[1] + v2[1]) / 4;
                const tz = (v0[2] + v1[2] + v2[2]) / 4;
                
                cx += vol * tx;
                cy += vol * ty;
                cz += vol * tz;
                totalVolume += vol;
            }
            
            if (abs(totalVolume) < EPSILON) return [0, 0, 0];
            return [cx / totalVolume, cy / totalVolume, cz / totalVolume];
        }
    }

    // ============================================================================
    // MESH PRIMITIVES FOR TESTING
    // ============================================================================

    const Primitives = {
        /**
         * Unit sphere (icosahedron subdivision)
         */
        icosphere(subdivisions = 2) {
            const t = (1 + sqrt(5)) / 2;
            
            let vertices = [
                [-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0],
                [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
                [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]
            ].map(v => {
                const len = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
                return [v[0] / len, v[1] / len, v[2] / len];
            });
            
            let faces = [
                [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
                [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
                [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
                [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]
            ];
            
            for (let s = 0; s < subdivisions; s++) {
                const newFaces = [];
                const midpointCache = new Map();
                
                const getMidpoint = (i1, i2) => {
                    const key = i1 < i2 ? `${i1}-${i2}` : `${i2}-${i1}`;
                    if (midpointCache.has(key)) return midpointCache.get(key);
                    
                    const v1 = vertices[i1], v2 = vertices[i2];
                    const mid = [(v1[0] + v2[0]) / 2, (v1[1] + v2[1]) / 2, (v1[2] + v2[2]) / 2];
                    const len = sqrt(mid[0] * mid[0] + mid[1] * mid[1] + mid[2] * mid[2]);
                    vertices.push([mid[0] / len, mid[1] / len, mid[2] / len]);
                    const idx = vertices.length - 1;
                    midpointCache.set(key, idx);
                    return idx;
                };
                
                for (const f of faces) {
                    const a = getMidpoint(f[0], f[1]);
                    const b = getMidpoint(f[1], f[2]);
                    const c = getMidpoint(f[2], f[0]);
                    newFaces.push([f[0], a, c], [f[1], b, a], [f[2], c, b], [a, b, c]);
                }
                faces = newFaces;
            }
            
            return { vertices, faces };
        },

        /**
         * Cube [-1, 1]³
         */
        cube() {
            const vertices = [
                [-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]
            ];
            
            const faces = [
                [0, 3, 2], [0, 2, 1],  // back
                [4, 5, 6], [4, 6, 7],  // front
                [0, 1, 5], [0, 5, 4],  // bottom
                [2, 3, 7], [2, 7, 6],  // top
                [0, 4, 7], [0, 7, 3],  // left
                [1, 2, 6], [1, 6, 5]   // right
            ];
            
            return { vertices, faces };
        },

        /**
         * Torus (major radius R, minor radius r)
         */
        torus(R = 2, r = 0.5, nu = 24, nv = 12) {
            const vertices = [];
            const faces = [];
            
            for (let i = 0; i < nu; i++) {
                const u = 2 * PI * i / nu;
                for (let j = 0; j < nv; j++) {
                    const v = 2 * PI * j / nv;
                    vertices.push([
                        (R + r * cos(v)) * cos(u),
                        (R + r * cos(v)) * sin(u),
                        r * sin(v)
                    ]);
                }
            }
            
            for (let i = 0; i < nu; i++) {
                for (let j = 0; j < nv; j++) {
                    const i1 = (i + 1) % nu;
                    const j1 = (j + 1) % nv;
                    const v0 = i * nv + j;
                    const v1 = i1 * nv + j;
                    const v2 = i1 * nv + j1;
                    const v3 = i * nv + j1;
                    faces.push([v0, v1, v2], [v0, v2, v3]);
                }
            }
            
            return { vertices, faces };
        },

        /**
         * Tetrahedron (regular)
         */
        tetrahedron() {
            const vertices = [
                [1, 1, 1],
                [1, -1, -1],
                [-1, 1, -1],
                [-1, -1, 1]
            ].map(v => {
                const len = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
                return [v[0] / len, v[1] / len, v[2] / len];
            });
            
            const faces = [
                [0, 1, 2],
                [0, 2, 3],
                [0, 3, 1],
                [1, 3, 2]
            ];
            
            return { vertices, faces };
        },

        /**
         * Octahedron (regular)
         */
        octahedron() {
            const vertices = [
                [1, 0, 0], [-1, 0, 0],
                [0, 1, 0], [0, -1, 0],
                [0, 0, 1], [0, 0, -1]
            ];
            
            const faces = [
                [0, 2, 4], [0, 4, 3], [0, 3, 5], [0, 5, 2],
                [1, 4, 2], [1, 3, 4], [1, 5, 3], [1, 2, 5]
            ];
            
            return { vertices, faces };
        }
    };

    // ============================================================================
    // EXPORTS
    // ============================================================================

    const DDG = {
        WedgeDeterminant,
        DiscreteCurvature,
        CGACurvature,
        PGACurvature,
        Primitives,
        EPSILON,
        PI
    };

    if (typeof module !== 'undefined' && module.exports) {
        module.exports = DDG;
    } else if (typeof define === 'function' && define.amd) {
        define([], () => DDG);
    } else {
        global.DDG = DDG;
    }

})(typeof window !== 'undefined' ? window : global);
