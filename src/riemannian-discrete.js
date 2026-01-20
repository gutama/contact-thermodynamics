/**
 * Discrete Riemannian Geometry on Triangle Meshes
 * 
 * Implements coordinate-free Riemannian geometry for triangle meshes:
 *   - Connection bivector ω (dihedral rotation at edges)
 *   - Curvature 2-form Ω (angle defect at vertices)
 *   - Parallel transport across edges
 *   - Bianchi identity verification (Gauss-Bonnet theorem)
 * 
 * Key correspondences:
 *   Continuous          Discrete                 Mesh Element
 *   ωᵢ                  Dihedral rotation        Edge
 *   Ω                   Angle defect (2π - Σθ)   Vertex
 *   K (Gaussian)        2π - Σ interior angles   Vertex
 *   ∇ ∧ Ω = 0          Σ K_v = 2π χ(M)          Global (Gauss-Bonnet)
 * 
 * @module riemannian-discrete
 * @license MIT
 */

(function (global) {
    'use strict';

    const EPSILON = 1e-10;
    const { abs, sqrt, sin, cos, PI, acos, atan2 } = Math;

    // ============================================================================
    // VECTOR/BIVECTOR UTILITIES
    // ============================================================================

    function dot(u, v) {
        return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    }

    function cross(u, v) {
        return [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0]
        ];
    }

    function norm(v) {
        return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }

    function normalize(v) {
        const n = norm(v);
        return n > EPSILON ? [v[0] / n, v[1] / n, v[2] / n] : [0, 0, 0];
    }

    function vecAdd(u, v) {
        return [u[0] + v[0], u[1] + v[1], u[2] + v[2]];
    }

    function vecSub(u, v) {
        return [u[0] - v[0], u[1] - v[1], u[2] - v[2]];
    }

    function vecScale(v, s) {
        return [v[0] * s, v[1] * s, v[2] * s];
    }

    /**
     * Bivector in 3D (three components: e₂₃, e₃₁, e₁₂)
     * Represents oriented planes / rotations.
     */
    class Bivector3D {
        constructor(e23 = 0, e31 = 0, e12 = 0) {
            this.e23 = e23;
            this.e31 = e31;
            this.e12 = e12;
        }

        static fromWedge(u, v) {
            // u ∧ v = (u₂v₃ - u₃v₂) e₂₃ + (u₃v₁ - u₁v₃) e₃₁ + (u₁v₂ - u₂v₁) e₁₂
            return new Bivector3D(
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0]
            );
        }

        /**
         * Create bivector from axis-angle representation.
         * B = θ * (axis as 2-blade dual)
         * For rotation around axis n by angle θ: B = θ * I·n where I is pseudoscalar
         */
        static fromAxisAngle(axis, angle) {
            // Dual of axis vector is the rotation plane bivector
            // In 3D: dual(e₁) = e₂₃, dual(e₂) = e₃₁, dual(e₃) = e₁₂
            return new Bivector3D(
                angle * axis[0],  // e₂₃ component
                angle * axis[1],  // e₃₁ component
                angle * axis[2]   // e₁₂ component
            );
        }

        toArray() {
            return [this.e23, this.e31, this.e12];
        }

        /**
         * Get the axis of rotation (dual of bivector).
         */
        axis() {
            return normalize([this.e23, this.e31, this.e12]);
        }

        /**
         * Get the rotation angle (magnitude of bivector).
         */
        angle() {
            return sqrt(this.e23 * this.e23 + this.e31 * this.e31 + this.e12 * this.e12);
        }

        add(other) {
            return new Bivector3D(
                this.e23 + other.e23,
                this.e31 + other.e31,
                this.e12 + other.e12
            );
        }

        scale(s) {
            return new Bivector3D(this.e23 * s, this.e31 * s, this.e12 * s);
        }

        neg() {
            return new Bivector3D(-this.e23, -this.e31, -this.e12);
        }

        normSquared() {
            return this.e23 * this.e23 + this.e31 * this.e31 + this.e12 * this.e12;
        }

        norm() {
            return sqrt(this.normSquared());
        }

        /**
         * Compute B × v = ½(Bv - vB) in 3D.
         * The bivector acts on vectors by rotation in its plane.
         * Equivalent to cross product with the dual vector.
         */
        commutatorWithVector(v) {
            // Dual of bivector is axial vector: B* = (e₂₃, e₃₁, e₁₂)
            return cross([this.e23, this.e31, this.e12], v);
        }

        toString() {
            return `Bivector3D(${this.e23.toFixed(4)} e₂₃ + ${this.e31.toFixed(4)} e₃₁ + ${this.e12.toFixed(4)} e₁₂)`;
        }
    }

    // ============================================================================
    // MESH CONNECTION BIVECTOR
    // ============================================================================

    /**
     * Discrete connection bivector on a triangle mesh.
     * 
     * The connection bivector at each edge encodes the rotation needed to
     * parallel transport a tangent vector from one face to the adjacent face.
     * 
     * For edge e with adjacent faces f₁, f₂:
     *   ω_e = θ_e · B_e
     * where θ_e is the dihedral angle and B_e is the rotation plane (edge ∧ avg_normal).
     */
    class MeshConnectionBivector {
        /**
         * @param {TriangleMesh} mesh 
         */
        constructor(mesh) {
            this.mesh = mesh;
            this._connectionBivectors = null;
        }

        /**
         * Compute connection bivectors at all edges.
         * @returns {Bivector3D[]} Array of bivectors, one per edge
         */
        compute() {
            if (this._connectionBivectors) return this._connectionBivectors;

            const mesh = this.mesh;
            const nE = mesh.nEdges;
            const dihedrals = mesh.dihedralAngles();
            const edgeVecs = mesh.edgeVectors();

            const bivectors = new Array(nE);

            for (let e = 0; e < nE; e++) {
                // Dihedral angle (deviation from flat = π)
                const theta = dihedrals[e] - PI;  // 0 for flat, positive for convex

                // Edge direction as rotation axis
                const axis = [
                    edgeVecs[e * 3],
                    edgeVecs[e * 3 + 1],
                    edgeVecs[e * 3 + 2]
                ];

                // Connection bivector: encodes rotation around edge by dihedral angle
                bivectors[e] = Bivector3D.fromAxisAngle(axis, theta);
            }

            this._connectionBivectors = bivectors;
            return bivectors;
        }

        /**
         * Get connection bivector at a specific edge.
         * @param {number} edgeIdx 
         * @returns {Bivector3D}
         */
        at(edgeIdx) {
            return this.compute()[edgeIdx];
        }

        /**
         * Compute the rotation operator (rotor) for an edge.
         * R = exp(-ω/2) = cos(θ/2) - sin(θ/2) * B̂
         * 
         * @param {number} edgeIdx 
         * @returns {Object} {scalar, bivector} components of rotor
         */
        rotor(edgeIdx) {
            const omega = this.at(edgeIdx);
            const angle = omega.angle();

            if (angle < EPSILON) {
                return { scalar: 1, bivector: new Bivector3D(0, 0, 0) };
            }

            const halfAngle = angle / 2;
            const cosHalf = cos(halfAngle);
            const sinHalf = sin(halfAngle);

            // Normalized bivector
            const Bhat = omega.scale(1 / angle);

            return {
                scalar: cosHalf,
                bivector: Bhat.scale(-sinHalf)
            };
        }

        /**
         * Invalidate cache.
         */
        invalidate() {
            this._connectionBivectors = null;
        }
    }

    // ============================================================================
    // MESH CURVATURE 2-FORM
    // ============================================================================

    /**
     * Discrete curvature 2-form on a triangle mesh.
     * 
     * Curvature is measured as angle defect at vertices:
     *   K_v = 2π - Σⱼ θⱼ
     * where θⱼ are the interior angles of triangles at vertex v.
     * 
     * This is the discrete analog of Ω = dω + ω ∧ ω.
     */
    class MeshCurvature2Form {
        /**
         * @param {TriangleMesh} mesh 
         */
        constructor(mesh) {
            this.mesh = mesh;
            this._angleDefects = null;
            this._interiorAngles = null;
        }

        /**
         * Compute interior angles at each vertex of each face.
         * Returns array of size nFaces * 3.
         * @returns {Float64Array}
         */
        _computeInteriorAngles() {
            if (this._interiorAngles) return this._interiorAngles;

            const mesh = this.mesh;
            const nF = mesh.nFaces;
            const angles = new Float64Array(nF * 3);
            const verts = mesh.vertices;
            const faces = mesh.faces;

            for (let t = 0; t < nF; t++) {
                const i = faces[t * 3];
                const j = faces[t * 3 + 1];
                const k = faces[t * 3 + 2];

                // Vertex positions
                const pi = [verts[i * 3], verts[i * 3 + 1], verts[i * 3 + 2]];
                const pj = [verts[j * 3], verts[j * 3 + 1], verts[j * 3 + 2]];
                const pk = [verts[k * 3], verts[k * 3 + 1], verts[k * 3 + 2]];

                // Edges from each vertex
                const eij = vecSub(pj, pi);
                const eik = vecSub(pk, pi);
                const eji = vecSub(pi, pj);
                const ejk = vecSub(pk, pj);
                const eki = vecSub(pi, pk);
                const ekj = vecSub(pj, pk);

                // Interior angles via dot product
                // cos(θ) = (a · b) / (|a| |b|)
                const cosI = dot(eij, eik) / (norm(eij) * norm(eik) + EPSILON);
                const cosJ = dot(eji, ejk) / (norm(eji) * norm(ejk) + EPSILON);
                const cosK = dot(eki, ekj) / (norm(eki) * norm(ekj) + EPSILON);

                // Clamp to [-1, 1] for numerical stability
                angles[t * 3] = acos(Math.max(-1, Math.min(1, cosI)));
                angles[t * 3 + 1] = acos(Math.max(-1, Math.min(1, cosJ)));
                angles[t * 3 + 2] = acos(Math.max(-1, Math.min(1, cosK)));
            }

            this._interiorAngles = angles;
            return angles;
        }

        /**
         * Compute angle defect (discrete Gaussian curvature) at each vertex.
         * K_v = 2π - Σ(interior angles at v)
         * For boundary vertices: K_v = π - Σ(interior angles at v)
         * 
         * @returns {Float64Array} Angle defect per vertex (in radians)
         */
        angleDefects() {
            if (this._angleDefects) return this._angleDefects;

            const mesh = this.mesh;
            const nV = mesh.nVertices;
            const nF = mesh.nFaces;
            const faces = mesh.faces;
            const angles = this._computeInteriorAngles();
            const defects = new Float64Array(nV);

            // Start with 2π for interior vertices, π for boundary
            for (let v = 0; v < nV; v++) {
                defects[v] = mesh.boundaryVertices[v] ? PI : 2 * PI;
            }

            // Subtract interior angles
            for (let t = 0; t < nF; t++) {
                const i = faces[t * 3];
                const j = faces[t * 3 + 1];
                const k = faces[t * 3 + 2];

                defects[i] -= angles[t * 3];
                defects[j] -= angles[t * 3 + 1];
                defects[k] -= angles[t * 3 + 2];
            }

            this._angleDefects = defects;
            return defects;
        }

        /**
         * Get Gaussian curvature at a vertex.
         * This is the angle defect (pointwise, not averaged over area).
         * 
         * @param {number} vertexIdx 
         * @returns {number} Angle defect in radians
         */
        gaussianCurvature(vertexIdx) {
            return this.angleDefects()[vertexIdx];
        }

        /**
         * Compute total Gaussian curvature (sum of all angle defects).
         * By Gauss-Bonnet: Σ K_v = 2π χ(M)
         * 
         * @returns {number} Total curvature in radians
         */
        totalCurvature() {
            const defects = this.angleDefects();
            let total = 0;
            for (let i = 0; i < defects.length; i++) {
                total += defects[i];
            }
            return total;
        }

        /**
         * Compute Euler characteristic from mesh topology.
         * χ = V - E + F
         * 
         * @returns {number}
         */
        eulerCharacteristic() {
            const mesh = this.mesh;
            return mesh.nVertices - mesh.nEdges + mesh.nFaces;
        }

        /**
         * Get curvature as a bivector at a vertex.
         * The curvature bivector lies in the tangent plane.
         * 
         * @param {number} vertexIdx 
         * @returns {Bivector3D} Curvature bivector
         */
        curvatureBivector(vertexIdx) {
            const K = this.gaussianCurvature(vertexIdx);

            // Compute average normal at vertex (area-weighted)
            const mesh = this.mesh;
            const vertTris = mesh.vertexTriangles[vertexIdx];
            const normals = mesh.faceNormals();
            const areas = mesh.faceAreas();

            let nx = 0, ny = 0, nz = 0;
            let totalArea = 0;

            for (let i = 0; i < vertTris.length; i++) {
                const t = vertTris[i];
                const a = areas[t];
                nx += normals[t * 3] * a;
                ny += normals[t * 3 + 1] * a;
                nz += normals[t * 3 + 2] * a;
                totalArea += a;
            }

            if (totalArea > EPSILON) {
                nx /= totalArea;
                ny /= totalArea;
                nz /= totalArea;
                const len = sqrt(nx * nx + ny * ny + nz * nz);
                if (len > EPSILON) {
                    nx /= len;
                    ny /= len;
                    nz /= len;
                }
            }

            // Curvature bivector is K times the tangent pseudoscalar
            // In 3D embedded case, we use the dual of the normal: I·n
            return new Bivector3D(K * nx, K * ny, K * nz);
        }

        /**
         * Invalidate cache.
         */
        invalidate() {
            this._angleDefects = null;
            this._interiorAngles = null;
        }
    }

    // ============================================================================
    // MESH PARALLEL TRANSPORT
    // ============================================================================

    /**
     * Parallel transport of tangent vectors on a triangle mesh.
     * 
     * Uses the connection bivector (dihedral rotation) to transport vectors
     * across edges from one face to an adjacent face.
     */
    class MeshParallelTransport {
        /**
         * @param {TriangleMesh} mesh 
         */
        constructor(mesh) {
            this.mesh = mesh;
            this.connection = new MeshConnectionBivector(mesh);
        }

        /**
         * Transport a tangent vector across an edge from face t0 to face t1.
         * 
         * The vector is rotated by the dihedral angle around the edge axis.
         * 
         * @param {number[]} vector - Tangent vector in 3D (will be projected)
         * @param {number} edgeIdx - Edge to cross
         * @param {number} fromFace - Source face index
         * @returns {number[]} Transported vector
         */
        transportAcrossEdge(vector, edgeIdx, fromFace) {
            const mesh = this.mesh;
            const t0 = mesh.edgeTriangles[edgeIdx * 2];
            const t1 = mesh.edgeTriangles[edgeIdx * 2 + 1];

            // Boundary edge - no transport
            if (t1 < 0) return vector.slice();

            // Determine direction of transport
            const toFace = (fromFace === t0) ? t1 : t0;
            const sign = (fromFace === t0) ? 1 : -1;

            // Get edge vector (rotation axis)
            const edgeVecs = mesh.edgeVectors();
            const axis = [
                edgeVecs[edgeIdx * 3],
                edgeVecs[edgeIdx * 3 + 1],
                edgeVecs[edgeIdx * 3 + 2]
            ];

            // Get dihedral angle
            const dihedrals = mesh.dihedralAngles();
            const theta = (dihedrals[edgeIdx] - PI) * sign;

            // Rotate vector around axis by angle theta
            return this._rotateAroundAxis(vector, axis, theta);
        }

        /**
         * Rotate a vector around an axis by a given angle.
         * Uses Rodrigues' rotation formula:
         *   v' = v cos(θ) + (k × v) sin(θ) + k (k · v)(1 - cos(θ))
         * 
         * @param {number[]} v - Vector to rotate
         * @param {number[]} k - Unit rotation axis
         * @param {number} theta - Rotation angle
         * @returns {number[]} Rotated vector
         */
        _rotateAroundAxis(v, k, theta) {
            const cosT = cos(theta);
            const sinT = sin(theta);
            const kdotv = dot(k, v);
            const kcrossv = cross(k, v);

            return [
                v[0] * cosT + kcrossv[0] * sinT + k[0] * kdotv * (1 - cosT),
                v[1] * cosT + kcrossv[1] * sinT + k[1] * kdotv * (1 - cosT),
                v[2] * cosT + kcrossv[2] * sinT + k[2] * kdotv * (1 - cosT)
            ];
        }

        /**
         * Transport a vector along a path of edges.
         * 
         * @param {number[]} vector - Initial tangent vector
         * @param {Array<{edge: number, fromFace: number}>} path - Sequence of edge crossings
         * @returns {number[]} Final transported vector
         */
        transportAlongPath(vector, path) {
            let v = vector.slice();
            for (const step of path) {
                v = this.transportAcrossEdge(v, step.edge, step.fromFace);
            }
            return v;
        }

        /**
         * Compute holonomy around a vertex by transporting a vector around its 1-ring.
         * 
         * The holonomy angle equals the angle defect (Gaussian curvature at vertex).
         * 
         * @param {number} vertexIdx - Vertex to go around
         * @param {number[]} vector - Initial tangent vector
         * @returns {Object} {finalVector, angle, defect}
         */
        holonomyAroundVertex(vertexIdx, vector) {
            const mesh = this.mesh;
            const vertTris = mesh.vertexTriangles[vertexIdx];

            if (vertTris.length === 0) {
                return { finalVector: vector.slice(), angle: 0, defect: 0 };
            }

            // Build path around vertex (sequence of face transitions)
            const path = [];

            // Find ordered ring of faces around vertex
            const orderedFaces = this._orderFacesAroundVertex(vertexIdx);

            // For each consecutive pair of faces, find the shared edge
            for (let i = 0; i < orderedFaces.length; i++) {
                const f0 = orderedFaces[i];
                const f1 = orderedFaces[(i + 1) % orderedFaces.length];

                // Find shared edge between f0 and f1
                const sharedEdge = this._findSharedEdge(f0, f1);
                if (sharedEdge >= 0) {
                    path.push({ edge: sharedEdge, fromFace: f0 });
                }
            }

            // Transport vector around the loop
            const finalVector = this.transportAlongPath(vector, path);

            // Compute rotation angle between initial and final vectors
            const v0 = normalize(vector);
            const v1 = normalize(finalVector);
            const dotProduct = dot(v0, v1);
            const angle = acos(Math.max(-1, Math.min(1, dotProduct)));

            // Get theoretical defect
            const curvature = new MeshCurvature2Form(mesh);
            const defect = curvature.gaussianCurvature(vertexIdx);

            return { finalVector, angle, defect };
        }

        /**
         * Order faces around a vertex in consistent winding.
         * @private
         */
        _orderFacesAroundVertex(vertexIdx) {
            const mesh = this.mesh;
            const vertTris = mesh.vertexTriangles[vertexIdx];

            if (vertTris.length <= 1) return Array.from(vertTris);

            // Start with first face
            const ordered = [vertTris[0]];
            const used = new Set([vertTris[0]]);

            // Greedily add adjacent faces
            while (ordered.length < vertTris.length) {
                const lastFace = ordered[ordered.length - 1];
                let foundNext = false;

                for (let i = 0; i < vertTris.length; i++) {
                    const candidate = vertTris[i];
                    if (used.has(candidate)) continue;

                    // Check if candidate shares an edge with lastFace (not at vertexIdx)
                    const sharedEdge = this._findSharedEdge(lastFace, candidate);
                    if (sharedEdge >= 0) {
                        ordered.push(candidate);
                        used.add(candidate);
                        foundNext = true;
                        break;
                    }
                }

                if (!foundNext) break;  // Boundary vertex or non-manifold
            }

            return ordered;
        }

        /**
         * Find the edge shared between two faces.
         * @private
         */
        _findSharedEdge(face0, face1) {
            const mesh = this.mesh;
            const faces = mesh.faces;

            // Get vertices of each face
            const verts0 = [faces[face0 * 3], faces[face0 * 3 + 1], faces[face0 * 3 + 2]];
            const verts1 = [faces[face1 * 3], faces[face1 * 3 + 1], faces[face1 * 3 + 2]];

            // Find common vertices
            const common = verts0.filter(v => verts1.includes(v));

            if (common.length === 2) {
                return mesh.getEdgeIndex(common[0], common[1]);
            }

            return -1;
        }
    }

    // ============================================================================
    // BIANCHI IDENTITY VERIFIER
    // ============================================================================

    /**
     * Verifies the Bianchi identity ∇ ∧ Ω = 0 on triangle meshes.
     * 
     * The discrete Bianchi identity is the Gauss-Bonnet theorem:
     *   Σᵥ K_v = 2π χ(M)
     * 
     * where K_v is the angle defect at vertex v and χ(M) is the Euler characteristic.
     */
    class BianchiIdentityVerifier {
        /**
         * @param {TriangleMesh} mesh 
         */
        constructor(mesh) {
            this.mesh = mesh;
            this.curvature = new MeshCurvature2Form(mesh);
        }

        /**
         * Verify the Bianchi identity (Gauss-Bonnet theorem).
         * 
         * Checks that: Σ K_v ≈ 2π χ(M)
         * 
         * @returns {Object} Verification result:
         *   - passed: boolean - Whether identity holds within tolerance
         *   - totalCurvature: number - Sum of angle defects
         *   - expectedCurvature: number - 2π χ(M)
         *   - error: number - Absolute difference
         *   - relativeError: number - Relative difference
         *   - eulerCharacteristic: number - χ = V - E + F
         */
        verify(tolerance = 1e-6) {
            const totalK = this.curvature.totalCurvature();
            const chi = this.curvature.eulerCharacteristic();
            const expected = 2 * PI * chi;
            const error = abs(totalK - expected);
            const relativeError = abs(expected) > EPSILON ? error / abs(expected) : error;

            return {
                passed: error < tolerance,
                totalCurvature: totalK,
                expectedCurvature: expected,
                error: error,
                relativeError: relativeError,
                eulerCharacteristic: chi
            };
        }

        /**
         * Detailed verification with per-vertex breakdown.
         * 
         * @returns {Object} Detailed results including per-vertex curvatures
         */
        verifyDetailed() {
            const mesh = this.mesh;
            const defects = this.curvature.angleDefects();
            const basic = this.verify();

            // Statistics
            let minK = Infinity, maxK = -Infinity;
            let positiveCount = 0, negativeCount = 0, zeroCount = 0;

            for (let v = 0; v < mesh.nVertices; v++) {
                const K = defects[v];
                minK = Math.min(minK, K);
                maxK = Math.max(maxK, K);
                if (abs(K) < EPSILON) zeroCount++;
                else if (K > 0) positiveCount++;
                else negativeCount++;
            }

            return {
                ...basic,
                minCurvature: minK,
                maxCurvature: maxK,
                positiveVertices: positiveCount,
                negativeVertices: negativeCount,
                flatVertices: zeroCount,
                topology: {
                    vertices: mesh.nVertices,
                    edges: mesh.nEdges,
                    faces: mesh.nFaces,
                    euler: basic.eulerCharacteristic,
                    genus: (2 - basic.eulerCharacteristic) / 2  // For closed surfaces
                }
            };
        }
    }

    // ============================================================================
    // EXPORTS
    // ============================================================================

    const RiemannianDiscrete = {
        // Bivector
        Bivector3D,

        // Core classes
        MeshConnectionBivector,
        MeshCurvature2Form,
        MeshParallelTransport,
        BianchiIdentityVerifier,

        // Utilities
        EPSILON,
        dot,
        cross,
        norm,
        normalize,
        vecAdd,
        vecSub,
        vecScale
    };

    // Export for different module systems
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = RiemannianDiscrete;
    }
    if (typeof global !== 'undefined') {
        global.RiemannianDiscrete = RiemannianDiscrete;
    }
    if (typeof window !== 'undefined') {
        window.RiemannianDiscrete = RiemannianDiscrete;
    }

})(typeof globalThis !== 'undefined' ? globalThis : this);
