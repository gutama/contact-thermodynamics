/**
 * Triangle Mesh Data Structure for FTGC (Geometric Calculus on Meshes)
 * 
 * Provides a typed-array-based mesh representation with precomputed topology
 * and geometric quantities needed for discrete geometric calculus.
 * 
 * Key features:
 * - Typed arrays (Float64Array, Uint32Array) for performance
 * - Auto-detected boundary vertices
 * - Mixed Voronoi dual areas (robust for obtuse triangles)
 * - Cotan weights for metric encoding
 * 
 * @module mesh
 * @license MIT
 */

(function (global) {
    'use strict';

    // ============================================================================
    // IMPORT SHARED UTILITIES
    // ============================================================================

    let Utils;
    if (typeof require !== 'undefined') {
        try {
            Utils = require('./utils.js');
        } catch (e) {
            Utils = global.ContactThermoUtils || {};
        }
    } else {
        Utils = global.ContactThermoUtils || {};
    }

    // Use shared EPSILON for consistency, fallback to local if utils not loaded
    const EPSILON = Utils.EPSILON || 1e-10;
    const { abs, sqrt, PI, acos, min, max } = Math;

    // ============================================================================
    // VECTOR UTILITIES (3D)
    // ============================================================================

    function vec3Sub(a, aOffset, b, bOffset, out, outOffset) {
        out[outOffset] = a[aOffset] - b[bOffset];
        out[outOffset + 1] = a[aOffset + 1] - b[bOffset + 1];
        out[outOffset + 2] = a[aOffset + 2] - b[bOffset + 2];
    }

    function vec3Dot(a, aOffset, b, bOffset) {
        return a[aOffset] * b[bOffset] +
            a[aOffset + 1] * b[bOffset + 1] +
            a[aOffset + 2] * b[bOffset + 2];
    }

    function vec3Cross(a, aOffset, b, bOffset, out, outOffset) {
        const ax = a[aOffset], ay = a[aOffset + 1], az = a[aOffset + 2];
        const bx = b[bOffset], by = b[bOffset + 1], bz = b[bOffset + 2];
        out[outOffset] = ay * bz - az * by;
        out[outOffset + 1] = az * bx - ax * bz;
        out[outOffset + 2] = ax * by - ay * bx;
    }

    function vec3Norm(a, offset) {
        const x = a[offset], y = a[offset + 1], z = a[offset + 2];
        return sqrt(x * x + y * y + z * z);
    }

    function vec3NormSq(a, offset) {
        const x = a[offset], y = a[offset + 1], z = a[offset + 2];
        return x * x + y * y + z * z;
    }

    // ============================================================================
    // TRIANGLE MESH
    // ============================================================================

    /**
     * Triangle mesh with typed array storage and precomputed topology.
     * 
     * Vertices are stored as Float64Array with stride 3 (x, y, z per vertex).
     * Faces are stored as Uint32Array with stride 3 (i, j, k per triangle).
     */
    class TriangleMesh {
        /**
         * @param {Float64Array} vertices - Vertex positions (n_vertices * 3)
         * @param {Uint32Array} faces - Triangle indices (n_faces * 3)
         * @param {Object} opts - Options
         * @param {boolean} opts.buildTopology - Auto-build topology (default: true)
         */
        constructor(vertices, faces, opts = {}) {
            const { buildTopology = true } = opts;

            // Store vertex and face data
            this.vertices = vertices;
            this.faces = faces;
            this.nVertices = vertices.length / 3 | 0;
            this.nFaces = faces.length / 3 | 0;

            // Topology (built lazily or on request)
            this.edges = null;           // Uint32Array (n_edges * 2)
            this.nEdges = 0;
            this.edgeToIdx = null;       // Map<string, number>
            this.vertexTriangles = null; // Array<Uint32Array>
            this.vertexEdges = null;     // Array<Uint32Array>
            this.edgeTriangles = null;   // Uint32Array (n_edges * 2), -1 if boundary

            // Boundary
            this.boundaryVertices = null;  // Uint8Array mask (1 = boundary)
            this.boundaryEdges = null;     // Uint8Array mask (1 = boundary)

            // Geometry cache
            this._faceAreas = null;
            this._vertexDualAreas = null;
            this._cotanWeights = null;
            this._edgeLengths = null;
            this._faceNormals = null;
            this._edgeVectors = null;
            this._dihedralAngles = null;

            if (buildTopology) {
                this.buildTopology();
            }
        }

        // ========================================================================
        // TOPOLOGY CONSTRUCTION
        // ========================================================================

        /**
         * Build edge list and adjacency structures.
         */
        buildTopology() {
            const nV = this.nVertices;
            const nF = this.nFaces;
            const faces = this.faces;

            // Temporary structures
            const edgeMap = new Map();
            const vertexTriLists = [];
            for (let i = 0; i < nV; i++) vertexTriLists.push([]);

            // Collect edges and vertex-triangle adjacency
            for (let t = 0; t < nF; t++) {
                const i = faces[t * 3];
                const j = faces[t * 3 + 1];
                const k = faces[t * 3 + 2];

                vertexTriLists[i].push(t);
                vertexTriLists[j].push(t);
                vertexTriLists[k].push(t);

                const triEdges = [[i, j], [j, k], [k, i]];
                for (const [a, b] of triEdges) {
                    const key = a < b ? `${a},${b}` : `${b},${a}`;
                    if (!edgeMap.has(key)) {
                        edgeMap.set(key, { idx: edgeMap.size, v0: min(a, b), v1: max(a, b), tris: [] });
                    }
                    edgeMap.get(key).tris.push(t);
                }
            }

            // Build edge array
            const nE = edgeMap.size;
            this.nEdges = nE;
            this.edges = new Uint32Array(nE * 2);
            this.edgeToIdx = new Map();
            this.edgeTriangles = new Int32Array(nE * 2).fill(-1);

            for (const [key, info] of edgeMap) {
                const eIdx = info.idx;
                this.edges[eIdx * 2] = info.v0;
                this.edges[eIdx * 2 + 1] = info.v1;
                this.edgeToIdx.set(key, eIdx);
                for (let i = 0; i < info.tris.length && i < 2; i++) {
                    this.edgeTriangles[eIdx * 2 + i] = info.tris[i];
                }
            }

            // Convert vertex-triangle lists to typed arrays
            this.vertexTriangles = vertexTriLists.map(arr => new Uint32Array(arr));

            // Build vertex-edge adjacency
            const vertexEdgeLists = [];
            for (let i = 0; i < nV; i++) vertexEdgeLists.push([]);
            for (let e = 0; e < nE; e++) {
                const v0 = this.edges[e * 2];
                const v1 = this.edges[e * 2 + 1];
                vertexEdgeLists[v0].push(e);
                vertexEdgeLists[v1].push(e);
            }
            this.vertexEdges = vertexEdgeLists.map(arr => new Uint32Array(arr));

            // Detect boundary edges and vertices
            this.boundaryEdges = new Uint8Array(nE);
            this.boundaryVertices = new Uint8Array(nV);

            for (let e = 0; e < nE; e++) {
                if (this.edgeTriangles[e * 2 + 1] < 0) {
                    this.boundaryEdges[e] = 1;
                    this.boundaryVertices[this.edges[e * 2]] = 1;
                    this.boundaryVertices[this.edges[e * 2 + 1]] = 1;
                }
            }
        }

        /**
         * Get edge index for vertex pair (order-independent).
         * @param {number} v0 
         * @param {number} v1 
         * @returns {number} Edge index, or -1 if not found
         */
        getEdgeIndex(v0, v1) {
            const key = v0 < v1 ? `${v0},${v1}` : `${v1},${v0}`;
            const idx = this.edgeToIdx.get(key);
            return idx !== undefined ? idx : -1;
        }

        // ========================================================================
        // GEOMETRY: AREAS, LENGTHS, ANGLES
        // ========================================================================

        /**
         * Get vertex position as [x, y, z].
         */
        getVertex(idx) {
            const offset = idx * 3;
            return [this.vertices[offset], this.vertices[offset + 1], this.vertices[offset + 2]];
        }

        /**
         * Compute triangle areas (cached).
         * @returns {Float64Array} Area per face
         */
        faceAreas() {
            if (this._faceAreas) return this._faceAreas;

            const nF = this.nFaces;
            const areas = new Float64Array(nF);
            const v = this.vertices;
            const f = this.faces;
            const tmp1 = new Float64Array(3);
            const tmp2 = new Float64Array(3);
            const cross = new Float64Array(3);

            for (let t = 0; t < nF; t++) {
                const i = f[t * 3], j = f[t * 3 + 1], k = f[t * 3 + 2];
                vec3Sub(v, j * 3, v, i * 3, tmp1, 0);
                vec3Sub(v, k * 3, v, i * 3, tmp2, 0);
                vec3Cross(tmp1, 0, tmp2, 0, cross, 0);
                areas[t] = 0.5 * vec3Norm(cross, 0);
            }

            this._faceAreas = areas;
            return areas;
        }

        /**
         * Compute edge lengths (cached).
         * @returns {Float64Array} Length per edge
         */
        edgeLengths() {
            if (this._edgeLengths) return this._edgeLengths;

            const nE = this.nEdges;
            const lengths = new Float64Array(nE);
            const v = this.vertices;
            const e = this.edges;
            const tmp = new Float64Array(3);

            for (let i = 0; i < nE; i++) {
                const v0 = e[i * 2], v1 = e[i * 2 + 1];
                vec3Sub(v, v1 * 3, v, v0 * 3, tmp, 0);
                lengths[i] = vec3Norm(tmp, 0);
            }

            this._edgeLengths = lengths;
            return lengths;
        }

        /**
         * Compute cotangent of angle at vertex v0 in triangle (v0, v1, v2).
         * cot(θ) = (a·b) / |a×b|
         */
        _computeCotan(v0, v1, v2) {
            const verts = this.vertices;
            const tmp1 = new Float64Array(3);
            const tmp2 = new Float64Array(3);
            const cross = new Float64Array(3);

            vec3Sub(verts, v1 * 3, verts, v0 * 3, tmp1, 0);
            vec3Sub(verts, v2 * 3, verts, v0 * 3, tmp2, 0);
            vec3Cross(tmp1, 0, tmp2, 0, cross, 0);

            const crossNorm = vec3Norm(cross, 0);
            if (crossNorm < EPSILON) return 0;

            return vec3Dot(tmp1, 0, tmp2, 0) / crossNorm;
        }

        /**
         * Compute cotan weights per edge (cached).
         * Each edge gets sum of cotangents from its two adjacent triangles.
         * @returns {Float64Array} Weight per edge
         */
        cotanWeights() {
            if (this._cotanWeights) return this._cotanWeights;

            const nE = this.nEdges;
            const weights = new Float64Array(nE);
            const f = this.faces;

            for (let t = 0; t < this.nFaces; t++) {
                const i = f[t * 3], j = f[t * 3 + 1], k = f[t * 3 + 2];

                // Cotan at i → opposite edge is (j, k)
                const cotI = this._computeCotan(i, j, k);
                // Cotan at j → opposite edge is (k, i)
                const cotJ = this._computeCotan(j, k, i);
                // Cotan at k → opposite edge is (i, j)
                const cotK = this._computeCotan(k, i, j);

                // Add to corresponding edges
                const eJK = this.getEdgeIndex(j, k);
                const eKI = this.getEdgeIndex(k, i);
                const eIJ = this.getEdgeIndex(i, j);

                if (eJK >= 0) weights[eJK] += cotI * 0.5;
                if (eKI >= 0) weights[eKI] += cotJ * 0.5;
                if (eIJ >= 0) weights[eIJ] += cotK * 0.5;
            }

            this._cotanWeights = weights;
            return weights;
        }

        /**
         * Compute mixed Voronoi dual areas per vertex (cached).
         * Uses circumcentric for non-obtuse triangles, barycentric fallback for obtuse.
         * @returns {Float64Array} Dual area per vertex
         */
        vertexDualAreas() {
            if (this._vertexDualAreas) return this._vertexDualAreas;

            const nV = this.nVertices;
            const nF = this.nFaces;
            const areas = new Float64Array(nV);
            const verts = this.vertices;
            const f = this.faces;
            const faceAreas = this.faceAreas();

            const tmp1 = new Float64Array(3);
            const tmp2 = new Float64Array(3);

            for (let t = 0; t < nF; t++) {
                const vi = f[t * 3], vj = f[t * 3 + 1], vk = f[t * 3 + 2];
                const triArea = faceAreas[t];

                // Compute squared edge lengths (opposite to each vertex)
                // Edge (vj, vk) is opposite to vi
                vec3Sub(verts, vk * 3, verts, vj * 3, tmp1, 0);
                const len_jk_sq = vec3NormSq(tmp1, 0);
                // Edge (vk, vi) is opposite to vj
                vec3Sub(verts, vi * 3, verts, vk * 3, tmp1, 0);
                const len_ki_sq = vec3NormSq(tmp1, 0);
                // Edge (vi, vj) is opposite to vk
                vec3Sub(verts, vj * 3, verts, vi * 3, tmp1, 0);
                const len_ij_sq = vec3NormSq(tmp1, 0);

                // Compute angles via dot products
                vec3Sub(verts, vj * 3, verts, vi * 3, tmp1, 0);
                vec3Sub(verts, vk * 3, verts, vi * 3, tmp2, 0);
                const dotI = vec3Dot(tmp1, 0, tmp2, 0);

                vec3Sub(verts, vi * 3, verts, vj * 3, tmp1, 0);
                vec3Sub(verts, vk * 3, verts, vj * 3, tmp2, 0);
                const dotJ = vec3Dot(tmp1, 0, tmp2, 0);

                vec3Sub(verts, vi * 3, verts, vk * 3, tmp1, 0);
                vec3Sub(verts, vj * 3, verts, vk * 3, tmp2, 0);
                const dotK = vec3Dot(tmp1, 0, tmp2, 0);

                // Check for obtuse angles (dot < 0 means angle > 90°)
                const obtuseI = dotI < 0;
                const obtuseJ = dotJ < 0;
                const obtuseK = dotK < 0;
                const isObtuse = obtuseI || obtuseJ || obtuseK;

                if (isObtuse) {
                    // Use barycentric fallback for obtuse triangles
                    if (obtuseI) {
                        // Obtuse at i: i gets half, j and k get quarter each
                        areas[vi] += triArea * 0.5;
                        areas[vj] += triArea * 0.25;
                        areas[vk] += triArea * 0.25;
                    } else if (obtuseJ) {
                        areas[vi] += triArea * 0.25;
                        areas[vj] += triArea * 0.5;
                        areas[vk] += triArea * 0.25;
                    } else {
                        areas[vi] += triArea * 0.25;
                        areas[vj] += triArea * 0.25;
                        areas[vk] += triArea * 0.5;
                    }
                } else {
                    // Non-obtuse: use Voronoi (cotan formula)
                    // A_v = (1/8) * Σ (cot α) * |e|²
                    // where α is the angle opposite to edge e
                    const cotI = this._computeCotan(vi, vj, vk);
                    const cotJ = this._computeCotan(vj, vk, vi);
                    const cotK = this._computeCotan(vk, vi, vj);

                    // Contribution to vi: edges (vi,vj) and (vi,vk)
                    // Edge (vi,vj) has opposite angle at k (cotK), length len_ij_sq
                    // Edge (vi,vk) has opposite angle at j (cotJ), length len_ki_sq
                    const areaI = (cotK * len_ij_sq + cotJ * len_ki_sq) / 8;

                    // Contribution to vj: edges (vj,vi) and (vj,vk)
                    // Edge (vj,vi) has opposite angle at k (cotK), length len_ij_sq
                    // Edge (vj,vk) has opposite angle at i (cotI), length len_jk_sq
                    const areaJ = (cotK * len_ij_sq + cotI * len_jk_sq) / 8;

                    // Contribution to vk: edges (vk,vi) and (vk,vj)
                    // Edge (vk,vi) has opposite angle at j (cotJ), length len_ki_sq
                    // Edge (vk,vj) has opposite angle at i (cotI), length len_jk_sq
                    const areaK = (cotJ * len_ki_sq + cotI * len_jk_sq) / 8;

                    // Clamp to positive (numerical safety)
                    areas[vi] += max(0, areaI);
                    areas[vj] += max(0, areaJ);
                    areas[vk] += max(0, areaK);
                }
            }

            // Ensure minimum area for numerical stability
            for (let i = 0; i < nV; i++) {
                if (areas[i] < EPSILON) areas[i] = EPSILON;
            }

            this._vertexDualAreas = areas;
            return areas;
        }

        /**
         * Compute unit face normals (cached).
         * Normal direction follows right-hand rule from triangle winding.
         * @returns {Float64Array} Normal vectors (nFaces * 3)
         */
        faceNormals() {
            if (this._faceNormals) return this._faceNormals;

            const nF = this.nFaces;
            const normals = new Float64Array(nF * 3);
            const v = this.vertices;
            const f = this.faces;
            const tmp1 = new Float64Array(3);
            const tmp2 = new Float64Array(3);
            const cross = new Float64Array(3);

            for (let t = 0; t < nF; t++) {
                const i = f[t * 3], j = f[t * 3 + 1], k = f[t * 3 + 2];
                vec3Sub(v, j * 3, v, i * 3, tmp1, 0);
                vec3Sub(v, k * 3, v, i * 3, tmp2, 0);
                vec3Cross(tmp1, 0, tmp2, 0, cross, 0);
                const len = vec3Norm(cross, 0);
                if (len > EPSILON) {
                    normals[t * 3] = cross[0] / len;
                    normals[t * 3 + 1] = cross[1] / len;
                    normals[t * 3 + 2] = cross[2] / len;
                } else {
                    // Degenerate triangle, use z-up
                    normals[t * 3 + 2] = 1;
                }
            }

            this._faceNormals = normals;
            return normals;
        }

        /**
         * Get face normal as [nx, ny, nz].
         * @param {number} faceIdx
         * @returns {number[]}
         */
        getFaceNormal(faceIdx) {
            const normals = this.faceNormals();
            const offset = faceIdx * 3;
            return [normals[offset], normals[offset + 1], normals[offset + 2]];
        }

        /**
         * Compute unit edge vectors (cached).
         * Edge direction goes from lower-index vertex to higher-index vertex.
         * @returns {Float64Array} Edge vectors (nEdges * 3)
         */
        edgeVectors() {
            if (this._edgeVectors) return this._edgeVectors;

            const nE = this.nEdges;
            const edgeVecs = new Float64Array(nE * 3);
            const v = this.vertices;
            const e = this.edges;
            const tmp = new Float64Array(3);

            for (let i = 0; i < nE; i++) {
                const v0 = e[i * 2], v1 = e[i * 2 + 1];
                vec3Sub(v, v1 * 3, v, v0 * 3, tmp, 0);
                const len = vec3Norm(tmp, 0);
                if (len > EPSILON) {
                    edgeVecs[i * 3] = tmp[0] / len;
                    edgeVecs[i * 3 + 1] = tmp[1] / len;
                    edgeVecs[i * 3 + 2] = tmp[2] / len;
                }
            }

            this._edgeVectors = edgeVecs;
            return edgeVecs;
        }

        /**
         * Get edge vector as [ex, ey, ez].
         * @param {number} edgeIdx
         * @returns {number[]}
         */
        getEdgeVector(edgeIdx) {
            const edgeVecs = this.edgeVectors();
            const offset = edgeIdx * 3;
            return [edgeVecs[offset], edgeVecs[offset + 1], edgeVecs[offset + 2]];
        }

        /**
         * Compute dihedral angles at each edge (cached).
         * Dihedral angle = angle between face normals at shared edge.
         * Convention: angle in [0, 2π), measured as rotation from n1 to n2 around edge.
         * Boundary edges have angle = π (flat).
         * @returns {Float64Array} Dihedral angle per edge (in radians)
         */
        dihedralAngles() {
            if (this._dihedralAngles) return this._dihedralAngles;

            const nE = this.nEdges;
            const angles = new Float64Array(nE);
            const normals = this.faceNormals();
            const edgeVecs = this.edgeVectors();

            for (let e = 0; e < nE; e++) {
                const t0 = this.edgeTriangles[e * 2];
                const t1 = this.edgeTriangles[e * 2 + 1];

                if (t1 < 0) {
                    // Boundary edge: treat as flat (π)
                    angles[e] = PI;
                    continue;
                }

                // Get normals of adjacent faces
                const n0 = [normals[t0 * 3], normals[t0 * 3 + 1], normals[t0 * 3 + 2]];
                const n1 = [normals[t1 * 3], normals[t1 * 3 + 1], normals[t1 * 3 + 2]];

                // Edge vector (for signed angle)
                const ev = [edgeVecs[e * 3], edgeVecs[e * 3 + 1], edgeVecs[e * 3 + 2]];

                // Compute angle between normals
                // cos(θ) = n0 · n1
                const cosTheta = n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2];

                // sin(θ) = (n0 × n1) · edge / |edge|
                // Cross product n0 × n1
                const cross = [
                    n0[1] * n1[2] - n0[2] * n1[1],
                    n0[2] * n1[0] - n0[0] * n1[2],
                    n0[0] * n1[1] - n0[1] * n1[0]
                ];
                const sinTheta = cross[0] * ev[0] + cross[1] * ev[1] + cross[2] * ev[2];

                // Dihedral angle using atan2 for full range
                // Convention: 0 = coplanar same direction, π = coplanar opposite
                angles[e] = PI - Math.atan2(sinTheta, cosTheta);
            }

            this._dihedralAngles = angles;
            return angles;
        }

        /**
         * Invalidate cached geometry (call after modifying vertices).
         */
        invalidateGeometryCache() {
            this._faceAreas = null;
            this._vertexDualAreas = null;
            this._cotanWeights = null;
            this._edgeLengths = null;
            this._faceNormals = null;
            this._edgeVectors = null;
            this._dihedralAngles = null;
        }

        // ========================================================================
        // MESH GENERATORS (Static)
        // ========================================================================

        /**
         * Create a regular grid mesh (nx × ny vertices, triangulated).
         * @param {number} nx - Vertices in x direction
         * @param {number} ny - Vertices in y direction
         * @param {number} Lx - Domain size in x
         * @param {number} Ly - Domain size in y
         * @returns {TriangleMesh}
         */
        static createGrid(nx, ny, Lx = 1, Ly = 1) {
            const nV = nx * ny;
            const nF = 2 * (nx - 1) * (ny - 1);

            const vertices = new Float64Array(nV * 3);
            const faces = new Uint32Array(nF * 3);

            const dx = Lx / (nx - 1);
            const dy = Ly / (ny - 1);

            // Create vertices
            for (let j = 0; j < ny; j++) {
                for (let i = 0; i < nx; i++) {
                    const idx = j * nx + i;
                    vertices[idx * 3] = i * dx;
                    vertices[idx * 3 + 1] = j * dy;
                    vertices[idx * 3 + 2] = 0;
                }
            }

            // Create triangles
            let fIdx = 0;
            for (let j = 0; j < ny - 1; j++) {
                for (let i = 0; i < nx - 1; i++) {
                    const v00 = j * nx + i;
                    const v10 = j * nx + (i + 1);
                    const v01 = (j + 1) * nx + i;
                    const v11 = (j + 1) * nx + (i + 1);

                    // Lower triangle
                    faces[fIdx * 3] = v00;
                    faces[fIdx * 3 + 1] = v10;
                    faces[fIdx * 3 + 2] = v11;
                    fIdx++;

                    // Upper triangle
                    faces[fIdx * 3] = v00;
                    faces[fIdx * 3 + 1] = v11;
                    faces[fIdx * 3 + 2] = v01;
                    fIdx++;
                }
            }

            return new TriangleMesh(vertices, faces);
        }

        /**
         * Create a triangulated disk mesh.
         * @param {number} nRadial - Number of radial divisions
         * @param {number} nAngular - Number of angular divisions
         * @param {number} radius - Disk radius
         * @returns {TriangleMesh}
         */
        static createDisk(nRadial, nAngular, radius = 1) {
            const nV = 1 + nRadial * nAngular; // center + rings
            const nF = nAngular + 2 * nAngular * (nRadial - 1);

            const vertices = new Float64Array(nV * 3);
            const faces = new Uint32Array(nF * 3);

            // Center vertex
            vertices[0] = 0;
            vertices[1] = 0;
            vertices[2] = 0;

            // Ring vertices
            for (let r = 1; r <= nRadial; r++) {
                const rFrac = r / nRadial;
                for (let a = 0; a < nAngular; a++) {
                    const theta = 2 * PI * a / nAngular;
                    const idx = 1 + (r - 1) * nAngular + a;
                    vertices[idx * 3] = radius * rFrac * Math.cos(theta);
                    vertices[idx * 3 + 1] = radius * rFrac * Math.sin(theta);
                    vertices[idx * 3 + 2] = 0;
                }
            }

            // Center triangles (fan)
            let fIdx = 0;
            for (let a = 0; a < nAngular; a++) {
                const a1 = 1 + a;
                const a2 = 1 + ((a + 1) % nAngular);
                faces[fIdx * 3] = 0;
                faces[fIdx * 3 + 1] = a1;
                faces[fIdx * 3 + 2] = a2;
                fIdx++;
            }

            // Ring triangles
            for (let r = 1; r < nRadial; r++) {
                for (let a = 0; a < nAngular; a++) {
                    const a1 = (a + 1) % nAngular;
                    const v00 = 1 + (r - 1) * nAngular + a;
                    const v01 = 1 + (r - 1) * nAngular + a1;
                    const v10 = 1 + r * nAngular + a;
                    const v11 = 1 + r * nAngular + a1;

                    faces[fIdx * 3] = v00;
                    faces[fIdx * 3 + 1] = v10;
                    faces[fIdx * 3 + 2] = v11;
                    fIdx++;

                    faces[fIdx * 3] = v00;
                    faces[fIdx * 3 + 1] = v11;
                    faces[fIdx * 3 + 2] = v01;
                    fIdx++;
                }
            }

            return new TriangleMesh(vertices, faces);
        }

        /**
         * Create an icosahedron (closed mesh, 12 vertices, 20 faces).
         * @param {number} radius - Radius of circumscribed sphere
         * @returns {TriangleMesh}
         */
        static createIcosahedron(radius = 1) {
            const phi = (1 + sqrt(5)) / 2; // Golden ratio
            const a = radius / sqrt(1 + phi * phi);
            const b = a * phi;

            const vertices = new Float64Array([
                -a, b, 0, a, b, 0, -a, -b, 0, a, -b, 0,
                0, -a, b, 0, a, b, 0, -a, -b, 0, a, -b,
                b, 0, -a, b, 0, a, -b, 0, -a, -b, 0, a
            ]);

            const faces = new Uint32Array([
                0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11,
                1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
                3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9,
                4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1
            ]);

            return new TriangleMesh(vertices, faces);
        }
    }

    // ============================================================================
    // EXPORTS
    // ============================================================================

    const MeshModule = {
        TriangleMesh
    };

    // CommonJS
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = MeshModule;
    }
    // AMD
    else if (typeof define === 'function' && define.amd) {
        define(function () { return MeshModule; });
    }
    // Browser global
    else {
        global.ContactThermo = global.ContactThermo || {};
        Object.assign(global.ContactThermo, MeshModule);
    }

})(typeof globalThis !== 'undefined' ? globalThis : typeof window !== 'undefined' ? window : this);
