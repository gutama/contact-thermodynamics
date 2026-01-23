/**
 * Discrete Entropic Flow on Triangle Meshes
 * 
 * Implements the dynamical evolution of particles driven by the relative entropy
 * gradient on a discrete manifold (triangle mesh).
 */

const { TriangleMesh } = require('./mesh.js');
const { MeshGeometricDerivative } = require('./mesh-ftgc.js');

const { sqrt, exp, log, max, abs } = Math;

class EntropicMeshflow {
    /**
     * @param {TriangleMesh} mesh - The geometry
     * @param {Object} options
     * @param {number} options.mass - Particle mass (inertia)
     * @param {number} options.coupling - Entropic coupling constant (alpha)
     * @param {number} options.damping - Friction/damping (gamma)
     */
    constructor(mesh, options = {}) {
        this.mesh = mesh;
        this.mass = options.mass || 1.0;
        this.coupling = options.coupling || 0.1;
        this.damping = options.damping || 0.0; // Dissipative flow if > 0

        // Fields
        this.phi = new Float64Array(mesh.nVertices); // Scalar field
        this.entropyDensity = new Float64Array(mesh.nFaces); // S per face
        this.entropicPotential = new Float64Array(mesh.nVertices); // V per vertex

        // Cache for optimization
        this.g_lengths = mesh.edgeLengths();
        this.G_lengths = new Float64Array(mesh.nEdges);
        this.areas_g = mesh.faceAreas();
        this.vertexDualAreas = mesh.vertexDualAreas();

        // Derivative operator
        this.nabla = new MeshGeometricDerivative(mesh);
    }

    /**
     * Update the scalar field phi(x) and recompute induced metric G and Entropy S.
     * @param {Function|Float64Array} fieldSource 
     */
    updateField(fieldSource) {
        // 1. Update Phi
        if (fieldSource instanceof Float64Array) {
            this.phi.set(fieldSource);
        } else if (typeof fieldSource === 'function') {
            for (let i = 0; i < this.mesh.nVertices; i++) {
                const p = this.mesh.getVertex(i);
                this.phi[i] = fieldSource(p, i);
            }
        }

        // 2. Update G lengths: L² = l² + (dφ)²
        for (let e = 0; e < this.mesh.nEdges; e++) {
            const v0 = this.mesh.edges[2 * e];
            const v1 = this.mesh.edges[2 * e + 1];
            const dphi = this.phi[v1] - this.phi[v0];
            const l = this.g_lengths[e];
            this.G_lengths[e] = sqrt(l * l + dphi * dphi);
        }

        // 3. Update Entropy Density on Faces
        const areas_G = new Float64Array(this.mesh.nFaces);

        for (let f = 0; f < this.mesh.nFaces; f++) {
            const v = [
                this.mesh.faces[3 * f],
                this.mesh.faces[3 * f + 1],
                this.mesh.faces[3 * f + 2]
            ];

            // Get edge indices (naive search for now, optimization needed for huge meshes)
            const e0 = this.mesh.getEdgeIndex(v[0], v[1]);
            const e1 = this.mesh.getEdgeIndex(v[1], v[2]);
            const e2 = this.mesh.getEdgeIndex(v[2], v[0]);

            const L0 = this.G_lengths[e0];
            const L1 = this.G_lengths[e1];
            const L2 = this.G_lengths[e2];

            // Heron's formula for Area_G
            const s = (L0 + L1 + L2) / 2;
            areas_G[f] = sqrt(max(0, s * (s - L0) * (s - L1) * (s - L2)));

            // Relative Entropy Density s = R log R - (R - 1)
            const area_g = this.areas_g[f];
            const R = areas_G[f] / (area_g + 1e-12);

            this.entropyDensity[f] = (R * log(R + 1e-12)) - (R - 1);
        }

        // 4. Compute Entropic Potential V at vertices (Vertex averaging)
        // V_i = (1/A_i) * sum_{faces f adj i} (area_f/3) * s_f
        this.entropicPotential.fill(0);

        // Accumulate weighted entropy
        for (let f = 0; f < this.mesh.nFaces; f++) {
            const s_val = this.entropyDensity[f];
            const area = this.areas_g[f];
            const contribution = s_val * area / 3;

            for (let k = 0; k < 3; k++) {
                const v = this.mesh.faces[3 * f + k];
                this.entropicPotential[v] += contribution;
            }
        }

        // Normalize by vertex dual area
        for (let i = 0; i < this.mesh.nVertices; i++) {
            if (this.vertexDualAreas[i] > 1e-12) {
                this.entropicPotential[i] /= this.vertexDualAreas[i];
            }
        }
    }

    /**
     * Compute gradient of entropic potential at a barycentric point.
     * @param {Object} point { face, u, v, w } (barycentric)
     * @returns {Float64Array} Gradient vector (tangent 3D vector)
     */
    getForce(point) {
        // F = - α ∇V
        // For barycentric interpolation, gradient is constant per face.
        // We can interpolate vertex gradients or just use the face gradient.
        // Using face gradient of V (from vertex values) is standard FEM.

        // grad V on face f: sum_i V_i * grad(phi_i)
        // grad basis function phi_i is a vector in the triangle plane.

        const f = point.face;
        const v_indices = [
            this.mesh.faces[3 * f],
            this.mesh.faces[3 * f + 1],
            this.mesh.faces[3 * f + 2]
        ];

        const V = [
            this.entropicPotential[v_indices[0]],
            this.entropicPotential[v_indices[1]],
            this.entropicPotential[v_indices[2]]
        ];

        // Vertices positions
        const p0 = this.mesh.getVertex(v_indices[0]);
        const p1 = this.mesh.getVertex(v_indices[1]);
        const p2 = this.mesh.getVertex(v_indices[2]);

        // Edges
        const u = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        const v = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];

        // Normal N
        const N = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0]
        ];
        const doubleArea = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
        const n = [N[0] / doubleArea, N[1] / doubleArea, N[2] / doubleArea];

        // Basis gradients
        // grad(lambda_0) = (N x e23) / 2A ...

        // Standard formula: grad f = (1/2A) sum_i f_i (N x e_jk)

        const e12 = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]]; // opposing v0
        const e20 = [p0[0] - p2[0], p0[1] - p2[1], p0[2] - p2[2]]; // opposing v1
        const e01 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]]; // opposing v2

        const edges = [e12, e20, e01];

        const gradV = [0, 0, 0];

        for (let i = 0; i < 3; i++) {
            // N x e_opposing
            const cross = [
                n[1] * edges[i][2] - n[2] * edges[i][1],
                n[2] * edges[i][0] - n[0] * edges[i][2],
                n[0] * edges[i][1] - n[1] * edges[i][0]
            ];

            gradV[0] += V[i] * cross[0];
            gradV[1] += V[i] * cross[1];
            gradV[2] += V[i] * cross[2];
        }

        gradV[0] /= doubleArea;
        gradV[1] /= doubleArea;
        gradV[2] /= doubleArea;

        // Force is -alpha * gradV
        return [-this.coupling * gradV[0], -this.coupling * gradV[1], -this.coupling * gradV[2]];
    }
}

module.exports = { EntropicMeshflow };
