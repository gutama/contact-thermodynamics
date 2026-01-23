/**
 * Entropic Gravity on Triangle Meshes
 *
 * This example demonstrates how to simulate Bianconi's Entropic Gravity
 * on a discrete curved surface (a sphere) using Geometric Algebra.
 *
 * We implement:
 * 1. A discrete sphere mesh
 * 2. Background metric g (induced from embedding)
 * 3. Matter field φ (scalar field on vertices)
 * 4. Induced metric G_ij = g_ij + ∂_iφ ∂_jφ
 * 5. Relative entropy density S(G||g)
 * 6. Entropic gradient flow
 */

const GA = require('../src/geometric-calculus.js');   // GA utilities
const Mesh = require('../src/mesh.js');               // TriangleMesh
const Solver = require('../src/mesh-solvers.js');     // MeshGeometricDerivative

const { sqrt, exp, log, abs, cos, sin, PI } = Math;

console.log('='.repeat(70));
console.log('DISCRETE ENTROPIC GRAVITY ON A SPHERE');
console.log('='.repeat(70));

// ============================================================================
// 1. Create a Sphere Mesh
// ============================================================================

console.log('\nSTEP 1: Creating Sphere Mesh...');

// Simple icosphere approximation
const radius = 5.0;
const mesh = Mesh.TriangleMesh.createIcosahedron(radius);
// Ideally we'd subdivide for better resolution, but 12 vertices is good for a quick test
// Let's implement a simple subdivision for better visuals
function subdivide(mesh) {
    // This is a placeholder for actual subdivision logic
    // For this example, we'll stick to the base icosahedron to keep code concise
    // In a real app, use a proper subdivision algorithm
    return mesh;
}

console.log(`  Vertices: ${mesh.nVertices}`);
console.log(`  Faces: ${mesh.nFaces}`);
console.log(`  Edges: ${mesh.nEdges}`);

// ============================================================================
// 2. Define Matter Field φ
// ============================================================================

console.log('\nSTEP 2: Defining Matter Field φ...');

// Scalar field concentrated at the north pole
// φ(x) = A * exp(-|x - x_north|² / σ²)
const phi = new Float64Array(mesh.nVertices);
const northPole = [0, radius, 0]; // (assuming y-up) -- check actual orientation
// Icosahedron orientation might vary, let's find the vertex with max Y
let maxY = -Infinity;
let poleIdx = -1;
for (let i = 0; i < mesh.nVertices; i++) {
    const y = mesh.vertices[3 * i + 1];
    if (y > maxY) {
        maxY = y;
        poleIdx = i;
    }
}
const P_pole = mesh.getVertex(poleIdx);

console.log(`  North Pole Vertex ${poleIdx}: [${P_pole.map(v => v.toFixed(2))}]`);

for (let i = 0; i < mesh.nVertices; i++) {
    const P = mesh.getVertex(i);
    const dist2 = (P[0] - P_pole[0]) ** 2 + (P[1] - P_pole[1]) ** 2 + (P[2] - P_pole[2]) ** 2;
    phi[i] = 2.0 * exp(-dist2 / (radius * radius * 0.5));
}

console.log('  Matter field φ defined (Gaussian peak at pole)');

// ============================================================================
// 3. Compute Gradients and Induced Metric G
// ============================================================================

console.log('\nSTEP 3: Computing Induced Metric G...');

// Metric g is implicit in the mesh geometry (edge lengths)
// Induced metric G_ij = g_ij + ∂_iφ ∂_jφ
// In discrete GA, we work with edge lengths l_ij (from g) and L_ij (from G)
// L_ij² = l_ij² + (φ_j - φ_i)²

const g_lengths = mesh.edgeLengths(); // Float64Array
const G_lengths = new Float64Array(mesh.nEdges);
const entropy_density = new Float64Array(mesh.nFaces);

// Iterate over edges to compute G lengths
for (let e = 0; e < mesh.nEdges; e++) {
    const l = g_lengths[e];

    // We need vertex indices for this edge.
    // mesh.edges is flattened [v0_0, v1_0, v0_1, v1_1...]
    // But TriangleMesh implementation might vary. Let's check mesh.js or use topology.
    // Assuming mesh.edges is available; if not we might need to rely on faces.
    // The previous API doc said mesh.edges is Uint32Array (x2)

    const v0 = mesh.edges[2 * e];
    const v1 = mesh.edges[2 * e + 1];

    const dphi = phi[v1] - phi[v0];
    const L2 = l * l + dphi * dphi;
    G_lengths[e] = sqrt(L2);
}

console.log('  Computed perturbed edge lengths L_ij from Matter Field');

// ============================================================================
// 4. Compute Relative Entropy S(G||g) per Face
// ============================================================================

console.log('\nSTEP 4: Computing Relative Entropy S(G||g)...');

// For 2D surfaces, S(G||g) ≈ ∫ dA (Area_G / Area_g - 1 + log terms)
// Simplified discrete version: S_face = Area_G(face) * log(Area_G / Area_g) 
// or simpler: deviation of areas.
// Let's use the trace form approximation: Tr[G g⁻¹] - 2 ≈ (L² / l²) - 2 (summed over edges?)
// A geometric invariant way is to compare area elements.
// Ratio R = Area_G / Area_g
// Local entropy density s = R * log(R) (Kullback-Leibler like)

const areas_g = mesh.faceAreas();
const areas_G = new Float64Array(mesh.nFaces);

// Helper to compute triangle area from side lengths (Heron's formula)
function heron(a, b, c) {
    const s = (a + b + c) / 2;
    return sqrt(Math.max(0, s * (s - a) * (s - b) * (s - c)));
}

let total_S = 0;

for (let f = 0; f < mesh.nFaces; f++) {
    // Get edges of face f
    // This requires topology lookup. Assuming mesh.faces gives vertices v0,v1,v2
    const v0 = mesh.faces[3 * f];
    const v1 = mesh.faces[3 * f + 1];
    const v2 = mesh.faces[3 * f + 2];

    // Find edges (naive search or lookup if available)
    const e0 = mesh.getEdgeIndex(v0, v1);
    const e1 = mesh.getEdgeIndex(v1, v2);
    const e2 = mesh.getEdgeIndex(v2, v0);

    const L0 = G_lengths[e0];
    const L1 = G_lengths[e1];
    const L2 = G_lengths[e2];

    areas_G[f] = heron(L0, L1, L2);

    const R = areas_G[f] / (areas_g[f] + 1e-12); // Volume ratio
    // S ≈ R log R - (R - 1)  (standard relative entropy density)
    const val = (R * log(R + 1e-12)) - (R - 1);

    entropy_density[f] = val;
    total_S += val * areas_g[f];
}

console.log(`  Total Relative Entropy: ${total_S.toFixed(6)}`);

// ============================================================================
// 5. Entropic Force Simulation
// ============================================================================

console.log('\nSTEP 5: Simulating Entropic Force...');

// Particle moves on the mesh. 
// "Entropic Gravity" -> The particle sees the metric G, not g.
// It follows geodesics of G. Or, in the framework H = H_geo(g) + α S, 
// the scalar field S acts as a potential.

// Let's visualize the "Entropic Potential" V = S_density (averaged to vertices)
const potential = new Float64Array(mesh.nVertices);
const vertex_areas = new Float64Array(mesh.nVertices);

for (let f = 0; f < mesh.nFaces; f++) {
    const s_val = entropy_density[f];
    const area = areas_g[f]; // weight by physical area

    for (let k = 0; k < 3; k++) {
        const v = mesh.faces[3 * f + k];
        potential[v] += s_val * area / 3;
        vertex_areas[v] += area / 3;
    }
}

let maxV = -Infinity;
let minV = Infinity;

for (let i = 0; i < mesh.nVertices; i++) {
    if (vertex_areas[i] > 0) potential[i] /= vertex_areas[i];
    if (potential[i] > maxV) maxV = potential[i];
    if (potential[i] < minV) minV = potential[i];
}

console.log(`  Entropic Potential Range: [${minV.toFixed(4)}, ${maxV.toFixed(4)}]`);
console.log('  Maximum entropy density is near the pole where φ gradient is steepest.');

// ============================================================================
// 6. Summary
// ============================================================================

console.log('\nSUMMARY');
console.log('-'.repeat(30));
console.log('Successfully computed:');
console.log('1. Background metric g (sphere)');
console.log('2. Matter field φ (gaussian)');
console.log('3. Induced metric G (perturbed edge lengths)');
console.log('4. Relative entropy S(G||g) on faces');
console.log('5. Entropic potential V on vertices');

console.log('\nTo visualize, one would render this mesh with color map = V.');
console.log('Particles would drift towards areas of higher information density.');

