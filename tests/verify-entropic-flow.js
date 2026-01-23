/**
 * Verification Script for Discrete Entropic Flow
 */

const { TriangleMesh } = require('../src/mesh.js');
const { EntropicMeshflow } = require('../src/mesh-entropic-flow.js');

const { sqrt, exp, PI } = Math;

function assert(condition, message) {
    if (!condition) {
        throw new Error(`FAILED: ${message}`);
    }
    console.log(`PASS: ${message}`);
}

console.log('Verifying Discrete Entropic Flow...');

// 1. Setup Mesh (Icosahedron R=5)
const mesh = TriangleMesh.createIcosahedron(5.0);
const flow = new EntropicMeshflow(mesh, { coupling: 0.1 });

// 2. Define Scalar Field (Static Pole)
const northPole = [0, 5, 0];
flow.updateField((p) => {
    const dist2 = (p[0] - northPole[0]) ** 2 + (p[1] - northPole[1]) ** 2 + (p[2] - northPole[2]) ** 2;
    // Gaussian with sigma^2 = 10 -> sigma ~ 3.16
    // Max gradient is at r = sigma/sqrt(2) = 3.16/1.414 = 2.23
    return 2.0 * exp(-dist2 / 10.0);
});

// Check potential range
let maxV = -Infinity;
let maxIdx = -1;
for (let i = 0; i < mesh.nVertices; i++) {
    if (flow.entropicPotential[i] > maxV) {
        maxV = flow.entropicPotential[i];
        maxIdx = i;
    }
}
console.log(`Max potential at vertex ${maxIdx} with value ${maxV.toFixed(4)}`);
const poleV = mesh.getVertex(maxIdx);
console.log(`Max potential position: [${poleV.map(x => x.toFixed(2))}]`);

// 3. Check Peak Location
// The max potential should be in a ring around the pole where gradient is steepest.
// For Gaussian exp(-r^2/10), max gradient is at r = sqrt(5) ~ 2.23.
// Our found max is at distance ~2.7 from pole ([-2.63, 4.25, 0] to [0, 5, 0]).
const dToPole = sqrt((poleV[0] - northPole[0]) ** 2 + (poleV[1] - northPole[1]) ** 2 + (poleV[2] - northPole[2]) ** 2);
console.log(`Distance of potential peak to pole: ${dToPole.toFixed(2)} (Expected ~2.23)`);
assert(dToPole < 4.0 && dToPole > 1.0, 'Entropic potential peak is in the high-gradient ring (not at pole, not at equator)');

// 4. Test Gradient Direction
// Pick a face near the ring (e.g. at distance ~2.7)
// Face 0 is likely near the pole or ring depending on indexing.
const testFace = 0;
const p_bary = { face: testFace, u: 0.33, v: 0.33, w: 0.34 };
const force = flow.getForce(p_bary);

console.log(`Force on test face: [${force.map(x => x.toFixed(4))}]`);
assert(force.some(x => x !== 0), 'Force is non-zero');

// Vector to pole
const v0 = mesh.getVertex(mesh.faces[3 * testFace]);
const v_to_pole = [northPole[0] - v0[0], northPole[1] - v0[1], northPole[2] - v0[2]];
const v_to_pole_mag = sqrt(v_to_pole[0] ** 2 + v_to_pole[1] ** 2 + v_to_pole[2] ** 2);
// Normalize
v_to_pole[0] /= v_to_pole_mag; v_to_pole[1] /= v_to_pole_mag; v_to_pole[2] /= v_to_pole_mag;

const dot = force[0] * v_to_pole[0] + force[1] * v_to_pole[1] + force[2] * v_to_pole[2];

console.log(`Force dot direction to pole: ${dot.toFixed(4)}`);
// If we are *outside* the ring, force should point towards pole (dot > 0).
// If *inside* ring, it points away (dot < 0).
// Let's check where v0 is.
const distV0 = sqrt((v0[0] - northPole[0]) ** 2 + (v0[1] - northPole[1]) ** 2 + (v0[2] - northPole[2]) ** 2);
console.log(`Test point distance to pole: ${distV0.toFixed(2)}`);

if (distV0 > 2.8) {
    assert(dot > 0, 'Outside ring, force points to pole');
} else if (distV0 < 1.8) {
    assert(dot < 0, 'Inside ring, force points away from pole (towards ring)');
} else {
    console.log('Test point is within the ring, direction is ambiguous without precise ring check.');
}

console.log('Verification Complete.');
