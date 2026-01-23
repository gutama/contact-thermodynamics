/**
 * Verification for Spacetime Geometric Algebra
 * 
 * Tests the SpacetimeManifoldGA class on:
 * 1. Minkowski Space (Flat, should have zero connection)
 * 2. Schwarzschild Metric (Curved, should have non-zero connection)
 */

const { SpacetimeManifoldGA } = require('../src/riemannian-spacetime.js');
const { abs, sqrt, sin, cos, PI } = Math;

function assert(condition, msg) {
    if (!condition) {
        console.error(`FAILED: ${msg}`);
        process.exit(1);
    }
    console.log(`PASS: ${msg}`);
}

function approx(a, b, tol = 1e-6) {
    return abs(a - b) < tol;
}

console.log('='.repeat(50));
console.log('VERIFYING SPACETIME GEOMETRIC ALGEBRA');
console.log('='.repeat(50));

// ============================================================================
// 1. Minkowski Space (Flat)
// ============================================================================
console.log('\n1. MINKOWSKI SPACE');

const minkowski = new SpacetimeManifoldGA(x => [
    [1, 0, 0, 0],
    [0, -1, 0, 0],
    [0, 0, -1, 0],
    [0, 0, 0, -1]
]);

// Point check
const p0 = [0, 0, 0, 0];
const tetradFlat = minkowski.computeTetrad(p0);

// Check tetrad orthogonality e^a . e^b = eta^ab (for inverse vectors)
// Our computeTetrad returns E^a_u (forms) and e_a^u (vectors)
// We check e_a^u
const e_vec = tetradFlat.vectors; // e_a^u (rows are a, cols are u)

// e_0 = [1,0,0,0], e_1 = [0,1,0,0] ...
console.log('Tetrad Vectors (Minkowski):');
console.log(e_vec);

assert(approx(e_vec[0][0], 1), 'e_0^0 = 1');
assert(approx(e_vec[1][1], 1), 'e_1^1 = 1');

// Connection should be zero
const omegaFlat = minkowski.connectionBivector(p0);
let maxMag = 0;
omegaFlat.forEach(B => {
    // Check magnitude of bivector coefficients
    for (const k in B.coeffs) {
        maxMag = Math.max(maxMag, abs(B.coeffs[k]));
    }
});
console.log(`Max connection component (Flat): ${maxMag.toFixed(9)}`);
assert(maxMag < 1e-9, 'Connection vanishes in flat space');


// ============================================================================
// 2. Schwarzschild Metric (Curved)
// ============================================================================
console.log('\n2. SCHWARZSCHILD METRIC (M=1)');

const M = 1.0;
const schwarzschild = new SpacetimeManifoldGA(x => {
    // x = [t, r, theta, phi]
    const r = x[1];
    const theta = x[2];
    const f = 1.0 - 2.0 * M / r;

    // Diagonal metric
    // g_tt = f, g_rr = -1/f, g_thth = -r^2, g_phph = -r^2 sin^2 th
    // Avoid horizon singularity for test
    return [
        [f, 0, 0, 0],
        [0, -1.0 / f, 0, 0],
        [0, 0, -r * r, 0],
        [0, 0, 0, -r * r * sin(theta) * sin(theta)]
    ];
});

// Test point outside horizon
const p1 = [0, 5.0, PI / 2, 0]; // r=5, equatorial plane
console.log(`Point: r=${p1[1]}, theta=${p1[2]}`);

const omegas = schwarzschild.connectionBivector(p1);

// We expect specific non-zero components for Schwarzschild in standard tetrad
// Non-zero gamma coeff: gamma_010 = M/r^2 (acceleration due to gravity)
// This corresponds to omega_t (a=0) having a e^1^e^0 term? Or omega_0 having 01 term? 
// Spin connection has indices omega^ab_mu.
// The object returned is omega_mu (bivector).
// omega_t should contain e_1 ^ e_0 term proportional to M/r^2.

const omega_t = omegas[0]; // Time component
console.log('omega_t (time component of connection):');
console.log(omega_t.toString());

// Check for e1^e0 term (or e0^e1)
// In our basis, e0=1, e1=2 (bitmap 1, 2). e0^e1 is bitmap 3.
// e1^e0 = -e0^e1.
// Expected value: approx M/r^2 = 1/25 = 0.04
// Let's inspect coefficients.

// Note: Connection coefficients depend on the tetrad choice.
// computeTetrad uses diagonal sqrt.
// e^0 = sqrt(f) dt, e^1 = 1/sqrt(f) dr, e^2 = r dth, e^3 = r sin th dph
// This is the standard static observer frame.
// Known result: omega^0_1 = (M/r^2) e^0  (1-form) ?? No.
// omega^0_1 = (M/r^2) dt ??
// omega_t has component e0^e1 of magnitude related to gravitational force.

// Let's just assert it is non-zero
let magT = omega_t.norm();
console.log(`omega_t norm: ${magT.toFixed(6)}`);
assert(magT > 0.03 && magT < 0.05, 'Time connection has magnitude ~ M/r^2');

// omega_phi (a=3) should have term related to cos(theta) and sin(theta)
// On equator (theta=pi/2), cot(theta)=0, so some terms vanish
// r=5, theta=pi/2.
const omega_phi = omegas[3];
console.log('omega_phi (azimuthal component):');
console.log(omega_phi.toString());

console.log('\nVerification Complete.');
