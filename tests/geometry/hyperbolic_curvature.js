/**
 * Hyperbolic Plane Curvature Test
 * 
 * Verifies constant negative curvature K = -1 for the hyperbolic plane (H²).
 * Uses the upper half-plane model with metric ds² = (dx² + dy²) / y².
 */

const RiemannianGA = require('../../src/riemannian-ga.js');

const { HyperbolicPlane, Curvature2Form, EPSILON } = RiemannianGA;
const { abs, random } = Math;

function assertEqual(actual, expected, tolerance, message) {
    const diff = abs(actual - expected);
    if (diff > tolerance) {
        console.error(`FAIL: ${message}`);
        console.error(`  Expected: ${expected}`);
        console.error(`  Actual: ${actual}`);
        console.error(`  Difference: ${diff}`);
        return false;
    }
    console.log(`PASS: ${message}`);
    return true;
}

function testHyperbolicCurvature() {
    console.log('=== Hyperbolic Plane Curvature Test: K = -1 ===\n');

    let passed = true;
    const H2 = new HyperbolicPlane();
    const expectedK = -1;

    console.log('Testing at various points in upper half-plane:');

    // Test at multiple random points (y > 0)
    const nTests = 10;
    for (let i = 0; i < nTests; i++) {
        const x = (random() - 0.5) * 10;  // x ∈ [-5, 5]
        const y = 0.5 + random() * 5;     // y ∈ [0.5, 5.5]
        const coords = [x, y];

        const curvature = new Curvature2Form(H2);
        const K = curvature.gaussianCurvature(coords);

        passed &= assertEqual(K, expectedK, 0.01,
            `K at (${x.toFixed(2)}, ${y.toFixed(2)}) = -1`);
    }

    // Test at specific points
    console.log('\n--- Specific Points ---');

    const testPoints = [
        [0, 1],
        [0, 2],
        [1, 1],
        [-1, 0.5],
        [3, 3]
    ];

    for (const [x, y] of testPoints) {
        const curvature = new Curvature2Form(H2);
        const K = curvature.gaussianCurvature([x, y]);
        passed &= assertEqual(K, expectedK, 0.01, `K at (${x}, ${y}) = -1`);
    }

    // Verify metric scaling with y
    console.log('\n--- Metric Scaling ---');
    const frame1 = H2.frame([0, 1]);
    const frame2 = H2.frame([0, 2]);

    const g1 = frame1.metric();
    const g2 = frame2.metric();

    // At y=1: g = I, at y=2: g = I/4
    passed &= assertEqual(g1[0][0], 1, 0.01, 'g_xx at y=1 is 1');
    passed &= assertEqual(g2[0][0], 0.25, 0.01, 'g_xx at y=2 is 0.25');

    return passed;
}

// Run test
const passed = testHyperbolicCurvature();
console.log(`\n${passed ? 'ALL TESTS PASSED' : 'SOME TESTS FAILED'}`);
process.exit(passed ? 0 : 1);
