/**
 * Sphere Curvature Test
 * 
 * Verifies Gaussian curvature K = 1/R² for 2D sphere (S²) embedded in R³.
 * Tests at multiple random points to ensure curvature is constant.
 */

const RiemannianGA = require('../../src/riemannian-ga.js');

const { Sphere2D, Curvature2Form, EPSILON } = RiemannianGA;
const { abs, PI, random } = Math;

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

function testSphereCurvature() {
    console.log('=== Sphere Curvature Test: K = 1/R² ===\n');

    let passed = true;
    const radii = [0.5, 1.0, 2.0, 5.0];

    for (const R of radii) {
        const sphere = new Sphere2D(R);
        const expectedK = 1 / (R * R);

        console.log(`\nTesting R = ${R}, expected K = ${expectedK.toFixed(6)}`);

        // Test at multiple random points (avoid poles: θ ∈ (0.2, π-0.2))
        const nTests = 5;
        for (let i = 0; i < nTests; i++) {
            const theta = 0.2 + random() * (PI - 0.4);
            const phi = random() * 2 * PI;
            const coords = [theta, phi];

            const curvature = new Curvature2Form(sphere);
            const K = curvature.gaussianCurvature(coords);

            passed &= assertEqual(K, expectedK, 0.01,
                `K at (θ=${theta.toFixed(2)}, φ=${phi.toFixed(2)}) = 1/R²`);
        }

        // Test at specific known point: equator
        const K_equator = new Curvature2Form(sphere).gaussianCurvature([PI / 2, 0]);
        passed &= assertEqual(K_equator, expectedK, 0.01, `K at equator = 1/R²`);
    }

    return passed;
}

// Run test
const passed = testSphereCurvature();
console.log(`\n${passed ? 'ALL TESTS PASSED' : 'SOME TESTS FAILED'}`);
process.exit(passed ? 0 : 1);
