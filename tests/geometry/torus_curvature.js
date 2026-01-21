/**
 * Torus Curvature Test
 * 
 * Verifies Gaussian curvature formula K = cos(θ) / (r(R + r cos(θ)))
 * where R is major radius, r is minor radius.
 * 
 * Key facts:
 * - K > 0 on outer part (θ near 0)
 * - K < 0 on inner part (θ near π)
 * - K = 0 on top/bottom (θ = ±π/2)
 * - Total curvature ∫∫K dA = 0 (Gauss-Bonnet for torus)
 */

const RiemannianGA = require('../../src/riemannian-ga.js');

const { Torus2D, Curvature2Form, EPSILON } = RiemannianGA;
const { abs, PI, cos, random } = Math;

function assertEqual(actual, expected, tolerance, message) {
    const diff = abs(actual - expected);
    if (diff > tolerance) {
        console.error(`FAIL: ${message}`);
        console.error(`  Expected: ${expected.toFixed(6)}`);
        console.error(`  Actual: ${actual.toFixed(6)}`);
        return false;
    }
    console.log(`PASS: ${message}`);
    return true;
}

function assertTrue(condition, message) {
    if (!condition) {
        console.error(`FAIL: ${message}`);
        return false;
    }
    console.log(`PASS: ${message}`);
    return true;
}

/**
 * Analytical Gaussian curvature for torus.
 */
function torusCurvature(R, r, theta) {
    return cos(theta) / (r * (R + r * cos(theta)));
}

function testTorusCurvature() {
    console.log('=== Torus Curvature Test: K = cos(θ) / (r(R + r cos(θ))) ===\n');

    let passed = true;
    const R = 3.0;  // Major radius
    const r = 1.0;  // Minor radius
    const torus = new Torus2D(R, r);

    console.log(`Torus parameters: R = ${R}, r = ${r}`);

    // Test at specific angles
    const testCases = [
        { theta: 0, desc: 'Outer equator' },
        { theta: PI / 2, desc: 'Top' },
        { theta: PI, desc: 'Inner equator' },
        { theta: 3 * PI / 2, desc: 'Bottom' }
    ];

    for (const { theta, desc } of testCases) {
        const phi = 0;  // Fix azimuthal angle
        const expectedK = torusCurvature(R, r, theta);

        const curvature = new Curvature2Form(torus);
        const K = curvature.gaussianCurvature([theta, phi]);

        console.log(`\n${desc} (θ = ${(theta * 180 / PI).toFixed(0)}°):`);
        console.log(`  Expected K: ${expectedK.toFixed(6)}`);
        console.log(`  Computed K: ${K.toFixed(6)}`);

        passed &= assertEqual(K, expectedK, 0.05, `Curvature at ${desc}`);
    }

    // Verify curvature signs
    console.log('\n--- Curvature Sign Tests ---');

    const curvature = new Curvature2Form(torus);

    // Outer part: K > 0
    const K_outer = curvature.gaussianCurvature([0.1, 0]);
    passed &= assertTrue(K_outer > 0, 'K > 0 on outer part');

    // Inner part: K < 0
    const K_inner = curvature.gaussianCurvature([PI, 0]);
    passed &= assertTrue(K_inner < 0, 'K < 0 on inner part');

    // Top/bottom: K ≈ 0
    const K_top = curvature.gaussianCurvature([PI / 2, 0]);
    passed &= assertEqual(K_top, 0, 0.01, 'K ≈ 0 at top');

    // Test Gauss-Bonnet: ∫∫K dA = 2πχ = 0 for torus
    console.log('\n--- Gauss-Bonnet Integral ---');
    console.log('For torus, ∫∫K dA should equal 0');

    let totalCurvature = 0;
    const nTheta = 50;
    const nPhi = 50;
    const dTheta = 2 * PI / nTheta;
    const dPhi = 2 * PI / nPhi;

    for (let i = 0; i < nTheta; i++) {
        const theta = i * dTheta;
        for (let j = 0; j < nPhi; j++) {
            const phi = j * dPhi;
            const K = curvature.gaussianCurvature([theta, phi]);
            // Area element: dA = r(R + r cos θ) dθ dφ
            const dA = r * (R + r * cos(theta)) * dTheta * dPhi;
            totalCurvature += K * dA;
        }
    }

    console.log(`  Total curvature: ${totalCurvature.toFixed(6)}`);
    passed &= assertEqual(totalCurvature, 0, 0.1, 'Gauss-Bonnet: ∫∫K dA = 0');

    return passed;
}

// Run test
const passed = testTorusCurvature();
console.log(`\n${passed ? 'ALL TESTS PASSED' : 'SOME TESTS FAILED'}`);
process.exit(passed ? 0 : 1);
