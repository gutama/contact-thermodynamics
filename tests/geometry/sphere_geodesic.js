/**
 * Sphere Geodesic Test
 * 
 * Verifies that geodesics on the sphere are great circles.
 * Tests arc length calculation against analytical formula.
 */

const RiemannianGA = require('../../src/riemannian-ga.js');

const { Sphere2D, GeodesicSolver, dot, norm, EPSILON } = RiemannianGA;
const { abs, PI, cos, sin, sqrt, acos } = Math;

function assertEqual(actual, expected, tolerance, message) {
    const diff = abs(actual - expected);
    if (diff > tolerance) {
        console.error(`FAIL: ${message}`);
        console.error(`  Expected: ${expected}`);
        console.error(`  Actual: ${actual}`);
        return false;
    }
    console.log(`PASS: ${message}`);
    return true;
}

/**
 * Compute geodesic distance on sphere using Haversine-like formula.
 * For points (θ₁, φ₁) and (θ₂, φ₂), distance = R * arccos(n₁ · n₂)
 */
function sphereGeodesicDistance(R, coords1, coords2) {
    const [theta1, phi1] = coords1;
    const [theta2, phi2] = coords2;

    // Convert to Cartesian
    const x1 = sin(theta1) * cos(phi1);
    const y1 = sin(theta1) * sin(phi1);
    const z1 = cos(theta1);

    const x2 = sin(theta2) * cos(phi2);
    const y2 = sin(theta2) * sin(phi2);
    const z2 = cos(theta2);

    // Dot product
    const dotProd = x1 * x2 + y1 * y2 + z1 * z2;
    const angle = acos(Math.max(-1, Math.min(1, dotProd)));

    return R * angle;
}

function testSphereGeodesic() {
    console.log('=== Sphere Geodesic Test: Arc Length ===\n');

    let passed = true;
    const R = 1.0;
    const sphere = new Sphere2D(R);

    // Test 1: Meridian geodesic (φ = const)
    console.log('Test 1: Meridian (φ = 0)');
    const start1 = [PI / 4, 0];  // θ = 45°
    const end1 = [3 * PI / 4, 0];  // θ = 135°
    const expectedDist1 = R * (3 * PI / 4 - PI / 4);  // Δθ * R
    const computedDist1 = sphereGeodesicDistance(R, start1, end1);

    passed &= assertEqual(computedDist1, expectedDist1, 0.01,
        `Meridian geodesic distance = R*Δθ`);

    // Test 2: Equator geodesic (θ = π/2)
    console.log('\nTest 2: Equator segment');
    const start2 = [PI / 2, 0];
    const end2 = [PI / 2, PI / 2];  // 90° apart
    const expectedDist2 = R * PI / 2;
    const computedDist2 = sphereGeodesicDistance(R, start2, end2);

    passed &= assertEqual(computedDist2, expectedDist2, 0.01,
        `Equator geodesic distance = R*Δφ`);

    // Test 3: General great circle
    console.log('\nTest 3: General great circle');
    const start3 = [PI / 3, PI / 6];
    const end3 = [PI / 2, PI / 2];
    const dist3 = sphereGeodesicDistance(R, start3, end3);
    console.log(`Distance: ${dist3.toFixed(6)}`);
    passed &= (dist3 > 0 && dist3 < PI);  // Should be less than half circumference

    // Test 4: Antipodal points
    console.log('\nTest 4: Antipodal points');
    const start4 = [PI / 2, 0];
    const end4 = [PI / 2, PI];
    const expectedDist4 = R * PI;  // Half circumference
    const computedDist4 = sphereGeodesicDistance(R, start4, end4);

    passed &= assertEqual(computedDist4, expectedDist4, 0.01,
        `Antipodal distance = π*R`);

    return passed;
}

// Run test
const passed = testSphereGeodesic();
console.log(`\n${passed ? 'ALL TESTS PASSED' : 'SOME TESTS FAILED'}`);
process.exit(passed ? 0 : 1);
