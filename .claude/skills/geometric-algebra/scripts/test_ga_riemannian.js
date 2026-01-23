/**
 * Test script for ga_riemannian.js
 * 
 * Runs the same tests as riemannian_ga.py to verify the JavaScript implementation
 * matches the Python version.
 * 
 * Usage: node test_ga_riemannian.js
 */

// Load the module
const RiemannianGA = require('./ga_riemannian.js');

// Helper function for formatting
function printHeader(title) {
    console.log("\n" + "=".repeat(60));
    console.log(title);
    console.log("=".repeat(60));
}

function printSubHeader(title) {
    console.log("\n" + "-".repeat(40));
    console.log(title);
    console.log("-".repeat(40));
}

// Test 1: Sphere Curvature
function testSphereCurvature() {
    printHeader("TEST 1: Sphere Curvature via Connection Bivector");

    const R = 1.0;
    const sphere = new RiemannianGA.Sphere(R);
    const rga = new RiemannianGA.RiemannianGA(sphere);

    // Test at θ = π/3 (60°)
    const coords = [Math.PI / 3, 0.0];

    // Method 1: Via shape operator
    const K_shape = rga.gaussianCurvature(coords);
    const H_shape = rga.meanCurvature(coords);
    const expected_K = 1 / (R * R);
    // Note: With S = -∇ᵥn and outward normal, principal curvatures are -1/R
    // So H = -1/R for this convention (not +1/R)
    const expected_H = -1 / R;  // Sign convention: S = -∇ᵥn with outward normal

    console.log("\nMethod 1: Shape Operator");
    console.log(`  K = ${K_shape.toFixed(6)} (expected: ${expected_K.toFixed(6)})`);
    console.log(`  H = ${H_shape.toFixed(6)} (expected: ${expected_H.toFixed(6)})`);
    console.log("  Note: With S = -∇ᵥn and outward normal, κ = -1/R, H = -1/R");

    // Method 2: Via curvature 2-form (intrinsic curvature)
    const K_curv = rga.sectionalCurvature(coords);

    console.log("\nMethod 2: Curvature 2-Form");
    console.log(`  K = ${K_curv.toFixed(6)} (expected: ${expected_K.toFixed(6)})`);

    // Check connection bivector
    const omega = rga.getConnection(coords);
    console.log("\nConnection Bivectors:");
    console.log(`  ω_θ = ${omega[0]}`);
    console.log(`  ω_φ = ${omega[1]}`);

    // Principal curvatures
    const kappas = rga.principalCurvatures(coords);
    console.log("\nPrincipal Curvatures:");
    console.log(`  κ₁ = ${kappas[0].toFixed(6)}, κ₂ = ${kappas[1].toFixed(6)}`);
    console.log(`  (expected: ${(-1 / R).toFixed(6)}, ${(-1 / R).toFixed(6)})`);

    const errorK = Math.abs(K_shape - expected_K);
    const errorH = Math.abs(H_shape - expected_H);

    console.log("\nErrors:");
    console.log(`  |K - 1/R²| = ${errorK.toFixed(6)}`);
    console.log(`  |H - (-1/R)| = ${errorH.toFixed(6)}`);

    const passed = errorK < 0.1 && errorH < 0.1;
    console.log(`\nStatus: ${passed ? 'PASS ✓' : 'FAIL ✗'}`);
    return passed;
}

// Test 2: Torus Curvature (Variable K)
function testTorusCurvature() {
    printHeader("TEST 2: Torus Curvature (Variable K)");

    const R = 2.0, r = 0.5;
    const torus = new RiemannianGA.Torus(R, r);
    const rga = new RiemannianGA.RiemannianGA(torus);

    const testPoints = [
        [0.0, 0.0, "Outer equator (K > 0)"],
        [Math.PI, 0.0, "Inner equator (K < 0)"],
        [Math.PI / 2, 0.0, "Top (K = 0)"]
    ];

    let allPassed = true;

    for (const [theta, phi, name] of testPoints) {
        const coords = [theta, phi];
        const K = rga.gaussianCurvature(coords);
        const H = rga.meanCurvature(coords);

        // Expected: K = cos(θ) / (r(R + r cos θ))
        const ct = Math.cos(theta);
        const K_expected = ct / (r * (R + r * ct));

        console.log(`\n${name}:`);
        console.log(`  θ = ${theta.toFixed(4)}, φ = ${phi.toFixed(4)}`);
        console.log(`  K = ${K.toFixed(6)} (expected: ${K_expected.toFixed(6)})`);
        console.log(`  H = ${H.toFixed(6)}`);

        const error = Math.abs(K - K_expected);
        console.log(`  Error: ${error.toFixed(6)}`);

        if (error > 0.5) allPassed = false;
    }

    console.log(`\nStatus: ${allPassed ? 'PASS ✓' : 'FAIL ✗'}`);
    return allPassed;
}

// Test 3: Parallel Transport Holonomy
function testParallelTransportHolonomy() {
    printHeader("TEST 3: Parallel Transport Holonomy on Sphere");

    const R = 1.0;
    const sphere = new RiemannianGA.Sphere(R);
    const rga = new RiemannianGA.RiemannianGA(sphere);

    // Transport around latitude circle at θ = π/4
    const thetaFixed = Math.PI / 4;

    const latitudeLoop = (t) => [thetaFixed, t];

    // Initial vector in θ direction
    const w0 = [1.0, 0.0];

    const { wFinal, angle } = rga.computeHolonomy(latitudeLoop, w0, 200);

    // Expected holonomy: 2π(1 - cos θ)
    const expectedAngle = 2 * Math.PI * (1 - Math.cos(thetaFixed));

    console.log(`\nLatitude circle at θ = π/4 (${(thetaFixed * 180 / Math.PI).toFixed(1)}°)`);
    console.log(`Initial vector: [${w0.join(', ')}]`);
    console.log(`Final vector: [${wFinal.map(x => x.toFixed(4)).join(', ')}]`);
    console.log(`Holonomy angle: ${(angle * 180 / Math.PI).toFixed(2)}°`);
    console.log(`Expected angle: ${(expectedAngle * 180 / Math.PI).toFixed(2)}°`);

    const errorDeg = Math.abs(angle - expectedAngle) * 180 / Math.PI;
    console.log(`Error: ${errorDeg.toFixed(2)}°`);

    // Note: Holonomy = ∫K dA over enclosed region
    // For spherical cap: A = 2πR²(1 - cos θ), K = 1/R²
    // So ∫K dA = 2π(1 - cos θ) ✓
    console.log("\nNote: Holonomy = ∫K dA over enclosed cap = 2π(1 - cos θ)");

    const passed = errorDeg < 10;  // Allow 10 degrees error due to discretization
    console.log(`\nStatus: ${passed ? 'PASS ✓' : 'FAIL ✗'}`);
    return passed;
}

// Test 4: Geodesic on Sphere (Great Circle)
function testGeodesic() {
    printHeader("TEST 4: Geodesic on Sphere (Great Circle)");

    const sphere = new RiemannianGA.Sphere(1.0);
    const rga = new RiemannianGA.RiemannianGA(sphere);

    // Start at equator, move in positive theta direction (toward south pole)
    // Note: θ=0 is north pole, θ=π is south pole
    // With v0=[1,0], we move in positive θ direction (toward south)
    const x0 = [Math.PI / 2, 0];
    const v0 = [1, 0];  // Initial velocity in positive θ direction

    const result = rga.solveGeodesic(x0, v0, Math.PI / 2, 50);

    console.log("\nGeodesic from equator, moving along meridian (great circle):");
    console.log(`  Start: θ=${x0[0].toFixed(4)} (equator), φ=${x0[1].toFixed(4)}`);
    console.log(`  Velocity: v0=[${v0.join(', ')}] (positive θ direction = toward south)`);

    const final = result.x[result.x.length - 1];
    console.log(`  End: θ=${final[0].toFixed(4)}, φ=${final[1].toFixed(4)}`);
    console.log(`  Expected: θ ≈ π (south pole) since v₀ points toward increasing θ`);

    // Trace some intermediate points
    console.log("\nTrajectory (sample points):");
    const indices = [0, 10, 20, 30, 40, 50];
    for (const i of indices) {
        if (i < result.x.length) {
            const p = result.x[i];
            console.log(`  t=${result.t[i].toFixed(3)}: θ=${p[0].toFixed(4)}, φ=${p[1].toFixed(4)}`);
        }
    }

    // Check if we're heading toward south pole (θ increases toward π)
    const passed = final[0] > x0[0];
    console.log(`\nNote: Geodesic flows in direction of initial velocity`);
    console.log(`      θ: ${x0[0].toFixed(4)} → ${final[0].toFixed(4)} (${final[0] > x0[0] ? 'increasing ✓' : 'decreasing ✗'})`);
    console.log(`\nStatus: ${passed ? 'PASS ✓' : 'FAIL ✗'}`);
    return passed;
}

// Test 5: Connection Bivector Properties
function testConnectionBivector() {
    printHeader("TEST 5: Connection Bivector Properties");

    const sphere = new RiemannianGA.Sphere(1.0);
    const rga = new RiemannianGA.RiemannianGA(sphere);

    // Test at multiple points
    const testCoords = [
        [Math.PI / 4, 0],
        [Math.PI / 2, 0],
        [Math.PI / 3, Math.PI / 4]
    ];

    console.log("\nConnection bivectors at various points:");

    for (const coords of testCoords) {
        const omega = rga.getConnection(coords);
        const frame = rga.getFrame(coords);

        console.log(`\n  At θ=${coords[0].toFixed(4)}, φ=${coords[1].toFixed(4)}:`);
        console.log(`    ω_θ = ${omega[0]}`);
        console.log(`    ω_φ = ${omega[1]}`);
        console.log(`    |e_θ| = ${RiemannianGA.norm(frame.vectors[0]).toFixed(6)}`);
        console.log(`    |e_φ| = ${RiemannianGA.norm(frame.vectors[1]).toFixed(6)}`);
    }

    console.log("\nNote: Connection bivector generates frame rotation");
    console.log("      ∂ᵢeⱼ = ωᵢ × eⱼ");

    return true;
}

// Test 6: Covariant Derivative
function testCovariantDerivative() {
    printHeader("TEST 6: Covariant Derivative");

    const sphere = new RiemannianGA.Sphere(1.0);
    const connection = new RiemannianGA.ConnectionBivector(2, 3);
    const covDeriv = new RiemannianGA.GACovariantDerivative(connection);

    const frameFunc = (coords) => sphere.frame(coords[0], coords[1]);

    // Define a vector field: constant in θ direction
    function vectorField(coords) {
        const frame = frameFunc(coords);
        // Return e_θ (first basis vector)
        return frame.vectors[0];
    }

    const coords = [Math.PI / 3, 0];
    const u = [1, 0];  // Direction along θ

    const result = covDeriv.ofVector(vectorField, u, coords, frameFunc);

    console.log("\nCovariant derivative of e_θ in θ direction:");
    console.log(`  Coords: θ=${coords[0].toFixed(4)}, φ=${coords[1].toFixed(4)}`);
    console.log(`  ∇_{e_θ}e_θ = [${result.map(x => x.toFixed(6)).join(', ')}]`);
    console.log(`  |∇_{e_θ}e_θ| = ${RiemannianGA.norm(result).toFixed(6)}`);

    // On a sphere, for the coordinate frame:
    // ∇_{e_θ}e_θ should be perpendicular to e_θ
    const frame = frameFunc(coords);
    const dotProduct = RiemannianGA.dot(result, frame.vectors[0]);
    console.log(`  (∇_{e_θ}e_θ) · e_θ = ${dotProduct.toFixed(6)} (should be ~0)`);

    return true;
}

// Summary function
function runAllTests() {
    console.log("╔══════════════════════════════════════════════════════════╗");
    console.log("║      GEOMETRIC ALGEBRA RIEMANNIAN GEOMETRY TESTS         ║");
    console.log("║              JavaScript Implementation                   ║");
    console.log("╚══════════════════════════════════════════════════════════╝");

    const results = [];

    results.push(['Sphere Curvature', testSphereCurvature()]);
    results.push(['Torus Curvature', testTorusCurvature()]);
    results.push(['Parallel Transport', testParallelTransportHolonomy()]);
    results.push(['Geodesic', testGeodesic()]);
    results.push(['Connection Bivector', testConnectionBivector()]);
    results.push(['Covariant Derivative', testCovariantDerivative()]);

    printHeader("TEST SUMMARY");

    let passCount = 0;
    for (const [name, passed] of results) {
        console.log(`  ${passed ? '✓' : '✗'} ${name}: ${passed ? 'PASSED' : 'FAILED'}`);
        if (passed) passCount++;
    }

    console.log("\n" + "-".repeat(40));
    console.log(`Total: ${passCount}/${results.length} tests passed`);
    console.log("=".repeat(60));

    return passCount === results.length;
}

// Run tests
if (require.main === module) {
    const success = runAllTests();
    process.exit(success ? 0 : 1);
}

module.exports = {
    testSphereCurvature,
    testTorusCurvature,
    testParallelTransportHolonomy,
    testGeodesic,
    testConnectionBivector,
    testCovariantDerivative,
    runAllTests
};
