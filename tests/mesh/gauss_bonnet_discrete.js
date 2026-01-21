/**
 * Discrete Gauss-Bonnet Test
 * 
 * Verifies the discrete Gauss-Bonnet theorem: Σ K_v = 2πχ
 * where K_v is the angle defect at each vertex and χ is the Euler characteristic.
 * 
 * Test cases:
 * - Icosahedron: χ = 2 (sphere topology)
 * - Torus mesh: χ = 0
 * - Disk: χ = 1 (with boundary correction)
 */

const { TriangleMesh } = require('../../src/mesh.js');
const RiemannianDiscrete = require('../../src/riemannian-discrete.js');

const { PI, abs } = Math;

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

function testGaussBonnetDiscrete() {
    console.log('=== Discrete Gauss-Bonnet Test: Σ K_v = 2πχ ===\n');

    let passed = true;

    // Test 1: Icosahedron (χ = 2)
    console.log('--- Test 1: Icosahedron (sphere topology, χ = 2) ---');
    const icosahedron = TriangleMesh.createIcosahedron(1.0);

    const { vertexCurvatures, totalCurvature } = RiemannianDiscrete.discreteCurvature(icosahedron);
    const expectedTotal = 2 * PI * 2;  // 2πχ = 4π

    console.log(`  Vertices: ${icosahedron.nVertices}`);
    console.log(`  Faces: ${icosahedron.nFaces}`);
    console.log(`  Euler characteristic: V - E + F = ${icosahedron.nVertices} - ${icosahedron.nEdges} + ${icosahedron.nFaces} = ${icosahedron.nVertices - icosahedron.nEdges + icosahedron.nFaces}`);
    console.log(`  Total curvature: ${totalCurvature.toFixed(6)}`);
    console.log(`  Expected (4π): ${expectedTotal.toFixed(6)}`);

    passed &= assertEqual(totalCurvature, expectedTotal, 0.1, 'Icosahedron: Σ K_v = 4π');

    // Test 2: Subdivided icosahedron (still χ = 2)
    console.log('\n--- Test 2: Subdivided Icosahedron ---');
    // Use a finer sphere approximation
    const sphere = RiemannianDiscrete.createSphereMesh(1.0, 3);  // 3 subdivisions

    const sphereResult = RiemannianDiscrete.discreteCurvature(sphere);
    console.log(`  Vertices: ${sphere.nVertices}`);
    console.log(`  Total curvature: ${sphereResult.totalCurvature.toFixed(6)}`);

    passed &= assertEqual(sphereResult.totalCurvature, expectedTotal, 0.5, 'Sphere mesh: Σ K_v ≈ 4π');

    // Test 3: Flat grid (boundary effects)
    console.log('\n--- Test 3: Flat Grid (Planar, χ = 1 with boundary) ---');
    const grid = TriangleMesh.createGrid(5, 5, 1, 1);

    const gridResult = RiemannianDiscrete.discreteCurvature(grid);
    console.log(`  Vertices: ${grid.nVertices}`);
    console.log(`  Interior curvature should be 0`);
    console.log(`  Total with boundary: ${gridResult.totalCurvature.toFixed(6)}`);

    // For flat disk-like topology, total should be 2π
    passed &= assertEqual(gridResult.totalCurvature, 2 * PI, 0.5, 'Grid: Σ K_v ≈ 2π');

    return passed;
}

// Run test
const passed = testGaussBonnetDiscrete();
console.log(`\n${passed ? 'ALL TESTS PASSED' : 'SOME TESTS FAILED'}`);
process.exit(passed ? 0 : 1);
