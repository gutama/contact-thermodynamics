/**
 * Discrete Gauss-Bonnet Test
 * 
 * Verifies the discrete Gauss-Bonnet theorem: Σ K_v = 2πχ
 * where K_v is the angle defect at each vertex and χ is the Euler characteristic.
 * 
 * Test cases:
 * - Icosahedron: χ = 2 (sphere topology)
 * - Flat grid: χ = 1 (disk topology with boundary)
 */

const { TriangleMesh } = require('../../src/mesh.js');
const RiemannianDiscrete = require('../../src/riemannian-discrete.js');

const { MeshCurvature2Form } = RiemannianDiscrete;
const { PI, abs } = Math;

function assertEqual(actual, expected, tolerance, message) {
    const diff = abs(actual - expected);
    if (diff > tolerance) {
        console.error(`FAIL: ${message}`);
        console.error(`  Expected: ${expected.toFixed(6)}`);
        console.error(`  Actual: ${actual.toFixed(6)}`);
        console.error(`  Difference: ${diff.toFixed(6)}`);
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

    const curvature = new MeshCurvature2Form(icosahedron);
    const totalK = curvature.totalCurvature();
    const chi = curvature.eulerCharacteristic();
    const expectedTotal = 2 * PI * chi;  // 2πχ

    console.log(`  Vertices: ${icosahedron.nVertices}`);
    console.log(`  Edges: ${icosahedron.nEdges}`);
    console.log(`  Faces: ${icosahedron.nFaces}`);
    console.log(`  Euler characteristic χ = V - E + F = ${chi}`);
    console.log(`  Total curvature Σ K_v: ${totalK.toFixed(6)}`);
    console.log(`  Expected (2πχ = ${(2 * chi).toFixed(0)}π): ${expectedTotal.toFixed(6)}`);

    passed &= assertEqual(totalK, expectedTotal, 0.01, 'Icosahedron: Σ K_v = 4π');
    passed &= assertEqual(chi, 2, 0.001, 'Icosahedron Euler characteristic χ = 2');

    // Test 2: Flat grid (disk topology, χ = 1)
    console.log('\n--- Test 2: Flat Grid (disk topology, χ = 1) ---');
    const grid = TriangleMesh.createGrid(5, 5, 1, 1);

    const gridCurvature = new MeshCurvature2Form(grid);
    const gridTotalK = gridCurvature.totalCurvature();
    const gridChi = gridCurvature.eulerCharacteristic();
    const gridExpected = 2 * PI * gridChi;

    console.log(`  Vertices: ${grid.nVertices}`);
    console.log(`  Edges: ${grid.nEdges}`);
    console.log(`  Faces: ${grid.nFaces}`);
    console.log(`  Euler characteristic χ = V - E + F = ${gridChi}`);
    console.log(`  Total curvature Σ K_v: ${gridTotalK.toFixed(6)}`);
    console.log(`  Expected (2πχ): ${gridExpected.toFixed(6)}`);

    passed &= assertEqual(gridTotalK, gridExpected, 0.5, 'Flat grid: Σ K_v = 2πχ');

    // Test 3: Verify angle defects at individual vertices
    console.log('\n--- Test 3: Angle Defects Distribution ---');
    const defects = curvature.angleDefects();

    // For icosahedron, all vertices have same curvature (regular polyhedron)
    const firstDefect = defects[0];
    let allSame = true;
    for (let i = 1; i < defects.length; i++) {
        if (abs(defects[i] - firstDefect) > 0.001) {
            allSame = false;
            break;
        }
    }

    console.log(`  First vertex defect: ${firstDefect.toFixed(6)} rad`);
    console.log(`  Expected per vertex: ${(4 * PI / 12).toFixed(6)} rad (= π/3)`);

    passed &= assertEqual(firstDefect, PI / 3, 0.01, 'Icosahedron vertex defect = π/3');
    if (allSame) console.log('  PASS: All vertices have equal curvature (regular)');

    return passed;
}

// Run test
const passed = testGaussBonnetDiscrete();
console.log(`\n${passed ? 'ALL TESTS PASSED' : 'SOME TESTS FAILED'}`);
process.exit(passed ? 0 : 1);
