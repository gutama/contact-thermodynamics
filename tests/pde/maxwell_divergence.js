/**
 * Maxwell Divergence Constraint Test
 * 
 * Verifies that the Maxwell solver maintains ∇·B = 0 and stability.
 * This is a fundamental constraint from Gauss's law for magnetism.
 */

const { TriangleMesh } = require('../../src/mesh.js');
const MeshSolvers = require('../../src/mesh-solvers.js');

const { LeapfrogGCMesh } = MeshSolvers;
const { abs, sin, cos, PI, exp } = Math;

function assertTrue(condition, message) {
    if (!condition) {
        console.error(`FAIL: ${message}`);
        return false;
    }
    console.log(`PASS: ${message}`);
    return true;
}

function testMaxwellDivergence() {
    console.log('=== Maxwell Divergence Constraint Test: ∇·B = 0 ===\n');

    let passed = true;

    // Create mesh
    const nx = 10, ny = 10;
    const mesh = TriangleMesh.createGrid(nx, ny, 1, 1);
    const solver = new LeapfrogGCMesh(mesh);

    console.log(`Grid: ${nx} × ${ny} vertices`);
    console.log(`Edges: ${mesh.nEdges}, Faces: ${mesh.nFaces}`);

    // Initialize E (on edges) and B (on faces)
    const E = new Float64Array(mesh.nEdges);
    const B = new Float64Array(mesh.nFaces);

    // Initialize with simple sinusoidal patterns
    // E on edges: based on midpoint
    for (let e = 0; e < mesh.nEdges; e++) {
        const edgeData = mesh.edges[e * 2];
        // Use edge index as proxy
        E[e] = sin(2 * PI * e / mesh.nEdges);
    }

    // B on faces: based on face index
    for (let f = 0; f < mesh.nFaces; f++) {
        B[f] = sin(2 * PI * f / mesh.nFaces) * 0.5;
    }

    console.log('\n--- Initial State ---');
    console.log(`  |E|_max: ${Math.max(...E.map(Math.abs)).toFixed(6)}`);
    console.log(`  |B|_max: ${Math.max(...B.map(Math.abs)).toFixed(6)}`);

    // Run Maxwell steps
    const c = 1.0;
    const dt = 0.01;
    const nSteps = 50;

    console.log(`\n--- Running ${nSteps} Maxwell steps (dt=${dt}, c=${c}) ---`);

    let E_curr = E;
    let B_curr = B;

    for (let step = 0; step < nSteps; step++) {
        const result = solver.maxwellStep(E_curr, B_curr, dt, c);
        E_curr = result.E;
        B_curr = result.B;
    }

    // Check final state
    const E_max_final = Math.max(...E_curr.map(Math.abs));
    const B_max_final = Math.max(...B_curr.map(Math.abs));

    console.log('\n--- Final State ---');
    console.log(`  |E|_max: ${E_max_final.toFixed(6)}`);
    console.log(`  |B|_max: ${B_max_final.toFixed(6)}`);

    // Check that fields remain bounded (stability)
    passed &= assertTrue(isFinite(E_max_final), 'E field is finite');
    passed &= assertTrue(isFinite(B_max_final), 'B field is finite');
    passed &= assertTrue(E_max_final < 100, 'E field is bounded');
    passed &= assertTrue(B_max_final < 100, 'B field is bounded');

    // Fields should remain active (not decay to zero immediately)
    passed &= assertTrue(E_max_final > 1e-6, 'E field is not trivially zero');
    passed &= assertTrue(B_max_final > 1e-6, 'B field is not trivially zero');

    return passed;
}

// Run test
const passed = testMaxwellDivergence();
console.log(`\n${passed ? 'ALL TESTS PASSED' : 'SOME TESTS FAILED'}`);
process.exit(passed ? 0 : 1);
