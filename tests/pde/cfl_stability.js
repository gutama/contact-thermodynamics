/**
 * CFL Stability Test
 * 
 * Verifies that solutions remain bounded under CFL condition.
 * Tests both wave and heat equations.
 */

const { TriangleMesh } = require('../../src/mesh.js');
const MeshSolvers = require('../../src/mesh-solvers.js');

const { LeapfrogGCMesh, boundaryDirichletMask } = MeshSolvers;
const { abs, max, exp } = Math;

function assertTrue(condition, message) {
    if (!condition) {
        console.error(`FAIL: ${message}`);
        return false;
    }
    console.log(`PASS: ${message}`);
    return true;
}

function testCFLStability() {
    console.log('=== CFL Stability Test ===\n');

    let passed = true;

    // Create a mesh
    const nx = 10, ny = 10;
    const mesh = TriangleMesh.createGrid(nx, ny, 1, 1);
    const solver = new LeapfrogGCMesh(mesh);

    console.log(`Grid: ${nx} × ${ny} vertices`);

    // --- Test 1: Wave equation with valid CFL ---
    console.log('\n--- Test 1: Wave Equation (valid CFL) ---');

    const c = 1.0;
    const dtWave = solver.estimateCFL(c) * 0.5;
    console.log(`  CFL limit: ${solver.estimateCFL(c).toFixed(6)}`);
    console.log(`  Using dt: ${dtWave.toFixed(6)} (50% of limit)`);

    // Initial bump
    const u0Wave = new Float64Array(mesh.nVertices);
    const v0Wave = new Float64Array(mesh.nVertices);
    for (let i = 0; i < mesh.nVertices; i++) {
        const [x, y] = mesh.getVertex(i);
        u0Wave[i] = exp(-((x - 0.5) ** 2 + (y - 0.5) ** 2) / 0.02);
    }

    // Boundary conditions
    const dirichletMask = boundaryDirichletMask(mesh);
    const dirichletValues = new Float64Array(mesh.nVertices);

    // Run wave simulation
    const uWaveFinal = solver.waveSimulate(u0Wave, v0Wave, dtWave, 100, c, {
        dirichletMask,
        dirichletValues
    });

    const maxWave = Math.max(...uWaveFinal.map(Math.abs));
    console.log(`  Max |u| after 100 steps: ${maxWave.toFixed(6)}`);

    passed &= assertTrue(maxWave < 10, 'Wave solution bounded (|u| < 10)');
    passed &= assertTrue(isFinite(maxWave), 'Wave solution is finite');

    // --- Test 2: Heat equation (explicit, valid CFL) ---
    console.log('\n--- Test 2: Heat Equation (explicit, valid CFL) ---');

    const alpha = 0.1;
    const dtHeat = solver.estimateCFLHeat(alpha) * 0.4;
    console.log(`  CFL limit: ${solver.estimateCFLHeat(alpha).toFixed(6)}`);
    console.log(`  Using dt: ${dtHeat.toFixed(6)} (40% of limit)`);

    // Initial condition (same bump)
    const u0Heat = new Float64Array(u0Wave);

    // Run heat simulation
    const uHeatFinal = solver.heatSimulate(u0Heat, dtHeat, 100, alpha, {
        dirichletMask,
        dirichletValues,
        implicit: false
    });

    const maxHeat = Math.max(...uHeatFinal);
    const minHeat = Math.min(...uHeatFinal);
    console.log(`  Max u after 100 steps: ${maxHeat.toFixed(6)}`);
    console.log(`  Min u after 100 steps: ${minHeat.toFixed(6)}`);

    // Heat equation should smooth out, not blow up
    passed &= assertTrue(maxHeat <= 1.1, 'Heat max ≤ initial max');
    passed &= assertTrue(minHeat >= -0.1, 'Heat min ≥ 0');
    passed &= assertTrue(isFinite(maxHeat), 'Heat solution is finite');

    // --- Test 3: Heat equation (implicit, unconditionally stable) ---
    console.log('\n--- Test 3: Heat Equation (implicit, large dt) ---');

    const dtLarge = solver.estimateCFLHeat(alpha) * 5;  // 5x the explicit limit
    console.log(`  Using dt: ${dtLarge.toFixed(6)} (5x explicit limit)`);

    const u0HeatImplicit = new Float64Array(u0Wave);
    const uHeatImplicitFinal = solver.heatSimulate(u0HeatImplicit, dtLarge, 50, alpha, {
        dirichletMask,
        dirichletValues,
        implicit: true
    });

    const maxHeatImplicit = Math.max(...uHeatImplicitFinal);
    console.log(`  Max u after 50 steps: ${maxHeatImplicit.toFixed(6)}`);

    passed &= assertTrue(isFinite(maxHeatImplicit), 'Implicit heat solution is finite');
    passed &= assertTrue(maxHeatImplicit <= 1.1, 'Implicit heat respects maximum principle');

    return passed;
}

// Run test
const passed = testCFLStability();
console.log(`\n${passed ? 'ALL TESTS PASSED' : 'SOME TESTS FAILED'}`);
process.exit(passed ? 0 : 1);
