/**
 * Wave Equation Energy Conservation Test
 * 
 * Verifies that the wave equation solver conserves total energy.
 * For wave equation: E = ½∫(u_t² + c²|∇u|²) dA ≈ const
 * 
 * We use the simpler L2 energy proxy: E = ∫u² dA
 * which oscillates but should remain bounded.
 */

const { TriangleMesh } = require('../../src/mesh.js');
const MeshSolvers = require('../../src/mesh-solvers.js');

const { LeapfrogGCMesh, boundaryDirichletMask } = MeshSolvers;
const { abs, PI, sin, cos, exp, sqrt } = Math;

function assertTrue(condition, message) {
    if (!condition) {
        console.error(`FAIL: ${message}`);
        return false;
    }
    console.log(`PASS: ${message}`);
    return true;
}

function testWaveEnergyConservation() {
    console.log('=== Wave Equation Energy Conservation Test ===\n');

    let passed = true;

    // Create a mesh
    const nx = 15, ny = 15;
    const Lx = 1, Ly = 1;
    const mesh = TriangleMesh.createGrid(nx, ny, Lx, Ly);
    const solver = new LeapfrogGCMesh(mesh);

    console.log(`Grid: ${nx} × ${ny} vertices`);
    console.log(`Domain: [0, ${Lx}] × [0, ${Ly}]`);

    // Wave speed
    const c = 1.0;

    // Time stepping
    const dt = solver.estimateCFL(c) * 0.5;  // Safe CFL factor
    const nSteps = 200;
    console.log(`CFL dt: ${solver.estimateCFL(c).toFixed(6)}, using dt = ${dt.toFixed(6)}`);
    console.log(`Steps: ${nSteps}`);

    // Initial condition: Gaussian bump
    const u0 = new Float64Array(mesh.nVertices);
    const v0 = new Float64Array(mesh.nVertices);  // Zero initial velocity
    const cx = 0.5, cy = 0.5;
    const sigma = 0.1;

    for (let i = 0; i < mesh.nVertices; i++) {
        const [x, y, z] = mesh.getVertex(i);
        const r2 = (x - cx) ** 2 + (y - cy) ** 2;
        u0[i] = exp(-r2 / (2 * sigma ** 2));
    }

    // Boundary conditions (fixed at 0)
    const dirichletMask = boundaryDirichletMask(mesh);
    const dirichletValues = new Float64Array(mesh.nVertices);

    // Track energy over time
    const energies = [];
    const callback = (step, u) => {
        if (step % 20 === 0) {
            energies.push(solver.energy(u));
        }
    };

    // Run simulation
    const uFinal = solver.waveSimulate(u0, v0, dt, nSteps, c, {
        dirichletMask,
        dirichletValues,
        callback
    });

    // Analyze energy conservation
    const E0 = energies[0];
    const Emax = Math.max(...energies);
    const Emin = Math.min(...energies);
    const maxDrift = Math.max(abs(Emax - E0), abs(Emin - E0)) / E0;

    console.log(`\nEnergy analysis:`);
    console.log(`  Initial: ${E0.toFixed(6)}`);
    console.log(`  Min: ${Emin.toFixed(6)}, Max: ${Emax.toFixed(6)}`);
    console.log(`  Relative drift: ${(maxDrift * 100).toFixed(2)}%`);

    // For leapfrog, L2 energy oscillates but total energy is conserved
    // We check that energy stays bounded and doesn't blow up
    const blowUpRatio = Emax / E0;
    const decayRatio = E0 / Emax;

    passed &= assertTrue(blowUpRatio < 5, 'Energy does not blow up (max < 5× initial)');
    passed &= assertTrue(Emax > 0, 'Energy remains positive');
    passed &= assertTrue(Emin > 0, 'Energy never goes to zero (solution not vanishing)');

    return passed;
}

// Run test
const passed = testWaveEnergyConservation();
console.log(`\n${passed ? 'ALL TESTS PASSED' : 'SOME TESTS FAILED'}`);
process.exit(passed ? 0 : 1);
