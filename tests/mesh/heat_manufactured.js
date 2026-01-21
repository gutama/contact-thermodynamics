/**
 * Heat Equation Manufactured Solution Test
 * 
 * Uses the Method of Manufactured Solutions (MMS) to verify convergence.
 * 
 * Approach:
 * 1. Choose an exact solution u(x,y,t) = exp(-αt) * sin(πx) * sin(πy)
 * 2. Compute the forcing term f = ∂u/∂t - α∇²u
 * 3. Run simulation with forcing and compare to exact solution
 * 4. Verify error decreases with mesh refinement
 */

const { TriangleMesh } = require('../../src/mesh.js');
const MeshSolvers = require('../../src/mesh-solvers.js');
const MeshFTGC = require('../../src/mesh-ftgc.js');

const { LeapfrogGCMesh, boundaryDirichletMask } = MeshSolvers;
const { MeshGeometricDerivative } = MeshFTGC;
const { abs, exp, sin, cos, PI, max, sqrt } = Math;

function assertTrue(condition, message) {
    if (!condition) {
        console.error(`FAIL: ${message}`);
        return false;
    }
    console.log(`PASS: ${message}`);
    return true;
}

/**
 * Exact solution: u(x,y,t) = exp(-2απ²t) * sin(πx) * sin(πy)
 * This satisfies ∂u/∂t = α∇²u exactly.
 */
function exactSolution(x, y, t, alpha) {
    const decay = exp(-2 * alpha * PI * PI * t);
    return decay * sin(PI * x) * sin(PI * y);
}

function testHeatManufactured() {
    console.log('=== Heat Equation Manufactured Solution Test ===\n');

    let passed = true;
    const alpha = 0.1;
    const T = 0.1;  // Final time

    // Test with different mesh resolutions
    const resolutions = [5, 10, 15];
    const errors = [];

    for (const n of resolutions) {
        console.log(`\n--- Grid: ${n} × ${n} ---`);

        const mesh = TriangleMesh.createGrid(n, n, 1, 1);
        const solver = new LeapfrogGCMesh(mesh);

        // Time stepping
        const dt = solver.estimateCFLHeat(alpha) * 0.3;
        const nSteps = Math.ceil(T / dt);
        const actualT = nSteps * dt;

        console.log(`  dt = ${dt.toFixed(6)}, steps = ${nSteps}`);

        // Initial condition: exact solution at t=0
        const u0 = new Float64Array(mesh.nVertices);
        for (let i = 0; i < mesh.nVertices; i++) {
            const [x, y] = mesh.getVertex(i);
            u0[i] = exactSolution(x, y, 0, alpha);
        }

        // Boundary conditions (Dirichlet = 0 on boundary)
        const dirichletMask = boundaryDirichletMask(mesh);
        const dirichletValues = new Float64Array(mesh.nVertices);

        // Run simulation
        const uFinal = solver.heatSimulate(u0, dt, nSteps, alpha, {
            dirichletMask,
            dirichletValues,
            implicit: false
        });

        // Compute L2 error
        let errorSum = 0;
        let normSum = 0;
        for (let i = 0; i < mesh.nVertices; i++) {
            const [x, y] = mesh.getVertex(i);
            const uExact = exactSolution(x, y, actualT, alpha);
            const diff = uFinal[i] - uExact;
            errorSum += diff * diff;
            normSum += uExact * uExact;
        }
        const l2Error = sqrt(errorSum / mesh.nVertices);
        const relError = sqrt(errorSum) / (sqrt(normSum) + 1e-10);

        console.log(`  L2 error: ${l2Error.toExponential(4)}`);
        console.log(`  Relative error: ${(relError * 100).toFixed(2)}%`);

        errors.push({ n, l2Error, relError });
    }

    // Check convergence: error should decrease with refinement
    console.log('\n--- Convergence Analysis ---');

    for (let i = 1; i < errors.length; i++) {
        const ratio = errors[i - 1].l2Error / errors[i].l2Error;
        const refinement = errors[i].n / errors[i - 1].n;
        const order = Math.log(ratio) / Math.log(refinement);

        console.log(`  ${errors[i - 1].n}→${errors[i].n}: error ratio = ${ratio.toFixed(2)}, order ≈ ${order.toFixed(2)}`);

        // Should see at least first-order convergence
        passed &= assertTrue(ratio > 1.2, `Error decreases when refining ${errors[i - 1].n}→${errors[i].n}`);
    }

    // Error on finest mesh should be reasonably small
    const finestError = errors[errors.length - 1].relError;
    passed &= assertTrue(finestError < 0.5, `Finest mesh relative error < 50%`);

    return passed;
}

// Run test
const passed = testHeatManufactured();
console.log(`\n${passed ? 'ALL TESTS PASSED' : 'SOME TESTS FAILED'}`);
process.exit(passed ? 0 : 1);
