/**
 * Laplacian Linear Test
 * 
 * Verifies that the discrete Laplacian of a linear field is zero.
 * Δ(ax + by + c) = 0 at interior vertices.
 */

const { TriangleMesh } = require('../../src/mesh.js');
const MeshFTGC = require('../../src/mesh-ftgc.js');

const { MeshGeometricDerivative } = MeshFTGC;
const { abs, max } = Math;

function testLaplacianLinear() {
    console.log('=== Laplacian Linear Test: Δ(ax + by + c) = 0 ===\n');

    let passed = true;

    // Create a flat grid mesh
    const nx = 10, ny = 10;
    const Lx = 1, Ly = 1;
    const mesh = TriangleMesh.createGrid(nx, ny, Lx, Ly);
    const derivative = new MeshGeometricDerivative(mesh);

    console.log(`Grid: ${nx} × ${ny} vertices`);
    console.log(`Domain: [0, ${Lx}] × [0, ${Ly}]`);

    // Test with different linear functions
    const testCases = [
        { a: 1, b: 0, c: 0, desc: 'f = x' },
        { a: 0, b: 1, c: 0, desc: 'f = y' },
        { a: 1, b: 1, c: 0, desc: 'f = x + y' },
        { a: 2, b: -1, c: 3, desc: 'f = 2x - y + 3' },
        { a: -1, b: 2, c: -5, desc: 'f = -x + 2y - 5' }
    ];

    for (const { a, b, c, desc } of testCases) {
        // Create linear field: f(x, y) = ax + by + c
        const f = new Float64Array(mesh.nVertices);
        for (let i = 0; i < mesh.nVertices; i++) {
            const [x, y, z] = mesh.getVertex(i);
            f[i] = a * x + b * y + c;
        }

        // Compute Laplacian
        const lapF = derivative.laplacian(f);

        // Find max absolute error at interior vertices
        let maxError = 0;
        let interiorCount = 0;
        for (let i = 0; i < mesh.nVertices; i++) {
            if (!mesh.boundaryVertices[i]) {
                maxError = max(maxError, abs(lapF[i]));
                interiorCount++;
            }
        }

        console.log(`\n  ${desc}:`);
        console.log(`    Interior vertices: ${interiorCount}`);
        console.log(`    Max |Δf| at interior: ${maxError.toExponential(4)}`);

        const tolerance = 1e-8;
        if (maxError < tolerance) {
            console.log(`    PASS: Δ(${desc}) ≈ 0`);
        } else {
            console.log(`    FAIL: Δ(${desc}) = ${maxError.toExponential(4)} > ${tolerance}`);
            passed = false;
        }
    }

    return passed;
}

// Run test
const passed = testLaplacianLinear();
console.log(`\n${passed ? 'ALL TESTS PASSED' : 'SOME TESTS FAILED'}`);
process.exit(passed ? 0 : 1);
