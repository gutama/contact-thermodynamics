/**
 * Laplacian Patch Test
 * 
 * Verifies that the discrete Laplacian of a constant field is zero.
 * Δ(const) = 0 at interior vertices.
 */

const { TriangleMesh } = require('../../src/mesh.js');
const MeshFTGC = require('../../src/mesh-ftgc.js');

const { MeshGeometricDerivative } = MeshFTGC;
const { abs, max } = Math;

function testLaplacianConstant() {
    console.log('=== Laplacian Patch Test: Δ(const) = 0 ===\n');

    let passed = true;

    // Create a flat grid mesh
    const nx = 10, ny = 10;
    const mesh = TriangleMesh.createGrid(nx, ny, 1, 1);
    const derivative = new MeshGeometricDerivative(mesh);

    console.log(`Grid: ${nx} × ${ny} vertices`);
    console.log(`Vertices: ${mesh.nVertices}, Edges: ${mesh.nEdges}, Faces: ${mesh.nFaces}`);

    // Test with different constant values
    const constants = [0, 1, 5, -3, 100];

    for (const c of constants) {
        // Create constant field
        const f = new Float64Array(mesh.nVertices).fill(c);

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

        console.log(`\n  Constant c = ${c}:`);
        console.log(`    Interior vertices: ${interiorCount}`);
        console.log(`    Max |Δf| at interior: ${maxError.toExponential(4)}`);

        const tolerance = 1e-10;
        if (maxError < tolerance) {
            console.log(`    PASS: Δ(${c}) ≈ 0`);
        } else {
            console.log(`    FAIL: Δ(${c}) = ${maxError.toExponential(4)} > ${tolerance}`);
            passed = false;
        }
    }

    return passed;
}

// Run test
const passed = testLaplacianConstant();
console.log(`\n${passed ? 'ALL TESTS PASSED' : 'SOME TESTS FAILED'}`);
process.exit(passed ? 0 : 1);
