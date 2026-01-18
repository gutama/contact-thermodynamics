/**
 * FTGC Mesh Heat Diffusion Example
 * 
 * Demonstrates entropy diffusion on a triangulated grid using
 * discrete geometric calculus (FTGC) with Dirichlet boundary conditions.
 * 
 * Physical interpretation:
 *   - Scalar field u = entropy density S
 *   - Heat equation ∂S/∂t = α∇²S models entropy diffusion
 *   - Variance decay shows approach to thermal equilibrium
 * 
 * @example
 *   node examples/mesh-heat-ftgc.js
 */

const {
    TriangleMesh,
    MeshGeometricDerivative,
    LeapfrogGCMesh,
    boundaryDirichletMask
} = require('../src/index.js');

// ============================================================================
// SETUP: Create triangulated grid mesh
// ============================================================================

const nx = 16, ny = 16;
const Lx = 1.0, Ly = 1.0;

console.log('='.repeat(60));
console.log('FTGC Mesh Heat Diffusion (Entropy Relaxation)');
console.log('='.repeat(60));
console.log(`Grid: ${nx}×${ny} vertices → ${2*(nx-1)*(ny-1)} triangles`);

const mesh = TriangleMesh.createGrid(nx, ny, Lx, Ly);

console.log(`Vertices: ${mesh.nVertices}`);
console.log(`Edges: ${mesh.nEdges}`);
console.log(`Faces: ${mesh.nFaces}`);
console.log(`Boundary vertices: ${mesh.boundaryVertices.reduce((a,b) => a+b, 0)}`);

// ============================================================================
// VERIFY FTGC IDENTITIES
// ============================================================================

console.log('\n--- FTGC Verification ---');

const nabla = new MeshGeometricDerivative(mesh);

// Test 1: ∇∧(∇∧f) = 0 (curl of gradient is zero)
const testField = new Float64Array(mesh.nVertices);
for (let i = 0; i < mesh.nVertices; i++) {
    const x = mesh.vertices[i * 3];
    const y = mesh.vertices[i * 3 + 1];
    testField[i] = Math.sin(Math.PI * x) * Math.cos(Math.PI * y);
}
const curlGradErr = nabla.verifyCurlGradZero(testField);
console.log(`∇∧(∇∧f) = 0: max error = ${curlGradErr.toExponential(3)}`);

// Test 2: ∇²(constant) = 0
const lapConstErr = nabla.verifyLaplacianConstantZero(1.0);
console.log(`∇²(const) = 0: max error = ${lapConstErr.toExponential(3)}`);

// ============================================================================
// INITIAL CONDITION: Gaussian bump (entropy concentration)
// ============================================================================

const u0 = new Float64Array(mesh.nVertices);
const cx = 0.5, cy = 0.5, sigma = 0.1;

for (let i = 0; i < mesh.nVertices; i++) {
    const x = mesh.vertices[i * 3];
    const y = mesh.vertices[i * 3 + 1];
    const r2 = (x - cx) ** 2 + (y - cy) ** 2;
    u0[i] = Math.exp(-r2 / (2 * sigma * sigma));
}

// ============================================================================
// DIRICHLET BOUNDARY CONDITIONS: u = 0 on boundary
// ============================================================================

const dirichletMask = boundaryDirichletMask(mesh);
const dirichletValues = new Float64Array(mesh.nVertices); // zeros

// ============================================================================
// TIME STEPPING: Heat equation via explicit Euler
// ============================================================================

console.log('\n--- Heat Equation Simulation ---');

const solver = new LeapfrogGCMesh(mesh);
const alpha = 0.1;  // Diffusion coefficient
const dt = solver.estimateCFLHeat(alpha);
const nSteps = 100;

console.log(`Diffusion coefficient α = ${alpha}`);
console.log(`Time step dt = ${dt.toExponential(3)} (CFL estimate)`);
console.log(`Total steps = ${nSteps}`);

// Track variance decay (measure of equilibration)
const initialVariance = solver.variance(u0);
const initialMass = solver.mass(u0);

console.log(`\nInitial mass = ${initialMass.toFixed(6)}`);
console.log(`Initial variance = ${initialVariance.toExponential(4)}`);

const u = solver.heatSimulate(u0, dt, nSteps, alpha, {
    dirichletMask,
    dirichletValues,
    callback: (step, u) => {
        if (step % 25 === 0 || step === nSteps - 1) {
            const var_t = solver.variance(u);
            const mass_t = solver.mass(u);
            const relDecay = (1 - var_t / initialVariance) * 100;
            console.log(`  Step ${step.toString().padStart(3)}: variance = ${var_t.toExponential(4)}, decay = ${relDecay.toFixed(1)}%`);
        }
    }
});

const finalVariance = solver.variance(u);
const finalMass = solver.mass(u);

console.log(`\nFinal mass = ${finalMass.toFixed(6)} (conservation: ${(finalMass/initialMass*100).toFixed(2)}%)`);
console.log(`Final variance = ${finalVariance.toExponential(4)}`);
console.log(`Variance decay = ${((1 - finalVariance/initialVariance)*100).toFixed(1)}%`);

// ============================================================================
// VERIFY COTAN LAPLACIAN STRUCTURE
// ============================================================================

console.log('\n--- Cotan Laplacian Structure ---');

const L = nabla.laplacianMatrix();
const diagL = L.diagonal();

// Check row sums (should be ~0 for closed meshes / interior vertices)
let maxRowSum = 0;
for (let i = 0; i < mesh.nVertices; i++) {
    if (mesh.boundaryVertices[i]) continue;
    
    // Sum row i
    let rowSum = 0;
    for (let j = 0; j < mesh.nVertices; j++) {
        rowSum += L.get(i, j);
    }
    maxRowSum = Math.max(maxRowSum, Math.abs(rowSum));
}
console.log(`Interior row sums ≈ 0: max = ${maxRowSum.toExponential(3)}`);

// Diagonal should be negative (negative semi-definite)
let diagNeg = 0, diagPos = 0;
for (let i = 0; i < diagL.length; i++) {
    if (diagL[i] < 0) diagNeg++;
    else if (diagL[i] > 0) diagPos++;
}
console.log(`Diagonal signs: ${diagNeg} negative, ${diagPos} positive (expect all negative)`);

console.log('\n' + '='.repeat(60));
console.log('FTGC mesh heat diffusion example completed successfully!');
console.log('='.repeat(60));
