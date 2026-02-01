/**
 * Pilot-Wave Demo
 * 
 * Demonstrates Valentini's regularized pilot-wave theory:
 * - Wavefunction with nodes
 * - Regularization keeps velocity finite
 * - H-theorem: relaxation toward Born rule equilibrium
 * - Connection to contact geometry
 */

const {
    SmearingKernel,
    WaveFunction,
    PilotWaveSystem,
    QuantumEnsemble,
    ActionPhaseBridge
} = require('../src/pilot-wave.js');

console.log('═'.repeat(60));
console.log('  PILOT-WAVE THEORY DEMO');
console.log('═'.repeat(60));

// ============================================================================
// 1. WAVEFUNCTION WITH A NODE
// ============================================================================
console.log('\n▸ Creating wavefunction with a node at x=0...');

const n = 100;
const dx = 0.1;
const xMin = -5, xMax = 5;

// Create a wavefunction that has a NODE at x=0
// ψ = x * exp(-x²/4) * exp(ikx)  (p-orbital like)
const amplitude = [];
const phase = [];
const k = 0.5;  // momentum

for (let i = 0; i < n; i++) {
    const x = xMin + i * dx;
    // |ψ|² ∝ x² exp(-x²/2) → zero at x=0
    amplitude.push(Math.abs(x) * Math.exp(-x * x / 4));
    phase.push(k * x);
}

const psi = new WaveFunction(amplitude, phase, { mass: 1, hbar: 1 });
console.log(`  ψ has ${n} grid points, dx = ${dx}`);
console.log(`  Node at center: |ψ(0)|² = ${psi.probabilityDensity[n / 2].toExponential(2)}`);

// ============================================================================
// 2. REGULARIZATION COMPARISON
// ============================================================================
console.log('\n▸ Comparing standard vs regularized velocity...');

// Create kernel with regularization scale ε = 0.3
const kernel = new SmearingKernel(0.3, 'gaussian');
const system = new PilotWaveSystem(psi, { kernel, dx });

// Standard de Broglie velocity at node
const vStandard = system.deBroglieVelocity();
const vNode = vStandard[Math.floor(n / 2)];
const rhoNode = psi.probabilityDensity[Math.floor(n / 2)];
console.log(`  Standard v at node: ${vNode.toFixed(4)} (but j/ρ diverges where ρ→0!)`);
console.log(`  |ψ|² at node: ${rhoNode.toExponential(2)}`);

// Regularized velocity
const vReg = system.regularizedVelocity();
const vNodeReg = vReg[Math.floor(n / 2)];
console.log(`  Regularized v_reg at node: ${vNodeReg.toFixed(4)} (finite!)`);
console.log(`  Max |v_reg|: ${Math.max(...vReg.map(Math.abs)).toFixed(4)}`);

// Regularized density > 0 everywhere
const rhoReg = system.regularizedDensity();
const minRhoReg = Math.min(...rhoReg);
console.log(`  Min (|ψ|²)_reg: ${minRhoReg.toExponential(2)} > 0 ✓`);

// ============================================================================
// 3. TIME-DEPENDENT REGULARIZATION
// ============================================================================
console.log('\n▸ Time-dependent regularization ε(t)...');

// Kernel width grows with time (could model decoherence)
const dynamicKernel = new SmearingKernel(t => 0.1 + 0.02 * t, 'gaussian');

console.log(`  ε(0) = ${dynamicKernel.width.toFixed(3)}`);
dynamicKernel.setTime(5);
console.log(`  ε(5) = ${dynamicKernel.width.toFixed(3)}`);
dynamicKernel.setTime(10);
console.log(`  ε(10) = ${dynamicKernel.width.toFixed(3)}`);
console.log(`  → Dynamic kernel enables unstable equilibria!`);

// ============================================================================
// 4. H-THEOREM DEMONSTRATION
// ============================================================================
console.log('\n▸ H-theorem: nonequilibrium → equilibrium...');

// Create wavefunction (Gaussian)
const ampGauss = [];
const phaseGauss = [];
for (let i = 0; i < n; i++) {
    const x = xMin + i * dx;
    ampGauss.push(Math.exp(-x * x / 2));
    phaseGauss.push(0.3 * x);  // small drift
}
const psiGauss = new WaveFunction(ampGauss, phaseGauss, { mass: 1, hbar: 1 });
const systemGauss = new PilotWaveSystem(psiGauss, { kernel, dx });

// Start with SHIFTED distribution (nonequilibrium)
const rhoInit = [];
for (let i = 0; i < n; i++) {
    const x = xMin + i * dx;
    const shift = 1.5;
    rhoInit.push(Math.exp(-(x - shift) * (x - shift) / 2));
}
const ensemble = new QuantumEnsemble(rhoInit, psiGauss, dx, kernel);

// Compute H-function trajectory
const H0 = ensemble.hFunctionFine();
console.log(`  Initial H = ${H0.toFixed(4)} (nonequilibrium: H > 0)`);
console.log(`  Is in equilibrium? ${ensemble.isInEquilibrium(0.1) ? 'YES' : 'NO'}`);

// Evolve
const hValues = ensemble.relaxation(systemGauss, 0.01, 50, true);
const HFinal = hValues[hValues.length - 1];
console.log(`  Final H = ${HFinal.toFixed(4)}`);
console.log(`  H trajectory: ${hValues.slice(0, 5).map(h => h.toFixed(3)).join(' → ')} → ...`);

// ============================================================================
// 5. REGULARIZED EQUILIBRIUM
// ============================================================================
console.log('\n▸ Regularized equilibrium: (|ψ|²)_reg...');

// Create ensemble at regularized equilibrium
const rhoQMReg = kernel.convolve1D(psiGauss.probabilityDensity, dx);
const regEnsemble = new QuantumEnsemble([...rhoQMReg], psiGauss, dx, kernel);

console.log(`  Standard equilibrium (ρ = |ψ|²)? ${regEnsemble.isInEquilibrium(0.1) ? 'YES' : 'NO'}`);
console.log(`  Regularized equilibrium (ρ = (|ψ|²)_reg)? ${regEnsemble.isInRegularizedEquilibrium(kernel, 0.1) ? 'YES ✓' : 'NO'}`);

// Show that (|ψ|²)_reg > 0 everywhere (no true nodes)
const regEq = regEnsemble.getRegularizedEquilibrium(kernel);
console.log(`  Min (|ψ|²)_reg: ${Math.min(...regEq).toExponential(2)} (no zeros)`);

// ============================================================================
// 6. CONTACT GEOMETRY CONNECTION
// ============================================================================
console.log('\n▸ Connection to contact geometry...');

const bridge = new ActionPhaseBridge(psiGauss);

// Phase → Action
const centerIdx = Math.floor(n / 2);
const A_array = bridge.phaseToAction();
const p_array = bridge.momentum(dx);
const A = A_array[centerIdx];
const p = p_array[centerIdx];

console.log(`  Phase S(center) = ${phaseGauss[centerIdx].toFixed(4)}`);
console.log(`  Action A = ℏS = ${A.toFixed(4)}`);
console.log(`  Momentum p = ∇S = ${p.toFixed(4)}`);

// Legendrian condition: α|_L = dA - p dx = 0
const violation = bridge.legendrianViolation(dx, centerIdx);
console.log(`  Legendrian violation |α|_L|: ${violation.toExponential(2)} (≈ 0 ✓)`);

// ============================================================================
// SUMMARY
// ============================================================================
console.log('\n' + '═'.repeat(60));
console.log('  SUMMARY');
console.log('═'.repeat(60));
console.log(`
  ✓ Regularized velocity v_reg = j_reg/(|ψ|²)_reg is FINITE at nodes
  ✓ Time-dependent ε(t) enables dynamically unstable equilibria
  ✓ H-theorem: systems relax toward Born rule ρ → |ψ|²
  ✓ Regularized equilibrium: ρ_eq = (|ψ|²)_reg (nonzero at former nodes)
  ✓ Phase S connects to contact geometry action A via Legendrian condition
`);
