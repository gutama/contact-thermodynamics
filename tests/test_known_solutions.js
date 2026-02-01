/**
 * Known Quantum Solutions - Pilot-Wave Comparison
 * 
 * Tests pilot-wave regularization against analytically known wavefunctions:
 * 1. Hydrogen atom (2p state with nodal plane at z=0)
 * 2. Quantum harmonic oscillator (n=1, n=2 with Hermite polynomial nodes)
 * 
 * Verifies:
 * - Regularized velocity stays finite at nodes
 * - Velocity matches expected de Broglie prediction away from nodes
 * - Born rule equilibrium is achieved
 */

const PilotWave = require('../src/pilot-wave.js');
const { SmearingKernel, WaveFunction, PilotWaveSystem, QuantumEnsemble } = PilotWave;

console.log('========================================');
console.log('  KNOWN SOLUTIONS COMPARISON');
console.log('========================================');

let passed = 0;
let failed = 0;

function assert(condition, name) {
    if (condition) {
        console.log(`  ✓ ${name}`);
        passed++;
    } else {
        console.log(`  ✗ ${name}`);
        failed++;
    }
}

function assertApprox(got, expected, tol, name) {
    const ok = Math.abs(got - expected) < tol;
    assert(ok, `${name} (got ${got.toFixed(6)}, expected ${expected.toFixed(6)})`);
}

// ============================================================================
// 1. HYDROGEN ATOM - 2p_z ORBITAL
// ============================================================================
//
// ψ_2pz ∝ r·exp(-r/2a₀)·cos(θ)
//
// Has nodal plane at z=0 (where cos(θ) = 0, i.e., θ = π/2)
// The de Broglie velocity would diverge at this plane in standard theory.

console.log('\n--- Hydrogen 2p_z Orbital ---\n');

(function testHydrogen2p() {
    const a0 = 1.0;  // Bohr radius (dimensionless units)
    const n = 101;   // Grid points along z-axis
    const dz = 0.1;  // Grid spacing

    // Sample along z-axis (θ = 0 means z > 0, θ = π means z < 0)
    // At θ = π/2 (z = 0), wavefunction vanishes
    const amplitude = [];
    const phase = [];

    for (let i = 0; i < n; i++) {
        const z = (i - 50) * dz;  // z from -5 to +5
        const r = Math.abs(z);

        // ψ_2pz ∝ z·exp(-|z|/2a₀) along z-axis
        // (since along z-axis, cos(θ) = sign(z) and r = |z|)
        const R = Math.abs(z) * Math.exp(-r / (2 * a0));
        amplitude.push(R);

        // Phase is 0 for z > 0, π for z < 0 (sign change)
        phase.push(z >= 0 ? 0 : Math.PI);
    }

    const psi = new WaveFunction(amplitude, phase, { mass: 1, hbar: 1 });

    // Check: node at z=0
    const centerIdx = 50;
    assert(amplitude[centerIdx] < 0.01, 'Node present at z=0 (ψ_2pz)');
    assert(amplitude[30] > 0.1, 'Non-zero amplitude away from node');

    // Regularized system
    const kernel = new SmearingKernel(0.3, 'gaussian');
    const system = new PilotWaveSystem(psi, { kernel, dx: dz });

    // Standard velocity would have issues at node
    const vReg = system.regularizedVelocity();

    // Check velocity is finite everywhere
    let allFinite = true;
    for (const v of vReg) {
        if (!isFinite(v)) allFinite = false;
    }
    assert(allFinite, 'Regularized velocity finite everywhere (2p_z)');

    // Check velocity at node is bounded
    const vAtNode = vReg[centerIdx];
    assert(Math.abs(vAtNode) < 10, `Velocity bounded at node (v=${vAtNode.toFixed(4)})`);

    // Check regularized density > 0 at node
    const rhoReg = system.regularizedDensity();
    assert(rhoReg[centerIdx] > 0, `ρ_reg > 0 at node (=${rhoReg[centerIdx].toExponential(2)})`);

    console.log('\n  2p_z profile summary:');
    console.log(`    Max |ψ|²: ${Math.max(...psi.probabilityDensity).toFixed(4)}`);
    console.log(`    |ψ|² at node: ${psi.probabilityDensity[centerIdx].toExponential(2)}`);
    console.log(`    ρ_reg at node: ${rhoReg[centerIdx].toExponential(2)}`);
})();

// ============================================================================
// 2. HARMONIC OSCILLATOR - FIRST EXCITED STATE (n=1)
// ============================================================================
//
// ψ₁(x) ∝ x·exp(-x²/2)  (Hermite H₁(x) = 2x)
//
// Has node at x=0

console.log('\n--- Harmonic Oscillator n=1 ---\n');

(function testOscillator_n1() {
    const n = 101;
    const dx = 0.1;

    const amplitude = [];
    const phase = [];

    for (let i = 0; i < n; i++) {
        const x = (i - 50) * dx;  // x from -5 to +5

        // ψ₁ ∝ x·exp(-x²/2)
        amplitude.push(Math.abs(x) * Math.exp(-x * x / 2));
        phase.push(x >= 0 ? 0 : Math.PI);  // Sign change at origin
    }

    const psi = new WaveFunction(amplitude, phase, { mass: 1, hbar: 1 });
    const kernel = new SmearingKernel(0.25, 'gaussian');
    const system = new PilotWaveSystem(psi, { kernel, dx });

    // Check node at x=0
    const centerIdx = 50;
    assert(amplitude[centerIdx] < 0.01, 'Node present at x=0 (n=1 oscillator)');

    // Regularized quantities
    const vReg = system.regularizedVelocity();
    const rhoReg = system.regularizedDensity();

    assert(isFinite(vReg[centerIdx]), 'v_reg finite at node (n=1)');
    assert(rhoReg[centerIdx] > 0, `ρ_reg > 0 at node (=${rhoReg[centerIdx].toExponential(2)})`);

    // Check that far from node, velocity is bounded
    const maxV = Math.max(...vReg.map(Math.abs));
    assert(maxV < 20, `Max velocity bounded (=${maxV.toFixed(2)})`);

    // Equilibrium check: for noded wavefunctions, ρ = |ψ|² differs slightly 
    // from (|ψ|²)_reg due to convolution, so H is not exactly 0.
    // H should still be small (< 0.2) indicating near-equilibrium
    const eqEnsemble = QuantumEnsemble.equilibrium(psi, dx);
    const H = eqEnsemble.hFunction();
    assert(H < 0.2, `H-function small at near-equilibrium (H=${H.toFixed(4)})`);
})();

// ============================================================================
// 3. HARMONIC OSCILLATOR - SECOND EXCITED STATE (n=2)
// ============================================================================
//
// ψ₂(x) ∝ (4x² - 2)·exp(-x²/2)  (Hermite H₂(x) = 4x² - 2)
//
// Has nodes at x = ±1/√2 ≈ ±0.707

console.log('\n--- Harmonic Oscillator n=2 ---\n');

(function testOscillator_n2() {
    const n = 121;
    const dx = 0.1;

    const amplitude = [];
    const phase = [];

    // Node positions: 4x² - 2 = 0 → x² = 0.5 → x = ±0.707
    const nodePos = Math.sqrt(0.5);

    for (let i = 0; i < n; i++) {
        const x = (i - 60) * dx;  // x from -6 to +6

        // H₂(x) = 4x² - 2
        const H2 = 4 * x * x - 2;
        const R = Math.abs(H2) * Math.exp(-x * x / 2);
        amplitude.push(R);

        // Phase: 0 where H₂ > 0 (|x| > 0.707), π where H₂ < 0 (|x| < 0.707)
        phase.push(H2 >= 0 ? 0 : Math.PI);
    }

    const psi = new WaveFunction(amplitude, phase, { mass: 1, hbar: 1 });
    const kernel = new SmearingKernel(0.2, 'gaussian');
    const system = new PilotWaveSystem(psi, { kernel, dx });

    // Find indices near nodes (x ≈ ±0.707)
    const nodeIdx1 = Math.round(60 + nodePos / dx);  // ~67
    const nodeIdx2 = Math.round(60 - nodePos / dx);  // ~53

    // Check nodes are present
    assert(amplitude[nodeIdx1] < 0.1, `Node near x=+${nodePos.toFixed(3)} (n=2)`);
    assert(amplitude[nodeIdx2] < 0.1, `Node near x=-${nodePos.toFixed(3)} (n=2)`);

    // Regularized quantities
    const vReg = system.regularizedVelocity();
    const rhoReg = system.regularizedDensity();

    // Check velocity finite at both nodes
    assert(isFinite(vReg[nodeIdx1]) && isFinite(vReg[nodeIdx2]), 'v_reg finite at both nodes (n=2)');

    // Check density > 0 at nodes
    assert(rhoReg[nodeIdx1] > 0 && rhoReg[nodeIdx2] > 0, 'ρ_reg > 0 at both nodes (n=2)');

    console.log('\n  n=2 oscillator profile:');
    console.log(`    Node positions: x = ±${nodePos.toFixed(4)}`);
    console.log(`    |ψ|² near node1: ${psi.probabilityDensity[nodeIdx1].toExponential(2)}`);
    console.log(`    ρ_reg near node1: ${rhoReg[nodeIdx1].toExponential(2)}`);
})();

// ============================================================================
// 4. VELOCITY FIELD COMPARISON (Plane Wave Limit)
// ============================================================================
//
// For a plane wave ψ = exp(ikx), the de Broglie velocity should be exactly v = ℏk/m
// This is the flat-space, node-free case where regularization shouldn't change much.

console.log('\n--- Plane Wave Velocity Check ---\n');

(function testPlaneWave() {
    const n = 100;
    const dx = 0.1;
    const k = 3.0;  // Wavenumber
    const hbar = 1.0;
    const mass = 1.0;

    const amplitude = [];
    const phase = [];

    for (let i = 0; i < n; i++) {
        const x = i * dx;
        amplitude.push(1.0);  // Unit amplitude (no nodes)
        phase.push(k * x);    // Linear phase
    }

    const psi = new WaveFunction(amplitude, phase, { mass, hbar });
    const kernel = new SmearingKernel(0.3, 'gaussian');
    const system = new PilotWaveSystem(psi, { kernel, dx });

    // Expected de Broglie velocity: v = ℏk/m
    const expectedV = hbar * k / mass;

    // Standard velocity
    const vStd = system.deBroglieVelocity();
    // Regularized velocity (should be similar for uniform amplitude)
    const vReg = system.regularizedVelocity();

    // Check velocity at middle of grid
    const midIdx = 50;
    assertApprox(vStd[midIdx], expectedV, 0.01, 'Standard v = ℏk/m (plane wave)');
    assertApprox(vReg[midIdx], expectedV, 0.1, 'Regularized v ≈ ℏk/m (plane wave)');
})();

// ============================================================================
// 5. GROUND STATE GAUSSIAN (No Nodes)
// ============================================================================
//
// ψ₀(x) ∝ exp(-x²/2) — the harmonic oscillator ground state has NO nodes.
// This tests that regularization doesn't distort a well-behaved wavefunction.

console.log('\n--- Ground State (No Nodes) ---\n');

(function testGroundState() {
    const n = 101;
    const dx = 0.1;

    const amplitude = [];
    const phase = [];

    for (let i = 0; i < n; i++) {
        const x = (i - 50) * dx;
        amplitude.push(Math.exp(-x * x / 2));  // Gaussian (no nodes)
        phase.push(0);  // Stationary
    }

    const psi = new WaveFunction(amplitude, phase, { mass: 1, hbar: 1 });
    const kernel = new SmearingKernel(0.3, 'gaussian');
    const system = new PilotWaveSystem(psi, { kernel, dx });

    // No nodes → standard velocity should be well-defined
    const vStd = system.deBroglieVelocity();
    const vReg = system.regularizedVelocity();

    // For stationary state with phase = 0, velocity should be ~0
    const midIdx = 50;
    assertApprox(vStd[midIdx], 0, 0.01, 'Standard v ≈ 0 for ground state');
    assertApprox(vReg[midIdx], 0, 0.01, 'Regularized v ≈ 0 for ground state');

    // Density should be unchanged by regularization (within smoothing)
    const rho = psi.probabilityDensity[midIdx];
    const rhoReg = system.regularizedDensity()[midIdx];

    // rho_reg should be close to rho for smooth distributions
    const relDiff = Math.abs(rhoReg - rho) / rho;
    assert(relDiff < 0.3, `ρ_reg ≈ ρ for smooth ground state (diff=${(relDiff * 100).toFixed(1)}%)`);
})();

// ============================================================================
// SUMMARY
// ============================================================================

console.log('\n========================================');
console.log(`  RESULTS: ${passed} passed, ${failed} failed`);
console.log('========================================');

if (failed > 0) {
    process.exit(1);
}
