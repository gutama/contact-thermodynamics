/**
 * Contact Hamiltonian Dynamics Examples
 * 
 * Demonstrates:
 * - Creating Hamiltonians
 * - Computing vector fields
 * - Integrating trajectories
 * - Energy evolution
 */

const CT = require('../src/index.js');

console.log('═══════════════════════════════════════════════════════════════');
console.log('  Contact Hamiltonian Dynamics Examples');
console.log('═══════════════════════════════════════════════════════════════\n');

const M13 = CT.grandManifold();

// ============================================================================
// 1. Simple Hamiltonian: H = ω (frequency as energy)
// ============================================================================

console.log('1. SIMPLE HAMILTONIAN: H = ω\n');

const H_simple = coords => coords.omega;
const ham_simple = new CT.ContactHamiltonian(M13, H_simple);

const pt1 = M13.physicalPoint(
    0, 0, 0, 0, 0, 1,
    0.5, 0, 0, 1, 0, 1,
    0
);

console.log('Initial point:');
console.log('   ω =', pt1.get('omega'));
console.log('   t =', pt1.get('t'));
console.log();

console.log('Hamiltonian value: H =', ham_simple.evaluate(pt1));
console.log();

// Vector field
const X1 = ham_simple.vectorField(pt1);
console.log('Contact vector field X_H:');
console.log('   ṫ = ∂H/∂ω =', X1.t.toFixed(4), '(expect 1)');
console.log('   ω̇ = -∂H/∂t - ω·RH =', X1.omega.toFixed(4), '(expect 0)');
console.log('   Ȧ = p·∂H/∂p - H =', X1.A.toFixed(4));
console.log();

// ============================================================================
// 2. Dispersion Relation: H = ω - c|k|
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('2. DISPERSION RELATION: H = ω - c|k| (massless)\n');

const ham_disp = CT.ThermodynamicHamiltonian.dispersionRelation(M13, 1, 0);

// Point on the mass shell (H = 0)
const pt_onshell = M13.physicalPoint(
    0, 0, 0, 0, 0, 0,
    1, 0, 0,  // |k| = 1
    1,        // ω = 1 (so H = 1 - 1 = 0)
    0, 0,
    0
);

console.log('On-shell point: |k| = 1, ω = 1');
console.log('   H = ω - |k| =', ham_disp.evaluate(pt_onshell).toFixed(6), '(expect 0)');
console.log();

// Off-shell point
const pt_offshell = M13.physicalPoint(
    0, 0, 0, 0, 0, 0,
    1, 0, 0,
    2,        // ω = 2 ≠ |k|
    0, 0,
    0
);

console.log('Off-shell point: |k| = 1, ω = 2');
console.log('   H = ω - |k| =', ham_disp.evaluate(pt_offshell).toFixed(6), '(expect 1)');
console.log();

// ============================================================================
// 3. Massive Dispersion: H = ω - √(c²|k|² + m²c⁴)
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('3. MASSIVE DISPERSION: H = ω - √(|k|² + m²)\n');

const m = 1;
const ham_massive = CT.ThermodynamicHamiltonian.dispersionRelation(M13, 1, m);

// At rest (k = 0): ω = m
const pt_rest = M13.physicalPoint(
    0, 0, 0, 0, 0, 0,
    0, 0, 0,  // k = 0
    m,        // ω = m
    0, 0,
    0
);

console.log('Particle at rest: k = 0');
console.log('   ω = m =', m);
console.log('   H =', ham_massive.evaluate(pt_rest).toFixed(6), '(expect 0)');
console.log();

// Moving (k ≠ 0): ω = √(|k|² + m²)
const k_mag = 1;
const omega_moving = Math.sqrt(k_mag * k_mag + m * m);
const pt_moving = M13.physicalPoint(
    0, 0, 0, 0, 0, 0,
    k_mag, 0, 0,
    omega_moving,
    0, 0,
    0
);

console.log('Moving particle: |k| =', k_mag);
console.log('   ω = √(|k|² + m²) =', omega_moving.toFixed(4));
console.log('   H =', ham_massive.evaluate(pt_moving).toFixed(6), '(expect 0)');
console.log();

// ============================================================================
// 4. Free Particle: H = ½|k|²
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('4. FREE PARTICLE: H = ½|k|²\n');

const ham_free = new CT.ContactHamiltonian(M13, coords => {
    const k1 = coords.k1 || 0;
    const k2 = coords.k2 || 0;
    const k3 = coords.k3 || 0;
    return 0.5 * (k1*k1 + k2*k2 + k3*k3);
});

const pt_free = M13.physicalPoint(
    0, 0, 0, 0, 0, 0,
    1, 0, 0,  // k = (1, 0, 0)
    0, 0, 0,
    0
);

console.log('Initial: q = (0, 0, 0), k = (1, 0, 0)');
console.log('   H = ½|k|² =', ham_free.evaluate(pt_free).toFixed(4));
console.log();

// Integrate
const dt = 0.1;
const steps = 10;
const trajectory = ham_free.flow(pt_free, dt, steps);

console.log(`Flow for ${steps} steps with dt = ${dt}:`);
console.log();
console.log('   Step    q¹          k₁          A           H');
console.log('   ─────────────────────────────────────────────────────');

for (let i = 0; i <= steps; i += 2) {
    const p = trajectory[i];
    const H = ham_free.evaluate(p);
    console.log(`   ${i.toString().padStart(3)}     ${p.get('q1').toFixed(4).padStart(8)}    ${p.get('k1').toFixed(4).padStart(8)}    ${p.get('A').toFixed(4).padStart(8)}    ${H.toFixed(4).padStart(8)}`);
}
console.log();

// Verify momentum conservation
const k1_init = trajectory[0].get('k1');
const k1_final = trajectory[steps].get('k1');
console.log('Momentum conservation:');
console.log('   k₁(0) =', k1_init.toFixed(6));
console.log('   k₁(T) =', k1_final.toFixed(6));
console.log('   Δk₁ =', Math.abs(k1_final - k1_init).toExponential(2));
console.log();

// ============================================================================
// 5. Thermodynamic Hamiltonian
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('5. THERMODYNAMIC HAMILTONIAN: H = ½T - TS\n');

const ham_thermo = CT.ThermodynamicHamiltonian.equationOfState(M13, 'ideal');

const pt_thermo = M13.physicalPoint(
    0, 0, 0, 0, 0, 2,  // S = 2
    0, 0, 0,
    0, 0, 1,           // T = 1
    0
);

console.log('Initial: S = 2, T = 1');
console.log('   H =', ham_thermo.evaluate(pt_thermo).toFixed(4));
console.log();

// Flow
const thermo_traj = ham_thermo.flow(pt_thermo, 0.05, 20);
const H_values = ham_thermo.hamiltonianEvolution(thermo_traj);

console.log('Thermodynamic flow (20 steps):');
console.log('   H(0) =', H_values[0].toFixed(4));
console.log('   H(T) =', H_values[20].toFixed(4));
console.log('   ΔH =', (H_values[20] - H_values[0]).toFixed(4));
console.log();

// ============================================================================
// 6. Custom Hamiltonian with Gradient
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('6. CUSTOM HAMILTONIAN WITH ANALYTICAL GRADIENT\n');

// Harmonic oscillator: H = ½(k² + q²)
const ham_harmonic = new CT.ContactHamiltonian(
    M13,
    // Hamiltonian
    coords => 0.5 * ((coords.k1 || 0)**2 + (coords.q1 || 0)**2),
    // Analytical gradient
    coords => ({
        q1: coords.q1 || 0,
        q2: 0, q3: 0,
        t: 0, ell: 0, S: 0,
        k1: coords.k1 || 0,
        k2: 0, k3: 0,
        omega: 0, Delta: 0, T: 0,
        A: 0
    })
);

const pt_osc = M13.physicalPoint(
    1, 0, 0, 0, 0, 0,  // q = (1, 0, 0)
    0, 0, 0,           // k = (0, 0, 0)
    0, 0, 0,
    0
);

console.log('Harmonic oscillator: H = ½(k₁² + q₁²)');
console.log('Initial: q₁ = 1, k₁ = 0');
console.log('   H =', ham_harmonic.evaluate(pt_osc).toFixed(4));
console.log();

// Integrate one period (approximately)
const osc_traj = ham_harmonic.flow(pt_osc, 0.05, 126);

console.log('Oscillation (one period ≈ 2π):');
console.log('   t = 0:   q₁ =', osc_traj[0].get('q1').toFixed(4), ', k₁ =', osc_traj[0].get('k1').toFixed(4));
console.log('   t = π/2: q₁ =', osc_traj[31].get('q1').toFixed(4), ', k₁ =', osc_traj[31].get('k1').toFixed(4));
console.log('   t = π:   q₁ =', osc_traj[63].get('q1').toFixed(4), ', k₁ =', osc_traj[63].get('k1').toFixed(4));
console.log('   t = 2π:  q₁ =', osc_traj[126].get('q1').toFixed(4), ', k₁ =', osc_traj[126].get('k1').toFixed(4));
console.log();

console.log('═══════════════════════════════════════════════════════════════');
console.log('  Dynamics examples complete!');
console.log('═══════════════════════════════════════════════════════════════');
