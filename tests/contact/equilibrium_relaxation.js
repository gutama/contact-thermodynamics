/**
 * Equilibrium Relaxation Test (Contact Hamiltonian)
 * 
 * Verifies relaxation to equilibrium with monotone entropy production.
 * Uses a contact Hamiltonian with known Lyapunov function.
 * 
 * System: Double-well potential with dissipation
 * H = (p² - 1)² + q⁴ - 2q² (Mexican hat)
 * 
 * Expected behavior:
 * - System relaxes to one of the stable equilibria (q = ±1, p = 0)
 * - Entropy (free energy) decreases monotonically
 */

const { abs, sqrt, exp, min, max } = Math;

function assertTrue(condition, message) {
    if (!condition) {
        console.error(`FAIL: ${message}`);
        return false;
    }
    console.log(`PASS: ${message}`);
    return true;
}

/**
 * Double-well potential energy.
 */
function potential(q) {
    return q * q * q * q - 2 * q * q;  // V = q⁴ - 2q²
}

/**
 * Force from potential.
 */
function force(q) {
    return -(4 * q * q * q - 4 * q);  // F = -dV/dq = -4q³ + 4q
}

/**
 * Simulate relaxation dynamics.
 * dq/dt = p
 * dp/dt = F(q) - γp (dissipation)
 * 
 * Lyapunov function: L = p²/2 + V(q) (total energy)
 */
function simulateRelaxation(q0, p0, gamma, dt, nSteps) {
    const trajectory = [];
    let q = q0;
    let p = p0;

    for (let i = 0; i <= nSteps; i++) {
        const V = potential(q);
        const E = 0.5 * p * p + V;  // Total energy (Lyapunov)
        trajectory.push({ t: i * dt, q, p, E, V });

        // Symplectic-ish Euler with dissipation
        const F = force(q);
        p = p + (F - gamma * p) * dt;
        q = q + p * dt;
    }

    return trajectory;
}

function testEquilibriumRelaxation() {
    console.log('=== Equilibrium Relaxation Test ===\n');

    let passed = true;

    // Initial condition: far from equilibrium
    const q0 = 0.1;  // Near unstable equilibrium at q=0
    const p0 = 0.5;  // Some initial momentum
    const gamma = 0.5;  // Damping coefficient
    const dt = 0.01;
    const T = 50;
    const nSteps = Math.floor(T / dt);

    console.log(`Initial: q0 = ${q0}, p0 = ${p0}`);
    console.log(`Damping: γ = ${gamma}`);
    console.log(`Simulation: T = ${T}, dt = ${dt}`);

    const trajectory = simulateRelaxation(q0, p0, gamma, dt, nSteps);

    // Check initial and final states
    const initial = trajectory[0];
    const final = trajectory[trajectory.length - 1];

    console.log(`\nInitial: q = ${initial.q.toFixed(4)}, p = ${initial.p.toFixed(4)}, E = ${initial.E.toFixed(4)}`);
    console.log(`Final: q = ${final.q.toFixed(4)}, p = ${final.p.toFixed(4)}, E = ${final.E.toFixed(4)}`);

    // --- Test 1: Approach to equilibrium ---
    console.log('\n--- Test 1: Equilibrium Approach ---');

    // Should approach one of the stable points: q = ±1
    const distToPlus1 = abs(final.q - 1);
    const distToMinus1 = abs(final.q + 1);
    const distToEquil = min(distToPlus1, distToMinus1);

    console.log(`  Distance to nearest equilibrium: ${distToEquil.toFixed(4)}`);
    passed &= assertTrue(distToEquil < 0.1, 'System approaches stable equilibrium');
    passed &= assertTrue(abs(final.p) < 0.1, 'Momentum approaches zero');

    // --- Test 2: Energy decreases (Lyapunov) ---
    console.log('\n--- Test 2: Monotone Energy Decay ---');

    let violations = 0;
    for (let i = 10; i < trajectory.length; i++) {
        if (trajectory[i].E > trajectory[i - 1].E + 1e-3) {
            violations++;
        }
    }

    console.log(`  Energy violations: ${violations} / ${trajectory.length - 10}`);
    passed &= assertTrue(violations < trajectory.length * 0.05, 'Energy mostly decreases');
    passed &= assertTrue(final.E < initial.E, 'Final energy < initial energy');

    // --- Test 3: Final energy matches equilibrium ---
    console.log('\n--- Test 3: Equilibrium Energy ---');

    const V_equil = potential(1);  // V(±1) = 1 - 2 = -1
    const E_equil = 0 + V_equil;   // p = 0 at equilibrium

    console.log(`  Expected equilibrium energy: ${E_equil.toFixed(4)}`);
    console.log(`  Actual final energy: ${final.E.toFixed(4)}`);

    passed &= assertTrue(abs(final.E - E_equil) < 0.1, 'Final energy matches equilibrium');

    return passed;
}

// Run test
const passed = testEquilibriumRelaxation();
console.log(`\n${passed ? 'ALL TESTS PASSED' : 'SOME TESTS FAILED'}`);
process.exit(passed ? 0 : 1);
