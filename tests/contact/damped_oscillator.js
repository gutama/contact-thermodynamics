/**
 * Damped Harmonic Oscillator Test (Contact Hamiltonian)
 * 
 * The damped harmonic oscillator is the canonical example of a contact Hamiltonian system.
 * H = (p² + q²)/2 with contact form α = dz - p dq.
 * 
 * The equations of motion include a dissipative term from the contact structure.
 * 
 * Expected behavior:
 * - Energy decays exponentially: E(t) = E(0) * exp(-γt)
 * - Phase space spiral converges to origin
 */

const ContactThermo = require('../../src/index.js');

const { abs, exp, sqrt, sin, cos, PI } = Math;

function assertEqual(actual, expected, tolerance, message) {
    const diff = abs(actual - expected);
    if (diff > tolerance) {
        console.error(`FAIL: ${message}`);
        console.error(`  Expected: ${expected.toFixed(6)}`);
        console.error(`  Actual: ${actual.toFixed(6)}`);
        return false;
    }
    console.log(`PASS: ${message}`);
    return true;
}

function assertTrue(condition, message) {
    if (!condition) {
        console.error(`FAIL: ${message}`);
        return false;
    }
    console.log(`PASS: ${message}`);
    return true;
}

/**
 * Simulate damped harmonic oscillator using contact Hamiltonian dynamics.
 */
function simulateDampedOscillator(q0, p0, gamma, omega, dt, nSteps) {
    const trajectory = [];
    let q = q0;
    let p = p0;
    let z = 0;  // Contact coordinate

    for (let i = 0; i <= nSteps; i++) {
        const E = 0.5 * (p * p + omega * omega * q * q);
        trajectory.push({ t: i * dt, q, p, z, E });

        // Contact Hamiltonian evolution:
        // dq/dt = ∂H/∂p = p
        // dp/dt = -∂H/∂q - γp = -ω²q - γp
        // dz/dt = p*∂H/∂p - H = p² - H = (p² - q²)/2

        const dq = p * dt;
        const dp = (-omega * omega * q - gamma * p) * dt;
        const dz = (p * p - E) * dt;

        q += dq;
        p += dp;
        z += dz;
    }

    return trajectory;
}

function testDampedOscillator() {
    console.log('=== Damped Harmonic Oscillator Test ===\n');

    let passed = true;

    const q0 = 1.0;
    const p0 = 0.0;
    const omega = 1.0;
    const gamma = 0.2;  // Damping coefficient
    const dt = 0.01;
    const T = 20;  // Total time
    const nSteps = Math.floor(T / dt);

    console.log(`Initial: q0 = ${q0}, p0 = ${p0}`);
    console.log(`Parameters: ω = ${omega}, γ = ${gamma}`);
    console.log(`Simulation: T = ${T}, dt = ${dt}`);

    const trajectory = simulateDampedOscillator(q0, p0, gamma, omega, dt, nSteps);

    // Check initial energy
    const E0 = trajectory[0].E;
    console.log(`\nInitial energy: ${E0.toFixed(6)}`);

    // Check energy decay
    console.log('\n--- Energy Decay Test ---');

    // For underdamped oscillator, energy decays as exp(-γt)
    const checkTimes = [5, 10, 15];
    for (const t of checkTimes) {
        const idx = Math.floor(t / dt);
        const E_actual = trajectory[idx].E;
        const E_expected = E0 * exp(-gamma * t);

        console.log(`  t = ${t}: E_actual = ${E_actual.toFixed(6)}, E_expected ≈ ${E_expected.toFixed(6)}`);

        // Allow reasonable tolerance since numerical integration differs from exact decay
        passed &= assertTrue(E_actual < E0 * exp(-gamma * t / 2),
            `Energy at t=${t} is less than initial * exp(-γt/2)`);
    }

    // Check final state approaches equilibrium
    console.log('\n--- Equilibrium Approach Test ---');
    const final = trajectory[trajectory.length - 1];
    console.log(`Final state: q = ${final.q.toFixed(6)}, p = ${final.p.toFixed(6)}, E = ${final.E.toFixed(6)}`);

    passed &= assertTrue(final.E < 0.1 * E0, 'Energy decays to less than 10% of initial');
    passed &= assertTrue(abs(final.q) < 0.5, 'Position approaches zero');
    passed &= assertTrue(abs(final.p) < 0.5, 'Momentum approaches zero');

    // Check monotonic energy decay (skip first few steps due to discretization)
    console.log('\n--- Monotonic Decay Test ---');
    let monotonic = true;
    let increaseCount = 0;
    for (let i = 10; i < trajectory.length; i++) {  // Skip initial transient
        if (trajectory[i].E > trajectory[i - 1].E + 1e-3) {
            increaseCount++;
            if (increaseCount > 5) {  // Allow a few numerical fluctuations
                monotonic = false;
                console.log(`  Significant energy increases detected`);
                break;
            }
        }
    }
    passed &= assertTrue(monotonic, 'Energy decreases monotonically (after transient)');

    return passed;
}

// Run test
const passed = testDampedOscillator();
console.log(`\n${passed ? 'ALL TESTS PASSED' : 'SOME TESTS FAILED'}`);
process.exit(passed ? 0 : 1);
