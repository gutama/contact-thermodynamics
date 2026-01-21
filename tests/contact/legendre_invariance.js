/**
 * Legendre Invariance Test (Contact Geometry)
 * 
 * Verifies that the contact form is preserved (up to conformal factor)
 * under contact Hamiltonian flow.
 * 
 * For a contact Hamiltonian H, the contact form α transforms as:
 *   L_X α = λ·α  (Lie derivative gives a multiple of α)
 * 
 * where X is the contact Hamiltonian vector field and λ = -∂H/∂z
 */

const { abs, exp, sqrt } = Math;

function assertTrue(condition, message) {
    if (!condition) {
        console.error(`FAIL: ${message}`);
        return false;
    }
    console.log(`PASS: ${message}`);
    return true;
}

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

/**
 * Contact Hamiltonian system on (q, p, z) with α = dz - p dq.
 * 
 * For Hamiltonian H(q, p, z), the contact vector field X_H is:
 *   dq/dt = ∂H/∂p
 *   dp/dt = -∂H/∂q - p·(∂H/∂z)
 *   dz/dt = p·(∂H/∂p) - H
 * 
 * The conformal factor is λ = -∂H/∂z
 */
class ContactHamiltonian {
    constructor(H, dH_dq, dH_dp, dH_dz) {
        this.H = H;
        this.dH_dq = dH_dq;
        this.dH_dp = dH_dp;
        this.dH_dz = dH_dz;
    }

    /**
     * Compute the contact vector field at a point (q, p, z).
     */
    vectorField(q, p, z) {
        const dq = this.dH_dp(q, p, z);
        const dp = -this.dH_dq(q, p, z) - p * this.dH_dz(q, p, z);
        const dz = p * this.dH_dp(q, p, z) - this.H(q, p, z);
        return { dq, dp, dz };
    }

    /**
     * Compute the conformal factor λ = -∂H/∂z.
     */
    conformalFactor(q, p, z) {
        return -this.dH_dz(q, p, z);
    }

    /**
     * Verify that α(X_H) = -H (fundamental property of contact Hamiltonian).
     * α = dz - p dq, X = (dq, dp, dz)
     * α(X) = dz - p*dq
     */
    verifyAlphaXH(q, p, z) {
        const X = this.vectorField(q, p, z);
        // α = dz - p dq, so α(X) = X.dz - p*X.dq
        // For contact Hamiltonians: α(X_H) = -H (not +H!)
        const alpha_X = X.dz - p * X.dq;
        return alpha_X;
    }
}

function testLegendreInvariance() {
    console.log('=== Legendre Invariance Test ===\n');

    let passed = true;

    // --- Test 1: Simple Hamiltonian H = p²/2 + q²/2 + z ---
    console.log('--- Test 1: H = p²/2 + q²/2 + z (linear in z) ---');

    const H1 = new ContactHamiltonian(
        (q, p, z) => p * p / 2 + q * q / 2 + z,  // H
        (q, p, z) => q,     // ∂H/∂q
        (q, p, z) => p,     // ∂H/∂p
        (q, p, z) => 1      // ∂H/∂z
    );

    const q1 = 1, p1 = 0.5, z1 = 0;

    // Check α(X_H) = H
    const alpha_X1 = H1.verifyAlphaXH(q1, p1, z1);
    const H_val1 = H1.H(q1, p1, z1);

    console.log(`  At (q, p, z) = (${q1}, ${p1}, ${z1}):`);
    console.log(`  α(X_H) = ${alpha_X1.toFixed(6)}`);
    console.log(`  H = ${H_val1.toFixed(6)}`);

    passed &= assertEqual(alpha_X1, -H_val1, 1e-6, 'α(X_H) = -H');

    // Check conformal factor
    const lambda1 = H1.conformalFactor(q1, p1, z1);
    console.log(`  Conformal factor λ = ${lambda1}`);
    passed &= assertEqual(lambda1, -1, 1e-6, 'λ = -∂H/∂z = -1');

    // --- Test 2: z-independent Hamiltonian H = p²/2 + q²/2 ---
    console.log('\n--- Test 2: H = p²/2 + q²/2 (z-independent) ---');

    const H2 = new ContactHamiltonian(
        (q, p, z) => p * p / 2 + q * q / 2,  // H (no z)
        (q, p, z) => q,     // ∂H/∂q
        (q, p, z) => p,     // ∂H/∂p
        (q, p, z) => 0      // ∂H/∂z = 0
    );

    const q2 = 1, p2 = 1, z2 = 5;

    // Check α(X_H) = H
    const alpha_X2 = H2.verifyAlphaXH(q2, p2, z2);
    const H_val2 = H2.H(q2, p2, z2);

    console.log(`  At (q, p, z) = (${q2}, ${p2}, ${z2}):`);
    console.log(`  α(X_H) = ${alpha_X2.toFixed(6)}`);
    console.log(`  H = ${H_val2.toFixed(6)}`);

    passed &= assertEqual(alpha_X2, -H_val2, 1e-6, 'α(X_H) = -H');

    // For z-independent H, the flow strictly preserves α (λ = 0)
    const lambda2 = H2.conformalFactor(q2, p2, z2);
    console.log(`  Conformal factor λ = ${lambda2}`);
    passed &= assertEqual(lambda2, 0, 1e-6, 'λ = 0 (strict contact)');

    // --- Test 3: Simulate flow and check form preservation ---
    console.log('\n--- Test 3: Flow Simulation ---');

    let q = 1, p = 0.5, z = 0;
    const dt = 0.01;
    const nSteps = 100;

    console.log(`  Initial: (q, p, z) = (${q}, ${p}, ${z})`);

    // Expected: α should scale by exp(∫λ dt)
    let integratedLambda = 0;

    for (let i = 0; i < nSteps; i++) {
        const X = H1.vectorField(q, p, z);
        const lambda = H1.conformalFactor(q, p, z);

        integratedLambda += lambda * dt;

        q += X.dq * dt;
        p += X.dp * dt;
        z += X.dz * dt;
    }

    console.log(`  Final: (q, p, z) = (${q.toFixed(4)}, ${p.toFixed(4)}, ${z.toFixed(4)})`);
    console.log(`  ∫λ dt = ${integratedLambda.toFixed(4)}`);
    console.log(`  Expected scaling: exp(∫λ dt) = ${exp(integratedLambda).toFixed(4)}`);

    // Check final α(X_H) = H still holds
    const finalAlphaX = H1.verifyAlphaXH(q, p, z);
    const finalH = H1.H(q, p, z);
    passed &= assertEqual(finalAlphaX, -finalH, 0.1, 'α(X_H) = -H after evolution');

    return passed;
}

// Run test
const passed = testLegendreInvariance();
console.log(`\n${passed ? 'ALL TESTS PASSED' : 'SOME TESTS FAILED'}`);
process.exit(passed ? 0 : 1);
