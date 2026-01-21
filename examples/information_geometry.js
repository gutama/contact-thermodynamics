/**
 * Information Geometry Verification Script
 * 
 * Verifies the contact geometry implementation for probability distributions.
 * Based on Baez's "Information Geometry" Part 18 & 19.
 */

const { ProbabilityManifold } = require('../src/information-geometry');

function runVerification() {
    console.log('--- Information Geometry Verification ---');

    const n = 2; // 2 microstates (coin flip)
    console.log(`\nInitializing Probability Manifold for N=${n} microstates...`);
    const manifold = new ProbabilityManifold(n);

    // 1. Define a probability distribution
    // q = [0.3, 0.7]
    const q = [0.3, 0.7];
    console.log(`\nDistribution q: [${q.join(', ')}]`);

    // 2. Compute Entropy
    const S = manifold.entropy(q);
    console.log(`Shannon Entropy S(q): ${S.toFixed(6)}`);
    // Expected: -(0.3 ln 0.3 + 0.7 ln 0.7)
    const expectedS = -(0.3 * Math.log(0.3) + 0.7 * Math.log(0.7));
    console.log(`Expected S:       ${expectedS.toFixed(6)}`);

    // 3. Compute Surprisal (Conjugate Variables)
    const p = manifold.surprisal(q);
    console.log(`\nSurprisal p (Conjugate Momentum):`);
    console.log(`p_1: ${p[0].toFixed(6)} (expected: ${(-Math.log(0.3) - 1).toFixed(6)})`);
    console.log(`p_2: ${p[1].toFixed(6)} (expected: ${(-Math.log(0.7) - 1).toFixed(6)})`);

    // 4. Verify Contact Structure (Legendrian Condition)
    console.log(`\nVerifying Legendrian Condition (alpha restricted to submanifold = 0)...`);
    const check = manifold.checkLegendrianCondition(q);

    let allPassed = true;
    check.forEach(res => {
        console.log(`Direction q_${res.k}: Contraction = ${res.contraction.toExponential(4)} [${res.passed ? 'PASS' : 'FAIL'}]`);
        if (!res.passed) allPassed = false;
    });

    if (allPassed) {
        console.log('SUCCESS: The submanifold of probability distributions is Legendrian.');
    } else {
        console.log('FAILURE: The submanifold is NOT Legendrian.');
    }

    // 5. Contact Volume element
    console.log(`\nChecking Contact Volume Form (alpha ^ (d alpha)^n)...`);
    try {
        const volume = manifold.volumeForm(q, p);
        console.log('Volume Form computed successfully.');
        console.log('Magnitude:', volume.norm());
        if (volume.norm() > 1e-10) {
            console.log('SUCCESS: Contact volume form is non-vanishing (Contact Manifold structure verified).');
        } else {
            console.log('WARNING: Contact volume form vanishes (Degenerate contact structure).');
        }
    } catch (e) {
        console.log('Error computing volume form:', e.message);
    }
}

runVerification();
