/**
 * GA-Native Contact Form Tests (Phase 1)
 *
 * Validates that the GA-native representation of the contact structure
 * (src/contact/ga-contact-form.js) reproduces the conventional
 * ContactManifold / ContactHamiltonian results exactly, that the genuine
 * α ∧ (dα)^n non-degeneracy check works (positive and negative cases), and
 * that the closed-form bivector exponentials (gap G3) and rotor sign
 * convention behave as documented.
 */

const {
    createTestTracker, assert, assertApprox, section, summary
} = require('../test-utils.js');

const { ContactManifold } = require('../../src/contact/manifold.js');
const { ContactHamiltonian } = require('../../src/contact/hamiltonian.js');
const { GAContactForm } = require('../../src/contact/ga-contact-form.js');
const RGA = require('../../src/geometry/riemannian-ga.js');
const GA = require('../../src/algebra/multivector.js');

const t = createTestTracker();
const { abs } = Math;

// ---------------------------------------------------------------------------
section('I. Contact 1-form α and its exterior derivative dα');
// ---------------------------------------------------------------------------

{
    const M = new ContactManifold(['q'], ['p'], 'A');
    const ga = new GAContactForm(M);
    const pt = M.point({ q: 0.7, p: 1.3, A: 0.2 });

    const alpha = ga.contactForm(pt);
    // α = dA − p dq  ⇒  coeff(dA)=1, coeff(dq)=−p
    assertApprox(t, ga.coeff(alpha, 'A'), 1, 1e-12, 'α has coeff +1 on dA');
    assertApprox(t, ga.coeff(alpha, 'q'), -1.3, 1e-12, 'α has coeff −p on dq');

    const dAlpha = ga.contactFormDerivative();
    // dα = dq ∧ dp is grade 2 and coordinate-independent
    assert(t, dAlpha.grade(2).norm() > 0, 'dα is a nonzero grade-2 bivector');
    assert(t, dAlpha.grade(1).isZero() && dAlpha.grade(0).isZero(),
        'dα is pure grade-2');
}

// ---------------------------------------------------------------------------
section('II. Genuine non-degeneracy α ∧ (dα)^n (replaces the n! stub)');
// ---------------------------------------------------------------------------

{
    // n = 1, 2, 3: the wedge coefficient must have magnitude n! and be nonzero.
    const cases = [
        { base: ['q1'], mom: ['p1'], nFact: 1 },
        { base: ['q1', 'q2'], mom: ['p1', 'p2'], nFact: 2 },
        { base: ['q1', 'q2', 'q3'], mom: ['p1', 'p2', 'p3'], nFact: 6 }
    ];
    for (const c of cases) {
        const M = new ContactManifold(c.base, c.mom, 'A');
        const ga = new GAContactForm(M);
        const pt = M.origin;
        const nd = ga.nonDegeneracy(pt);
        assert(t, nd.isNonDegenerate,
            `n=${c.base.length}: α∧(dα)^n ≠ 0 (contact condition holds)`);
        assertApprox(t, nd.magnitude, c.nFact, 1e-9,
            `n=${c.base.length}: |α∧(dα)^n| = n! = ${c.nFact}`);
    }
}

{
    // Cross-check: the genuine GA magnitude equals the conventional n! stub for
    // the canonical form — confirming consistency, not just agreement of stubs.
    const M = new ContactManifold(['q1', 'q2', 'q3'], ['p1', 'p2', 'p3'], 'A');
    const ga = new GAContactForm(M);
    assertApprox(t, ga.nonDegeneracy(M.origin).magnitude,
        M.verifyContactCondition(M.origin), 1e-9,
        'GA non-degeneracy magnitude matches conventional verifyContactCondition (n!)');
}

// ---------------------------------------------------------------------------
section('III. Negative test: a broken (degenerate) structure IS detected');
// ---------------------------------------------------------------------------

{
    // Build a deliberately degenerate 2-form by dropping one symplectic pair
    // from dα. Then α ∧ (dα)^n must vanish (the top blade is unreachable),
    // demonstrating the check genuinely detects degeneracy.
    const M = new ContactManifold(['q1', 'q2'], ['p1', 'p2'], 'A');
    const ga = new GAContactForm(M);
    const pt = M.origin;

    const alpha = ga.contactForm(pt);
    // Broken dα: only the (q1,p1) pair — rank-deficient, radical now 3-dim.
    const brokenDAlpha = ga._e('q1').wedge(ga._e('p1'));
    const top = ga.wedgeAlphaDAlphaPower(alpha, brokenDAlpha);
    const psBitmap = ga.algebra.size - 1;
    const coeff = top.coeffs[psBitmap] || 0;
    assert(t, abs(coeff) < 1e-12,
        'degenerate dα ⇒ α∧(dα)^n = 0 (degeneracy correctly detected)');

    // Sanity: the correct dα for the same manifold is non-degenerate.
    assert(t, ga.nonDegeneracy(pt).isNonDegenerate,
        'correct dα for the same manifold IS non-degenerate');
}

// ---------------------------------------------------------------------------
section('IV. Reeb field R: α(R) = 1, ι_R dα = 0 (via GA)');
// ---------------------------------------------------------------------------

{
    const M = new ContactManifold(['q1', 'q2'], ['p1', 'p2'], 'A');
    const ga = new GAContactForm(M);
    // Check at a nontrivial point (α depends on p; the Reeb conditions must
    // hold regardless).
    const pt = M.point({ q1: 0.3, q2: -0.5, p1: 0.9, p2: 1.1, A: 0.4 });
    const v = ga.verifyReeb(pt);
    assertApprox(t, v.alphaR, 1, 1e-12, 'α(R) = 1');
    assert(t, v.iRdAlphaIsZero, 'ι_R dα = 0');
}

// ---------------------------------------------------------------------------
section('V. Contact Hamiltonian vector field X_H reproduces Bravetti EOM');
// ---------------------------------------------------------------------------

function compareVectorFields(label, M, Hf, pt) {
    const ga = new GAContactForm(M);
    const H = new ContactHamiltonian(M, Hf);
    const conv = H.vectorField(pt);
    const { components } = ga.hamiltonianVectorField(H, pt);

    let maxDiff = 0;
    for (const c of M.allCoords) {
        maxDiff = Math.max(maxDiff, abs((components[c] || 0) - (conv[c] || 0)));
    }
    assert(t, maxDiff < 1e-6,
        `${label}: GA X_H matches conventional vectorField (max diff ${maxDiff.toExponential(2)})`);

    // GA defining equations are satisfied to floating-point.
    const ver = ga.verifyHamiltonianVectorField(H, pt);
    assert(t, ver.contactResidual < 1e-9,
        `${label}: α(X_H) = −H residual ~0 (${ver.contactResidual.toExponential(2)})`);
    assert(t, ver.symplecticResidual < 1e-9,
        `${label}: ι_{X_H}dα = dH − R(H)α residual ~0 (${ver.symplecticResidual.toExponential(2)})`);
}

{
    // (a) Damped-oscillator-like H (no explicit s-dependence).
    const M1 = new ContactManifold(['q'], ['p'], 'A');
    compareVectorFields('damped osc H=(p²+q²)/2', M1,
        c => 0.5 * (c.p * c.p + c.q * c.q),
        M1.point({ q: 0.7, p: 1.3, A: 0.2 }));

    // (b) Ideal-gas-like / kinetic H, n=2.
    const M2 = new ContactManifold(['q1', 'q2'], ['p1', 'p2'], 'A');
    compareVectorFields('kinetic H=½|p|²+V(q)', M2,
        c => 0.5 * (c.p1 * c.p1 + c.p2 * c.p2) + Math.sin(c.q1) - 0.3 * c.q2,
        M2.point({ q1: 0.3, q2: -0.5, p1: 0.9, p2: 1.1, A: 0.4 }));

    // (c) Explicit s-dependence — exercises the dissipative −p_a ∂H/∂s term.
    const M3 = new ContactManifold(['q'], ['p'], 'A');
    compareVectorFields('dissipative H=½p²+q+γA', M3,
        c => 0.5 * c.p * c.p + c.q + 0.3 * c.A,
        M3.point({ q: 0.7, p: 1.3, A: 0.2 }));

    // (d) s-dependence coupling to both q and p, n=2.
    const M4 = new ContactManifold(['q1', 'q2'], ['p1', 'p2'], 'A');
    compareVectorFields('H=½|p|²+0.2A·q1+0.1A²', M4,
        c => 0.5 * (c.p1 * c.p1 + c.p2 * c.p2) + 0.2 * c.A * c.q1 + 0.1 * c.A * c.A,
        M4.point({ q1: 0.3, q2: -0.5, p1: 0.9, p2: 1.1, A: 0.4 }));
}

// ---------------------------------------------------------------------------
section('VI. Non-degeneracy is invariant along the contact flow');
// ---------------------------------------------------------------------------

{
    // Physical consistency: the contact condition holds all along a trajectory.
    const M = new ContactManifold(['q'], ['p'], 'A');
    const ga = new GAContactForm(M);
    const H = new ContactHamiltonian(M, c => 0.5 * (c.p * c.p + c.q * c.q));
    const traj = H.flow(M.point({ q: 1, p: 0, A: 0 }), 0.05, 20);
    let allNonDeg = true;
    for (const pt of traj) {
        if (!ga.nonDegeneracy(pt).isNonDegenerate) allNonDeg = false;
    }
    assert(t, allNonDeg, 'α∧(dα)^n ≠ 0 at every point along the flow');
}

// ---------------------------------------------------------------------------
section('VII. Bivector exponentials / rotors (gap G3)');
// ---------------------------------------------------------------------------

{
    const A3 = new GA.Algebra(3, 0, 0);
    const e1 = A3.e(1), e2 = A3.e(2);

    // Bivector3D.rotor(θ) uses Hestenes R = exp(−θB̂/2): sandwich rotates CCW.
    const Bxy = new RGA.Bivector3D(0, 0, 1);           // e12 plane
    const R = Bxy.rotor(Math.PI / 2);
    const rotated = R.sandwich(e1);                    // expect +e2
    assertApprox(t, rotated.coeffs[2] || 0, 1, 1e-10,
        'Bivector3D.rotor(90°): e₁ ↦ +e₂ (Hestenes CCW convention)');
    assertApprox(t, rotated.coeffs[1] || 0, 0, 1e-10,
        'Bivector3D.rotor(90°): no residual e₁ component');

    // exp(B) is a unit rotor: R R̃ = 1.
    const Bgen = new RGA.Bivector3D(0.3, -0.5, 0.8);
    const Rg = Bgen.exp();
    assertApprox(t, Rg.mul(Rg.reverse()).scalar(), 1, 1e-10,
        'Bivector3D.exp(B) is a unit rotor (R R̃ = 1)');

    // Documented sign relation to Algebra.rotor: Algebra.rotor(B,θ) = exp(+θB̂/2)
    // sends e₁ ↦ −e₂, i.e. the opposite orientation.
    const codeRotor = A3.rotor(e1.wedge(e2), Math.PI / 2);
    assertApprox(t, codeRotor.sandwich(e1).coeffs[2] || 0, -1, 1e-10,
        'Algebra.rotor(e₁₂,90°): e₁ ↦ −e₂ (opposite orientation, as documented)');
}

{
    // Bivector4D: simple bivector uses closed form; non-simple uses series.
    const simple = new RGA.Bivector4D(1, 0, 0, 0, 0, 0);  // e12 only → simple
    assert(t, simple.isSimple(), 'Bivector4D e₁₂ is simple');
    const Rs = simple.exp();
    assertApprox(t, Rs.mul(Rs.reverse()).scalar(), 1, 1e-10,
        'Bivector4D.exp (simple) is a unit rotor');

    const nonSimple = new RGA.Bivector4D(1, 0, 0, 0, 0, 1); // e12 + e34 → non-simple
    assert(t, !nonSimple.isSimple(), 'Bivector4D e₁₂+e₃₄ is non-simple');
    const Rn = nonSimple.exp();
    // exp(e12+e34) = (cos1 + sin1 e12)(cos1 + sin1 e34); scalar part = cos²(1).
    assertApprox(t, Rn.scalar(), Math.cos(1) * Math.cos(1), 1e-8,
        'Bivector4D.exp (non-simple) scalar part = cos²(1)');
    assertApprox(t, Rn.mul(Rn.reverse()).scalar(), 1, 1e-8,
        'Bivector4D.exp (non-simple) is a unit rotor');
}

// ---------------------------------------------------------------------------
section('VIII. Analytic gradient omitting ∂H/∂s ⇒ R(H) defaults to 0 (no NaN)');
// ---------------------------------------------------------------------------

{
    // A common case: ∂H/∂s = 0, and the analytic gradient function simply omits
    // the fiber component. R(H) must default to 0 rather than injecting NaN.
    const M = new ContactManifold(['q'], ['p'], 'A');
    const ga = new GAContactForm(M);
    // dH deliberately omits 'A' (the fiber coord) → grad['A'] is undefined.
    const H = new ContactHamiltonian(
        M,
        c => 0.5 * c.p * c.p,
        c => ({ q: 0, p: c.p })
    );
    const pt = M.point({ q: 0.7, p: 1.3, A: 0.2 });

    const { components } = ga.hamiltonianVectorField(H, pt);
    const finite = Object.values(components).every(Number.isFinite);
    assert(t, finite,
        'X_H components are all finite (missing ∂H/∂s ⇒ R(H)=0, not NaN)');

    // With R(H)=0 the field reduces to the standard symplectic form: q̇ = p,
    // ṗ = 0, ṡ = p·q̇ − H = p² − ½p².
    assertApprox(t, components.q, 1.3, 1e-12, 'q̇ = ∂H/∂p = p');
    assertApprox(t, components.p, 0, 1e-12, 'ṗ = −∂H/∂q = 0');
    assertApprox(t, components.A, 1.3 * 1.3 - 0.5 * 1.3 * 1.3, 1e-12,
        'ṡ = p·q̇ − H is finite');

    const ver = ga.verifyHamiltonianVectorField(H, pt);
    assert(t, Number.isFinite(ver.contactResidual) && ver.contactResidual < 1e-9,
        'α(X_H) = −H residual finite and ~0 with missing ∂H/∂s');
    assert(t, Number.isFinite(ver.symplecticResidual) && ver.symplecticResidual < 1e-9,
        'ι_{X_H}dα residual finite and ~0 with missing ∂H/∂s');
}

process.exit(summary(t));
