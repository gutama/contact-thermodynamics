/**
 * GA-Native Contact Versor Tests (Phase 1b)
 *
 * Tests whether contactomorphisms (φ*α = f·α) can be represented as versors —
 * exponentials of bivectors applied by a sandwich — in the degenerate carrier
 * Cl(n,n,1), where dα = Σ eq_a∧ep_a is a genuine metric bivector.
 *
 * The result is MIXED and every assertion below checks the *mathematical
 * definition* of a contactomorphism via an independent pullback residual
 * (`contactResidual`), never by construction:
 *
 *   POSITIVE  — the symplectic squeeze (q,p,s) ↦ (e^θ q, e^{−θ} p, s) IS a
 *               versor sandwich and IS a strict contactomorphism (f = 1).
 *   NEGATIVE  — the Reeb flow (q,p,s) ↦ (q,p,s+t) is a translation (affine); a
 *               linear sandwich fixes the origin, so no versor reproduces it,
 *               and the null-generator "translator" versors give non-contact
 *               shears.
 *   NEGATIVE  — the conformal dilation (q,p,s) ↦ (q,λp,λs) IS a contactomorphism
 *               (f = λ) but requires det = λ² ≠ ±1; versor sandwiches are
 *               isometries (|det| = 1), so no versor reproduces it (f ≠ 1
 *               is unreachable).
 *
 * See docs/PHASE1B_VERSOR_CONTACTOMORPHISMS.md.
 */

const {
    createTestTracker, assert, assertApprox, section, summary
} = require('../test-utils.js');

const { GAContactVersor } = require('../../src/contact/ga-contact-versor.js');

const t = createTestTracker();
const { abs, log, exp } = Math;

// ---------------------------------------------------------------------------
section('I. Carrier Cl(n,n,1): dα is a genuine (hyperbolic) metric bivector');
// ---------------------------------------------------------------------------

{
    const V = new GAContactVersor(1);
    assert(t, V.algebra.p === 1 && V.algebra.q === 1 && V.algebra.r === 1,
        'n=1 carrier is Cl(1,1,1) (one +1, one −1, one null generator)');

    // Signature of the generators, as the construction requires.
    assertApprox(t, V.eq(0).mul(V.eq(0)).scalar(), 1, 1e-12, 'eq² = +1 (q-direction)');
    assertApprox(t, V.ep(0).mul(V.ep(0)).scalar(), -1, 1e-12, 'ep² = −1 (p-direction)');
    assertApprox(t, V.es().mul(V.es()).scalar(), 0, 1e-12, 'es² = 0 (null / Reeb axis)');

    // The Phase-1 contrast: here B = dα squares to +1 (metric does work),
    // whereas in the Euclidean Phase-1 carrier the signature was irrelevant.
    const B = V.metricBivector();
    assert(t, B.grade(2).norm() > 0 && B.grade(0).isZero(),
        'dα = eq∧ep is a pure grade-2 bivector');
    assertApprox(t, B.mul(B).scalar(), 1, 1e-12,
        'B² = +1: dα is a genuine hyperbolic metric bivector (unlike Phase-1)');

    // Light-cone identification pairing ⟨nq, np⟩ = 1, both null.
    const sym = (a, b) => a.mul(b).add(b.mul(a)).scale(0.5).scalar();
    assertApprox(t, sym(V.nq(0), V.nq(0)), 0, 1e-12, '⟨nq,nq⟩ = 0 (null)');
    assertApprox(t, sym(V.np(0), V.np(0)), 0, 1e-12, '⟨np,np⟩ = 0 (null)');
    assertApprox(t, sym(V.nq(0), V.np(0)), 1, 1e-12, '⟨nq,np⟩ = 1 (Darboux pairing)');
}

// ---------------------------------------------------------------------------
section('II. POSITIVE: the squeeze versor IS a strict contactomorphism (f=1)');
// ---------------------------------------------------------------------------

{
    const V = new GAContactVersor(1);
    const theta = log(2);                 // λ = e^θ = 2
    const S = V.squeezeVersor(theta);

    // It is a genuine unit versor (rotor): S S̃ = 1.
    assertApprox(t, S.mul(S.reverse()).scalar(), 1, 1e-10,
        'squeeze versor is a unit versor (S S̃ = 1)');

    // The sandwich acts as the diagonal squeeze in Darboux coordinates.
    const M = V.inducedMatrix(S);
    assertApprox(t, M[0][0], exp(theta), 1e-9, 'q ↦ e^{θ} q  (M_qq = λ)');
    assertApprox(t, M[1][1], exp(-theta), 1e-9, 'p ↦ e^{−θ} p (M_pp = 1/λ)');
    assertApprox(t, M[2][2], 1, 1e-9, 's ↦ s (fiber fixed)');
    // no off-diagonal leakage
    let offDiag = 0;
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) if (i !== j) offDiag = Math.max(offDiag, abs(M[i][j]));
    assert(t, offDiag < 1e-9, 'squeeze induced map is diagonal (no coordinate leakage)');

    // THE genuine, non-tautological check: pull α back and confirm φ*α = 1·α.
    const chk = V.versorContactResidual(S);
    assert(t, chk.isContacto, 'squeeze sandwich satisfies φ*α = f·α (pullback residual ~0)');
    assertApprox(t, chk.f, 1, 1e-9,
        'squeeze conformal factor f = 1 (a STRICT contactomorphism / symmetry)');

    // Versors are isometries: |det| = 1 (contrast with the dilation below).
    assertApprox(t, abs(GAContactVersor.det(M)), 1, 1e-9,
        'squeeze is an isometry (|det| = 1)');
}

{
    // Positive result generalises to n = 2 (independent rapidities per plane).
    const V = new GAContactVersor(2);
    const S = V.squeezeVersor([0.5, -0.3]);
    assertApprox(t, S.mul(S.reverse()).scalar(), 1, 1e-9, 'n=2 squeeze is a unit versor');
    const chk = V.versorContactResidual(S);
    assert(t, chk.isContacto, 'n=2 squeeze satisfies φ*α = f·α');
    assertApprox(t, chk.f, 1, 1e-9, 'n=2 squeeze is strict (f = 1)');
}

// ---------------------------------------------------------------------------
section('III. NEGATIVE: the Reeb flow is NOT a versor sandwich');
// ---------------------------------------------------------------------------

{
    const V = new GAContactVersor(1);

    // The Reeb flow (q,p,s) ↦ (q,p,s+t) is affine (linear part = I, translation
    // along es). It IS a contactomorphism (f = 1) as a coordinate map...
    const reeb = V.contactResidual([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 0, 0.7]);
    assert(t, reeb.isContacto && abs(reeb.f - 1) < 1e-12,
        'Reeb flow as an AFFINE map is a contactomorphism with f = 1 (baseline)');

    // ...but a versor sandwich is LINEAR and fixes the origin, so it cannot
    // produce the translation part. The candidate null-generator "translator"
    // versors instead produce a coordinate SHEAR of s, which is NOT a
    // contactomorphism.
    for (const dir of ['q', 'p']) {
        const Tw = V.nullTranslatorVersor(0, dir, 0.7);
        // The origin is fixed (hallmark of a linear sandwich) — cannot be a
        // nonzero translation.
        const originImg = V.readCoords(V.apply(Tw, V.pointVector(0, 0, 0)));
        assert(t, abs(originImg.q[0]) < 1e-12 && abs(originImg.p[0]) < 1e-12 && abs(originImg.s) < 1e-12,
            `es∧e${dir} translator fixes the origin (linear sandwich ⇒ cannot translate)`);
        const chk = V.versorContactResidual(Tw);
        assert(t, !chk.isContacto,
            `es∧e${dir} translator sandwich is NOT a contactomorphism (produces a non-contact shear)`);
    }
}

// ---------------------------------------------------------------------------
section('IV. NEGATIVE: the conformal dilation (f=λ≠1) is NOT a versor sandwich');
// ---------------------------------------------------------------------------

{
    const V = new GAContactVersor(1);
    const lambda = 2;

    // The dilation (q,p,s) ↦ (q, λp, λs) IS a genuine contactomorphism with a
    // NON-unit conformal factor f = λ...
    const dil = [[1, 0, 0], [0, lambda, 0], [0, 0, lambda]];
    const chk = V.contactResidual(dil);
    assert(t, chk.isContacto, 'dilation diag(1,λ,λ) is a contactomorphism (pullback residual ~0)');
    assertApprox(t, chk.f, lambda, 1e-12, 'dilation conformal factor f = λ ≠ 1');

    // ...but it has det = λ² ≠ ±1, so it is not an isometry and thus cannot be a
    // versor sandwich (which always preserves the metric).
    assertApprox(t, GAContactVersor.det(dil), lambda * lambda, 1e-12,
        'dilation det = λ² ≠ ±1 (a genuine non-isometric scaling)');

    // Concretely: EVERY versor sandwich here is an isometry (|det| = 1), so no
    // choice of rapidity can hit the dilation. Sample the squeeze family.
    let allIsometries = true;
    for (const th of [0.2, 0.5, 1.0, -0.7, 1.3]) {
        const M = V.inducedMatrix(V.squeezeVersor(th));
        if (abs(abs(GAContactVersor.det(M)) - 1) > 1e-9) allIsometries = false;
        // and none of them is a non-unit-f contactomorphism
        const c = V.versorContactResidual(V.squeezeVersor(th));
        if (c.isContacto && abs(c.f - 1) > 1e-9) allIsometries = false;
    }
    assert(t, allIsometries,
        'all squeeze versors are isometries with f = 1 — none yields the f = λ dilation');
}

// ---------------------------------------------------------------------------
section('V. HYBRID: versor∘dilation reproduces a general conformal contactomorphism');
// ---------------------------------------------------------------------------

{
    // Honest positive-with-caveat: since the squeeze (f=1) is a versor and the
    // scalar dilation supplies f=λ, their COMPOSITION reproduces a conformal
    // contactomorphism — but the dilation step is an explicit scaling, NOT a
    // versor operation. This is the "versor + dilation hybrid" conclusion.
    const V = new GAContactVersor(1);
    const lambda = 1.5, theta = 0.4;

    // matrix multiply helper
    const mul = (A, B) => A.map((row, i) => B[0].map((_, j) =>
        row.reduce((acc, _v, k) => acc + A[i][k] * B[k][j], 0)));

    const squeezeM = V.inducedMatrix(V.squeezeVersor(theta));         // f = 1
    const dilM = [[1, 0, 0], [0, lambda, 0], [0, 0, lambda]];         // f = λ
    const hybrid = mul(dilM, squeezeM);                              // dilation ∘ squeeze

    const chk = V.contactResidual(hybrid);
    assert(t, chk.isContacto,
        'versor∘dilation hybrid IS a contactomorphism (φ*α = f·α)');
    assertApprox(t, chk.f, lambda, 1e-9,
        'hybrid conformal factor f = λ (from the dilation factor, not the versor)');
    assert(t, abs(abs(GAContactVersor.det(hybrid)) - lambda * lambda) < 1e-9,
        'hybrid has det = λ² — confirming a non-versor scaling was required');
}

process.exit(summary(t));
