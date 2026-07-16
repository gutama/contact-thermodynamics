# Phase 1b — Versor Contactomorphisms in the Degenerate Carrier Cl(n,n,1)

**A focused follow-up to Phase 1** (`docs/PHASE1_GA_CONTACT_FORM.md`),
prompted by an external reviewer. It answers one precise question and reports a
**mixed result**: one construction genuinely works, two provably do not.

New code is **additive**: `src/contact/ga-contact-versor.js` (class
`GAContactVersor`), exported through `src/contact/index.js`, tested in
`tests/contact/ga_contact_versor.js`. Nothing in Phase 1 is changed.

---

## 0. Why this phase exists

Phase 1 represented the contact structure (`α`, `dα`, Reeb `R`, `X_H`) in a
**Euclidean** carrier `Cl(2n+1, 0, 0)`, using only the outer product for
`α`/`dα`/non-degeneracy and the left contraction for the pairing. A reviewer
pointed out that this construction is **signature-invariant**: running the exact
same `α ∧ (dα)^n` (and the Reeb / `X_H` machinery) in `Cl(5,0,0)`, `Cl(0,5,0)`,
`Cl(3,2,0)`, `Cl(2,2,1)` gives **numerically identical** results. That is a
proof that the geometric product — the Clifford/metric part — does **no work**
in Phase 1: it is exterior (Grassmann) algebra with a Clifford wrapper that adds
nothing.

The reviewer's follow-up: the Phase-1 §3.3 *rejected* carrier `Cl(n,n,1)` — with
opposite signs on the q- and p-generators and a genuine null generator for `ds`
— makes `dα` a real **metric bivector**. In that carrier, can
**contactomorphisms** (diffeomorphisms `φ` with `φ*α = f·α`, `f` a nonvanishing
conformal factor) be represented as **versors**: exponentials of bivectors
applied by a sandwich `v ↦ V v V⁻¹`, exactly as ordinary rotors in `Cl(p,q)`
generate rotations/boosts?

This is genuinely open. The answer here is **partly yes, mostly no**, stated
honestly below.

---

## 1. The carrier `Cl(n,n,1)`

```
eq_1..eq_n   e² = +1     (q-directions)
ep_1..ep_n   e² = −1     (p-directions)
es           e² =  0     (null generator = the ds / Reeb axis)
```

The contact bivector is

```
B = dα = Σ_a eq_a ∧ ep_a ,      (eq_a ep_a)² = −eq_a² ep_a² = −(+1)(−1) = +1
```

so **`B² = +1` per pair — a genuine hyperbolic metric bivector**. This is the
key contrast with Phase 1: the signature now *matters*. (Verified:
`tests/contact/ga_contact_versor.js` §I.)

### 1.1 Darboux ↔ light-cone identification

A hyperbolic versor `V = exp(θ B/2)` acts on the orthonormal `(eq, ep)` basis as
a *boost* (mixing eq and ep). To recover the physically meaningful **diagonal
squeeze**, identify the Darboux coordinates `(q_a, p_a)` with the **null
(light-cone) basis** of each hyperbolic plane:

```
nq_a = eq_a − ep_a          ⟨nq_a, nq_a⟩ = 0
np_a = (eq_a + ep_a)/2      ⟨np_a, np_a⟩ = 0 ,  ⟨nq_a, np_a⟩ = 1
```

A point `(q, p, s)` is carried by the grade-1 vector

```
X = Σ_a ( q_a nq_a + p_a np_a ) + s·es .
```

In this basis the versor `V = exp(θ B/2)` acts **diagonally**:
`nq_a ↦ e^{θ} nq_a`, `np_a ↦ e^{−θ} np_a`, i.e. the symplectic squeeze
`q_a ↦ e^{θ} q_a`, `p_a ↦ e^{−θ} p_a`.

---

## 2. The honest test (`contactResidual`)

Every claim below is checked against the **definition** of a contactomorphism,
not by construction. Given an induced (affine) coordinate map `φ(X) = M·X + b`,
we pull `α = ds − Σ_a p_a dq^a` back:

```
(φ*α)_X(v) = (Mv)_s − Σ_a (M X + b)_{p_a} (Mv)_{q_a}
```

and require it to equal `f·α_X(v) = f·(v_s − Σ_a X_{p_a} v_{q_a})` for a
**constant** `f`, sampling many random `(X, v)`. The worst residual over all
coefficients and samples is reported. A tautology-free check: it fails for
non-contact maps (§4 shear residual ≈ 0.7–6) and passes only for genuine
contactomorphisms (residual ≈ 1e-15).

---

## 3. POSITIVE — the squeeze IS a versor-generated strict contactomorphism

The symplectic squeeze

```
φ_θ : (q, p, s) ↦ (e^{θ} q, e^{−θ} p, s)
```

is realized **exactly** by the sandwich with `V = exp(θ·eq∧ep / 2)`:

- `V` is a genuine unit versor (`V Ṽ = 1`), because `B² = +1` (hyperbolic rotor).
- Its induced Darboux map is `diag(e^{θ}, e^{−θ}, 1)` — no coordinate leakage.
- **`φ_θ*α = α` exactly (`f = 1`), pullback residual ~1e-15**: it is a *strict*
  contactomorphism (an actual contact-form-preserving symmetry, `p dq` is
  preserved since `(e^{−θ}p)·d(e^{θ}q) = p dq`).
- It is an isometry (`|det| = 1`), as any versor sandwich must be.
- Generalises cleanly to `n = 2` with independent rapidities per plane.

**This is a real geometric-algebra fact where the metric does genuine work.**
`B² = +1` (hyperbolic, from the `(+,−)` signature) is exactly what makes the
scaling `e^{±θ}`; in a Euclidean carrier the same bivector would give a rotation,
not a squeeze. Unlike Phase 1, the signature is *load-bearing* here. (Tests §II.)

---

## 4. NEGATIVE — the Reeb flow is NOT a versor sandwich

The Reeb flow `φ_t : (q, p, s) ↦ (q, p, s + t)` is the simplest
contactomorphism (`f = 1`). As an **affine** coordinate map (linear part `I`,
translation `t` along `es`) it passes the contact test — that is the baseline.

But a versor sandwich `v ↦ V v V⁻¹` is **linear** and therefore **fixes the
origin** (`V·0·V⁻¹ = 0`), whereas the Reeb flow sends `(0,0,0) ↦ (0,0,t)`. So no
versor sandwich in the direct vector representation can reproduce it.

The natural candidates — the **null-generator "translator" versors**
`exp(t·es∧eq / 2)` and `exp(t·es∧ep / 2)` (parabolic, since `(es∧e)² = 0`) — do
fix the origin and instead produce a **coordinate shear of `s`**
(`s ↦ s + c·q` or `s ↦ s + c·p`), which is **not a contactomorphism** at all
(`φ*α = ds − (p−c) dq ≠ f·α`; residual ≈ 0.7 and up). This is precisely the
Phase-1 §3.3 prediction that contractions/operations involving the null
generator behave badly for the operations we need. (Tests §III.)

Representing a constant translation as a versor would require a *separate*
homogeneous/projective generator (as in PGA), i.e. enlarging the algebra beyond
`Cl(n,n,1)` — which is exactly the extra structure this phase was testing
whether we could avoid.

---

## 5. NEGATIVE — the conformal dilation (f ≠ 1) is NOT a versor sandwich

The contact dilation

```
φ_λ : (q, p, s) ↦ (q, λ p, λ s) ,      φ_λ*α = λ·α   (f = λ ≠ 1)
```

is a genuine contactomorphism with a **non-unit conformal factor** (verified:
pullback gives `f = λ`, residual 0). But its linear part has **`det = λ² ≠ ±1`**:
it is a genuine **non-isometric scaling**. A versor sandwich `V x V⁻¹` is an
isometry (`|det| = 1`) *regardless of normalization* — the scalar `V Ṽ` cancels
between `V` and `V⁻¹` — so **no versor can produce a conformal factor `f ≠ 1`**.
We confirm numerically that every squeeze versor is an isometry with `f = 1`
(Tests §IV).

This is the actual obstruction anticipated in the task brief: **versors generate
isometries (`f = 1`); genuine conformal rescalings (`f ≠ 1`) are outside the
versor group.** In ordinary GA, dilations are versors only in the *conformal*
model `Cl(n+1,1)` (two extra generators) — again, structure beyond `Cl(n,n,1)`.

---

## 6. HYBRID — versor ∘ explicit dilation reproduces conformal contactomorphisms

Composing the squeeze versor (`f = 1`) with an **explicit scalar dilation** on
the `(p, s)` block (`f = λ`) reproduces a conformal contactomorphism with any
`f = λ` (Tests §V, `det = λ²` confirms a non-versor scaling was needed). But the
dilation step is **not** a versor operation. Honest conclusion:

> **Versors alone are insufficient for the contact group in `Cl(n,n,1)`. A
> versor-plus-explicit-dilation hybrid is required, and the versor part only
> ever supplies the isometric (`f = 1`) piece.**

---

## 7. Verdict

| Contactomorphism | `f` | Versor sandwich in `Cl(n,n,1)`? |
|---|---|---|
| Symplectic squeeze `(e^{θ}q, e^{−θ}p, s)` | 1 | **YES** — `exp(θ eq∧ep/2)`, exact, metric does real work |
| Reeb flow `(q, p, s+t)` | 1 | **NO** — affine translation; sandwich fixes the origin |
| Conformal dilation `(q, λp, λs)` | λ | **NO** — needs `det = λ²`; versors are isometries |
| General conformal | ≠1 | only as **versor ∘ dilation hybrid** (not a pure versor) |

**Is `Cl(n,n,1)` a genuinely useful GA representation for contact geometry?**
*Partially, and honestly scoped.* It is a real improvement over the Phase-1
Euclidean carrier in one specific sense: the metric is **not** inert — `dα` is a
true hyperbolic bivector and its versors generate the **linear symplectic
isometry subgroup** of contactomorphisms (the strict, `f = 1`, squeeze/boost
transformations) genuinely, with the signature doing load-bearing work.

But the full contactomorphism group is **not** a versor group here. Two of its
most basic elements — the Reeb translation and the conformal dilation — are
provably outside the reach of pure versor sandwiches (one because it is affine,
the other because it is non-isometric). Capturing them requires enlarging the
algebra (a projective generator for translations, conformal generators for
dilations) or falling back to an explicit hybrid — i.e. reintroducing exactly
the extra structure that would undercut the claim of a clean versor
representation.

**Bottom line for a paper:** this supports an *honest, positive-but-bounded*
contribution — "in `Cl(n,n,1)`, contact *symplectic isometries* (`f = 1`) are
genuine bivector-generated versors, with `dα` as the metric bivector" — and an
*honest negative* — "the Reeb flow and conformal dilations are **not** versors;
the contact group is not a versor group in this carrier." It does **not**
support a claim that general contactomorphisms are GA versors. Framed as the
former, it is a genuine (if modest) contribution; framed as the latter, it would
be an overclaim of the kind Phase 0/1 were careful to avoid.
