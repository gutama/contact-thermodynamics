# Phase 1 — GA-Native Contact Form

**Phase 1 of the research plan** to unify contact-geometric thermodynamics with
geometric algebra (GA). This phase makes the contact 1-form `α`, its exterior
derivative `dα`, the Reeb field `R`, and the contact Hamiltonian vector field
`X_H` into genuine `Multivector` objects, and replaces the stubbed
non-degeneracy check with a real `α ∧ (dα)^n` wedge computation. It targets gaps
**G1** and **G2** (with **G3** addressed as far as it supports the bridge). See
`docs/RESEARCH_FOUNDATIONS.md` (Phase 0) for the notation table and gap list.

Everything here is **additive**: `ContactManifold` / `ContactHamiltonian` are
untouched, so existing behavior (and the tests that assert
`verifyContactCondition() === n!`) is preserved. The new code lives in
`src/contact/ga-contact-form.js` (class `GAContactForm`), exported through
`src/contact/index.js`.

---

## 1. The construction

### 1.1 Carrier algebra

The cotangent basis 1-forms `{dq^a, dp_a, ds}` of a canonical `(2n+1)`-dim
contact manifold are modeled as an **orthonormal grade-1 basis of a
Clifford/Grassmann algebra `Cl(2n+1, 0, 0)`**, ordered

```
index 0 .. n-1    →  dq^a   (base coordinates)
index n .. 2n-1   →  dp_a   (conjugate momenta)
index 2n          →  ds     (fiber / Reeb axis; codebase calls it 'A')
```

Two different GA products carry the two different geometric operations we need,
and crucially **only one of them touches the metric**:

- **Outer product** (`op` / `wedge`) builds `α`, `dα`, and the non-degeneracy
  wedge `α ∧ (dα)^n`. The outer product is metric-independent — it only keeps
  disjoint-blade terms with the Koszul reordering sign — so these results are
  pure Grassmann-algebra facts, valid for *any* signature. The Euclidean
  signature of the carrier plays no role here.

- **Left contraction** (`ip` / `dot`) realizes the natural pairing `α(v) = ⟨α,
  v⟩` and the interior product `ι_v ω`. Identifying the tangent basis `∂_i` with
  the cotangent basis `dx^i` under the Euclidean *coordinate identity metric*
  makes `ip(v, ω)` reproduce `⟨ω, v⟩` and `ι_v ω` **exactly in Darboux
  coordinates**. This is the coordinate musical isomorphism; it is faithful for
  the pairing/contraction operations the contact structure is defined by (which
  are themselves metric-free — the identity metric is just a computational
  carrier for the component-wise pairing `Σ α_i v^i`).

Verified numerically (`Cl(3)`): `ip(e₁, e₁∧e₂) = e₂`, `ip(e₂, e₁∧e₂) = −e₁` —
i.e. `ip` is the left contraction and reproduces `ι_v`.

### 1.2 The forms

```
α  = ds − p_a dq^a                 (grade-1, point-dependent through p_a)
dα = d(−p_a dq^a) = −dp_a ∧ dq^a = Σ_a dq^a ∧ dp_a   (grade-2, constant)
```

`dα` is exactly the canonical symplectic 2-form on the `(q,p)`-plane. It
**annihilates `ds`**: `ds` spans `ker(dα)`, the *radical* of the degenerate
2-form.

### 1.3 Where the "Heisenberg / null Reeb direction" lives

The foundations doc suggested a Heisenberg-type algebra with a distinguished
null direction for the Reeb axis. In this construction that null direction is
realized **intrinsically as `ker(dα)`**: `dα` is nondegenerate on the
symplectic `(q,p)`-block and identically zero along `ds`, so the Reeb line is
precisely the kernel/radical of `dα`, and `R = ∂/∂s` is its canonical generator
normalized by `α(R) = 1`. The presymplectic/degenerate structure is present
without introducing a separately-signed or null generator.

**Alternative considered and rejected.** A degenerate Clifford algebra
`Cl(n, n, 1)` with an explicit null generator for `ds` and opposite signs for
`q` vs `p` (so that `dα` would *be* the metric bivector). Rejected because:

1. The contact conditions `α(R)=1`, `ι_R dα=0`, `α(X_H)=−H`,
   `ι_{X_H} dα = dH − R(H)α` are all defined through the natural pairing and the
   interior product, neither of which needs a metric. A Euclidean-identity
   carrier is therefore sufficient *and* simpler.
2. A genuinely null generator (`e² = 0`) makes `ip`-based contractions into that
   direction ambiguous — which would undermine the very operations (`α(R)`,
   `ι_R dα`) the representation exists to compute.

The rejected route becomes attractive only if a later phase wants `dα` itself to
be a *metric* object (e.g. to build boosts/rotors that mix `q`, `p`, `s`); that
is flagged for Phase 2, not needed for G1/G2.

---

## 2. Reeb field and contact Hamiltonian vector field

`R = ∂/∂s` is the `ds` basis vector. The GA identities are checked directly:
`α(R) = ip(R, α) = 1` and `ι_R dα = ip(R, dα) = 0`, both confirmed at nontrivial
points (they are independent of `p`).

`X_H` is **solved from its GA defining equations**, not copied from the Bravetti
EOM. Writing the RHS `dH − R(H)α` as a grade-1 multivector and reading off
coefficients:

```
ι_X dα = Σ_a [ (X·dq^a) dp_a − (X·dp_a) dq^a ]        ⇒
    X^{q_a} =  coeff_{dp_a}(dH − R(H)α)
    X^{p_a} = −coeff_{dq_a}(dH − R(H)α)
α(X) = −H                                             ⇒
    X^{s}   = −H + Σ_a p_a X^{q_a}
```

Expanding the RHS coefficients recovers exactly the Bravetti equations of motion

```
q̇^a = ∂H/∂p_a,   ṗ_a = −∂H/∂q^a − p_a ∂H/∂s,   ṡ = p_a ∂H/∂p_a − H,
```

which is why the GA `X_H` matches `ContactHamiltonian.vectorField()` to floating
point (see §4). The module also exposes `verifyHamiltonianVectorField`, which
recomputes both defining equations with GA operations and returns residuals
(`α(X_H)+H` and `‖ι_{X_H}dα − (dH−R(H)α)‖`), both `~0`.

---

## 3. Rotor sign-convention resolution (G3)

**Question.** Hestenes writes `R = exp(−θB/2)`; the codebase's
`Algebra.rotor(B, θ)` returns `cos(θ/2) + sin(θ/2) B̂ = exp(+½ θ B̂)`. Are they
consistent?

**Resolution (verified numerically).** They are the same rotation up to the
**orientation of the unit bivector** `B̂` — equivalently, they rotate by opposite
signs of the angle. Rotating `e₁` in the `e₁e₂` plane by `+90°`:

| Construction | Formula | `sandwich(e₁)` |
|---|---|---|
| `Algebra.rotor(e₁₂, π/2)` | `exp(+½·(π/2)·e₁₂)` | `−e₂` (clockwise) |
| Hestenes `exp(−(π/2)e₁₂/2)` | `exp(−½·(π/2)·e₁₂)` | `+e₂` (counter-clockwise) |

So `Algebra.rotor(B, θ) = exp(−θ(−B̂)/2)`: it agrees with Hestenes after flipping
the orientation of `B̂` (or negating `θ`). To get the standard active CCW
rotation (`e₁ ↦ e₂` for `+θ`), either call `Algebra.rotor(B, −θ)`, negate the
bivector, or use `Algebra.exp(B̂.scale(−θ/2))` directly.

**Convention adopted for Phase 1 code.** The new `Bivector3D.rotor(angle)`
helper uses the **Hestenes convention** `R = exp(−(angle/2) B̂)`, so
`Bivector3D(0,0,1).rotor(π/2)` sends `e₁ ↦ +e₂` (tested). Any code mixing the two
rotor sources must pin this orientation.

### 3.1 Closed-form bivector exponentials (G3)

`Bivector3D` and `Bivector4D` previously had no `exp`/rotor. Per "prefer reuse",
they now **delegate to the generic `Algebra.exp`** via a `toMultivector(algebra)`
conversion, rather than reimplementing the case-split exponential:

- `Bivector3D.exp()` — 3D bivectors are always *simple* (`B²` is a negative
  scalar), so `Algebra.exp` is exact. Also adds `rotor(angle)` (Hestenes conv.).
- `Bivector4D.exp()` — general 4D bivectors are **not** simple (`B∧B ≠ 0`), where
  the scalar-`B²` closed form does not apply. `exp()` tests `isSimple()`
  (`B∧B = 0`): simple → `Algebra.exp` (closed form); non-simple → scaling-and-
  squaring of the Taylor series (`expSeriesMultivector`, built only from the
  generic geometric product/addition). Verified: `exp(e₁₂+e₃₄)` has scalar part
  `cos²(1)` and is a unit rotor.

---

## 4. Numerical validation results

From `tests/contact/ga_contact_form.js` (37 assertions, all pass):

- **α / dα**: `α` has the correct `+1·ds`, `−p·dq` coefficients; `dα` is pure
  grade-2 and coordinate-independent.
- **Genuine non-degeneracy**: `|α ∧ (dα)^n|` computed via GA wedge equals `n!`
  for `n = 1, 2, 3` (1, 2, 6) — and equals the conventional
  `verifyContactCondition` value, confirming the old stub was numerically right
  but is now *actually computed*.
- **Negative/degeneracy test**: dropping one symplectic pair from `dα` makes
  `α ∧ (dα)^n = 0` — the check genuinely detects degeneracy (it is not a
  constant).
- **Reeb**: `α(R) = 1`, `ι_R dα = 0` at nontrivial points.
- **X_H vs conventional**: GA `X_H` matches `ContactHamiltonian.vectorField()`
  to **max diff 0.0** for four Hamiltonians, including two with explicit
  `s`-dependence (exercising the dissipative `−p_a ∂H/∂s` term). Defining-
  equation residuals are `0`.
- **Flow invariance**: `α ∧ (dα)^n ≠ 0` at every point of a 20-step contact
  flow.
- **G3 rotors**: convention table above verified; `exp(B)` unit-rotor property
  `R R̃ = 1` holds for 3D and both simple/non-simple 4D bivectors.

Full pre-existing suite is unaffected: `npm test` → **79 passed / 0 failed**;
`node tests/run_all.js` → **7 suites passed** (the Contact suite now also runs
`ga_contact_form.js`). `eslint` on the new/changed files: **0 errors**.

---

## 5. Gap status after Phase 1

| Gap | Status | Notes |
|---|---|---|
| **G1** — GA object for `α`, `dα`, real non-degeneracy | ✅ **addressed** | `GAContactForm.contactForm/contactFormDerivative/nonDegeneracy`; genuine `α∧(dα)^n` wedge with positive+negative tests. |
| **G2** — bridge GA ↔ contact layers | ✅ **addressed** | `src/contact/ga-contact-form.js` imports `src/algebra/multivector.js`; exported via `src/contact/index.js`. `R` and `X_H` built and verified with GA ops, cross-checked against `hamiltonian.js`. |
| **G3** — closed-form rotor on `Bivector3D/4D` + rotor sign | ✅ **addressed** (as needed) | `exp()`/`rotor()` delegate to generic `Algebra.exp`; non-simple 4D handled by series; rotor-sign convention resolved and documented (§3). |
| **G4** — AdS extended-phase-space thermodynamics | ❌ deferred | Out of Phase-1 scope; G1 (done) is its prerequisite. |
| **G5** — horizon / Born-rule-breakdown coupling | ❌ deferred | Downstream of G1 and G4. |

### Remaining sub-gaps for Phase 2

1. **Large manifolds.** `GAContactForm` uses the generic `Algebra`, whose
   product-table precompute is `O(4^n)`. It is capped at `MAX_ALGEBRA_SIZE =
   4096` (dim ≤ 11, `n ≤ 5`), so the **Holographic M₇** (`n=3`) is fully
   supported but the **Grand M₁₃** (`n=6`, 8192 blades) throws a clear error. A
   sparse exterior-algebra backend (wedge/contraction over sorted blade lists,
   no full table) would lift this — a natural Phase-2 task.
2. **Non-canonical contact forms.** The construction assumes Darboux coordinates
   (`α = ds − p_a dq^a`). Contact transformations / general `α` with a nontrivial
   conformal factor are not yet represented as GA objects.
3. **Metric `dα` for boosts.** If a later phase needs rotors that mix `q`, `p`,
   `s` (treating `dα` as a metric bivector), the rejected `Cl(n,n,1)` route (§1.3)
   should be revisited.
