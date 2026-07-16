# Phase 2 — GA-Native Pilot-Wave Guidance & Node Regularization

**Phase 2 of the research plan** to unify contact-geometric thermodynamics with
geometric algebra (GA). This phase re-expresses the de Broglie–Bohm **guidance
equation** and Valentini's **node-regularization** using GA's rotor/bivector
representation of the wavefunction phase, and asks the concrete question **RQ2**:

> Does expressing the guidance equation and node-regularization in GA/STA
> multivector form (rather than standard vector calculus) produce a more natural
> or more numerically stable regularization near wavefunction nodes, and does it
> extend cleanly toward the relativistic/curved-spacetime case?

Everything is **additive**: the conventional `src/physics/pilot-wave.js`
(`SmearingKernel`, `WaveFunction`, `PilotWaveSystem`, `QuantumEnsemble`,
`CurvedSpacePilotWave`) is untouched. The new code lives in
`src/physics/ga-pilot-wave.js` and is exported through `src/physics/index.js`
alongside the existing module. It **reuses** `SmearingKernel` and
`QuantumEnsemble` verbatim, mirroring how Phase 1 wired GA in alongside the
existing contact-form code.

---

## 1. The GA construction

### 1.1 Wavefunction as an even multivector (spinor)

Following the Hestenes / Doran–Lasenby GA treatment of the Pauli/Schrödinger
equation, the complex phase factor is represented in the plane algebra
`Cl(2,0,0)` with the unit pseudoscalar playing the role of the imaginary unit:

```
i   ↦   I = e₁∧e₂          (I² = −1)
ψ   ↦   ψ_r + I ψ_i  =  R · U,     U = exp(I S)   (a unit rotor)
```

`U = ψ/|ψ|` is the **phase rotor**. This is the same "phase as a rotor/bivector"
idea the task calls for, and it parallels the `ConnectionBivector` pattern
already in `riemannian-ga.js` (a bivector that generates a rotation).

### 1.2 Guidance velocity as a rotor logarithmic derivative

The guidance velocity is the **grade-1 part of the phase rotor's logarithmic
derivative, de-rotated by `I⁻¹`**:

```
v  =  (ħ/m) ⟨ (∇U) U⁻¹ I⁻¹ ⟩₁ ,       U = exp(I S).
```

**Proof of equivalence (flat, non-relativistic).** With `∇ = Σ_k e_k ∂_k` the
ordinary vector derivative,

```
∇U = (∇S) I U        (∇S = Σ_k e_k ∂_k S, an ordinary gradient vector)
(∇U) U⁻¹ = (∇S) I
(∇U) U⁻¹ I⁻¹ = (∇S) I I⁻¹ = ∇S.
```

So `v = (ħ/m) ∇S` **exactly**. Analytically the GA form therefore carries **no
new predictive content** in flat space — it is the identical velocity field.
This is verified to machine precision (≤ 2.2e-16) in test §I for a plane wave, a
quadratic phase, and a phase vortex.

> **Note (grade-1 projection is a no-op in `Cl(2,0,0)`).** In the 2-D plane
> algebra `(∇U)U⁻¹I⁻¹` is *already* pure grade 1 for a unit rotor `U`, so the
> `⟨·⟩₁` projection in the formula above is a **redundant no-op here** — it
> discards nothing. We keep it written explicitly only because it becomes
> load-bearing in higher-dimensional / STA generalizations (where `(∇U)U⁻¹I⁻¹`
> can acquire grade-3 components that *must* be projected out). In this phase it
> is notation for the future case, not an active operation.

The interesting content of RQ2 is not the continuum identity but **how the two
forms discretize it**, which is where they diverge numerically (§3).

### 1.3 Alternative considered and rejected

**The "full-spinor" current** `j = (ħ/m)⟨ψ̃ ∇ψ I⁻¹⟩` using the *un-normalized*
spinor `ψ = ψ_r + I ψ_i` (a literal transcription of `Im(ψ*∇ψ)`).

Rejected because in the 2-D plane algebra a real vector `a` and a
pseudoscalar-times-vector `I·b` are **both grade 1**, so the amplitude-gradient
term `∇R` and the phase-gradient term `ρ∇S` in `ψ̃∇ψ = R∇R + ρ(∇S)I` do **not
separate by a clean grade projection**. Extracting the phase part requires
exactly the reversion/normalization bookkeeping that normalizing to the unit
rotor `U = ψ/|ψ|` performs up front. Normalizing first is simpler, removes the
`∇R` contamination by construction, and — the whole point for RQ2 — differences
the *smooth* rotor rather than the wrapped phase angle.

The un-normalized route offers no compensating advantage for guidance, but it is
the natural object once amplitude dynamics / the quantum potential enter (it
carries `∇R`). That is flagged for a later phase, not needed here. (This mirrors
Phase 1's rejected `Cl(n,n,1)` route, which likewise becomes attractive only for
a later, more demanding use.)

---

## 2. What was built

| Class / fn (`src/physics/ga-pilot-wave.js`) | Role |
|---|---|
| `phaseRotor(S)`, `unitRotorFromComponents(re,im)` | build `U = exp(IS)` |
| `guidanceFromRotor(U, dU, ħ/m)` | core operation `⟨(∇U)U⁻¹I⁻¹⟩₁` |
| `GAPilotWaveSystem1D` | 1-D grid guidance; **same interface** as `PilotWaveSystem` (`regularizedVelocity()`), so it drops into `QuantumEnsemble` unchanged |
| `GAPilotWaveField2D` | 2-D grid; exposes both `deBroglieVelocityGA()` (rotor) and `deBroglieVelocityPhase()` (conventional baseline) plus `regularizedVelocityGA()` |
| `GACurvedPilotWave` | curved guidance `v^i=(ħ/m)g^{ij}∂_jS`; composes with `ConnectionBivector` for phase-bivector parallel transport |

Node regularization is the identical Valentini/Bell construction
`v_reg = (j∗μ)/(ρ∗μ)` with the current `j = ρv` built from the GA rotor velocity
and the **reused** `SmearingKernel` convolution.

---

## 3. Near-node numerical comparison (the RQ2 result)

**Setup.** The phase vortex `ψ ∝ (x+iy) e^{−r²}` on a 41×41 grid over
`[−2,2]²` (`dx≈0.1`). Its exact guidance velocity is the azimuthal field
`v = (ħ/m)(−y, x)/r²`. This is the canonical case with a genuine phase
singularity (a node *and* a `2π` phase winding at the origin).

**Three methods compared.** To avoid a strawman we now compare **three**
discretizations of the *same* guidance velocity:

1. **`∇`atan2** — central-difference the wrapped phase angle
   `S = atan2(ψ_i, ψ_r)`. This is a *naive* implementation, **not** "the
   conventional guidance equation."
2. **Im-current** — the honest standard baseline: `v = (ħ/m) Im(ψ*∇ψ)/|ψ|²`
   central-differenced directly on `Re ψ`, `Im ψ` (never touching `atan2`).
   This is what textbooks and simulation codes actually use.
3. **GA rotor** — the `⟨(∇U)U⁻¹I⁻¹⟩` form of §1.2, differencing the unit rotor
   `U = ψ/|ψ|`.

**Metric.** `max|v|` and `max error vs analytic`, over all grid points
**excluding the physical node ball `r < 0.5`** — so what remains isolates the
*representational* artifact, not the real physical divergence at the node. A
second **far-field** column (`r > 1`) isolates the region where amplitude
gradients, not the phase singularity, dominate.

| Quantity (grid) | `∇`atan2 (naive) | Im-current (standard) | **GA rotor** |
|---|---|---|---|
| `max \|v\|` (`r>0.5`) | **30.9** | 1.98 | 1.98 |
| `max error vs analytic` (`r>0.5`) | **31.4** | 0.023 | 0.039 |
| `max error vs analytic` (far field `r>1`) | (huge) | 0.023 | **0.0050** |
| regularized `max \|v\|` | (—) | (—) | 1.54 (finite through node) |

**Reading — the branch cut is an atan2 artifact, not a property of the standard
equation.** The `∇`atan2 form differences the wrapped phase angle, which has a
`2π` branch cut along the whole negative-`x` axis. Central-differencing across
that cut produces a spurious `≈2π/dx` spike on an entire codimension-1 locus —
hence `max|v| = 30.9`, error `31.4`, even *far from* the node. **But the standard
current form has no branch cut at all**: it differences `Re ψ`, `Im ψ` (smooth
everywhere), so it tracks the analytic velocity to `0.023` and is bounded at
`1.98`, essentially the same as GA. The originally-reported ~3-orders-of-
magnitude "GA beats the standard equation" gap was **entirely the atan2
strawman** — it is not a property of the de Broglie guidance equation, which the
Im-current implements robustly.

**GA's genuine (modest) advantage.** GA and Im-current are the *same order of
magnitude* near the node. GA's real edge is a **constant-factor accuracy gain
away from the singularity**, because it differences the *amplitude-normalized*
rotor `U = ψ/|ψ|`, whose `O(dx²)` truncation error is **independent of `∇R`**,
whereas the Im-current differences the full `ψ` so its error couples to the
amplitude gradient. In the far field this shows as `0.0050` (GA) vs `0.023`
(Im-current) — a real but ~4–5× effect, not a categorical fix (§3.2).

### 3.1 The 1-D real node — GA and Im-current are both *correct*

For a real wavefunction `ψ ∝ x` (a sign-flip node), the true de Broglie velocity
is **zero everywhere** (`Im(∇ψ/ψ)=0`).

| Method | `max \|v\|` (true = 0) |
|---|---|
| `∇`atan2 (naive) | 15.7 raw / **0.22** regularized (spurious spike) |
| Im-current (standard) | **0.0** (exact) |
| GA rotor | **0.0** (exact) |

The `∇`atan2 form encodes the sign flip as a `0→π` phase jump and differencing
it spikes. Both the **Im-current** baseline (`Im ψ ≡ 0 ⇒ v ≡ 0`) and the **GA
rotor** (the jump lands in the internal bivector direction, off the physical
grade-1 axis) return the physically correct `v = 0`. So again: the spike is an
`atan2` artifact, and the honest standard form is already correct — GA matches
it, it does not uniquely fix it.

### 3.2 Honest scope of the improvement

- Against a **strawman** (`∇`atan2) GA looks like a categorical fix; against the
  **honest standard** (Im-current) the branch-cut/real-node pathologies vanish
  and GA is merely *comparable* near singularities. The originally-documented
  "31.4 → 0.039" headline was a strawman comparison and has been corrected.
- GA's **real, defensible** advantage is a **constant-factor** accuracy gain
  (~5× on a 1-D smooth phase with varying amplitude: test §II reports
  `Im/GA ≈ 4.9×` at `dx=0.10`, `5.0×` at `dx=0.05`; ~4–5× in the vortex far
  field). Both methods are `O(dx²)`; GA's error prefactor is smaller because it
  decouples the amplitude gradient. This is a genuine but *modest* numerical
  benefit, not a change in convergence order or a repair of a broken equation.
- Away from any phase singularity the rotor and phase discretizations agree to
  `O(dx²)` and converge to each other as `dx→0` (test §II:
  `dx=0.10→2.9e-2`, `dx=0.05→1.9e-3`).
- The GA rotor-difference has its own `O(dx²)` error (a `sinc(k·dx)` factor for a
  linear phase) and is *not* bit-identical to the standard grid at finite `dx`;
  the *operation* is exact (§1.2, §I), the *grid discretization* is merely
  convergent.

---

## 4. H-theorem convergence (quantum relaxation)

**The original experiment was buggy.** It drove the ensemble with a **static**
wavefunction (a fixed shifted Gaussian, linear phase). A static, simply-connected
phase gives only a constant *drift* velocity: it translates the ensemble but does
not *mix* it, and — critically — the target `|ψ|²` is **not stationary** under
that drift, so there is nothing for the coarse-grained `H` to relax toward. The
result was `H` *increasing* `1.2187 → 1.4194` (both fields), the opposite of the
Valentini H-theorem. That is a broken setup, not a physical anti-relaxation.

**The corrected experiment** uses the canonical relaxation scenario: a
**time-dependent** two-energy-eigenstate superposition in the infinite square
well `[0, L]` (`ℏ = m = 1`),

```
ψ(x,t) = (1/√2)[φ₁ e^{−iE₁t} + φ₂ e^{−iE₂t}],  φ_n = √(2/L) sin(nπx/L),
E_n = n²π²/(2L²).
```

The velocity field `v = j/|ψ|²` now oscillates and has a **blinking interior
node**, so it genuinely stirs the ensemble; `|ψ(x,t)|²` is the equivariant moving
target the H-function is measured against. Starting from a **uniform** (maximally
nonequilibrium) ensemble and advecting with the *same* `QuantumEnsemble` solver
driven by each velocity field (`L=π`, `n=101`, `dt=0.02`, 50 steps, Gaussian
smearing `ε=0.25`, coarse-grain cell 8):

| | initial `H` | final `H` (50 steps) | total drop | max step increase |
|---|---|---|---|---|
| standard (Im-J) | 0.9327 | 0.4504 | **51.7 %** | 9.0e-3 |
| **GA rotor** | 0.9327 | 0.4387 | **53.0 %** | 9.6e-3 |

Now `H` **decreases** monotonically (bar a small bump, below) for **both** fields
— the sign is corrected and relaxation toward quantum equilibrium is recovered.
As expected the two trajectories are nearly identical (the smooth-phase parts of
the field agree to `O(dx²)`); GA relaxes at the same rate as the honest standard
field.

**The residual ~1e-2 bump is honest solver noise, not anti-relaxation.** A single
step occasionally raises `H` by `~9e-3` (< 2.5 % of the ~0.48 total drop). It is
a **node-crossing artifact of the coarse first-order upwind advection**: it
*grows* as `dt→0` (ruling out a time-integration error) and is **identical for
the standard and GA fields** (ruling out any velocity-field cause). It is
spatial-discretization noise from the blinking node sweeping the grid, not a
property of either guidance form. A higher-order (flux-limited) advection scheme
would remove it; that is a solver upgrade orthogonal to RQ2.

---

## 5. Curved-space limit

`GACurvedPilotWave.curvedGuidanceVelocity` extracts the covariant phase gradient
`∂_j S` with the **same** rotor operation as the flat case and raises the index
with the manifold's inverse metric: `v^i = (ħ/m) g^{ij} ∂_j S`. It composes with
the existing `ConnectionBivector`/`parallelTransport` machinery in
`riemannian-ga.js` (phase-bivector transport solves `dB/dλ + ω(v)×B = 0` via
`commutatorWithBivector`, exactly as `CurvedSpacePilotWave` does).

Verified (test §V):

- **Flat reduction:** with `g^{ij}=δ^{ij}` the curved velocity equals the flat GA
  velocity to `0.0` (exact) — the required consistency check.
- **Index raising:** a diagonal `g^{ij}=diag(2, ½)` scales the components by
  exactly `(2, ½)` — the metric genuinely acts.
- **Zero connection:** on a flat identity-frame manifold the connection bivector
  is zero, so phase-bivector parallel transport leaves the bivector unchanged
  (drift `0.0`).

This is the natural bridge toward a future relativistic/curved-spacetime phase
(Phase 4): the guidance operation is already metric-covariant and rotor-based.

---

## 6. Test results

`tests/test_ga_pilot_wave.js` — **38 assertions, all pass**:

- §I GA guidance **operation** reproduces analytic `∇S` to ≤ 2.2e-16 (plane
  wave, quadratic phase, phase vortex).
- §II rotor-difference converges to phase-difference as `dx→0`; on a
  varying-amplitude smooth phase the GA rotor is ~5× more accurate than the
  Im-current baseline (`Im/GA ≈ 4.9×` at `dx=0.10`), both `O(dx²)`.
- §III phase-vortex stability, **3-way**: the `∇`atan2 strawman blows up
  (`max|v|=30.9`, error `31.4`); the honest **Im-current** standard baseline is
  accurate (`0.023`) with **no** branch cut; **GA** is comparable near the node
  (`0.039`) and ~4–5× better in the far field (`0.0050` vs `0.023`). Regularized
  GA finite through the node.
- §IV 1-D real node, **3-way**: both **Im-current** and **GA** give the correct
  `v=0` exactly; only the `∇`atan2 form spikes (reg `0.22`).
- §V curved reduces to flat exactly; index raising and zero-connection transport
  behave.
- §VI H-theorem: corrected time-dependent two-state relaxation — `H` **decreases**
  `0.9327 → 0.45` (~52 %) for **both** the standard and GA fields (was
  erroneously *increasing* before the fix).

Pre-existing suite unaffected: `npm test` → **79 passed / 0 failed**;
`node tests/run_all.js` → **7 suites passed** (the Pilot-Wave suite now also runs
`test_ga_pilot_wave.js`); `eslint` on the new/changed files → **0 errors**.

---

## 7. RQ2 verdict

| Sub-question | Finding |
|---|---|
| More **natural** formulation? | Yes — the guidance velocity is a single rotor logarithmic derivative, and phase-as-rotor reuses the `ConnectionBivector` pattern. But in flat space it is analytically **identical** to `∇S` (no new physics). |
| More **numerically stable** near nodes? | **Modestly, and only vs a strawman.** The dramatic gap (`0.039` vs `31.4`) is against a *naive* `∇`atan2 form and is an `atan2` branch-cut artifact, **not** a flaw of the guidance equation. The honest standard baseline (`Im(ψ*∇ψ)/|ψ|²`, central-differenced on `Re/Im ψ`) has no branch cut and is *comparable* to GA near nodes. GA's genuine benefit is a **constant-factor** (~4–5×) accuracy gain away from singularities, from differencing the amplitude-normalized rotor. |
| Extends cleanly to **curved** space? | Yes — the operation is metric-covariant, reduces exactly to flat when `g=I`/`ω=0`, and composes with the existing connection-bivector transport. Sets up Phase 4. |

**Bottom line, stated without overclaiming:** the GA reformulation adds **no new
predictive content** in the flat continuum (it is provably the same guidance
equation). Its numerical benefit is **not** a repair of a broken standard
equation — the honest current form `Im(ψ*∇ψ)/|ψ|²` is already free of the
branch-cut and real-node artifacts that only afflict a naive `∇`atan2
implementation. GA's defensible contribution is narrower: **a rotor-based
discretization of the non-relativistic guidance equation that avoids wrapped-
phase artifacts *by construction* and carries a modest ~4–5× constant-factor
accuracy edge** (from differencing `ψ/|ψ|`, whose truncation error decouples
from the amplitude gradient), extending unchanged to curved space.

### Relation to prior work (novelty scope)

GA/Clifford formulations of pilot-wave theory are **not new**: the relativistic
Dirac-current guidance in Spacetime Algebra is standard (Hestenes; Doran &
Lasenby), and quantum-relaxation studies for Dirac fermions already exist (e.g.
Colin's work on de Broglie–Bohm relaxation). This phase does **not** claim to
originate GA pilot-wave theory or the STA guidance law. The scoped contribution
here is the concrete, tested observation that a **rotor / unit-spinor
discretization of the *non-relativistic* guidance equation** avoids wrapped-phase
branch-cut artifacts by construction and gives a modest constant-factor accuracy
gain — plus the metric-covariant framing that sets up a future curved/relativistic
phase.

### Deferred / limitations

1. **2-D regularization is Euclidean-grid convolution.** The curved-space
   regularized average (`regularizedVelocityCurved` in the existing module) uses
   a geodesic-ball Monte-Carlo; a GA-native geodesic-ball version composing with
   the connection is left for Phase 4.
2. **Relativistic (STA) case not implemented.** The construction is set up to
   generalize `I = e₁∧e₂` to an STA bivector, but the Dirac-current guidance and
   its regularization are out of Phase-2 scope.
3. **Amplitude dynamics / quantum potential** are not modeled (only guidance);
   these need the un-normalized spinor rejected in §1.3.
