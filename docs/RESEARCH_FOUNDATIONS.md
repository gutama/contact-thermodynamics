# Notation & Foundations

**Phase 0 — Foundations & literature consolidation.**

This document is a single shared reference that reconciles the mathematical
notation of five load-bearing papers with the actual variable / class / method
names used in this codebase. It exists so that later phases (GA-native contact
form, G₂ remnant thermodynamics, Valentini regularization) can be implemented
against a consistent vocabulary.

> **Scope note.** Phase 0 is *documentation only*. Nothing in the source tree is
> changed by this document. Where a paper concept has no code equivalent yet,
> that gap is flagged explicitly (see Parts 3 and 4).

> **Verification caveat.** This sandbox has no live web access. The paper
> summaries below rely on the task briefing plus standard knowledge of these
> well-established formalisms. Claims that should be re-checked against the
> primary source before being treated as authoritative are marked
> `TODO(verify)`.

---

## Part 1 — Close-reading summaries of the five key papers

### 1. Bravetti, Cruz, Tapias — *Contact Hamiltonian Mechanics* (arXiv:1604.08266, 2016)

A contact manifold is a pair `(M, α)` where `M` has odd dimension `2n+1` and the
contact 1-form `α` satisfies the non-degeneracy condition `α ∧ (dα)^n ≠ 0`. In
Darboux coordinates `(q^a, p_a, s)` the canonical form is

```
α = ds − p_a dq^a.
```

The **Reeb vector field** `R` is the unique field defined by `α(R) = 1` and
`ι_R dα = 0`; for the canonical form `R = ∂/∂s`.

For a **contact Hamiltonian** `H(q, p, s)`, the contact Hamiltonian vector field
`X_H` is defined by

```
α(X_H) = −H,     ι_{X_H} dα = dH − R(H) α.
```

This yields the equations of motion

```
dq^a/dt = ∂H/∂p_a,
dp_a/dt = −∂H/∂q^a − p_a ∂H/∂s,
ds/dt   = p_a ∂H/∂p_a − H.
```

The distinctive feature versus symplectic mechanics is the `−p_a ∂H/∂s` term in
`ṗ` and the `−H` term in `ṡ`: the flow does **not** preserve a volume form. The
divergence of `X_H` is `−(n+1) ∂H/∂s` `TODO(verify: exact prefactor/sign
convention)`, so `∂H/∂s ≠ 0` drives phase-space contraction/expansion. This
non-conservation is exactly what encodes dissipation / entropy production, which
is why contact geometry is the natural home for thermodynamics rather than the
volume-preserving symplectic setting.

### 2. de León, Lainz — *A review on contact Hamiltonian and Lagrangian systems* (arXiv:2011.05579)

A modern, comprehensive treatment of contact geometry for mechanics. It develops
the same `(M, α, R, X_H)` structure as Bravetti et al. but organizes it around
the musical isomorphism `♭ : TM → T*M`, `v ↦ ι_v dα + α(v) α`, which sends
`X_H ↦ dH − (R(H) + H) α` `TODO(verify: precise form of the flat map image)`.

Key objects for us:

- **Legendrian submanifolds**: `n`-dimensional submanifolds `L ⊂ M^{2n+1}` on
  which `α|_L = 0` (maximal integral manifolds of the contact distribution
  `ker α`). These are the geometric encoding of *equations of state* /
  equilibrium surfaces in thermodynamics: a generating function `S(q)` produces
  `L = { p_a = ∂S/∂q^a, s = S(q) }`.
- **Legendre / contact symmetries**: diffeomorphisms preserving `ker α` (i.e.
  `φ*α = f α` for some non-vanishing `f`). These generalize canonical
  transformations and include the thermodynamic Legendre transforms between
  ensembles.
- Applications surveyed: dissipative mechanics, thermodynamics, cosmology, and
  control theory. The review is the standard entry point for the "modern"
  contact-geometry vocabulary used throughout the rest of this project.

### 3. Ghosh & Bhamidipati — *Contact geometry and thermodynamics of black holes in AdS spacetimes* (Phys. Rev. D, 2019)

Applies the Bravetti et al. contact-Hamiltonian machinery to **AdS black hole
thermodynamics** in the extended phase space. The thermodynamic phase space is
built from the black hole variables — mass `M` (interpreted as enthalpy in the
extended framework), entropy `S`, temperature `T`, pressure
`P = −Λ/(8π) = 3/(8π ℓ²)` (from the cosmological constant / AdS radius `ℓ`), and
thermodynamic volume `V` — treated as coordinates on a contact manifold.

- The **first law** `dM = T dS + V dP + …` (plus `Φ dQ`, `Ω dJ` for charge and
  angular momentum) is the statement that physical states lie on a **Legendrian
  submanifold** `α|_L = 0`, with `α = dM − T dS − V dP − …`.
- Concrete **equations of state** (e.g. from the RN-AdS or Kerr-AdS metric) fix
  `T(S,P)` and generate `L`. For RN-AdS this reproduces the `P–V` isotherms with
  a Van-der-Waals-like structure.
- **P–V criticality** and the small/large black hole phase transition appear as
  the breakdown / bifurcation of the Legendrian section under a Legendre
  symmetry; the critical point is where `∂P/∂V = ∂²P/∂V² = 0`
  `TODO(verify: exact criticality conditions and which paper Table they cite)`.

The takeaway for the codebase: black hole thermodynamics is "just" a Legendrian
submanifold of a contact manifold whose coordinates are the extended-phase-space
thermodynamic variables, with equations of state as generating data.

### 4. Hestenes — *Spacetime physics with geometric algebra* (Am. J. Phys., 2003)

Introduces **Spacetime Algebra (STA)** = the Clifford algebra `Cl(1,3)`
generated by an orthonormal frame `{γ_μ}`, `μ = 0,1,2,3`, with

```
γ_μ · γ_ν = η_μν = diag(+, −, −, −).
```

The **geometric product** of two vectors decomposes into symmetric and
antisymmetric parts:

```
ab = a·b + a∧b,   a·b = ½(ab + ba),   a∧b = ½(ab − ba).
```

- **Bivectors** (grade-2 elements such as `γ_1 γ_2`) form a 6-dimensional space
  spanning the Lorentz Lie algebra. A **rotor** `R = exp(−θB/2)` (with `B` a unit
  bivector) implements rotations and boosts; a vector transforms by the
  **sandwich** map `a ↦ R a R̃`, where `R̃` is the reverse of `R`. Rotors form a
  group under the geometric product (the spin group), and `R R̃ = 1`.
- The **unit pseudoscalar** `I = γ_0 γ_1 γ_2 γ_3` satisfies `I² = −1` and
  implements Hodge duality by multiplication (`*A = A I⁻¹` up to convention);
  it also supplies the algebra's complex structure without an imaginary unit.
- The formalism rewrites Dirac, Maxwell, and gravitation coordinate-free. Maxwell
  becomes the single equation `∇F = J`, with `F = E + I B` a bivector field.

STA is the target formalism for the GA-native rewrites in later phases of this
project.

### 5. Valentini — *Pilot-wave theory and the search for new physics* (arXiv:2411.10782, 2024)

Works in de Broglie–Bohm pilot-wave theory. Writing `ψ = |ψ| e^{iS/ħ}`, the
**de Broglie guidance equation** gives the configuration-space velocity

```
v = (ħ/m) ∇S  =  (ħ/m) Im(∇ψ/ψ)  =  j / |ψ|²,
```

where `j = (ħ/m)|ψ|² ∇S` is the quantum probability current. Because `S` is
ill-defined where `ψ = 0`, `v` **diverges at the nodes** of the wavefunction.

Two threads relevant to this project:

1. **Node regularization.** Following Bell's 1987 smearing proposal, one smooths
   the velocity field near nodes — e.g. via a regularized velocity
   `v_reg = (j ∗ μ)/(|ψ|² ∗ μ)` with a normalized kernel `μ` of width `ε`. This
   keeps trajectories finite and produces **short-distance corrections to
   Born-rule statistics**: the equilibrium distribution becomes the regularized
   `(|ψ|²)_reg` rather than `|ψ|²`, with deviations at `O(ε²)` near former nodes.
   Relaxation to (regularized) quantum equilibrium is governed by a subquantum
   H-theorem, `H̄ = ∫ ρ̄ ln(ρ̄/|ψ̄|²)`, `dH̄/dt ≤ 0`.
2. **Horizon physics.** Valentini argues that a **breakdown of the Born rule near
   black hole horizons** (in the "deep quantum gravity regime") could imprint
   observable anomalies on Hawking radiation — quantum non-equilibrium as a
   potential new-physics signal. `TODO(verify: exact claims and any quantitative
   estimates in §... of the 2024 paper)`.

---

## Part 2 — Codebase notation audit

All paths are relative to the repository root. Note the module layout: several
files under `src/*.js` are thin re-export shims for the "real" implementations
under `src/<area>/*.js` (verified by `tests/run_all.js`, which asserts the shims
re-export identical objects). Canonical locations are used below.

### 2.1 Contact geometry — `src/contact/`

#### `ContactManifold` — `src/contact/manifold.js`

The base contact-manifold class. Its coordinate convention differs from the
papers' `(q, p, s)`:

| Role | Paper | Constructor arg / field |
|---|---|---|
| base coordinates `q^a` | `q^a` | `baseCoords` (field `this.baseCoords`) |
| conjugate momenta `p_a` | `p_a` | `momentaCoords` (field `this.momentaCoords`) |
| fiber / action coordinate `s` | `s` | `fiberCoord` (default `'A'`) |

- `n` = `baseCoords.length` (= `dim Q`); `dim` = `2n + 1` (also `get dimension`).
- `get allCoords` — returns the coordinate names in canonical order
  `[...baseCoords, fiberCoord, ...momentaCoords]` (backed by private
  `_allCoords`).
- `hasCoord(coord)` — membership test against the private `_coordSet`.
- `point(coords)` → `ContactPoint`; `get origin` → all-zero point.
- `contactFormSymbolic()` — returns a **string** rendering of
  `α = d(fiber) − p_a d(base)` (e.g. `"dA -p1·dq1 ..."`). There is **no**
  Multivector / differential-form object; the contact form is symbolic only.
- `evaluateContactForm(pt, tangent)` — numerically evaluates `α(v) = du(v) −
  p_a dx^a(v)` for a tangent vector supplied as a `{coord: value}` dictionary.
- `verifyContactCondition(pt)` — returns `n!` (a **stub**: the true
  `α ∧ (dα)^n` wedge is not computed; it just reports the known non-degeneracy
  constant for the canonical form).
- `reebField(pt)` — returns the Reeb field `R = ∂/∂(fiber)` as a
  `{coord: value}` dictionary with a `1` in the fiber slot.

Subclasses: `GrandContactManifold` (the 13-D `M₁₃`, base
`(q1,q2,q3,t,ell,S)`, momenta `(k1,k2,k3,omega,Delta,T)`, fiber `A`) and
`HolographicContactManifold` (the 7-D `M₇`, base `(t,ell,S)`, momenta
`(omega,Delta,T)`). Both are defined in the same file. `GaugeExtendedManifold`
(in `src/index.js`) extends `GrandContactManifold` with a `(φ, I)` gauge pair.

#### `ContactPoint` — `src/contact/manifold.js`

- Fields: `manifold`, `coords` (a `{coord: value}` dictionary initialized to `0`
  for every coordinate of the manifold).
- `get(coord)` / `set(coord, value)` (set is guarded by `hasCoord`).
- `clone()`.
- `add(tangent, dt = 1)` — returns a new point advanced by `tangent * dt`
  (the primitive used by the flow integrator).

#### `ContactHamiltonian` — `src/contact/hamiltonian.js`

This is where the **contact Hamiltonian vector field** is computed.

- Constructor: `(manifold, H, dH = null)` — `H(coords) → number`, optional
  analytic gradient `dH`.
- `evaluate(pt)` = `H(pt.coords)`.
- `gradient(pt, h = 1e-7)` — analytic `dH` if provided, else central-difference
  numerical gradient over all coordinates (returns a `{coord: ∂H}` dictionary).
- `reebComponent(pt)` — `R(H) = ∂H/∂(fiber)` (the `RH` term).
- **`vectorField(pt)`** — the `X_H` implementation. Returns a `{coord: value}`
  dictionary with exactly the Bravetti equations of motion:
  - `X[x^a] = ∂H/∂p_a`
  - `X[p_a] = −∂H/∂x^a − p_a · RH`
  - `X[fiber] = Σ_a p_a ∂H/∂p_a − H`
- `flow(pt, dt, steps)` — RK4 integration of `X_H`, returns the trajectory as an
  array of `ContactPoint`s.

`ThermodynamicHamiltonian` (in `src/index.js`) extends `ContactHamiltonian` with
a concrete `H` (kinetic `½m|k|²`, frequency `−ω`, dilatation, thermal `−T·S`).

#### `LegendrianSubmanifold` — `src/contact/legendrian.js`

- Constructor: `(manifold, generatingFunc, dA = null)` — `generatingFunc` is the
  generating function `A(x)` (stored as `this.A`).
- `gradient(x, h)` — `∂_a A` (analytic or central difference).
- `lift(x)` — the Legendrian embedding: `fiber = A(x)`, `p_a = ∂_a A(x)`; returns
  a `ContactPoint`. This is the equation-of-state map.
- `verifyLegendrianCondition(x)` — returns `true` (automatic for the generating-
  function construction; not a numerical check).
- `hamiltonJacobiResidual(x, hamiltonian)` — evaluates `H` at the lifted point
  (the Hamilton–Jacobi constraint `H(x, A, ∂A) = 0`).

### 2.2 Geometric algebra — `src/algebra/multivector.js`

Exports the object `GA` (also mirrored by the shim `src/multivector.js`).

#### `Algebra` — Clifford algebra `Cl(p, q, r)`

- Constructor `(p, q, r = 0, basisNames = null)`; fields `p, q, r`, `n = p+q+r`,
  `size = 2^n`, `signature` (array of `+1 / −1 / 0`), `basisNames`.
- Precomputes `productSign[a][b]` and `productResult[a][b]` over all blade
  bitmaps (`_precomputeProducts` / `_computeProduct`).
- Blades are represented as **bitmaps** over basis vectors; grade = popcount.
- Factories: `scalar(s)`, `zero()`, `e(i)` (1-indexed basis vector),
  `vector(components)`, `bivector(components)`, `pseudoscalar()`.
- Products: `gp(a,b)` = **geometric product**, `op(a,b)` = **outer/wedge**,
  `ip(a,b)` = **inner (left contraction)**, `dual/undual`, `vee` (regressive).
- **`rotor(bivector, angle)`** — builds a rotor. Convention used:
  `R = cos(θ/2) + sin(θ/2) B̂` for the elliptic case (`B² < 0`),
  `cosh + sinh B̂` for hyperbolic (`B² > 0`), `1 + (θ/2)B` for parabolic
  (`B² = 0`). Note the sign: this is `exp(+½ θ B̂)`, whereas Hestenes writes
  `R = exp(−θB/2)`; the two differ by the orientation of `B` (see Part 3).
- **`exp(bivector)`** — case-split exponential (`cos/sin`, `cosh/sinh`, or
  `1+B`) keyed on `B²`.
- `log(rotor)` — inverse of `exp`.

#### `Multivector`

- Sparse `coeffs` keyed by blade bitmap; `algebra` back-reference.
- `grade(g)`, `scalar()` (grade-0 part), `reverse()` (reversion `Ã`),
  `involute()` (grade involution), `conjugate()` (Clifford conjugate).
- Arithmetic: `add`, `sub`, `scale`.
- Products: `mul` (= `gp`), `wedge` (= `op`), `dot` (= `ip`), `vee`,
  `commutator` `A×B = ½(AB−BA)`, `anticommutator` `½(AB+BA)`, and
  **`sandwich(other)` = `this * other * ~this`** — i.e. the rotor action
  `R a R̃`.
- `dual/undual`, `normSq`, `norm`, `normalized`, `inverse` (versor inverse),
  `isZero`.

#### `ContactAlgebra` factory + helpers

- `ContactAlgebra.create(baseDim, type)` — `euclidean` / `spacetime` / `pga` /
  `cga`.
- `ContactAlgebra.spacetime()` — **`Cl(1,3)`** with basis names `e0..e3` and
  signature `(+,−,−,−)`. This is the STA algebra used by
  `src/physics/spacetime.js`.
- `ContactAlgebra.euclidean3D()` — `Cl(3,0,0)`.
- Top-level shortcuts: `cl3()`, `pga3()`, `cga3()`, `sta()` (= `Cl(1,3)`).
- `classifyBivector(B)` — returns `{ type: 'elliptic'|'hyperbolic'|'parabolic',
  Bsq, description }` keyed on `B²`.

### 2.3 Riemannian GA — `src/geometry/riemannian-ga.js`

Exports `RiemannianGA`. This module supplies the concrete low-dimensional
bivector classes used by the pilot-wave connection.

- **`Bivector2D`** — single component `e12`; `commutatorWithVector(v)` = 90°
  rotation.
- **`Bivector3D`** — components `(e23, e31, e12)`. Methods: `fromArray`,
  `fromWedge(u,v)`, `toArray`, `add/sub/scale/neg`, `normSquared`, `norm`,
  `commutatorWithVector(v)` (= cross product with the dual axial vector),
  `commutatorWithBivector(other)` (Lie bracket via `cross3`),
  `innerWithBivector(other)`. **No `exp` / rotor method.**
- **`Bivector4D`** — six components `(e12,e13,e14,e23,e24,e34)`; same method set,
  with an explicit `so(4)` Lie bracket in `commutatorWithBivector`.
- `TangentFrame` — frame `{e_i}` + reciprocal `{e^i}`, `metric()`,
  `metricInverse()`, `lower`/`raise`.
- `RiemannianManifold` — abstract; `frame(coords)`, `metric`, `metricInverse`,
  `frameDerivatives`.
- **`ConnectionBivector`** — `computeAt(coords)` implements
  `ω_i = ½ e^j ∧ (∂_i e_j)` (returns an array of `Bivector3D`/`Bivector4D`),
  `along(coords, u)` = `u^i ω_i`. Replaces Christoffel symbols.
- `Curvature2Form` — `compute(coords)` implements `Ω = dω + ω∧ω`
  (`Ω_ij = ∂_i ω_j − ∂_j ω_i + ω_i × ω_j`); plus Gaussian/scalar curvature.
- `GACovariantDerivative` — `∇_u A = ∂_u A + ω(u) × A`.
- Concrete manifolds: `Sphere2D`, `Torus2D`, `HyperbolicPlane`, `Sphere3D`,
  `Torus3D`, `HyperbolicSpace3D`.

### 2.4 Pilot-wave dynamics — `src/physics/pilot-wave.js`

Exports `PilotWave`.

- **`SmearingKernel`** — the Valentini/Bell **node-regularization kernel** `μ`.
  Constructor `(width = 1e-10, type = 'gaussian')`; `width` may be a function
  `ε(t)` (time-dependent regularization). Types `gaussian`/`uniform`/
  `lorentzian`. `evaluate(r)` = `μ(r)`; `convolve1D(field, dx)` /
  `convolve2D(...)` perform normalized convolution `(f ∗ μ)`.
- **`WaveFunction`** — `ψ = R·exp(iS)`. Fields `amplitude` (`R`), `phase` (`S`,
  scalar / array / multivector), `mass`, `hbar`, `dimension`. `probabilityDensity`
  = `|ψ|² = R²`. `current(dx)` = `j = (ħ/m) R² ∇S`. Factories `fromCartesian`,
  `planeWave`, `gaussianPacket`.
- **`PilotWaveSystem`** — the guidance dynamics.
  - `deBroglieVelocity(position)` = **`v = (ħ/m) ∇S`** (comment explicitly warns
    it diverges at nodes).
  - `regularizedVelocity()` = **`v_reg = j_reg / ρ_reg = (j∗μ)/(|ψ|²∗μ)`** — the
    finite-at-nodes velocity.
  - `regularizedDensity()` = `(|ψ|²)_reg`; `evolveTrajectory(x0, dt, nSteps)`.
- **`QuantumEnsemble`** — tracks `ρ` vs `|ψ|²`. `hFunction(cellSize)` /
  `hFunctionFine()` = the subquantum **H-function** `∫ ρ̄ ln(ρ̄/|ψ̄|²)`;
  `isInEquilibrium`, `isInRegularizedEquilibrium`, `getRegularizedEquilibrium`,
  `evolve` (continuity equation), `relaxation`.
- **`ActionPhaseBridge`** — connects pilot-wave phase to contact geometry:
  `phaseToAction` (`A = ħS`), `momentum` (`p = ∇S`), `legendrianViolation`
  (checks `dA − p_a dx^a = 0`), `hamiltonJacobi`.
- **`CurvedSpacePilotWave`** — curved-space guidance
  `v^i = (ħ/m) g^{ij} ∂_j S` (`curvedDeBroglieVelocity`); lazily loads a
  `ConnectionBivector` (`get connection`); `parallelTransport` solves
  `dS/dλ + ω(v) × S = 0` using `Bivector*.commutatorWithBivector`;
  `regularizedVelocityCurved` uses a geodesic-ball Monte-Carlo average.

### 2.5 Spacetime GA (GR) — `src/physics/spacetime.js`

Exports `RiemannianSpacetime` with the single class **`SpacetimeManifoldGA`**.
Built on `GA.ContactAlgebra.spacetime()` = `Cl(1,3)`.

- `computeTetrad(x)` — vierbein `E^a_u` (forms) and inverse `e_a^u` (vectors).
- `frame(x)` — tetrad components as the local Lorentz frame `e_a ≈ γ_a`.
- `computeConnection(x)` — Ricci rotation coefficients `γ_abc`.
- `connectionBivector(x)` — the **spin connection** `ω_a = ½ γ_abc e^b∧e^c` as an
  array of `Multivector` bivectors.
- `curvature2Form(x)` — `Ω_ab = ∂_a ω_b − ∂_b ω_a + [ω_a, ω_b] − …`
  (Cartan structure equation, using `Multivector.commutator`).
- `riemannTensor`, `ricciTensor`, `ricciScalar`, `gaussianCurvature`,
  `curvatureNormSquared`.

### 2.6 G₂ black-hole remnant — `src/physics/g2-black-hole-remnant.js`

Relevant to the (asymptotically flat) black-hole thermodynamics side. Classes:
`G2Manifold`, `EinsteinCartanG2`, `KaluzaKleinReduction`, **`BlackHoleRemnant`**,
`G2RicciFlow`. `BlackHoleRemnant` provides `schwarzschildRadius`,
`hawkingTemperature` (`T_H = ħc³/(8πGMk_B)`), `hawkingLuminosity`, `massLossRate`,
`remnantMass/Radius`, `evaporation` simulation, `qubits`. This is
**Schwarzschild-type** (flat asymptotics), **not** the AdS extended-phase-space
setup of Ghosh & Bhamidipati (see gap G4).

---

## Part 3 — Notation correspondence table

Legend: ✅ implemented and wired in; ⚠️ present but partial / not wired into the
relevant module; ❌ no codebase equivalent yet (gap).

### Contact geometry (Bravetti; de León–Lainz; Ghosh–Bhamidipati)

| Paper concept | Paper notation | Codebase identifier | File | Status |
|---|---|---|---|---|
| Contact manifold `(M, α)` | `M^{2n+1}` | `ContactManifold` (`dim = 2n+1`) | `src/contact/manifold.js` | ✅ |
| Base coordinates | `q^a` | `baseCoords` | `src/contact/manifold.js` | ✅ |
| Conjugate momenta | `p_a` | `momentaCoords` | `src/contact/manifold.js` | ✅ |
| Fiber / action coordinate | `s` | `fiberCoord` (default `'A'`) | `src/contact/manifold.js` | ✅ (renamed `s→A/u`) |
| Contact 1-form | `α = ds − p_a dq^a` | `contactFormSymbolic()` (string only); `evaluateContactForm(pt, tangent)` | `src/contact/manifold.js` | ⚠️ symbolic/scalar; no form object |
| Non-degeneracy | `α ∧ (dα)^n ≠ 0` | `verifyContactCondition()` (returns `n!` stub) | `src/contact/manifold.js` | ⚠️ stubbed, not computed |
| Reeb field | `R`, `α(R)=1` | `reebField(pt)` (dict `∂/∂A`) | `src/contact/manifold.js` | ✅ (as dict) |
| Contact Hamiltonian | `H(q,p,s)` | `ContactHamiltonian` (`H`, `evaluate`) | `src/contact/hamiltonian.js` | ✅ |
| Reeb derivative | `R(H) = ∂H/∂s` | `reebComponent(pt)` | `src/contact/hamiltonian.js` | ✅ |
| Contact Hamiltonian vector field | `X_H` | `vectorField(pt)` | `src/contact/hamiltonian.js` | ✅ |
| EOM `q̇` | `∂H/∂p_a` | `X[x^a] = grad[p_a]` | `src/contact/hamiltonian.js` | ✅ |
| EOM `ṗ` | `−∂H/∂q^a − p_a ∂H/∂s` | `X[p_a] = −grad[x^a] − p·RH` | `src/contact/hamiltonian.js` | ✅ |
| EOM `ṡ` | `p_a ∂H/∂p_a − H` | `X[fiber] = Σ p·grad[p] − H` | `src/contact/hamiltonian.js` | ✅ |
| Contact flow | `exp(t X_H)` | `flow(pt, dt, steps)` (RK4) | `src/contact/hamiltonian.js` | ✅ |
| Legendrian submanifold | `L`, `α\|_L = 0` | `LegendrianSubmanifold` | `src/contact/legendrian.js` | ✅ |
| Generating function / EoS | `S(q)`, `p=∂S`, `s=S` | `generatingFunc` / `lift(x)` | `src/contact/legendrian.js` | ✅ |
| Hamilton–Jacobi constraint | `H(q, S, ∂S) = 0` | `hamiltonJacobiResidual(x, H)` | `src/contact/legendrian.js` | ✅ |
| BH first law as Legendrian | `dM = TdS + VdP + …` | `LegendrianSubmanifold` (generic) | `src/contact/legendrian.js` | ⚠️ generic only |
| Extended-phase-space pressure/volume | `P = −Λ/8π`, `V` | — | — | ❌ (gap G4) |
| P–V criticality / phase transition | `∂P/∂V = ∂²P/∂V² = 0` | — | — | ❌ (gap G4) |

### Geometric algebra / STA (Hestenes)

| Paper concept | Paper notation | Codebase identifier | File | Status |
|---|---|---|---|---|
| Spacetime algebra | `Cl(1,3)`, `{γ_μ}` | `ContactAlgebra.spacetime()` / `sta()` (basis `e0..e3`) | `src/algebra/multivector.js` | ✅ |
| Metric signature | `η = (+,−,−,−)` | `Algebra(1,3,0)` `signature` | `src/algebra/multivector.js` | ✅ |
| Frame vector | `γ_μ` | `Algebra.e(i)` (1-indexed) | `src/algebra/multivector.js` | ✅ |
| Geometric product | `ab = a·b + a∧b` | `gp` / `Multivector.mul` | `src/algebra/multivector.js` | ✅ |
| Inner product | `a·b` | `ip` / `Multivector.dot` (left contraction) | `src/algebra/multivector.js` | ✅ |
| Outer product | `a∧b` | `op` / `Multivector.wedge` | `src/algebra/multivector.js` | ✅ |
| Bivector | `γ_1γ_2` (grade 2) | `e(2).wedge(e(3))`, `Algebra.bivector([...])` | `src/algebra/multivector.js` | ✅ |
| Rotor / bivector exponential | `R = exp(−θB/2)` | `Algebra.rotor(B, θ)` (uses `+½θB̂`), `Algebra.exp(B)` | `src/algebra/multivector.js` | ⚠️ present; sign/orientation convention differs (see note) |
| Rotor application | `a ↦ R a R̃` | `Multivector.sandwich(a)` | `src/algebra/multivector.js` | ✅ |
| Reverse | `R̃` | `Multivector.reverse()` | `src/algebra/multivector.js` | ✅ |
| Pseudoscalar | `I = γ_0γ_1γ_2γ_3`, `I²=−1` | `Algebra.pseudoscalar()` | `src/algebra/multivector.js` | ✅ |
| Hodge duality | `*A = A I⁻¹` | `dual` / `undual` / `Multivector.dual()` | `src/algebra/multivector.js` | ✅ |
| Bivector exponential (low-dim classes) | `exp(B)` | — (`Bivector3D`/`4D` have only `commutator*`) | `src/geometry/riemannian-ga.js` | ❌ no closed-form rotor on `Bivector3D/4D` (gap G3) |
| Spin connection | `ω_μ` (bivector) | `SpacetimeManifoldGA.connectionBivector(x)` | `src/physics/spacetime.js` | ✅ |
| Curvature 2-form | `Ω = dω + ω∧ω` | `SpacetimeManifoldGA.curvature2Form(x)`; `Curvature2Form.compute` | `src/physics/spacetime.js`, `src/geometry/riemannian-ga.js` | ✅ |

> **Rotor sign-convention note.** Hestenes writes `R = exp(−θB/2)`. The
> codebase's `Algebra.rotor(B, θ)` returns `cos(θ/2) + sin(θ/2) B̂`, i.e.
> `exp(+½ θ B̂)`. These agree once the orientation of the unit bivector `B̂` is
> fixed consistently (`−B` in one convention = `+B` in the other). Any Phase-1
> code that mixes the two must pin down a single orientation convention.
> `TODO(verify: confirm the intended active-vs-passive rotation convention
> against Hestenes §... before relying on rotor signs.)`

### Pilot-wave (Valentini)

| Paper concept | Paper notation | Codebase identifier | File | Status |
|---|---|---|---|---|
| Wavefunction | `ψ = \|ψ\| e^{iS/ħ}` | `WaveFunction` (`amplitude` `R`, `phase` `S`) | `src/physics/pilot-wave.js` | ✅ |
| Probability density | `\|ψ\|² = R²` | `WaveFunction.probabilityDensity` | `src/physics/pilot-wave.js` | ✅ |
| Probability current | `j = (ħ/m)\|ψ\|²∇S` | `WaveFunction.current(dx)` | `src/physics/pilot-wave.js` | ✅ |
| De Broglie guidance velocity | `v = (ħ/m)∇S` | `PilotWaveSystem.deBroglieVelocity(pos)` | `src/physics/pilot-wave.js` | ✅ |
| Node smearing / regularization | Bell kernel `μ`, width `ε` | `SmearingKernel` | `src/physics/pilot-wave.js` | ✅ |
| Regularized velocity | `v_reg = (j∗μ)/(\|ψ\|²∗μ)` | `PilotWaveSystem.regularizedVelocity()` | `src/physics/pilot-wave.js` | ✅ |
| Regularized equilibrium | `(\|ψ\|²)_reg` | `regularizedDensity()`, `QuantumEnsemble.getRegularizedEquilibrium()` | `src/physics/pilot-wave.js` | ✅ |
| Subquantum H-function | `H̄ = ∫ρ̄ ln(ρ̄/\|ψ̄\|²)` | `QuantumEnsemble.hFunction()` / `hFunctionFine()` | `src/physics/pilot-wave.js` | ✅ |
| Curved-space guidance | `v^i = (ħ/m)g^{ij}∂_j S` | `CurvedSpacePilotWave.curvedDeBroglieVelocity()` | `src/physics/pilot-wave.js` | ✅ |
| Spinor phase transport | `dS/dλ + ω(v)×S = 0` | `CurvedSpacePilotWave.parallelTransport()` | `src/physics/pilot-wave.js` | ✅ |
| Phase ↔ action bridge | `A = ħS`, `p = ∇S` | `ActionPhaseBridge` | `src/physics/pilot-wave.js` | ✅ |
| Born-rule breakdown near horizon | (§ horizon) | — | — | ❌ (gap G5) |

---

## Part 4 — Identified gaps for Phase 1 (GA-native contact form)

Phase 1's goal is a **GA-native representation of the contact 1-form and the
contact Hamiltonian vector field** — i.e. expressing `α`, `dα`, `R`, and `X_H` as
first-class geometric-algebra objects rather than the current
dictionary/string-based scalars. The audit surfaces the following concrete gaps.

**G1 — No GA object for the contact 1-form `α` or its exterior derivative `dα`.**
`ContactManifold` represents `α` only as (a) a display string
(`contactFormSymbolic()`) and (b) a scalar evaluator on a `{coord: value}`
tangent dictionary (`evaluateContactForm`). There is no `Multivector`- or
differential-form-valued `α`, and `dα` is never constructed. `verifyContactCondition`
returns the constant `n!` instead of computing the wedge `α ∧ (dα)^n`. Phase 1
needs: a form/multivector encoding of `α = ds − p_a dq^a`, a real `dα`, and a
genuine non-degeneracy check.

**G2 — The GA layer and the contact layer are disconnected.** `src/algebra/`
(`Algebra`, `Multivector`, `ContactAlgebra`) is never imported by
`src/contact/`. The `ContactAlgebra` factory advertises contact-oriented signatures
but is unused. A GA-native contact form will need either a degenerate Clifford /
Heisenberg-type algebra (e.g. `Cl(n, n, 1)` with one null direction for the Reeb
axis) or an explicit graded module wired into `ContactManifold`. This bridge does
not exist yet and is the core Phase-1 construction.

**G3 — No closed-form bivector exponential / rotor on the low-dimensional
bivector classes.** `Bivector3D` / `Bivector4D` (`src/geometry/riemannian-ga.js`)
implement only `commutatorWithVector` / `commutatorWithBivector`. The general
`Algebra.rotor` / `Algebra.exp` exist in `src/algebra/multivector.js` but operate
on generic `Multivector`s, not on these classes, and are not used for parallel
transport (which currently uses first-order commutator stepping in
`CurvedSpacePilotWave.parallelTransport`). If Phase 1 wants exact rotor-based
transport of the Reeb frame or of spinor phases, a `Bivector*.exp()` (or a
conversion into the `Algebra`/`Multivector` world) must be added. Also resolve the
rotor **sign/orientation convention** (Hestenes `exp(−θB/2)` vs code
`cos(θ/2)+sin(θ/2)B̂`) before mixing the two.

**G4 — No AdS extended-phase-space thermodynamics (Ghosh–Bhamidipati).** The
codebase has generic Legendrian submanifolds and a Schwarzschild-type
`BlackHoleRemnant` (flat asymptotics, `T_H = ħc³/8πGMk_B`), but no
pressure/volume `P = −Λ/8π`, `V` pair, no RN-AdS/Kerr-AdS equation of state, and
no `P–V` criticality machinery. Encoding the AdS first law as a Legendrian
section is a later-phase task but is noted here because the contact-form
representation (G1) is its prerequisite.

**G5 — No horizon coupling for Born-rule breakdown (Valentini §horizon).** The
regularized pilot-wave machinery (`SmearingKernel`, `QuantumEnsemble`) and
curved-space guidance (`CurvedSpacePilotWave`) exist, but nothing ties the
regularization width `ε` or the equilibrium deviation to a horizon / Hawking-
radiation background. This is downstream of both G1 and G4.

**Minimal Phase-1 starting point.** The critical path is G1 + G2: introduce a GA
carrier for `α`, `dα`, `R`, and `X_H` inside `src/contact/`, backed by the
existing `Algebra`/`Multivector` machinery, with the canonical form
`α = d(fiber) − p_a d(base)` and a real non-degeneracy check replacing the
current `n!` stub. G3 (rotor exponentials on concrete bivectors) is the natural
companion once the Reeb axis is represented geometrically.

---

## Appendix — sanity check

Existing test suites were run during this audit (no source code was modified):

- `npm test` (`node tests/test.js`) → **79 passed, 0 failed**.
- `node tests/run_all.js` → **7 suites passed, 0 failed** (including the shim
  re-export identity checks that confirm `src/*.js` mirror `src/<area>/*.js`).
