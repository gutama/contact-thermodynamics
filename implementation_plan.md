# Comprehensive Test Bed Implementation Plan

A multi-layer test bed for validating the contact-thermodynamics library.

---

## 1. Analytic Geometry Tests (Continuous)

### A. Riemannian GA: Sphere/Torus/Hyperbolic

**Sphere S²(R):**
- Gaussian curvature K = 1/R² at random (θ, φ)
- Great-circle geodesics: verify arc length matches analytical
- Holonomy for latitude loop: angle = 2π(1 − cos θ)

**Torus:**
- Curvature varies with position; integrate numerically
- Sign/region sanity check vs known formula K = cos θ / (r(R + r cos θ))

**Hyperbolic Plane:**
- Constant negative curvature K = -1 (within tolerance)

### B. Parallel Transport / Holonomy Invariants

Use invariant-driven tests, not convention-driven:
- Compare holonomy to surface integral of curvature (Gauss-Bonnet)
- Verify ∮ω = ∫∫K dA for bounded loops

---

## 2. Contact-Thermodynamics Physics Tests

### A. Damped Harmonic Oscillator
- Standard pedagogical contact Hamiltonian example
- Verify energy decay rate matches analytic solution

### B. Relaxation to Equilibrium
- Contact Hamiltonian with known Lyapunov function
- Confirm monotone entropy production
- Verify approach to equilibrium state

### C. Legendre Invariance (if implemented)
- Contact form preserved (up to conformal factor)
- Observables transform consistently

---

## 3. Discrete Geometry + FTGC Mesh Validation

### A. Discrete Laplacian / FTGC Identities

| Test | Expected |
|------|----------|
| Constant field | Δf ≈ 0 (interior) |
| Linear field on plane | Δf ≈ 0 (boundary effects only) |
| Rectangle eigenvalues | Match first few analytical eigenpairs |

### B. Heat Equation (Manufactured Solutions)
- Use manufactured solution u(x,y,t), derive forcing term
- Measure error convergence under mesh refinement
- Track energy/variance decay monotonicity

### C. Maxwell/Wave Solver Sanity
- Divergence constraints: ∇·B ≈ 0, ∇·E ≈ ρ
- CFL stability: solutions bounded under standard IC
- Energy conservation (in conservative limits)

---

## 4. Demo-to-Test Harness

Convert interactive demos to regression tests:

1. Run simulation kernels headless (Node.js) for N steps
2. Compare against golden snapshots with tolerances
3. Keep rendering code separate from correctness tests

**Demos to convert:**
- EM Wave 3D
- Riemannian Geometry Demo
- Mesh Heat Diffusion
- Phase Space Demo

---

## 5. Proposed Directory Layout

```
tests/
├── geometry/
│   ├── sphere_curvature.js       # K = 1/R²
│   ├── sphere_geodesic.js        # Arc length
│   ├── sphere_holonomy.js        # 2π(1 - cos θ)
│   ├── torus_curvature.js        # Variable K
│   └── hyperbolic_curvature.js   # K = -1
├── contact/
│   ├── damped_oscillator.js      # Energy decay
│   ├── equilibrium_relaxation.js # Entropy production
│   └── legendre_invariance.js    # Form preservation
├── mesh/
│   ├── laplacian_patch.js        # Δ(const) = 0
│   ├── laplacian_linear.js       # Δ(linear) = 0
│   ├── heat_manufactured.js      # Convergence
│   └── gauss_bonnet_discrete.js  # Σ K_v = 2πχ
├── pde/
│   ├── wave_invariants.js        # Energy conservation
│   ├── maxwell_divergence.js     # ∇·B = 0
│   └── cfl_stability.js          # Boundedness
└── bench/
    ├── mesh_scaling.js           # Performance vs size
    └── solver_throughput.js      # Steps/sec
```

---

## Priority Order

1. **High**: `tests/geometry/` (validates core GA formulations)
2. **High**: `tests/mesh/gauss_bonnet_discrete.js` (already partial)
3. **Medium**: `tests/mesh/` Laplacian + heat tests
4. **Medium**: `tests/pde/` stability + invariants
5. **Lower**: `tests/contact/` (requires physics integration)
6. **Lower**: `bench/` (performance, not correctness)
