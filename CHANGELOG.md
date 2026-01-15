# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2025-01-15

### Added

#### Core Framework
- `ContactManifold` base class for 1-jet bundles J¹(Q)
- `ContactPoint` class for points on contact manifolds
- Canonical contact form α = du - p_a dx^a
- Reeb vector field computation
- Contact non-degeneracy verification

#### Contact Manifold Models
- `GrandContactManifold` (M₁₃): Full 13-dimensional extended phase space
  - Base coordinates: (q¹, q², q³, t, ℓ, S)
  - Conjugate momenta: (k₁, k₂, k₃, ω, Δ, T)
- `HolographicContactManifold` (M₇): Reduced 7-dimensional model with emergent space
  - Base coordinates: (t, ℓ, S)
  - Emergent spatial fields: q^i(t, ℓ, S)
- `GaugeExtendedManifold` (M₁₅): Grand model with gauge pair (φ, I)

#### Hamiltonian Dynamics
- `ContactHamiltonian` class for contact Hamiltonian systems
- Vector field computation (X_H)
- RK4 integration for trajectory evolution
- Hamilton's equations in contact form:
  - ẋ^a = ∂H/∂p_a
  - ṗ_a = -∂H/∂x^a - p_a · ∂H/∂u
  - u̇ = p_a · ∂H/∂p_a - H

#### Specialized Hamiltonians
- `ThermodynamicHamiltonian` with factory methods:
  - `dispersionRelation()`: H = ω - c|k| or H = ω - √(c²|k|² + m²)
  - `equationOfState()`: Ideal gas and Van der Waals

#### Legendrian Submanifolds
- `LegendrianSubmanifold` class for generating function construction
- Automatic lift from base to contact manifold
- Hamilton-Jacobi residual computation
- Sampling methods

#### Gravitational Extension
- `SpacetimeMetric` class for curved spacetime
- Built-in metrics:
  - Minkowski (flat)
  - Schwarzschild (black hole)
  - FLRW (cosmology)
- `RelativisticHamiltonian` for mass-shell constraint
- Geodesic integration

#### Documentation
- Comprehensive README with mathematical foundations
- Full API reference (docs/API.md)
- Mathematical theory guide (docs/THEORY.md)
- Five tutorials covering all major concepts
- Interactive browser demo

#### Developer Experience
- TypeScript type definitions
- 44 validation tests
- Four runnable examples
- GitHub Actions CI workflow

### Technical Notes

- Signature convention: (+, -, -, -) for spacetime metrics
- Coordinates follow physics conventions (ℓ = log λ for scale)
- RK4 integrator for numerical evolution
- Numerical differentiation with h = 10⁻⁷ default step

## [Unreleased]

### Planned
- Additional spacetime metrics (Kerr, Reissner-Nordström)
- Symplectic integrators for long-time evolution
- Browser bundle (webpack/rollup)
- Contact transformations
- Non-equilibrium thermodynamics extensions

---

## Version History Summary

| Version | Date | Highlights |
|---------|------|------------|
| 1.0.0 | 2025-01-15 | Initial release with full framework |
