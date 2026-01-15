# Tutorial 5: Gravitational Extension

This tutorial explains how to couple contact geometry with general relativity through spacetime metrics.

## The Key Principle

```
┌─────────────────────────────────────────────┐
│  Kinematics  ←→  Contact structure α        │
│  Curvature   ←→  Metric g_μν inside H       │
└─────────────────────────────────────────────┘
```

The contact structure provides the **canonical phase space scaffold** (locally flat via Darboux theorem). Spacetime **curvature** enters through the Hamiltonian constraint.

## Spacetime Metrics

The library provides several spacetime metrics:

### Minkowski (Flat Spacetime)

```javascript
const CT = require('contact-thermodynamics');

const mink = CT.SpacetimeMetric.minkowski();

// Get metric components
const g = mink.covariant([0, 0, 0, 0]);
// g = diag(+1, -1, -1, -1)
```

Signature convention: (+, -, -, -) (timelike positive).

### Schwarzschild (Black Hole)

```javascript
const M = 1;  // Mass (geometric units: G = c = 1)
const schw = CT.SpacetimeMetric.schwarzschild(M);

// Metric at r = 10M
const g = schw.covariant([0, 10, Math.PI/2, 0]);
// g_tt = 1 - 2M/r = 0.8
// g_rr = -(1 - 2M/r)^{-1} = -1.25
```

Coordinates: (t, r, θ, φ)

Key features:
- Event horizon at r = 2M
- Photon sphere at r = 3M
- Asymptotically flat as r → ∞

### FLRW (Cosmology)

```javascript
// Scale factor a(t)
const H0 = 0.1;  // Hubble parameter
const a = t => Math.exp(H0 * t);

const flrw = CT.SpacetimeMetric.flrw(a, 0);  // k=0 (flat universe)
```

Coordinates: (t, χ, θ, φ) — comoving

The parameter k determines spatial curvature:
- k = +1: closed (spherical)
- k = 0: flat
- k = -1: open (hyperbolic)

## The Relativistic Hamiltonian

The **mass-shell constraint** is:

```
H = ½g^{μν}(p_μ - qA_μ)(p_ν - qA_ν) - ½m² = 0
```

For a free particle (A_μ = 0):

```
H = ½g^{μν}p_μp_ν - ½m² = 0
```

### Creating a Relativistic Hamiltonian

```javascript
const metric = CT.SpacetimeMetric.schwarzschild(1);
const mass = 1;
const relH = new CT.RelativisticHamiltonian(metric, mass);
```

### With Electromagnetic Field

```javascript
// Uniform electric field in x-direction
const E = 0.1;
const gaugePotential = x => [E * x[1], 0, 0, 0];  // A_t = E·x
const charge = 1;

const relH_em = new CT.RelativisticHamiltonian(
    metric, mass, gaugePotential, charge
);
```

## Evaluating the Hamiltonian

```javascript
const x = [0, 10, Math.PI/2, 0];  // (t, r, θ, φ)
const p = [1.05, 0, 0, 0.2];      // (p_t, p_r, p_θ, p_φ)

const H = relH.evaluate(x, p);
console.log('H =', H);  // Should be ≈ 0 if on-shell
```

## Hamilton's Equations

```javascript
const { xDot, pDot } = relH.geodesicEquations(x, p);

// xDot = [ṫ, ṙ, θ̇, φ̇]
// pDot = [ṗ_t, ṗ_r, ṗ_θ, ṗ_φ]
```

## Integrating Geodesics

```javascript
const x0 = [0, 10, Math.PI/2, 0];
const p0 = [1.05, 0, 0, 0.2];  // E = 1.05, L = 0.2

const dtau = 0.1;   // proper time step
const steps = 100;

const geodesic = relH.integrateGeodesic(x0, p0, dtau, steps);

// geodesic[i] = { x: [...], p: [...], tau: i*dtau }
```

## Example: Geodesic in Schwarzschild

```javascript
const CT = require('contact-thermodynamics');

// Create Schwarzschild spacetime
const M = 1;
const schw = CT.SpacetimeMetric.schwarzschild(M);
const relH = new CT.RelativisticHamiltonian(schw, 1);

// Initial conditions for a bound orbit
const r0 = 10;
const E = 0.95;   // Energy < 1 for bound orbit
const L = 4.0;    // Angular momentum

// Compute initial p_r from H = 0
const f = 1 - 2*M/r0;
const pr_sq = (E*E/f - 1 - L*L/(r0*r0)) / f;
const pr = Math.sqrt(Math.max(0, pr_sq));

const x0 = [0, r0, Math.PI/2, 0];
const p0 = [E, pr, 0, L];

// Integrate
const geodesic = relH.integrateGeodesic(x0, p0, 1, 200);

// Check conserved quantities
console.log('E(0) =', p0[0], 'E(T) =', geodesic[200].p[0]);
console.log('L(0) =', p0[3], 'L(T) =', geodesic[200].p[3]);
```

## Conserved Quantities

In Schwarzschild spacetime:
- **E = p_t**: Energy (from time translation symmetry)
- **L = p_φ**: Angular momentum (from axial symmetry)

These remain constant along geodesics.

## Hamilton-Jacobi in Curved Spacetime

The relativistic Hamilton-Jacobi equation is:

```
½g^{μν}(∂_μS)(∂_νS) - ½m² = 0
```

Solutions S(x) are generating functions for Legendrian submanifolds.

For separable problems (like Schwarzschild), we can write:

```
S = -Et + L·φ + S_r(r) + S_θ(θ)
```

## Summary

1. **Contact structure** provides canonical kinematics
2. **Metric g_μν** encodes spacetime curvature inside H
3. **Mass-shell constraint** H = 0 determines geodesics
4. **Conserved quantities** from Killing vectors
5. **Hamilton-Jacobi** connects to Legendrian geometry

## Available Metrics

| Metric | Constructor | Coordinates | Physics |
|--------|-------------|-------------|---------|
| Minkowski | `minkowski()` | (t,x,y,z) | Flat spacetime |
| Schwarzschild | `schwarzschild(M)` | (t,r,θ,φ) | Black hole |
| FLRW | `flrw(a,k)` | (t,χ,θ,φ) | Cosmology |

## Custom Metrics

```javascript
// Create a custom metric
const customMetric = new CT.SpacetimeMetric(
    // Covariant metric g_μν(x)
    x => [
        [f(x), 0, 0, 0],
        [0, -1/f(x), 0, 0],
        [0, 0, -r(x)**2, 0],
        [0, 0, 0, -r(x)**2 * Math.sin(theta(x))**2]
    ],
    // Optional: contravariant metric g^μν(x)
    null  // Will be computed numerically if not provided
);
```

---

This completes the tutorial series! You now have the tools to:

1. Create contact manifolds
2. Choose appropriate models
3. Define and evolve Hamiltonians
4. Construct Legendrian submanifolds
5. Couple to curved spacetime
