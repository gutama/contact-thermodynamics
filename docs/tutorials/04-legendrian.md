# Tutorial 4: Legendrian Submanifolds and Hamilton-Jacobi

This tutorial explains how generating functions define Legendrian submanifolds and connect to Hamilton-Jacobi theory.

## What is a Legendrian Submanifold?

A **Legendrian submanifold** L ⊂ M is an n-dimensional submanifold where the contact form vanishes:

```
α|_L = 0
```

This is the maximum dimension for which this can hold (half the contact dimension).

## Construction via Generating Functions

The most important construction uses a **generating function** A(x):

```
u = A(x)              (value of the potential)
p_a = ∂A/∂x^a         (momenta = gradients)
```

This automatically satisfies the Legendrian condition:
```
α|_L = dA - (∂A/∂x^a)dx^a = 0  ✓
```

## Physical Interpretation

| Context | Generating Function | Interpretation |
|---------|---------------------|----------------|
| Mechanics | S(q,t) | Hamilton's principal function |
| Optics | φ(x) | Eikonal (optical path length) |
| Thermodynamics | F(T,V) or G(T,P) | Thermodynamic potential |
| Waves | Ψ(x) | Phase of wave function |

## Creating a Legendrian

```javascript
const CT = require('contact-thermodynamics');
const M13 = CT.grandManifold();

// Plane wave: A = k₀·q - ω₀·t
const k0 = [1, 0, 0];
const omega0 = 1;

const genFunc = x => {
    return k0[0]*x.q1 + k0[1]*x.q2 + k0[2]*x.q3 - omega0*x.t;
};

const L = new CT.LegendrianSubmanifold(M13, genFunc);
```

## Lifting Points to the Contact Manifold

The `lift` method takes base coordinates and produces a point on the Legendrian:

```javascript
const basePoint = { q1: 2, q2: 0, q3: 0, t: 1, ell: 0, S: 0 };
const liftedPoint = L.lift(basePoint);

// Now liftedPoint has:
// - A = A(base) = k₀·q - ω₀·t = 2 - 1 = 1
// - k₁ = ∂A/∂q¹ = k₀[0] = 1
// - ω = -∂A/∂t = -(-ω₀) = 1 (sign depends on convention)
```

## Hamilton-Jacobi Theory

A function A(x) generates a Legendrian that **solves** the dynamics if:

```
H(x, A(x), ∂A(x)) = 0
```

This is the **Hamilton-Jacobi equation**.

### Example: Dispersion Relation

For H = ω - c|k| = 0:
- The HJ equation is: -∂A/∂t - c|∇A| = 0
- Solution: A = k·q - c|k|·t (plane wave)

```javascript
const H_disp = CT.ThermodynamicHamiltonian.dispersionRelation(M13, 1, 0);

// Check HJ residual
const residual = L.hamiltonJacobiResidual(basePoint, H_disp);
console.log('HJ residual:', residual);  // Should be ≈ 0
```

## Spherical Waves

```javascript
// Inward spherical wave: A = -|q|
const sphericalGen = x => {
    const r2 = x.q1**2 + x.q2**2 + x.q3**2;
    return -Math.sqrt(r2 + 0.01);  // regularized at origin
};

const L_spherical = new CT.LegendrianSubmanifold(M13, sphericalGen);

// The momenta point radially inward
const pt = L_spherical.lift({ q1: 1, q2: 0, q3: 0, t: 0, ell: 0, S: 0 });
console.log('k₁ =', pt.get('k1'));  // ≈ -1 (radial, inward)
```

## Geometric Optics (Eikonal)

In geometric optics, the eikonal equation is:

```
|∇S|² = n²(x)
```

where n(x) is the refractive index.

```javascript
// Constant refractive index
const n = 1.5;

// Hamiltonian: H = ½(|k|² - n²)
const H_eikonal = new CT.ContactHamiltonian(M13, coords => {
    return 0.5 * (coords.k1**2 + coords.k2**2 + coords.k3**2 - n*n);
});

// Plane wave in medium: A = n·x
const eikonalGen = x => n * x.q1;
const L_eikonal = new CT.LegendrianSubmanifold(M13, eikonalGen);

// Verify HJ equation is satisfied
const res = L_eikonal.hamiltonJacobiResidual(
    { q1: 1, q2: 0, q3: 0, t: 0, ell: 0, S: 0 },
    H_eikonal
);
console.log('Eikonal HJ residual:', res);  // ≈ 0
```

## Thermodynamic Potentials

In thermodynamics, different potentials (U, F, G, H) define Legendrian submanifolds related by Legendre transforms.

```
Internal Energy:  U(S,V) → T = ∂U/∂S, P = -∂U/∂V
Helmholtz:        F(T,V) → S = -∂F/∂T, P = -∂F/∂V
Gibbs:            G(T,P) → S = -∂G/∂T, V = ∂G/∂P
```

The contact transformation between potentials preserves the contact structure.

## Sampling a Legendrian

```javascript
// Random samples on the spherical wave
const sampler = () => ({
    q1: (Math.random() - 0.5) * 2,
    q2: (Math.random() - 0.5) * 2,
    q3: (Math.random() - 0.5) * 2,
    t: 0, ell: 0, S: 0
});

const samples = L_spherical.sample(sampler, 100);

// All samples lie on the Legendrian
for (const pt of samples) {
    // Verify α|_L = 0 (implicitly true by construction)
}
```

## Summary

1. **Legendrian submanifolds** are where α = 0 (maximal non-integrability)
2. **Generating functions** A(x) automatically produce Legendrians
3. **Hamilton-Jacobi** connects Legendrians to dynamics
4. **Physical interpretation** includes action, eikonal, and thermodynamic potentials

## Next Steps

In [Tutorial 5](05-gravity.md), we'll explore the gravitational extension with curved spacetime metrics.
