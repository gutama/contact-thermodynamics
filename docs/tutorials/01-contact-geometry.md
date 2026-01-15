# Tutorial 1: Understanding Contact Geometry

This tutorial introduces the fundamental concepts of contact geometry and how they relate to thermodynamics and mechanics.

## What is Contact Geometry?

Contact geometry is the **odd-dimensional sibling** of symplectic geometry. While symplectic geometry describes phase spaces (position + momentum), contact geometry adds one more dimension—typically a potential or action.

### The Contact Condition

A contact manifold (M, α) consists of:
- A (2n+1)-dimensional manifold M
- A 1-form α satisfying the **non-degeneracy condition**:

```
α ∧ (dα)^n ≠ 0  everywhere
```

This means α cannot be integrated to give a codimension-1 foliation—it's "maximally non-integrable."

## The Canonical Example: 1-Jet Bundle

The most important contact manifold for physics is the **1-jet bundle** J¹(Q).

Given a configuration space Q with coordinates (x¹, ..., xⁿ), the 1-jet bundle has:
- Base coordinates: x^a
- Fiber coordinate: u (the "value" of a function)
- Momentum coordinates: p_a (the "derivatives")

The canonical contact form is:

```
α = du - p_a dx^a
```

**Dimension check:** n + 1 + n = 2n + 1 ✓

## Why Contact for Physics?

### Thermodynamics

In thermodynamics, we have:
- Extensive variables (V, N, S, ...)
- Intensive variables (P, μ, T, ...)
- Thermodynamic potentials (U, F, G, H, ...)

The number of variables is always **odd**: n extensive + n intensive + 1 potential = 2n + 1.

The fundamental relation `dU = TdS - PdV + μdN` is exactly the contact form!

### Mechanics

In classical mechanics:
- q^a = positions
- p_a = momenta  
- S = action

Hamilton's principal function S(q, t) satisfies `dS = p_a dq^a - H dt`, which is again the contact form structure.

## Hands-On Example

```javascript
const CT = require('contact-thermodynamics');

// Create a simple contact manifold with 2 base coordinates
// This will be a 5-dimensional manifold (2 + 1 + 2 = 5)
const M = new CT.ContactManifold(
    ['x', 'y'],      // base coordinates
    ['px', 'py'],    // momenta
    'u'              // fiber (potential)
);

console.log('Dimension:', M.dim);  // 5
console.log('Contact form:', M.contactFormSymbolic());
// Output: du - px·dx - py·dy

// Create a point
const pt = M.point({ x: 1, y: 2, px: 0.5, py: -0.3, u: 0 });

// The Reeb field is always ∂/∂u
const R = M.reebField(pt);
console.log('Reeb field:', R);  // { x: 0, y: 0, px: 0, py: 0, u: 1 }

// Verify α(R) = 1
console.log('α(R) =', M.evaluateContactForm(pt, R));  // 1
```

## Key Takeaways

1. Contact manifolds have dimension 2n+1
2. The canonical form is α = du - p_a dx^a
3. The Reeb field R = ∂/∂u satisfies α(R) = 1
4. Contact geometry naturally describes systems with potentials

## Next Steps

In [Tutorial 2](02-models.md), we'll explore the three specific models:
- Grand model M₁₃
- Holographic model M₇
- Gauge-extended model M₁₅
