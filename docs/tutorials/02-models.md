# Tutorial 2: The Grand vs Holographic Models

This tutorial explains the three scale models in Contact Thermodynamics and when to use each.

## Overview

| Model | Base Dim | Total Dim | Use Case |
|-------|----------|-----------|----------|
| Grand M₁₃ | 6 | 13 | Full phase space with explicit space |
| Holographic M₇ | 3 | 7 | Reduced description, emergent space |
| Gauge-Extended M₁₅ | 7 | 15 | Including gauge degrees of freedom |

## The Grand Model M₁₃

This is the "honest" full phase space where all physical variables are independent coordinates.

### Coordinates

**Base Configuration Q₆:**
```
x^a = (q¹, q², q³, t, ℓ, S)
```

Where:
- q^i = spatial position (i = 1, 2, 3)
- t = time
- ℓ = log(λ) = logarithmic scale factor
- S = entropy

**Conjugate Momenta:**
```
p_a = (k₁, k₂, k₃, ω, Δ, T)
```

Where:
- k_i = wavenumber (conjugate to position)
- ω = frequency (conjugate to time)
- Δ = dilatation/anomalous dimension (conjugate to scale)
- T = temperature (conjugate to entropy)

### Usage

```javascript
const CT = require('contact-thermodynamics');

const M13 = CT.grandManifold();

// Create a physical point
const pt = M13.physicalPoint(
    1, 0, 0,      // q¹, q², q³
    0,            // t
    0,            // ℓ
    1,            // S
    0.5, 0, 0,    // k₁, k₂, k₃
    1,            // ω
    0,            // Δ
    1,            // T
    0             // A (action)
);

console.log('Dimension:', M13.dim);  // 13
console.log('Position:', M13.spatialPosition(pt));  // [1, 0, 0]
console.log('Wave vector:', M13.waveVector(pt));    // [0.5, 0, 0]
```

### Contact Form

```
α = dA - k_i dq^i - ω dt - Δ dℓ - T dS
```

This encodes:
- dA = k·dq (wave mechanics)
- dA = -ω dt (time evolution)
- dA = T dS (thermodynamic identity)

## The Holographic Model M₇

In this reduced description, spatial coordinates become **dependent fields** rather than independent coordinates.

### Motivation

In many situations (thermodynamic limit, coarse-graining, holography), we don't need to track every spatial degree of freedom. Instead, space "emerges" from more fundamental variables.

### Coordinates

**Base Configuration Q₃:**
```
x^a = (t, ℓ, S)
```

**Conjugate Momenta:**
```
p_a = (ω, Δ, T)
```

**Emergent Space:**
```
q^i = q^i(t, ℓ, S)  — scalar fields on Q₃
```

### Usage

```javascript
const M7 = CT.holographicManifold();

// Create a point in the reduced space
const pt = M7.holographicPoint(
    0,    // t
    0.5,  // ℓ = log(scale)
    1,    // S
    1,    // ω
    0,    // Δ
    1,    // T
    0     // A
);

// Define how space emerges
const emergent = M7.createEmergentSpace(pt, (t, ell, S) => {
    const a = Math.exp(ell);  // scale factor
    return [
        a * Math.cos(t),  // q¹
        a * Math.sin(t),  // q²
        0                 // q³
    ];
});

console.log('Emergent position:', emergent);
```

### Physical Interpretation

The holographic model captures situations where:
- The thermodynamic limit has been taken
- Spatial degrees of freedom are coarse-grained
- We're interested in bulk properties, not microscopic details

## The Gauge-Extended Model M₁₅

When gauge symmetries are important, we add an additional canonical pair.

### Additional Coordinates

```
φ = gauge phase
I = gauge flux / current
```

### Usage

```javascript
const M15 = CT.gaugeExtended();

console.log('Dimension:', M15.dim);  // 15
console.log('Base coords:', M15.baseCoords);
// ['q1', 'q2', 'q3', 't', 'ell', 'S', 'phi']
```

### Contact Form

```
α = dA - k_i dq^i - ω dt - Δ dℓ - T dS - I dφ
```

## Choosing the Right Model

| Situation | Model |
|-----------|-------|
| Full wave mechanics in space | Grand M₁₃ |
| Thermodynamic calculations | Holographic M₇ |
| Statistical mechanics | Holographic M₇ |
| Particle physics / QFT | Grand or Gauge-Extended |
| Gauge theories | Gauge-Extended M₁₅ |

## Dimensional Counting

The formula **dim(M) = 2n + 1** always holds:

- Grand: 2×6 + 1 = 13
- Holographic: 2×3 + 1 = 7
- Gauge-Extended: 2×7 + 1 = 15

## Next Steps

In [Tutorial 3](03-dynamics.md), we'll learn how to evolve systems using Contact Hamiltonian dynamics.
