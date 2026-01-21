# Tutorial 7: Riemannian Geometry via Geometric Algebra

This tutorial covers coordinate-free Riemannian geometry using Geometric Algebra (GA).

## Overview

Traditional Riemannian geometry uses Christoffel symbols Γᵏᵢⱼ — 3-index tensors that are cumbersome and coordinate-dependent. In Geometric Calculus, we replace them with:

- **Connection bivector ω** — encodes parallel transport as rotation
- **Curvature 2-form Ω** — emerges naturally as dω + ω∧ω

This approach is:
- **Coordinate-free**: No index gymnastics
- **Geometrically intuitive**: Rotations instead of symbols
- **Dimension-agnostic**: Same formulas work in 2D, 3D, 4D

## Key Concepts

### 1. Tangent Frames

At each point on a manifold, we have a **tangent frame** {e₁, e₂, ...} — a basis for the tangent space.

```javascript
const RiemannianGA = require('./src/riemannian-ga.js');
const { Sphere2D, TangentFrame } = RiemannianGA;

// Create a sphere of radius 1
const sphere = new Sphere2D(1.0);

// Get tangent frame at (θ = π/4, φ = π/3)
const coords = [Math.PI/4, Math.PI/3];
const frame = sphere.frame(coords);

console.log(`Frame dimension: ${frame.dim}`);
console.log(`e₁ (∂θ): ${frame.frame[0]}`);
console.log(`e₂ (∂φ): ${frame.frame[1]}`);
```

### 2. The Metric

The metric tensor g captures distances and angles:

```javascript
// Get metric tensor gᵢⱼ = eᵢ · eⱼ
const g = frame.metric();
console.log('Metric tensor:');
console.log(`  g₁₁ = ${g[0][0].toFixed(4)}`);  // |e₁|²
console.log(`  g₁₂ = ${g[0][1].toFixed(4)}`);  // e₁ · e₂
console.log(`  g₂₂ = ${g[1][1].toFixed(4)}`);  // |e₂|²
```

### 3. Connection Bivector

The **connection bivector ω** replaces Christoffel symbols. It describes how the frame rotates as you move along the manifold:

```
ωᵢ = ½ eʲ ∧ (∂ᵢeⱼ)
```

In GA, parallel transport of a vector v along direction u is:

```
∇ᵤv = u · ∂v + ω(u) × v
```

```javascript
const { ConnectionBivector } = RiemannianGA;

const connection = new ConnectionBivector(sphere);
const omegas = connection.computeAt(coords);

console.log(`Connection bivectors at point:`);
console.log(`  ω_θ: ${omegas[0].toString()}`);
console.log(`  ω_φ: ${omegas[1].toString()}`);

// Connection along tangent vector u = ∂θ
const u = [1, 0];  // Pure θ direction
const omega_u = connection.along(coords, u);
console.log(`ω(∂θ): ${omega_u.toString()}`);
```

### 4. Curvature

**Gaussian curvature K** measures how much the surface deviates from flat space.

**Curvature 2-form Ω** is computed as:

```
Ω = dω + ω ∧ ω
```

For 2D surfaces, this reduces to a single number K.

```javascript
const { Curvature2Form } = RiemannianGA;

const curvature = new Curvature2Form(sphere);
const K = curvature.gaussianCurvature(coords);

console.log(`Gaussian curvature K: ${K.toFixed(6)}`);
console.log(`Expected (1/R²): ${(1/1).toFixed(6)}`);
```

## Example: Sphere Curvature

```javascript
const RiemannianGA = require('./src/riemannian-ga.js');
const { Sphere2D, Curvature2Form } = RiemannianGA;

// Sphere of radius R
const R = 2.0;
const sphere = new Sphere2D(R);

// Sample points
const points = [
    [Math.PI/4, 0],        // 45° latitude
    [Math.PI/2, Math.PI],  // Equator
    [Math.PI/3, Math.PI/6] // General point
];

for (const [theta, phi] of points) {
    const curvature = new Curvature2Form(sphere);
    const K = curvature.gaussianCurvature([theta, phi]);
    const expected = 1 / (R * R);
    
    console.log(`(θ=${(theta*180/Math.PI).toFixed(0)}°, φ=${(phi*180/Math.PI).toFixed(0)}°): K = ${K.toFixed(4)}, expected = ${expected.toFixed(4)}`);
}
```

Output:
```
(θ=45°, φ=0°): K = 0.2500, expected = 0.2500
(θ=90°, φ=180°): K = 0.2500, expected = 0.2500
(θ=60°, φ=30°): K = 0.2500, expected = 0.2500
```

## Example: Torus with Variable Curvature

The torus has curvature that varies with position:

```javascript
const { Torus2D, Curvature2Form } = RiemannianGA;

const R = 3.0;  // Major radius
const r = 1.0;  // Minor radius
const torus = new Torus2D(R, r);

// Curvature formula: K = cos(θ) / (r(R + r cos(θ)))
const testAngles = [0, Math.PI/2, Math.PI];

for (const theta of testAngles) {
    const curvature = new Curvature2Form(torus);
    const K = curvature.gaussianCurvature([theta, 0]);
    
    console.log(`θ = ${(theta*180/Math.PI).toFixed(0)}°: K = ${K.toFixed(4)}`);
}
```

Output:
```
θ = 0°: K = 0.2500    (outer, positive curvature)
θ = 90°: K = 0.0000   (top, zero curvature)
θ = 180°: K = -0.5000 (inner, negative curvature)
```

## Example: Hyperbolic Plane

The hyperbolic plane (upper half-plane model) has constant negative curvature:

```javascript
const { HyperbolicPlane, Curvature2Form } = RiemannianGA;

const H2 = new HyperbolicPlane();

// Test at various y values
const points = [[0, 1], [0, 2], [1, 3]];

for (const [x, y] of points) {
    const K = new Curvature2Form(H2).gaussianCurvature([x, y]);
    console.log(`(x=${x}, y=${y}): K = ${K.toFixed(4)}`);
}
```

Output:
```
(x=0, y=1): K = -1.0000
(x=0, y=2): K = -1.0000
(x=1, y=3): K = -1.0000
```

## Available Manifolds

| Manifold | Class | Curvature |
|----------|-------|-----------|
| Sphere S² | `Sphere2D(R)` | K = 1/R² (constant positive) |
| Torus T² | `Torus2D(R, r)` | K varies with position |
| Hyperbolic H² | `HyperbolicPlane()` | K = -1 (constant negative) |
| 3-Sphere S³ | `Sphere3D(R)` | K = 1/R² |
| 3-Torus T³ | `Torus3D(r1, r2, r3)` | K = 0 (flat) |
| Hyperbolic H³ | `HyperbolicSpace3D()` | K = -1 |

## The Gauss-Bonnet Theorem

A deep result connecting curvature to topology:

```
∫∫ K dA = 2π χ(M)
```

Where χ is the Euler characteristic:
- Sphere: χ = 2, so ∫∫K = 4π
- Torus: χ = 0, so ∫∫K = 0

See `tests/geometry/torus_curvature.js` for numerical verification.

## Holonomy

Parallel transporting a vector around a loop rotates it by the enclosed curvature:

```
Holonomy angle = ∫∫ K dA
```

For a latitude loop at colatitude θ on a sphere:

```
Holonomy = 2π(1 - cos θ)
```

See `tests/geometry/sphere_holonomy.js` for details.

## Running the Tests

```bash
# All geometry tests
node tests/run_geometry.js

# Individual tests
node tests/geometry/sphere_curvature.js
node tests/geometry/torus_curvature.js
node tests/geometry/hyperbolic_curvature.js
```

## Further Reading

- Hestenes & Sobczyk, *Clifford Algebra to Geometric Calculus*
- Doran & Lasenby, *Geometric Algebra for Physicists*
- `docs/THEORY.md` — Theoretical background
- `docs/API.md` — Full API reference
