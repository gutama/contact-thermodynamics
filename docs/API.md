# API Reference

Complete API documentation for the Contact Thermodynamics library.

## Table of Contents

- [Module Overview](#module-overview)
- [Algebra Module](#algebra-module)
- [Geometry Module](#geometry-module)
- [Calculus Module](#calculus-module)
- [Contact Module](#contact-module)
- [Physics Module](#physics-module)
- [Legacy Flat API](#legacy-flat-api)

---

## Module Overview

The library is organized into five namespaced modules accessible via the main entry point:

```javascript
const CT = require('contact-thermodynamics');

// Access modules
CT.algebra   // Geometric Algebra, Number Systems
CT.geometry  // Riemannian Geometry, Geodesics
CT.calculus  // Discrete Calculus, Mesh Operations
CT.contact   // Contact Manifolds (extracted)
CT.physics   // Physics Applications (extracted)
```

---

## Algebra Module

`CT.algebra` — Core Geometric Algebra and number systems.

### Algebra

Clifford Algebra Cl(p,q,r) factory.

```javascript
const { Algebra } = CT.algebra;

const sta = new Algebra(1, 3);  // Cl(1,3) Spacetime Algebra
const e0 = sta.e(1);            // e₀ (timelike)
const e1 = sta.e(2);            // e₁ (spacelike)
```

#### Constructor

```javascript
new Algebra(p, q, r = 0)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `p` | number | Positive signature dimensions |
| `q` | number | Negative signature dimensions |
| `r` | number | Zero signature dimensions (dual numbers) |

#### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `e(i)` | Multivector | Basis vector eᵢ |
| `scalar(s)` | Multivector | Scalar multivector |
| `pseudoscalar()` | Multivector | Highest-grade element |
| `rotor(B, angle)` | Multivector | exp(-B·angle/2) |

### Multivector

Multivector operations in Cl(p,q,r).

```javascript
const a = sta.e(1);
const b = sta.e(2);
const c = a.mul(b);        // Geometric product
const d = a.wedge(b);      // Outer product
const e = a.inner(b);      // Inner product
```

#### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `add(other)` | Multivector | Addition |
| `sub(other)` | Multivector | Subtraction |
| `mul(other)` | Multivector | Geometric product |
| `wedge(other)` | Multivector | Outer product ∧ |
| `inner(other)` | Multivector | Inner product · |
| `scale(s)` | Multivector | Scalar multiplication |
| `reverse()` | Multivector | Reversion |
| `grade(k)` | Multivector | Grade-k projection |
| `norm()` | number | Magnitude |

### Complex

Complex numbers (i² = -1) — elliptic type.

```javascript
const { Complex } = CT.algebra;

const z = new Complex(1, 2);  // 1 + 2i
z.add(new Complex(3, -1));    // 4 + 1i
z.mul(new Complex(3, -1));    // 5 + 5i
z.norm();                      // √5
Complex.exp(Math.PI);          // e^(iπ) = -1
```

#### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `add(z)` | Complex | Addition |
| `mul(z)` | Complex | Multiplication |
| `div(z)` | Complex | Division |
| `conjugate()` | Complex | Complex conjugate |
| `norm()` | number | Magnitude |
| `arg()` | number | Phase angle |
| `pow(n)` | Complex | Power |

#### Static Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `Complex.exp(θ)` | Complex | e^(iθ) = cos(θ) + i·sin(θ) |
| `Complex.unitCurve(t)` | Complex | e^(it) (unit circle) |

### Dual

Dual numbers (ε² = 0) — parabolic type, automatic differentiation.

```javascript
const { Dual } = CT.algebra;

// f(x) = x² → compute f(3) and f'(3) simultaneously
const x = new Dual(3, 1);     // x + ε
const fx = x.mul(x);          // x² + 2x·ε
console.log(fx.a, fx.b);      // 9, 6 (value, derivative)

Dual.exp(1);                   // 1 + ε (series terminates!)
```

#### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `add(d)` | Dual | Addition |
| `mul(d)` | Dual | Multiplication (ε² = 0) |
| `div(d)` | Dual | Division |
| `pow(n)` | Dual | Power with derivative |

#### Static Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `Dual.exp(θ)` | Dual | 1 + θε |
| `Dual.differentiate(f, x)` | {value, derivative} | Auto-diff f at x |

### Hyperbolic

Hyperbolic (split-complex) numbers (j² = +1) — hyperbolic type.

```javascript
const { Hyperbolic } = CT.algebra;

const h = new Hyperbolic(3, 1);  // 3 + 1j
h.normSq();                       // 3² - 1² = 8 (Minkowski)
h.isTimelike();                   // true
h.rapidity();                     // arctanh(1/3)

Hyperbolic.exp(1);                // cosh(1) + j·sinh(1)
Hyperbolic.boost(1, 0, 0.5);      // Lorentz boost
Hyperbolic.velocityAdd(0.5, 0.5); // Relativistic: 0.8
```

#### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `add(h)` | Hyperbolic | Addition |
| `mul(h)` | Hyperbolic | Multiplication (j² = +1) |
| `normSq()` | number | Minkowski norm a² - b² |
| `isTimelike()` | boolean | a² - b² > 0 |
| `isSpacelike()` | boolean | a² - b² < 0 |
| `isNull()` | boolean | a² - b² ≈ 0 |
| `rapidity()` | number | arctanh(b/a) |

#### Static Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `Hyperbolic.exp(φ)` | Hyperbolic | cosh(φ) + j·sinh(φ) |
| `Hyperbolic.boost(t, x, φ)` | {t, x} | Lorentz boost |
| `Hyperbolic.velocityAdd(v1, v2)` | number | Relativistic addition |

### Bivector Classification

```javascript
const { BivectorType, classifyBivectorSquare } = CT.algebra;

classifyBivectorSquare(-1);  // 'elliptic' (rotation)
classifyBivectorSquare(0);   // 'parabolic' (translation)
classifyBivectorSquare(1);   // 'hyperbolic' (boost)
```

---

## Geometry Module

`CT.geometry` — Riemannian geometry via Geometric Algebra.

### Sphere2D

Unit 2-sphere with Gaussian curvature K = 1/R².

```javascript
const { Sphere2D } = CT.geometry;
const sphere = new Sphere2D(1.0);  // R = 1
```

### Torus2D

Torus with major radius R, minor radius r.

```javascript
const { Torus2D } = CT.geometry;
const torus = new Torus2D(2.0, 0.5);  // R = 2, r = 0.5
```

### HyperbolicPlane

Hyperbolic plane with constant negative curvature K = -1.

```javascript
const { HyperbolicPlane } = CT.geometry;
const H2 = new HyperbolicPlane();
```

### ConnectionBivector

Compute connection bivectors ωᵢ from frame fields.

```javascript
const { ConnectionBivector } = CT.geometry;
const connection = new ConnectionBivector(manifold);
const omega = connection.compute([θ, φ], 0);  // ω₀ at point
```

### Curvature2Form

Gaussian curvature via Cartan structure equation.

```javascript
const { Curvature2Form } = CT.geometry;
const curvature = new Curvature2Form(manifold);
const K = curvature.gaussianCurvature([Math.PI/4, 0]);
```

### GAGeodesicSolver

Solve geodesic equation ∇ᵥv = 0.

```javascript
const { GAGeodesicSolver } = CT.geometry;
const solver = new GAGeodesicSolver(manifold);
const path = solver.solve(
    [Math.PI/4, 0],  // Initial position
    [1, 0],          // Initial velocity
    Math.PI          // Arc length
);
```

### GAParallelTransport

Parallel transport vectors along curves.

```javascript
const { GAParallelTransport } = CT.geometry;
const transport = new GAParallelTransport(manifold);
const holonomy = transport.holonomyAngle(
    t => [Math.PI/4, 2*Math.PI*t],  // Loop
    [1, 0],                          // Initial vector
    200                              // Steps
);
```

---

## Calculus Module

`CT.calculus` — Discrete geometric calculus on grids and meshes.

### ScalarField

Scalar field on a regular grid.

```javascript
const { ScalarField } = CT.calculus;
const field = new ScalarField([64, 64], 0.1);
field.set(32, 32, 1.0);
```

### VectorField

Vector field on a regular grid.

```javascript
const { VectorField } = CT.calculus;
const vfield = new VectorField([64, 64], 2, 0.1);
```

### SplitDifferentialOperator

Discrete ∇ = Σᵢ eⁱ ⊗ ∂/∂xᵢ on grids.

```javascript
const { SplitDifferentialOperator } = CT.calculus;
const nabla = new SplitDifferentialOperator(2, 0.1);
const gradF = nabla.grad(scalarField);
const divV = nabla.div(vectorField);
const lapF = nabla.laplacian(scalarField);
```

### TriangleMesh

Triangle mesh data structure with automatic topology.

```javascript
const { TriangleMesh } = CT.calculus;

// Create grid mesh
const mesh = TriangleMesh.createGrid(16, 16, 1.0, 1.0);

// Create disk mesh
const disk = TriangleMesh.createDisk(1.0, 32);

// Create closed icosahedron
const ico = TriangleMesh.createIcosahedron(1.0);
```

#### Properties

| Property | Type | Description |
|----------|------|-------------|
| `nVertices` | number | Number of vertices |
| `nFaces` | number | Number of faces |
| `nEdges` | number | Number of edges |
| `vertices` | Float64Array | Vertex positions [x,y,z,...] |
| `faces` | Uint32Array | Face indices [v0,v1,v2,...] |

### MeshGeometricDerivative

Discrete ∇ on triangle meshes (FTGC).

```javascript
const { MeshGeometricDerivative } = CT.calculus;
const nabla = new MeshGeometricDerivative(mesh);

const gradF = nabla.grad(scalarField);      // Scalar → Edge vectors
const divV = nabla.div(vectorField);        // Edge → Scalar
const lapF = nabla.laplacian(scalarField);  // ∇²f = ∇·∇f
```

### LeapfrogGCMesh

Leapfrog time-stepping for PDEs on meshes.

```javascript
const { LeapfrogGCMesh } = CT.calculus;
const solver = new LeapfrogGCMesh(mesh);

// Wave equation
const dt = solver.estimateCFLWave(1.0);
const u = solver.waveSimulate(u0, v0, dt, 100, 1.0);

// Heat equation
const dtHeat = solver.estimateCFLHeat(0.1);
const uHeat = solver.heatSimulate(u0, dtHeat, 100, 0.1);
```

---

## Contact Module

`CT.contact` — Contact geometry (extracted from main index).

### ContactManifold

Base class for contact manifolds.

```javascript
const { ContactManifold, ContactPoint } = CT.contact;
const M = new ContactManifold(['x', 'y'], ['px', 'py'], 'u');
const pt = M.point({ x: 1, y: 0, px: 0.5, py: 0, u: 0 });
```

### GrandContactManifold

13-dimensional grand model M₁₃ = J¹(Q₆).

```javascript
const { GrandContactManifold } = CT.contact;
const M13 = new GrandContactManifold();
const pt = M13.physicalPoint(1, 0, 0, 0, 0, 1, 0.5, 0, 0, 1, 0, 1, 0);
```

### HolographicContactManifold

7-dimensional holographic model M₇ = J¹(Q₃) with emergent space.

```javascript
const { HolographicContactManifold } = CT.contact;
const M7 = new HolographicContactManifold();
const pt = M7.holographicPoint(0, 0, 1, 1, 0, 1, 0);
```

### ContactHamiltonian

Contact Hamiltonian vector fields and dynamics.

```javascript
const { ContactHamiltonian } = CT.contact;
const H = coords => coords.omega - Math.sqrt(
    coords.k1**2 + coords.k2**2 + coords.k3**2
);
const hamiltonian = new ContactHamiltonian(M13, H);
const trajectory = hamiltonian.flow(pt, 0.1, 100);
```

### LegendrianSubmanifold

Legendrian submanifolds and Hamilton-Jacobi theory.

```javascript
const { LegendrianSubmanifold } = CT.contact;
const A = x => 0.5 * (x.q1**2 + x.q2**2);  // Generating function
const L = new LegendrianSubmanifold(M13, A);
const liftedPt = L.lift({ q1: 1, q2: 0, q3: 0, t: 0, ell: 0, S: 1 });
```

---

## Physics Module

`CT.physics` — Physics applications.

### PilotWaveSystem

Pilot-wave dynamics with Valentini regularization.

```javascript
const PilotWave = CT.physics;
// See physics/pilot-wave.js for full API
```

### EntropicGravity

Bianconi two-metric entropic gravity.

```javascript
const EntropicGravity = CT.physics;
// See physics/entropic-gravity.js for full API
```

---

## Legacy Flat API

All exports are also available directly on the main `CT` object for backward compatibility:

```javascript
const CT = require('contact-thermodynamics');

// These all still work:
const M13 = new CT.GrandContactManifold();
const mesh = CT.TriangleMesh.createGrid(16, 16, 1.0, 1.0);
const sphere = new CT.Sphere2D(1.0);
const nabla = new CT.MeshGeometricDerivative(mesh);
```

---

## Constants

| Constant | Value | Description |
|----------|-------|-------------|
| `CT.EPSILON` | 1e-12 | Floating-point tolerance |

---

## Factory Functions

| Function | Returns | Description |
|----------|---------|-------------|
| `CT.grandManifold()` | GrandContactManifold | Create M₁₃ |
| `CT.holographicManifold()` | HolographicContactManifold | Create M₇ |
| `CT.gaugeExtended()` | GaugeExtendedManifold | Create M₁₅ |
