# Contact Thermodynamics

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Node.js](https://img.shields.io/badge/Node.js-18%2B-green.svg)](https://nodejs.org/)
[![Tests](https://img.shields.io/badge/Tests-106%20passing-brightgreen.svg)](#testing)

A modular JavaScript framework for **Contact Geometry**, **Geometric Algebra**, **Discrete Geometric Calculus**, and **Pilot-Wave Theory**.

## 🌟 Features

### Core Modules

| Module | Description | Key Components |
|--------|-------------|----------------|
| **algebra/** | Geometric Algebra | Multivector, Complex, Dual, Hyperbolic numbers |
| **geometry/** | Riemannian Geometry | Connection bivectors, curvature, geodesics |
| **calculus/** | Discrete Calculus | Split ∇ operator, mesh derivatives, PDE solvers |
| **contact/** | Contact Geometry | Manifolds, Hamiltonians, Legendrian submanifolds |
| **physics/** | Physics Applications | Pilot-wave, entropic gravity, spacetime |

### Highlights

- **Contact Manifolds**: 1-jet bundles J¹(Q) with canonical contact form α = du - p_a dx^a
- **Geometric Algebra**: Cl(p,q,r) algebras, rotors, bivector classification
- **Number Systems**: Complex (i²=-1), Dual (ε²=0), Hyperbolic (j²=+1) with auto-diff
- **Riemannian GA**: Connection bivectors ωᵢ, curvature 2-form Ω, coordinate-free geodesics
- **FTGC Meshes**: Triangle mesh calculus with cotan weights, wave/heat/Maxwell solvers
- **Pilot-Wave Theory**: Valentini regularization, H-theorem, curved-space guidance

## 📁 Project Structure

```
src/
├── algebra/                # Core Geometric Algebra
│   ├── index.js
│   ├── multivector.js      # Cl(p,q,r) algebras, Multivector, rotors
│   └── number-systems.js   # Complex, Dual, Hyperbolic numbers
├── geometry/               # Riemannian Geometry
│   ├── index.js
│   ├── riemannian-ga.js    # Connection bivector, curvature
│   ├── riemannian-discrete.js
│   └── geodesic.js         # Geodesic solver, parallel transport
├── calculus/               # Discrete Geometric Calculus
│   ├── index.js
│   ├── grid.js             # Split ∇ on regular grids
│   ├── mesh.js             # Triangle mesh data structure
│   ├── mesh-derivative.js  # MeshGeometricDerivative
│   └── solvers.js          # Leapfrog, wave, heat solvers
├── contact/                # Contact Geometry
│   ├── index.js
│   ├── manifold.js         # ContactManifold, Grand, Holographic
│   ├── hamiltonian.js      # ContactHamiltonian dynamics
│   └── legendrian.js       # Legendrian submanifolds
├── physics/                # Physics Applications
│   ├── index.js
│   ├── pilot-wave.js       # Valentini regularization
│   ├── entropic-gravity.js # Bianconi framework
│   └── spacetime.js        # Spacetime GA
└── index.js                # Main entry (backward compatible)
```

## 🚀 Quick Start

### Installation

```bash
npm install contact-thermodynamics
```

### Usage Patterns

#### Namespaced Imports (Recommended)

```javascript
const CT = require('contact-thermodynamics');

// Algebra module
const { Complex, Dual, Hyperbolic } = CT.algebra;
const z = new Complex(1, 2);
console.log(z.mul(z));  // -3 + 4i

// Automatic differentiation with Dual numbers
const f = x => x.mul(x);  // f(x) = x²
const result = f(new Dual(3, 1));
console.log(result.a, result.b);  // 9, 6 (value, derivative)

// Calculus module
const { SplitDifferentialOperator, ScalarField } = CT.calculus;

// Geometry module
const { Sphere2D, ConnectionBivector, GAGeodesicSolver } = CT.geometry;
const sphere = new Sphere2D(1.0);
const solver = new GAGeodesicSolver(sphere);
```

#### Legacy Flat Imports (Still Works)

```javascript
const CT = require('contact-thermodynamics');

// Create 13D Grand Contact Manifold
const M13 = CT.grandManifold();
const pt = M13.physicalPoint(1, 0, 0, 0, 0, 1, 0.5, 0, 0, 1, 0, 1, 0);

// Geometric Algebra
const { Algebra } = CT;
const sta = new Algebra(1, 3);  // Cl(1,3) Spacetime Algebra
```

### Riemannian Geometry Example

```javascript
const { Sphere2D, Curvature2Form, GAGeodesicSolver } = CT.geometry;

// Create unit sphere
const sphere = new Sphere2D(1.0);

// Compute Gaussian curvature (no Christoffel symbols!)
const curvature = new Curvature2Form(sphere);
const K = curvature.gaussianCurvature([Math.PI/4, 0]);
console.log('K =', K);  // 1.0

// Solve geodesic equation ∇ᵥv = 0
const solver = new GAGeodesicSolver(sphere);
const path = solver.solve([Math.PI/4, 0], [1, 0], Math.PI);
```

### FTGC Mesh Example

```javascript
const { TriangleMesh, MeshGeometricDerivative, LeapfrogGCMesh } = CT.calculus;

// Create triangulated grid
const mesh = TriangleMesh.createGrid(16, 16, 1.0, 1.0);

// Build discrete geometric derivative ∇
const nabla = new MeshGeometricDerivative(mesh);

// Compute Laplacian
const f = new Float64Array(mesh.nVertices);
for (let i = 0; i < mesh.nVertices; i++) {
    f[i] = Math.sin(Math.PI * mesh.vertices[i * 3]);
}
const lapF = nabla.laplacian(f);

// Solve wave equation
const solver = new LeapfrogGCMesh(mesh);
const dt = solver.estimateCFLWave(1.0);
const u = solver.waveSimulate(f, new Float64Array(mesh.nVertices), dt, 100, 1.0);
```

## 📐 Mathematical Foundation

### Contact Form
$$\alpha = du - p_a \, dx^a$$

### Geometric Algebra
$$ab = a \cdot b + a \wedge b$$
$$\nabla = \sum_i e^i \partial_i$$

### Connection Bivector
$$\omega_i = \frac{1}{2} e^j \wedge (\partial_i e_j)$$

### Curvature 2-Form
$$\Omega = d\omega + \omega \wedge \omega$$

## 📊 Interactive Demos

- **[3D EM Wave](examples/em-wave-3d.html)** — Maxwell solver visualization
- **[Riemannian Geodesics](examples/riemannian-ga-demo.html)** — Sphere, Torus, Hyperbolic Plane
- **[Pilot-Wave Demo](examples/pilot-wave-demo.js)** — Quantum relaxation

## 🧪 Testing

```bash
npm test
```

Test coverage: 79 core tests + 27 number systems tests = **106 total**

## 📄 License

MIT License — see [LICENSE](LICENSE) for details.
