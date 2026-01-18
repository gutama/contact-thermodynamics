# Contact Thermodynamics

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Node.js](https://img.shields.io/badge/Node.js-18%2B-green.svg)](https://nodejs.org/)
[![Tests](https://img.shields.io/badge/Tests-65%20passing-brightgreen.svg)](#testing)

A JavaScript implementation of **Contact Geometry for Extended Thermodynamics** based on 1-jet bundles, now enhanced with **Geometric Algebra** and **Discrete Geometric Calculus**.

<p align="center">
  <img src="docs/assets/contact-manifold-diagram.svg" alt="Contact Manifold Structure" width="600">
</p>

## 🌟 Features

- **Contact Manifolds**: Full implementation of 1-jet bundles J¹(Q) with canonical contact form
- **Geometric Algebra**:
  - **Multivectors**: Proper handling of Cl(p,q,r) algebras, geometric products, and rotors
  - **Bivector Classification**: Elliptic (rotation), Parabolic (translation), Hyperbolic (boost)
  - **Discrete Calculus**: Split differential operator ∇, gradient, divergence, curl
- **Three Model Scales**:
  - **Grand Model (M₁₃)**: 13-dimensional full phase space
  - **Holographic Model (M₇)**: 7-dimensional with emergent space
  - **Gauge-Extended (M₁₅)**: 15-dimensional with gauge degrees of freedom
- **Contact Hamiltonian Dynamics**: Vector fields, Reeb fields, RK4 integration
- **General Relativity Extension**: 
  - **Christoffel Symbols**: Automatic computation from spacetime metrics
  - **Covariant Derivative**: ∇ operator for vectors and tensors
  - **Parallel Transport**: Connection-preserving transport along curves
- **Interactive Visualization**: Browser-based demos including **3D EM Wave Simulation**

## 📐 Mathematical Foundation

### The Canonical Structure

On a configuration manifold Q of dimension n, the theory is constructed on the **1-jet bundle**:

$$M := J^1(Q), \quad \dim M = 2n+1$$

with canonical **contact 1-form**:

$$\boxed{\alpha = du - p_a \, dx^a}$$

### Geometric Algebra Integration

The framework now integrates Clifford Algebra Cl(p,q) to handle geometric transformations and field theories:

- **Geometric Product**: $ab = a \cdot b + a \wedge b$
- **Rotors**: $R = \exp(-B\phi/2)$ generate rotations (elliptic B) and boosts (hyperbolic B)
- **Split Derivative**: $\nabla = \sum e^i \partial_i$

## 🚀 Quick Start

### Installation

```bash
npm install contact-thermodynamics
```

### Basic Usage

```javascript
const CT = require('contact-thermodynamics');

// Create the 13D Grand Contact Manifold
const M13 = CT.grandManifold();

// Create a physical point
const pt = M13.physicalPoint(
    1, 0, 0,      // q¹, q², q³ (spatial position)
    0,            // t (time)
    0,            // ℓ = log(λ) (scale)
    1,            // S (entropy)
    0.5, 0, 0,    // k₁, k₂, k₃ (wavenumber)
    1,            // ω (frequency)
    0,            // Δ (dilatation)
    1,            // T (temperature)
    0             // A (action)
);
```

### Geometric Algebra Usage

```javascript
const { Algebra } = require('contact-thermodynamics/src/multivector');

// Create Spacetime Algebra Cl(1,3)
const sta = new Algebra(1, 3); 
const e0 = sta.e(1), e1 = sta.e(2);

// Boost (hyperbolic rotation)
const B = e1.mul(e0); // e1∧e0
const R = sta.rotor(B, 1.0); // Boost by rapidity 1.0
```

### Differential Geometry

```javascript
const metric = CT.SpacetimeMetric.schwarzschild(1);
const cov = new CT.CovariantDerivative(metric);

// Compute Laplacian of a scalar field
const lap = cov.laplacian(f_func, coords);
```

## 📊 Interactive Demos

1. **[3D EM Wave Simulation](examples/em-wave-3d.html)** (NEW)
   - Real-time Maxwell solver using discrete geometric calculus
   - Interactive 3D visualization of electric fields
   - Add Gaussian pulses and plane waves

2. **Phase Space Demo** (`examples/demo.html`)
   - Interactive exploration of contact manifolds and flows

## 🧪 Testing

Run the extended test suite (65 tests):

```bash
npm test
```

## 🗂 Project Structure

```
contact-thermodynamics/
├── src/
│   ├── index.js              # Main library (Manifolds, GR, DiffGeo)
│   ├── multivector.js        # Geometric Algebra core
│   └── geometric-calculus.js # Discrete operators & fields
├── examples/
│   ├── em-wave-3d.html       # 3D visualization
│   ├── em-wave.js            # Node.js Maxwell solver
│   └── ...                   # Other examples
├── docs/                     # Documentation
├── tests/
│   └── test.js               # Test suite
└── README.md
```

## 📄 License

MIT License — see [LICENSE](LICENSE) for details.
