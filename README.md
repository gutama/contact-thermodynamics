# Contact Thermodynamics

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Node.js](https://img.shields.io/badge/Node.js-18%2B-green.svg)](https://nodejs.org/)
[![Tests](https://img.shields.io/badge/Tests-44%20passing-brightgreen.svg)](#testing)

A JavaScript implementation of **Contact Geometry for Extended Thermodynamics** based on 1-jet bundles. This framework provides a rigorous geometric foundation for thermodynamic systems, wave mechanics, and their gravitational extensions.

<p align="center">
  <img src="docs/assets/contact-manifold-diagram.svg" alt="Contact Manifold Structure" width="600">
</p>

## 🌟 Features

- **Contact Manifolds**: Full implementation of 1-jet bundles J¹(Q) with canonical contact form
- **Three Model Scales**:
  - **Grand Model (M₁₃)**: 13-dimensional full phase space
  - **Holographic Model (M₇)**: 7-dimensional with emergent space
  - **Gauge-Extended (M₁₅)**: 15-dimensional with gauge degrees of freedom
- **Contact Hamiltonian Dynamics**: Vector fields, Reeb fields, RK4 integration
- **Legendrian Submanifolds**: Hamilton-Jacobi theory via generating functions
- **General Relativity Extension**: Spacetime metrics coupled through Hamiltonians
- **Interactive Visualization**: Browser-based demos with phase space plotting

## 📐 Mathematical Foundation

### The Canonical Structure

On a configuration manifold Q of dimension n, the theory is constructed on the **1-jet bundle**:

$$M := J^1(Q), \quad \dim M = 2n+1$$

with canonical **contact 1-form**:

$$\boxed{\alpha = du - p_a \, dx^a}$$

where:
- $x^a$ — base configuration coordinates (dimension n)
- $u$ — fiber coordinate (generating potential / action)
- $p_a$ — conjugate momenta coordinates

### Contact Non-Degeneracy

A contact manifold satisfies the **non-degeneracy condition**:

$$\alpha \wedge (d\alpha)^n \neq 0 \quad \text{everywhere on } M$$

This ensures the contact structure is maximally non-integrable.

### Dimensionality Theorem

| Model | Base Dim (n) | Total Dim (2n+1) |
|-------|--------------|------------------|
| Grand | 6 | **13** |
| Holographic | 3 | **7** |
| Gauge-Extended | 7 | **15** |

## 🚀 Quick Start

### Installation

```bash
npm install contact-thermodynamics
```

Or clone directly:

```bash
git clone https://github.com/gutama/contact-thermodynamics.git
cd contact-thermodynamics
npm install
```

### Basic Usage

```javascript
const CT = require('contact-thermodynamics');

// Create the 13D Grand Contact Manifold
const M13 = CT.grandManifold();
console.log(M13.toString());
// → Grand Contact Manifold M₁₃ = J¹(Q₆): dim=13

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

// Verify contact non-degeneracy
console.log('α ∧ (dα)⁶ =', M13.verifyContactCondition(pt));
// → 720 (= 6!)
```

### Contact Hamiltonian Dynamics

```javascript
// Dispersion relation: H = ω - c|k| (massless waves)
const H = CT.ThermodynamicHamiltonian.dispersionRelation(M13, 1, 0);

// Evolve the system
const trajectory = H.flow(pt, 0.1, 100);

// Check energy evolution
const energies = H.hamiltonianEvolution(trajectory);
console.log('Initial H:', energies[0]);
console.log('Final H:', energies[100]);
```

## 📖 Documentation

### [Full API Reference →](docs/API.md)

### Core Classes

| Class | Description |
|-------|-------------|
| [`ContactManifold`](docs/API.md#contactmanifold) | Base class for J¹(Q) bundles |
| [`GrandContactManifold`](docs/API.md#grandcontactmanifold) | M₁₃ with full (q,t,ℓ,S) coordinates |
| [`HolographicContactManifold`](docs/API.md#holographiccontactmanifold) | M₇ with emergent spatial fields |
| [`ContactHamiltonian`](docs/API.md#contacthamiltonian) | Hamiltonian dynamics on contact manifolds |
| [`LegendrianSubmanifold`](docs/API.md#legendriansubmanifold) | n-dimensional submanifolds with α|_L = 0 |
| [`SpacetimeMetric`](docs/API.md#spacetimemetric) | Curved spacetime metrics (Minkowski, Schwarzschild, FLRW) |
| [`RelativisticHamiltonian`](docs/API.md#relativistinghamiltonian) | Mass-shell constraint for GR coupling |

### Tutorials

1. [Understanding Contact Geometry](docs/tutorials/01-contact-geometry.md)
2. [The Grand vs Holographic Models](docs/tutorials/02-models.md)
3. [Hamiltonian Dynamics on Contact Manifolds](docs/tutorials/03-dynamics.md)
4. [Legendrian Submanifolds and Hamilton-Jacobi](docs/tutorials/04-legendrian.md)
5. [Gravitational Extension](docs/tutorials/05-gravity.md)

## 🔬 The Three Models

### Grand Model M₁₃

The "honest" full phase space where all variables are independent base coordinates.

**Base Configuration Q₆:**
$$x^a = (q^1, q^2, q^3, t, \ell, S), \quad \ell := \log\lambda$$

**Conjugate Momenta:**
$$p_a = (k_1, k_2, k_3, \omega, \Delta, T)$$

**Contact Form:**
$$\alpha = d\mathcal{A} - k_i\,dq^i - \omega\,dt - \Delta\,d\ell - T\,dS$$

```javascript
const grand = CT.grandManifold();
// dim = 13, coordinates: q¹,q²,q³,t,ℓ,S,A,k₁,k₂,k₃,ω,Δ,T
```

### Holographic Model M₇

Space is demoted to dependent scalar fields on the reduced base.

**Base Configuration Q₃:**
$$x^a = (t, \ell, S)$$

**Emergent Space:**
$$q^i = q^i(t, \ell, S) \quad \text{(scalar fields on } Q_3\text{)}$$

```javascript
const holo = CT.holographicManifold();

// Define emergent spatial configuration
const emergent = holo.createEmergentSpace(pt, (t, ell, S) => {
    const scale = Math.exp(ell);
    return [scale * Math.cos(t), scale * Math.sin(t), 0];
});
```

### Gauge-Extended Model M₁₅

Adds an independent gauge canonical pair (φ, I).

```javascript
const gauge = CT.gaugeExtended();
// Adds: φ (gauge phase), I (gauge flux/current)
// dim = 15
```

## ⚡ Dynamics

### Contact Hamiltonian Vector Field

Given H: M → ℝ, the contact Hamiltonian vector field X_H satisfies:

$$\iota_{X_H}\alpha = -H, \quad \iota_{X_H}d\alpha = dH - (RH)\alpha$$

### Hamilton's Equations (Contact Form)

$$\dot{x}^a = \frac{\partial H}{\partial p_a}$$

$$\dot{p}_a = -\frac{\partial H}{\partial x^a} - p_a \cdot \frac{\partial H}{\partial u}$$

$$\dot{u} = p_a \cdot \frac{\partial H}{\partial p_a} - H$$

### Reeb Vector Field

$$\alpha(R) = 1, \quad d\alpha(R, \cdot) = 0$$

For canonical α: R = ∂/∂u

```javascript
const reeb = manifold.reebField(pt);
// reeb.A = 1, all other components = 0
```

## 🌌 Gravitational Extension

The contact structure α provides **kinematics** (locally Darboux-flat). Spacetime **curvature** enters through the Hamiltonian.

### Relativistic Mass-Shell Constraint

$$H = \frac{1}{2}g^{\mu\nu}(x)(p_\mu - qA_\mu)(p_\nu - qA_\nu) - \frac{1}{2}m^2 = 0$$

### Hamilton-Jacobi Master PDE

$$\frac{1}{2}g^{\mu\nu}(x)(\partial_\mu\mathcal{A} - qA_\mu)(\partial_\nu\mathcal{A} - qA_\nu) - \frac{1}{2}m^2 = 0$$

### Available Metrics

```javascript
// Minkowski (flat)
const mink = CT.SpacetimeMetric.minkowski();

// Schwarzschild (black hole)
const schw = CT.SpacetimeMetric.schwarzschild(M);

// FLRW (cosmology)
const flrw = CT.SpacetimeMetric.flrw(a, k);
```

### Geodesic Integration

```javascript
const metric = CT.SpacetimeMetric.schwarzschild(1);
const relH = new CT.RelativisticHamiltonian(metric, 1);

// Initial conditions
const x0 = [0, 10, Math.PI/2, 0];  // t, r, θ, φ
const p0 = [1.05, 0, 0, 0.2];      // E, p_r, p_θ, L

// Integrate geodesic
const geodesic = relH.integrateGeodesic(x0, p0, 0.1, 1000);
```

## 🧪 Testing

Run the test suite:

```bash
npm test
```

All 44 tests validate:
- Dimensionality theorem
- Coordinate structures
- Contact form properties
- Reeb vector fields
- Hamiltonian dynamics
- Legendrian lifts
- GR metric computations
- Flow conservation

## 📊 Interactive Demo

Open `examples/demo.html` in a browser for an interactive visualization featuring:

- Model selection (Grand/Holographic/Gauge-Extended)
- Hamiltonian selection (dispersion, massive, thermodynamic)
- Parameter adjustment (mass, speed, temperature)
- Phase space trajectory plotting
- Real-time dynamics simulation

## 🗂 Project Structure

```
contact-thermodynamics/
├── src/
│   └── index.js              # Main library
├── examples/
│   ├── demo.html             # Interactive browser demo
│   ├── basic-usage.js        # Getting started examples
│   ├── dynamics.js           # Hamiltonian flow examples
│   ├── legendrian.js         # Hamilton-Jacobi examples
│   └── geodesics.js          # GR geodesic examples
├── docs/
│   ├── API.md                # Full API reference
│   ├── THEORY.md             # Mathematical background
│   └── tutorials/            # Step-by-step guides
├── tests/
│   └── test.js               # Test suite
├── package.json
├── LICENSE
└── README.md
```

## 📚 References

### Primary Sources

1. **Arnold, V.I.** — *Mathematical Methods of Classical Mechanics* (Springer, 1989)
2. **Libermann, P. & Marle, C.-M.** — *Symplectic Geometry and Analytical Mechanics* (Reidel, 1987)
3. **Geiges, H.** — *An Introduction to Contact Topology* (Cambridge, 2008)

### Contact Hamiltonian Systems

4. **Bravetti, A.** — "Contact Hamiltonian Dynamics: The Concept and Its Use" (Entropy, 2017)
5. **de León, M. & Sardón, C.** — "Geometry of Contact Hamiltonian Systems" (J. Phys. A, 2017)

### Thermodynamic Geometry

6. **Mrugała, R.** — "Geometric formulation of equilibrium phenomenological thermodynamics" (Rep. Math. Phys., 1978)
7. **Ruppeiner, G.** — "Riemannian geometry in thermodynamic fluctuation theory" (Rev. Mod. Phys., 1995)

### Jet Bundles

8. **Saunders, D.J.** — *The Geometry of Jet Bundles* (Cambridge, 1989)

## 🤝 Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Development Setup

```bash
git clone https://github.com/gutama/contact-thermodynamics.git
cd contact-thermodynamics
npm install
npm test
```

## 📄 License

MIT License — see [LICENSE](LICENSE) for details.

## 🙏 Acknowledgments

This implementation draws on the geometric mechanics tradition from Arnold, Marsden, and Weinstein, combined with modern contact Hamiltonian approaches developed by Bravetti, de León, and collaborators.

---

<p align="center">
  <i>Kinematics from α; curvature from g<sub>μν</sub>(x) inside H</i>
</p>
