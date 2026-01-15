# Tutorial 3: Hamiltonian Dynamics on Contact Manifolds

This tutorial covers how to define Hamiltonians and evolve systems on contact manifolds.

## Contact Hamiltonian Systems

Given a smooth function H: M → ℝ, there exists a unique vector field X_H satisfying:

```
ι_{X_H} α = -H
ι_{X_H} dα = dH - (RH)α
```

where R is the Reeb field and RH = ∂H/∂u.

## Hamilton's Equations

In canonical coordinates (x^a, u, p_a), the equations of motion are:

```
ẋ^a = ∂H/∂p_a           (configuration evolution)
ṗ_a = -∂H/∂x^a - p_a·RH  (momentum evolution with friction)
u̇ = p_a·∂H/∂p_a - H     (action/potential evolution)
```

### Key Difference from Symplectic Mechanics

The extra term `-p_a·RH` in the momentum equation allows for **dissipation**. This is why contact geometry is natural for:
- Thermodynamics (entropy production)
- Damped systems
- Open systems

## Creating a Hamiltonian

### Method 1: Direct Function

```javascript
const CT = require('contact-thermodynamics');
const M13 = CT.grandManifold();

// Simple Hamiltonian: H = ω (frequency as energy)
const H_function = coords => coords.omega;

const hamiltonian = new CT.ContactHamiltonian(M13, H_function);
```

### Method 2: With Analytical Gradient

For better numerical accuracy, provide the gradient:

```javascript
// Harmonic oscillator: H = ½(k² + q²)
const H_osc = coords => 0.5 * (coords.k1**2 + coords.q1**2);

const dH_osc = coords => ({
    q1: coords.q1,
    k1: coords.k1,
    // All other partials are 0
    q2: 0, q3: 0, t: 0, ell: 0, S: 0,
    k2: 0, k3: 0, omega: 0, Delta: 0, T: 0,
    A: 0
});

const hamiltonian = new CT.ContactHamiltonian(M13, H_osc, dH_osc);
```

### Method 3: Built-in Hamiltonians

```javascript
// Dispersion relation: H = ω - c|k|
const H_disp = CT.ThermodynamicHamiltonian.dispersionRelation(M13, 1, 0);

// Massive dispersion: H = ω - √(c²|k|² + m²c⁴)
const H_massive = CT.ThermodynamicHamiltonian.dispersionRelation(M13, 1, 1);

// Thermodynamic equation of state
const H_thermo = CT.ThermodynamicHamiltonian.equationOfState(M13, 'ideal');
```

## Evaluating the Hamiltonian

```javascript
const pt = M13.physicalPoint(
    0, 0, 0, 0, 0, 0,
    1, 0, 0,  // k = (1, 0, 0)
    1,        // ω = 1
    0, 0,
    0
);

console.log('H =', hamiltonian.evaluate(pt));
```

## Computing the Vector Field

```javascript
const X = hamiltonian.vectorField(pt);

console.log('ẋ^a:', X.q1, X.q2, X.q3);
console.log('ṗ_a:', X.k1, X.k2, X.k3);
console.log('u̇:', X.A);
```

## Integrating Trajectories

The library uses RK4 (4th-order Runge-Kutta) integration:

```javascript
const dt = 0.1;     // time step
const steps = 100;  // number of steps

const trajectory = hamiltonian.flow(pt, dt, steps);

// trajectory is an array of ContactPoints
console.log('Initial:', trajectory[0].toString());
console.log('Final:', trajectory[100].toString());
```

## Checking Energy Evolution

In symplectic mechanics, H is conserved. In contact mechanics, H may change:

```javascript
const H_values = hamiltonian.hamiltonianEvolution(trajectory);

console.log('H(0) =', H_values[0]);
console.log('H(T) =', H_values[100]);
console.log('ΔH =', H_values[100] - H_values[0]);
```

## Example: Free Particle

```javascript
const CT = require('contact-thermodynamics');
const M13 = CT.grandManifold();

// H = ½|k|²
const H_free = new CT.ContactHamiltonian(M13, coords => {
    return 0.5 * (coords.k1**2 + coords.k2**2 + coords.k3**2);
});

const pt = M13.physicalPoint(
    0, 0, 0, 0, 0, 0,
    1, 0, 0,  // initial momentum
    0, 0, 0,
    0
);

const traj = H_free.flow(pt, 0.1, 50);

// Position should evolve linearly: q = q₀ + k·t
for (let i = 0; i <= 50; i += 10) {
    console.log(`t=${i*0.1}: q1=${traj[i].get('q1').toFixed(3)}`);
}
```

## Example: Dispersion Relation

```javascript
// Massless wave: H = ω - |k| = 0 on shell
const H_wave = CT.ThermodynamicHamiltonian.dispersionRelation(M13, 1, 0);

// Start on the mass shell
const pt_wave = M13.physicalPoint(
    0, 0, 0, 0, 0, 0,
    1, 0, 0,  // |k| = 1
    1,        // ω = 1, so H = 0
    0, 0,
    0
);

console.log('On shell:', H_wave.evaluate(pt_wave));  // ≈ 0

// Evolve
const traj_wave = H_wave.flow(pt_wave, 0.1, 100);
```

## Physical Interpretation

| Hamiltonian | Physical System |
|-------------|-----------------|
| H = ½\|k\|² | Free particle |
| H = ω - c\|k\| | Massless waves |
| H = ω - √(c²\|k\|² + m²) | Massive particles |
| H = ½T - TS | Ideal gas thermodynamics |

## Next Steps

In [Tutorial 4](04-legendrian.md), we'll explore Legendrian submanifolds and Hamilton-Jacobi theory.
