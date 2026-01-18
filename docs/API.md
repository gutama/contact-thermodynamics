# API Reference

Complete API documentation for the Contact Thermodynamics library.

## Table of Contents

- [ContactManifold](#contactmanifold)
- [ContactPoint](#contactpoint)
- [GrandContactManifold](#grandcontactmanifold)
- [HolographicContactManifold](#holographiccontactmanifold)
- [GaugeExtendedManifold](#gaugeextendedmanifold)
- [ContactHamiltonian](#contacthamiltonian)
- [LegendrianSubmanifold](#legendriansubmanifold)
- [ThermodynamicHamiltonian](#thermodynamichamiltonian)
- [SpacetimeMetric](#spacetimemetric)
- [RelativisticHamiltonian](#relativistinghamiltonian)
- [DifferentialForm](#differentialform)
- [TriangleMesh](#trianglemesh) (NEW)
- [MeshGeometricDerivative](#meshgeometricderivative) (NEW)
- [LeapfrogGCMesh](#leapfroggcmesh) (NEW)
- [Utility Functions](#utility-functions)

---

## ContactManifold

Base class for contact manifolds constructed as 1-jet bundles J¹(Q).

### Constructor

```javascript
new ContactManifold(baseCoords, momentaCoords, fiberCoord = 'A')
```

**Parameters:**
- `baseCoords` *(string[])* — Names of base configuration coordinates x^a
- `momentaCoords` *(string[])* — Names of conjugate momenta p_a
- `fiberCoord` *(string)* — Name of fiber coordinate u (default: 'A')

**Throws:** Error if baseCoords and momentaCoords have different lengths.

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `n` | number | Dimension of base manifold Q |
| `dim` | number | Total dimension (2n+1) |
| `baseCoords` | string[] | Base coordinate names |
| `momentaCoords` | string[] | Momentum coordinate names |
| `fiberCoord` | string | Fiber coordinate name |

### Methods

#### `point(coords)`

Create a point on the contact manifold.

```javascript
const pt = manifold.point({ q1: 1, k1: 0.5, A: 0 });
```

**Parameters:**
- `coords` *(Object)* — Coordinate values as key-value pairs

**Returns:** `ContactPoint`

#### `get origin`

Returns the origin point (all coordinates zero).

```javascript
const origin = manifold.origin;
```

**Returns:** `ContactPoint`

#### `contactFormSymbolic()`

Get symbolic representation of the contact form.

```javascript
console.log(manifold.contactFormSymbolic());
// → "dA - k1·dq1 - k2·dq2 - ..."
```

**Returns:** `string`

#### `evaluateContactForm(pt, tangent)`

Evaluate contact form α at a point with a tangent vector.

```javascript
const value = manifold.evaluateContactForm(pt, {
    q1: 1, k1: 0, A: 0.5  // tangent components
});
```

**Parameters:**
- `pt` *(ContactPoint)* — Point on manifold
- `tangent` *(Object)* — Tangent vector components

**Returns:** `number` — Value of α(v)

#### `verifyContactCondition(pt)`

Verify non-degeneracy α ∧ (dα)^n ≠ 0.

```javascript
const volume = manifold.verifyContactCondition(pt);
// Returns n! for canonical contact form
```

**Returns:** `number` — Volume element (n! for canonical form)

#### `reebField(pt)`

Get the Reeb vector field R at a point.

```javascript
const R = manifold.reebField(pt);
// R.A = 1, all others = 0 for canonical form
```

**Returns:** `Object` — Vector field components

---

## ContactPoint

Represents a point on a contact manifold.

### Constructor

```javascript
new ContactPoint(manifold, coords = {})
```

### Methods

#### `get(coord)`

Get coordinate value.

```javascript
const q1 = pt.get('q1');
```

#### `set(coord, value)`

Set coordinate value (returns self for chaining).

```javascript
pt.set('q1', 2.5).set('k1', 0.3);
```

#### `clone()`

Create a copy of the point.

```javascript
const pt2 = pt.clone();
```

#### `add(tangent, dt = 1)`

Add scaled tangent vector (for integration).

```javascript
const newPt = pt.add({ q1: 0.1, k1: -0.05 }, 0.01);
```

**Returns:** `ContactPoint` — New point

---

## GrandContactManifold

The "honest" 13-dimensional Grand Contact Manifold M₁₃ = J¹(Q₆).

**Extends:** `ContactManifold`

### Coordinates

| Base (x^a) | Momenta (p_a) | Interpretation |
|------------|---------------|----------------|
| q¹, q², q³ | k₁, k₂, k₃ | Spatial position ↔ Wavenumber |
| t | ω | Time ↔ Frequency |
| ℓ = log(λ) | Δ | Scale ↔ Dilatation |
| S | T | Entropy ↔ Temperature |

### Constructor

```javascript
const grand = new GrandContactManifold();
// or
const grand = CT.grandManifold();
```

### Methods

#### `physicalPoint(q1, q2, q3, t, ell, S, k1, k2, k3, omega, Delta, T, A)`

Create point with physical naming.

```javascript
const pt = grand.physicalPoint(
    1, 0, 0,      // spatial position
    0,            // time
    0,            // log(scale)
    1,            // entropy
    0.5, 0, 0,    // wavenumber
    1,            // frequency
    0,            // dilatation
    1,            // temperature
    0             // action
);
```

#### `spatialPosition(pt)`

Extract spatial position (q¹, q², q³).

```javascript
const [x, y, z] = grand.spatialPosition(pt);
```

**Returns:** `number[]`

#### `waveVector(pt)`

Extract wave vector (k₁, k₂, k₃).

```javascript
const [kx, ky, kz] = grand.waveVector(pt);
```

**Returns:** `number[]`

---

## HolographicContactManifold

The 7-dimensional Holographic Contact Manifold M₇ = J¹(Q₃) with emergent space.

**Extends:** `ContactManifold`

### Coordinates

| Base (x^a) | Momenta (p_a) |
|------------|---------------|
| t | ω |
| ℓ | Δ |
| S | T |

Space q^i becomes emergent: q^i = q^i(t, ℓ, S)

### Constructor

```javascript
const holo = new HolographicContactManifold();
// or
const holo = CT.holographicManifold();
```

### Methods

#### `holographicPoint(t, ell, S, omega, Delta, T, A)`

Create holographic point.

```javascript
const pt = holo.holographicPoint(0, 0, 1, 1, 0, 1, 0);
```

#### `createEmergentSpace(pt, fieldFunc)`

Define emergent spatial configuration.

```javascript
const emergent = holo.createEmergentSpace(pt, (t, ell, S) => {
    const a = Math.exp(ell);
    return [a * Math.cos(t), a * Math.sin(t), 0];
});
// emergent = { q1: ..., q2: ..., q3: ... }
```

**Parameters:**
- `pt` *(ContactPoint)* — Base point
- `fieldFunc` *(Function)* — (t, ell, S) → [q1, q2, q3]

**Returns:** `Object` — Emergent coordinates { q1, q2, q3 }

---

## GaugeExtendedManifold

Gauge-extended 15-dimensional manifold M₁₅ with additional (φ, I) pair.

**Extends:** `GrandContactManifold`

### Additional Coordinates

| Coordinate | Conjugate | Interpretation |
|------------|-----------|----------------|
| φ | I | Gauge phase ↔ Gauge flux |

### Constructor

```javascript
const gauge = new GaugeExtendedManifold();
// or
const gauge = CT.gaugeExtended();
```

---

## ContactHamiltonian

Hamiltonian dynamics on contact manifolds.

### Constructor

```javascript
new ContactHamiltonian(manifold, H, dH = null)
```

**Parameters:**
- `manifold` *(ContactManifold)* — The contact manifold
- `H` *(Function)* — Hamiltonian H(coords) → number
- `dH` *(Function, optional)* — Gradient of H

### Methods

#### `evaluate(pt)`

Evaluate Hamiltonian at a point.

```javascript
const E = hamiltonian.evaluate(pt);
```

#### `gradient(pt, h = 1e-7)`

Compute numerical gradient.

```javascript
const grad = hamiltonian.gradient(pt);
// grad = { q1: ∂H/∂q1, k1: ∂H/∂k1, ... }
```

#### `reebComponent(pt)`

Get RH = ∂H/∂u component.

```javascript
const RH = hamiltonian.reebComponent(pt);
```

#### `vectorField(pt)`

Compute contact Hamiltonian vector field X_H.

```javascript
const X = hamiltonian.vectorField(pt);
// X = { q1: q̇1, k1: k̇1, A: Ȧ, ... }
```

**The equations:**
- ẋ^a = ∂H/∂p_a
- ṗ_a = -∂H/∂x^a - p_a · ∂H/∂u
- u̇ = p_a · ∂H/∂p_a - H

#### `flow(pt, dt, steps = 1)`

Integrate dynamics using RK4.

```javascript
const trajectory = hamiltonian.flow(pt, 0.1, 100);
// Returns array of ContactPoints
```

**Parameters:**
- `pt` *(ContactPoint)* — Initial point
- `dt` *(number)* — Time step
- `steps` *(number)* — Number of steps

**Returns:** `ContactPoint[]` — Trajectory (length = steps + 1)

#### `hamiltonianEvolution(trajectory)`

Get Hamiltonian values along trajectory.

```javascript
const Hvals = hamiltonian.hamiltonianEvolution(trajectory);
```

**Returns:** `number[]`

---

## LegendrianSubmanifold

n-dimensional submanifold L ⊂ M with α|_L = 0.

### Constructor

```javascript
new LegendrianSubmanifold(manifold, generatingFunc, dA = null)
```

**Parameters:**
- `manifold` *(ContactManifold)* — The ambient contact manifold
- `generatingFunc` *(Function)* — A(x) generating function
- `dA` *(Function, optional)* — Gradient ∂_a A(x)

### Theory

Given generating function A(x), the Legendrian submanifold is defined by:
- u = A(x)
- p_a = ∂_a A(x)

This automatically satisfies α|_L = 0 since:
α = du - p_a dx^a = dA - ∂_a A dx^a = 0

### Methods

#### `gradient(x, h = 1e-7)`

Compute ∂_a A numerically.

```javascript
const grad = legendrian.gradient({ q1: 1, t: 0 });
```

#### `lift(x)`

Lift base coordinates to contact manifold.

```javascript
const pt = legendrian.lift({ q1: 1, q2: 0, q3: 0, t: 0, ell: 0, S: 0 });
// pt has A = A(x) and p_a = ∂_a A
```

#### `verifyLegendrianCondition(x)`

Verify α|_L = 0.

```javascript
const valid = legendrian.verifyLegendrianCondition(x);
// Always true for generating function construction
```

#### `hamiltonJacobiResidual(x, hamiltonian)`

Compute H(x, A(x), ∂A(x)).

```javascript
const residual = legendrian.hamiltonJacobiResidual(x, hamiltonian);
// Should be 0 on solutions of Hamilton-Jacobi equation
```

#### `sample(sampler, n)`

Sample points on the Legendrian.

```javascript
const points = legendrian.sample(() => ({
    q1: Math.random(), q2: 0, q3: 0, t: 0, ell: 0, S: 0
}), 100);
```

---

## ThermodynamicHamiltonian

Specialized Hamiltonians for thermodynamic systems.

**Extends:** `ContactHamiltonian`

### Static Methods

#### `ThermodynamicHamiltonian.dispersionRelation(manifold, c, mass)`

Create dispersion relation Hamiltonian.

```javascript
// Massless: H = ω - c|k|
const Hmassless = ThermodynamicHamiltonian.dispersionRelation(M13, 1, 0);

// Massive: H = ω - √(c²|k|² + m²c⁴)
const Hmassive = ThermodynamicHamiltonian.dispersionRelation(M13, 1, 1);
```

#### `ThermodynamicHamiltonian.equationOfState(manifold, type)`

Create equation of state Hamiltonian.

```javascript
// Ideal gas
const Hideal = ThermodynamicHamiltonian.equationOfState(M13, 'ideal');

// Van der Waals
const HvdW = ThermodynamicHamiltonian.equationOfState(M13, 'van_der_waals');
```

---

## SpacetimeMetric

Spacetime metric g_μν(x) for GR coupling.

### Constructor

```javascript
new SpacetimeMetric(metricFunc, inverseFunc = null)
```

**Parameters:**
- `metricFunc` *(Function)* — (x) → 4×4 covariant metric
- `inverseFunc` *(Function, optional)* — (x) → 4×4 contravariant metric

### Methods

#### `covariant(x)`

Get g_μν at point x.

```javascript
const g = metric.covariant([0, 10, Math.PI/2, 0]);
// g[0][0] = g_tt, g[1][1] = g_rr, etc.
```

#### `contravariant(x)`

Get g^μν at point x (inverts numerically if not provided).

```javascript
const gInv = metric.contravariant(x);
```

### Static Methods

#### `SpacetimeMetric.minkowski()`

Flat Minkowski metric η_μν = diag(+1, -1, -1, -1).

```javascript
const mink = SpacetimeMetric.minkowski();
```

#### `SpacetimeMetric.schwarzschild(M)`

Schwarzschild black hole metric.

```javascript
const schw = SpacetimeMetric.schwarzschild(1);
// Coordinates: (t, r, θ, φ)
// g_tt = 1 - 2M/r, g_rr = -(1 - 2M/r)^{-1}, etc.
```

#### `SpacetimeMetric.flrw(a, k)`

FLRW cosmological metric.

```javascript
const flrw = SpacetimeMetric.flrw(t => Math.exp(0.1 * t), 0);
// a(t) = scale factor, k = curvature (+1, 0, -1)
```

---

## RelativisticHamiltonian

Mass-shell constraint Hamiltonian for relativistic dynamics.

### Constructor

```javascript
new RelativisticHamiltonian(metric, mass = 1, gaugePotential = null, charge = 0)
```

**Parameters:**
- `metric` *(SpacetimeMetric)* — Spacetime metric
- `mass` *(number)* — Particle mass
- `gaugePotential` *(Function, optional)* — A_μ(x) for EM coupling
- `charge` *(number)* — Electric charge q

### The Hamiltonian

$$H = \frac{1}{2}g^{\mu\nu}(p_\mu - qA_\mu)(p_\nu - qA_\nu) - \frac{1}{2}m^2$$

On shell: H = 0.

### Methods

#### `evaluate(x, p)`

Evaluate Hamiltonian.

```javascript
const H = relH.evaluate([0, 10, Math.PI/2, 0], [1, 0, 0, 0.2]);
```

#### `geodesicEquations(x, p)`

Get Hamilton's equations.

```javascript
const { xDot, pDot } = relH.geodesicEquations(x, p);
// xDot = [ṫ, ṙ, θ̇, φ̇]
// pDot = [ṗ_t, ṗ_r, ṗ_θ, ṗ_φ]
```

#### `integrateGeodesic(x0, p0, dtau, steps)`

Integrate geodesic with RK4.

```javascript
const geodesic = relH.integrateGeodesic(
    [0, 10, Math.PI/2, 0],  // initial position
    [1.05, 0, 0, 0.2],      // initial momentum
    0.1,                     // proper time step
    1000                     // number of steps
);
// Returns: [{ x, p, tau }, ...]
```

---

## DifferentialForm

Symbolic differential form representation.

### Constructor

```javascript
new DifferentialForm(degree, terms = [])
```

### Static Methods

#### `DifferentialForm.oneForm(coord)`

Create 1-form dx^a.

```javascript
const dx = DifferentialForm.oneForm('q1');
```

### Methods

#### `wedge(other)`

Wedge product α ∧ β.

```javascript
const twoForm = dx.wedge(dy);
```

---

## TriangleMesh

Triangle mesh data structure with typed arrays and precomputed topology.

### Constructor

```javascript
new TriangleMesh(vertices, faces, opts = {})
```

**Parameters:**
- `vertices` *(Float64Array)* — Vertex positions (nVertices × 3)
- `faces` *(Uint32Array)* — Triangle indices (nFaces × 3)
- `opts.buildTopology` *(boolean)* — Auto-build topology (default: true)

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `nVertices` | number | Number of vertices |
| `nFaces` | number | Number of triangles |
| `nEdges` | number | Number of edges |
| `vertices` | Float64Array | Vertex positions (×3) |
| `faces` | Uint32Array | Triangle indices (×3) |
| `edges` | Uint32Array | Edge vertex pairs (×2) |
| `boundaryVertices` | Uint8Array | Boundary mask (1 = boundary) |
| `boundaryEdges` | Uint8Array | Boundary edge mask |

### Methods

#### `buildTopology()`

Build edge list and adjacency structures.

#### `getEdgeIndex(v0, v1)`

Get edge index for vertex pair (order-independent).

```javascript
const eIdx = mesh.getEdgeIndex(0, 5);
```

#### `getVertex(idx)`

Get vertex position as [x, y, z].

```javascript
const [x, y, z] = mesh.getVertex(10);
```

#### `faceAreas()`

Compute triangle areas (cached).

```javascript
const areas = mesh.faceAreas(); // Float64Array
```

#### `edgeLengths()`

Compute edge lengths (cached).

#### `cotanWeights()`

Compute cotan weights per edge (cached). These encode the mesh metric.

```javascript
const weights = mesh.cotanWeights();
```

#### `vertexDualAreas()`

Compute mixed Voronoi dual areas per vertex (cached). Uses circumcentric for non-obtuse triangles, barycentric fallback for obtuse.

```javascript
const dualAreas = mesh.vertexDualAreas();
```

### Static Factory Methods

#### `TriangleMesh.createGrid(nx, ny, Lx, Ly)`

Create a regular grid mesh.

```javascript
const mesh = TriangleMesh.createGrid(16, 16, 1.0, 1.0);
// 256 vertices, 450 triangles
```

#### `TriangleMesh.createDisk(nRadial, nAngular, radius)`

Create a triangulated disk mesh.

```javascript
const disk = TriangleMesh.createDisk(8, 24, 1.0);
```

#### `TriangleMesh.createIcosahedron(radius)`

Create an icosahedron (closed mesh, 12 vertices, 20 faces).

```javascript
const ico = TriangleMesh.createIcosahedron(1.0);
```

---

## MeshGeometricDerivative

Discrete geometric derivative ∇ on triangle meshes (FTGC).

Implements the split differential operator: ∇ = Σᵢ eⁱ ⊗ ∂/∂xᵢ

**Key insight**: Cotan weights encode the mesh metric (reciprocal basis).

### Constructor

```javascript
new MeshGeometricDerivative(mesh)
```

### Methods

#### `grad(f)`

Gradient of scalar field.

```javascript
const gradF = nabla.grad(f); // Float64Array on edges
```

**Parameters:**
- `f` *(Float64Array)* — Scalar field on vertices

**Returns:** Float64Array — Vector field on edges

#### `div(V)`

Divergence of vector field.

```javascript
const divV = nabla.div(V); // Float64Array on vertices
```

#### `curl(V)`

Curl of vector field.

```javascript
const curlV = nabla.curl(V); // Float64Array on faces
```

#### `laplacian(f, opts)`

Laplacian of scalar field (cotan Laplacian).

```javascript
const lapF = nabla.laplacian(f, {
    dirichletMask,    // Uint8Array (1 = fixed)
    dirichletValues   // Float64Array
});
```

**Returns:** Float64Array — Laplacian on vertices

#### `laplacianMatrix()`

Build the Laplacian as a sparse matrix.

```javascript
const L = nabla.laplacianMatrix();
const Lu = L.matvec(u);
```

#### `wedge(F)` / `inner(F)` / `apply(F)`

Apply outer/inner/full geometric derivative to a MeshMultivectorField.

#### `verifyCurlGradZero(f)`

Verify ∇∧(∇∧f) = 0 (returns max error).

#### `verifyLaplacianConstantZero(c)`

Verify ∇²(constant) ≈ 0 (returns max error).

---

## MeshMultivectorField

Multivector field on a mesh with staggered storage.

- Grade 0 (scalars) → Vertices
- Grade 1 (vectors) → Edges
- Grade 2 (bivectors) → Faces

### Constructor

```javascript
new MeshMultivectorField(mesh, grade0, grade1, grade2)
```

### Methods

- `add(other)`, `sub(other)`, `scale(s)` — Arithmetic
- `gradeSelect(k)` — Extract grade-k component
- `clone()` — Deep copy

### Static Factory Methods

```javascript
const scalarField = MeshMultivectorField.scalarField(mesh, values);
const vectorField = MeshMultivectorField.vectorField(mesh, values);
const bivectorField = MeshMultivectorField.bivectorField(mesh, values);
```

---

## LeapfrogGCMesh

Leapfrog time-stepping solver using FTGC on triangle meshes.

### Constructor

```javascript
new LeapfrogGCMesh(mesh)
```

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `mesh` | TriangleMesh | The underlying mesh |
| `nabla` | MeshGeometricDerivative | The ∇ operator |

### Methods

#### `estimateCFL(c)`

Estimate maximum stable time step for wave equation.

```javascript
const dt = solver.estimateCFL(1.0); // wave speed c = 1
```

#### `estimateCFLHeat(alpha)`

Estimate maximum stable time step for heat equation.

```javascript
const dt = solver.estimateCFLHeat(0.1); // diffusion α = 0.1
```

#### `waveStep(uPrev, uCurr, dt, c, opts)`

One leapfrog step for wave equation: u(t+dt) = 2u(t) - u(t-dt) + c²dt²∇²u(t)

```javascript
const uNext = solver.waveStep(uPrev, uCurr, dt, 1.0, {
    dirichletMask,
    dirichletValues
});
```

#### `waveSimulate(u0, v0, dt, nSteps, c, opts)`

Simulate wave equation from initial conditions.

```javascript
const uFinal = solver.waveSimulate(u0, v0, dt, 100, 1.0, {
    dirichletMask,
    callback: (step, u) => console.log(step)
});
```

#### `heatStepExplicit(u, dt, alpha, opts)`

Explicit Euler step for heat equation.

#### `heatStepImplicit(u, dt, alpha, opts)`

Implicit Euler step (unconditionally stable).

```javascript
const uNext = solver.heatStepImplicit(u, dt, 0.1, {
    dirichletMask,
    maxIter: 100,
    tol: 1e-8
});
```

#### `heatSimulate(u0, dt, nSteps, alpha, opts)`

Simulate heat equation.

```javascript
const uFinal = solver.heatSimulate(u0, dt, 100, 0.1, {
    implicit: true,
    dirichletMask,
    callback: (step, u) => { /* track variance */ }
});
```

#### `maxwellStep(E, B, dt, c)`

Yee-style leapfrog for Maxwell equations.

- E: vector field on edges
- B: bivector field on faces

```javascript
const { E: E_new, B: B_new } = solver.maxwellStep(E, B, dt, 1.0);
```

#### `energy(u)` / `mass(u)` / `variance(u)`

Compute integral quantities for a scalar field.

```javascript
const E = solver.energy(u);   // ∫ u² dA
const M = solver.mass(u);     // ∫ u dA
const V = solver.variance(u); // ∫ (u - mean)² dA
```

---

## Utility Functions

### Mesh Boundary Helpers

```javascript
// Apply Dirichlet BCs to a field (in-place)
CT.applyDirichlet(f, mask, values);

// Create Dirichlet mask from mesh boundary
const mask = CT.boundaryDirichletMask(mesh);
```

### Factory Functions

```javascript
// Create manifolds
const M13 = CT.grandManifold();
const M7 = CT.holographicManifold();
const M15 = CT.gaugeExtended();
```

### `summaryTable()`

Get theory summary table.

```javascript
const table = CT.summaryTable();
// { columns: [...], rows: [[...], [...], ...] }
```

### Constants

```javascript
CT.EPSILON  // Numerical tolerance (1e-12)
```

---

## TypeScript Support

Type definitions are available in `src/index.d.ts`. Example:

```typescript
import * as CT from 'contact-thermodynamics';

const manifold: CT.GrandContactManifold = CT.grandManifold();
const pt: CT.ContactPoint = manifold.physicalPoint(1, 0, 0, 0, 0, 1, 0.5, 0, 0, 1, 0, 1, 0);
```
