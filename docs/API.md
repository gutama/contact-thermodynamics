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
- [TriangleMesh](#trianglemesh)
- [MeshGeometricDerivative](#meshgeometricderivative)
- [LeapfrogGCMesh](#leapfroggcmesh)
- [Algebra](#algebra) (NEW)
- [Multivector](#multivector) (NEW)
- [RiemannianGA](#riemannianga) (NEW)
- [ProbabilityManifold](#probabilitymanifold) (NEW)
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

---

## Entropic Gravity

Implementation of Bianconi's "Gravity from Entropy" framework, integrated with Contact Thermodynamics.

### Logic

Gravity emerges from the interplay between two metrics:
- **g**: The "quantum" spacetime metric (operator)
- **G**: The "classical" matter-induced metric (state)
- **S(G||g)**: Relative entropy action generating dynamics

### MatterInducedMetric

Constructs the metric G induced by matter fields.

#### Constructor

```javascript
new EntropicGravity.MatterInducedMetric(options)
```

**Parameters:**
- `options` *(Object)*
    - `scalarField` *(Function)*: φ(x) → number
    - `vectorPotential` *(Function)*: A(x) → [A₀, A₁, A₂, A₃]
    - `twoFormPotential` *(Function)*: B(x) → 4×4 antisymmetric matrix
    - `ricciTensor` *(Function)*: R(x) → 4×4 matrix (Ricci curvature of g)

#### Methods

- `covariant(x, gInv)`: Returns G_μν
- `contravariant(x, gInv)`: Returns G^μν
- `fieldStrength(x)`: Computes F_μν = ∂_μA_ν - ∂_νA_μ

### TwoMetricSystem

Manages the pair (g, G).

#### Constructor

```javascript
new EntropicGravity.TwoMetricSystem(g, G)
```

**Parameters:**
- `g`: SpacetimeMetric instance
- `G`: MatterInducedMetric instance

#### Methods

- `spacetimeMetric(x)`: Returns g_μν
- `matterMetric(x)`: Returns G_μν
- `metricRatio(x)`: Returns G · g⁻¹

### RelativeEntropyAction

Computes S(G||g).

#### Methods

- `localDensity(x)`: returns Tr(Gg⁻¹ ln(Gg⁻¹) - (Gg⁻¹ - I))
- `integrand(x)`: returns local density * √|g|

### EntropicGravityHamiltonian

Contact Hamiltonian including the entropic potential.

$$ H = H_{geo} + \alpha S(G||g) $$

#### Constructor

```javascript
new EntropicGravity.EntropicGravityHamiltonian(manifold, system, options)
```

**Parameters:**
- `manifold`: ContactManifold
- `system`: TwoMetricSystem
- `options`: { mass, entropicCoupling }

#### Methods

- `evaluate(pt)`: Returns H(pt)
- `flow(pt, dt, steps)`: Integrates trajectory


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

---

## Algebra

Clifford Algebra Cl(p, q, r) implementation for Geometric Algebra operations.

**Module:** `src/multivector.js`

### Constructor

```javascript
new Algebra(p, q, r = 0, basisNames = null)
```

**Parameters:**
- `p` *(number)* — Number of positive-signature basis vectors
- `q` *(number)* — Number of negative-signature basis vectors
- `r` *(number)* — Number of zero-signature basis vectors (default: 0)
- `basisNames` *(string[], optional)* — Custom basis vector names

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `p`, `q`, `r` | number | Signature |
| `dim` | number | Space dimension (p + q + r) |
| `bladeCount` | number | Total blades (2^dim) |
| `basisNames` | string[] | Basis vector names |

### Methods

#### `e(i)`

Get basis vector eᵢ (1-indexed).

```javascript
const algebra = new Algebra(3, 0, 0);
const e1 = algebra.e(1);  // e₁
const e2 = algebra.e(2);  // e₂
```

#### `scalar(s)`

Create scalar multivector.

```javascript
const half = algebra.scalar(0.5);
```

#### `zero()`

Create zero multivector.

#### `vector(components)`

Create vector from components.

```javascript
const v = algebra.vector([1, 2, 3]);  // e₁ + 2e₂ + 3e₃
```

---

## Multivector

Multivector in a Clifford Algebra.

### Constructor

```javascript
new Multivector(algebra, coeffs = {})
```

**Parameters:**
- `algebra` *(Algebra)* — The parent algebra
- `coeffs` *(Object)* — Coefficient map: blade bitmap → value

### Methods

#### Arithmetic

```javascript
const sum = a.add(b);
const diff = a.sub(b);
const scaled = a.scale(2);
const neg = a.neg();
```

#### Products

```javascript
const geom = a.mul(b);    // Geometric product ab
const outer = a.wedge(b); // Outer product a∧b
const inner = a.dot(b);   // Inner product a·b
```

#### Grade Extraction

```javascript
const s = mv.scalar();      // Grade-0 (scalar part)
const v = mv.grade(1);      // Grade-1 (vector part)
const biv = mv.grade(2);    // Grade-2 (bivector part)
```

#### Norms and Duals

```javascript
const rev = mv.reverse();   // Reversion ~A
const mag = mv.norm();      // |A| = √⟨A~A⟩₀
const dual = mv.dual();     // Hodge dual
```

---

## RiemannianGA

Coordinate-free Riemannian geometry via Geometric Algebra.

**Module:** `src/riemannian-ga.js`

### Bivector Classes

#### `Bivector2D`

Bivector in 2D (single component e₁₂).

```javascript
const biv = new Bivector2D(1.5);  // 1.5 e₁∧e₂
const rotated = biv.commutatorWithVector([1, 0]);  // Rotate vector
```

#### `Bivector3D`

Bivector in 3D (three components: e₂₃, e₃₁, e₁₂).

```javascript
const biv = new Bivector3D(0, 0, 1);  // e₁∧e₂
const biv2 = Bivector3D.fromWedge([1, 0, 0], [0, 1, 0]);
```

**Methods:**
- `add(other)`, `sub(other)`, `scale(s)`, `neg()`
- `norm()`, `axis()`, `angle()`
- `commutatorWithVector(v)` — ω × v rotation

### Manifold Classes

| Class | Description | Curvature |
|-------|-------------|-----------|
| `Sphere2D(R)` | 2-sphere (embedded in R³) | K = 1/R² |
| `Torus2D(R, r)` | Torus (R: major, r: minor) | K = cos(θ)/(r(R+r·cos(θ))) |
| `HyperbolicPlane()` | Upper half-plane H² | K = −1 |
| `Sphere3D(R)` | 3-sphere (embedded in R⁴) | K = 1/R² |
| `HyperbolicSpace3D()` | Hyperbolic 3-space H³ | K = −1 |

#### Example: Creating a Manifold

```javascript
const { Sphere2D, Curvature2Form } = require('./src/riemannian-ga.js');

const sphere = new Sphere2D(1.0);  // Unit sphere
const frame = sphere.frame([Math.PI/4, Math.PI/3]);  // Tangent frame at (θ, φ)
const metric = frame.metric();  // Metric tensor gᵢⱼ
```

### ConnectionBivector

Computes connection bivectors ωᵢ = ½ eʲ ∧ (∂ᵢeⱼ).

```javascript
const { ConnectionBivector } = require('./src/riemannian-ga.js');

const connection = new ConnectionBivector(sphere);
const omegas = connection.computeAt([Math.PI/4, 0]);  // Array of bivectors
const omega_u = connection.along(coords, [1, 0]);    // ω(∂θ)
```

### Curvature2Form

Computes curvature from connection: Ω = dω + ω∧ω.

```javascript
const { Curvature2Form } = require('./src/riemannian-ga.js');

const curvature = new Curvature2Form(sphere);
const K = curvature.gaussianCurvature([Math.PI/4, 0]);
// K = 1.0 for unit sphere
```

**Methods:**
- `gaussianCurvature(coords)` — Gaussian curvature K
- `scalar()` — Scalar curvature R
- `sectional(u, v, coords)` — Sectional curvature K(u, v)

---

## ProbabilityManifold

Contact geometry of probability distributions (Information Geometry).

**Module:** `src/information-geometry.js`

Based on John Baez's Information Geometry series (Parts 18, 19).

### Constructor

```javascript
new ProbabilityManifold(n)
```

**Parameters:**
- `n` *(number)* — Number of microstates (dimension of probability vector)

### Coordinates

The extended phase space has dimension 2n + 1:
- q^i: Probabilities (analogous to position)
- p_i: Surprisals (conjugate variables, analogous to momentum)
- S: Shannon entropy (analogous to action)

### Methods

#### `entropy(q)`

Compute Shannon entropy S(q) = −Σ qᵢ ln(qᵢ).

```javascript
const manifold = new ProbabilityManifold(2);
const S = manifold.entropy([0.3, 0.7]);
// S ≈ 0.6109
```

#### `surprisal(q)`

Compute surprisal pᵢ = −ln(qᵢ) − 1.

```javascript
const p = manifold.surprisal([0.3, 0.7]);
// p = [0.204, -0.643]
```

#### `contactForm(q, p)`

Construct contact form α = dS − pᵢ dqⁱ.

```javascript
const alpha = manifold.contactForm(q, p);  // Returns Multivector
```

#### `volumeForm(q, p)`

Compute α ∧ (dα)ⁿ (non-zero verifies contact structure).

```javascript
const vol = manifold.volumeForm(q, p);
console.log(vol.norm());  // Non-zero = valid contact manifold
```

#### `checkLegendrianCondition(q)`

Verify that the probability distribution submanifold is Legendrian.

```javascript
const results = manifold.checkLegendrianCondition([0.3, 0.7]);
// results = [{ k: 0, contraction: 0, passed: true }, ...]
```

---

## Discrete Riemannian Geometry

**Module:** `src/riemannian-discrete.js`

### MeshCurvature2Form

Computes discrete curvature on triangle meshes via angle defects.

```javascript
const { MeshCurvature2Form } = require('./src/riemannian-discrete.js');
const { TriangleMesh } = require('./src/mesh.js');

const mesh = TriangleMesh.createIcosahedron(1.0);
const curvature = new MeshCurvature2Form(mesh);

const defects = curvature.angleDefects();  // K_v at each vertex
const totalK = curvature.totalCurvature(); // Σ K_v = 2πχ
const chi = curvature.eulerCharacteristic(); // χ = V - E + F
```

### MeshConnectionBivector

Computes connection bivectors at mesh edges (dihedral rotations).

```javascript
const { MeshConnectionBivector } = require('./src/riemannian-discrete.js');

const connection = new MeshConnectionBivector(mesh);
const omegas = connection.compute();  // Bivector3D[] per edge
const omega_e = connection.at(edgeIdx);
```

---

## EntropicGravity (NEW)

The `entropic-gravity` module implements Bianconi's "Gravity from Entropy" framework.

### MatterInducedMetric

Defines the metric $G_{\mu\nu}$ induced by topological matter fields.

#### Constructor

```javascript
new EntropicGravity.MatterInducedMetric(options)
```

**Parameters:**
- `options.scalarField` *(Function)* — $\phi(x) \to number$
- `options.vectorPotential` *(Function)* — $A_\mu(x) \to [A_0, A_1, A_2, A_3]$
- `options.twoFormPotential` *(Function)* — $B_{\mu\nu}(x) \to 4 \times 4$ matrix

#### Methods

- `covariant(x, gInv)`: Compute $G_{\mu\nu}$ at $x$.
- `contravariant(x)`: Compute $G^{\mu\nu}$ at $x$.

### TwoMetricSystem

Manages the spacetime metric $g$ and matter metric $G$.

#### Constructor

```javascript
new EntropicGravity.TwoMetricSystem(g, G)
```

**Parameters:**
- `g` *(SpacetimeMetric)* — Spacetime metric
- `G` *(MatterInducedMetric)* — Matter-induced metric

#### Methods

- `spacetimeMetric(x)`: Get $g_{\mu\nu}(x)$.
- `matterMetric(x)`: Get $G_{\mu\nu}(x)$.
- `metricDifference(x)`: Get $G_{\mu\nu} - g_{\mu\nu}$.
- `metricRatio(x)`: Get $G \cdot g^{-1}$.

### RelativeEntropyAction

Computes the quantum relative entropy functional $S(G||g)$.

$$ S(G||g) = \int d^4x \sqrt{|g|} \text{Tr}[G(\ln G - \ln g)] $$

#### Constructor

```javascript
new EntropicGravity.RelativeEntropyAction(twoMetricSystem)
```

#### Methods

- `localDensity(x)`: Compute integrand at point $x$.
- `totalAction(bounds)`: Integrate action over a region.
- `variationWrtMetric(x)`: Compute entropic stress-energy tensor $\delta S/\delta g^{\mu\nu}$.

### EmergentCosmologicalConstant

Computes the emergent $\Lambda_G$ field.

$$ \Lambda_G = \frac{1}{4} \text{Tr}(G \cdot g^{-1} - 4I) $$

#### Constructor

```javascript
new EntropicGravity.EmergentCosmologicalConstant(twoMetricSystem)
```

#### Methods

- `localValue(x)`: Compute $\Lambda_G(x)$.
- `classify(x)`: Returns 'de Sitter', 'anti-de Sitter', or 'Minkowski'.

### EntropicGravityHamiltonian

Combines GMET kinematics with Bianconi dynamics.

$$ H = H_{geo} + \alpha S(G||g) $$

#### Constructor

```javascript
new EntropicGravity.EntropicGravityHamiltonian(manifold, system, options)
```

**Parameters:**
- `manifold` *(GrandContactManifold)* — 13D manifold
- `system` *(TwoMetricSystem)* — The (g, G) system
- `options.mass` *(number)* — Particle mass
- `options.entropicCoupling` *(number)* — Strength $\alpha$

#### Methods

- `evaluate(coords)`: Compute H.
- `flow(pt, dt, steps)`: Integrate trajectory.

---

## Discrete Mesh Entropic Gravity

For 2D surfaces represented as triangle meshes, there is a simplified discrete formulation.

### Key Differences from Continuum Module

| Aspect | Continuum (`entropic-gravity.js`) | Discrete (`mesh-entropic-gravity.js`) |
|--------|-----------------------------------|---------------------------------------|
| Geometry | 4D spacetime | 2D surface mesh |
| Metric g | Schwarzschild / FLRW | Implicit in edge lengths |
| Metric G | Full 4×4 with φ, A, B fields | L² = l² + (Δφ)² per edge |
| Entropy | Matrix trace formula | Area ratio: S = R log R - (R-1) |
| Dynamics | Contact Hamiltonian flow | Gradient-based (planned) |

### Discrete Formula

Given a scalar matter field φ on vertices:

1. **Perturbed Edge Lengths**: $L_{ij}^2 = l_{ij}^2 + (\phi_j - \phi_i)^2$
2. **Perturbed Face Areas**: Computed via Heron's formula with $L$ lengths
3. **Relative Entropy Density per Face**:
   $$ s_f = R \ln R - (R - 1), \quad R = \frac{A_G^{(f)}}{A_g^{(f)}} $$
4. **Total Entropy**: $S = \sum_f s_f \cdot A_g^{(f)}$

### Usage Example

```javascript
const Mesh = require('../src/mesh.js');

// Create sphere mesh
const mesh = Mesh.TriangleMesh.createIcosahedron(5.0);

// Define scalar field (Gaussian at north pole)
const phi = new Float64Array(mesh.nVertices);
// ... populate phi ...

// Compute perturbed edge lengths
const g_lengths = mesh.edgeLengths();
const G_lengths = new Float64Array(mesh.nEdges);
for (let e = 0; e < mesh.nEdges; e++) {
    const v0 = mesh.edges[2*e], v1 = mesh.edges[2*e + 1];
    const dphi = phi[v1] - phi[v0];
    G_lengths[e] = Math.sqrt(g_lengths[e]**2 + dphi**2);
}

// Compute entropy density per face using Heron's formula
// ... (see examples/mesh-entropic-gravity.js)
```

