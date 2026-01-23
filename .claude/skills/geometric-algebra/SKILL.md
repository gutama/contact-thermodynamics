---
name: geometric-algebra
description: "Geometric Algebra / Clifford algebra help: multivectors, geometric/outer/inner products, grade projection & duality; rotors/motors (PGA) for 2D/3D points-lines-planes and rigid motions; CGA for circles/spheres and conformal transforms; Hybrid CGA versors and exponentials; Versor-style operations (Boosters for conformal bending, Tangent/Direction elements, extended Round/Flat ops); arbitrary metric algebra generation; DDG on triangle meshes (curvature, Steiner, flows) and discrete geodesics (exp/log, parallel transport, cut locus, Karcher mean); **discrete geometric calculus** (split differential operator, finite differences, leapfrog time-stepping, DEC for meshes, **FTGC on triangle meshes**); **Riemannian geometry in GA** (connection bivector ω, curvature 2-form Ω, coordinate-free geodesics, holonomy, Gauss-Bonnet). Includes JS/Python implementations and interactive demos."
---

# Geometric Algebra

Use this Skill when the user asks for **GA/Clifford computations**, **PGA/CGA modeling**, **rotors/motors**, **CGA versors/exponentials**, **mesh/geometry tasks** (DDG, curvature, geodesics, topology links), **discrete geometric calculus** (∇ operator, PDEs, wave equations, DEC, **FTGC**), or **Riemannian geometry** (connection, curvature, geodesics without Christoffel symbols).

## Workflow (do this in order)

1. **Identify the algebra** (and state it explicitly):
   - Euclidean rotors: `Cl(3,0,0)`
   - **PGA 2D/3D**: `Cl(2,0,1)` / `Cl(3,0,1)`
   - **CGA 2D/3D**: `Cl(3,1,0)` / `Cl(4,1,0)`
   - Spacetime: `Cl(1,3,0)`
   - **Riemannian manifold**: Frame-based (no fixed signature)
   - If unsure, ask one focused question (dimension + geometry type).

2. **Pick an implementation path** (optimize for reliability):
   - **Sandbox / offline / strict CSP** → use the bundled **vanilla JS** (`scripts/ga_vanilla.js`).
   - **Python numeric/symbolic mix** → prefer **kingdon**.
   - **Pure symbolic + LaTeX** → prefer **galgebra/SymPy**.
   - **Browser visualization** → use bundled HTML demos; ganja.js is optional and may fail without network.
   - **Discrete geometric calculus / PDEs** → use discrete implementation or **numga** (JAX).
   - **FTGC on triangle meshes** → use `scripts/leapfrog_gc_mesh.py` (unified ∇).
   - **Riemannian geometry (coordinate-free)** → use `scripts/riemannian_ga.py` (connection bivector).
   - **Riemannian geometry (coordinate-based)** → use `scripts/covariant_derivative.py` (Christoffel symbols).

3. **Solve + validate**:
   - Show the key identities used (not a full derivation unless asked).
   - Provide at least one sanity check (units, invariants, limiting case, symmetry).
   - If doing mesh algorithms, report assumptions (manifoldness, boundary, metric, approximation quality).
   - For Riemannian: verify Bianchi identity ∇∧Ω = 0, metric compatibility, geodesic arc-length preservation.

4. **Present**:
   - Default to concise math + a minimal working code snippet.
   - If the user wants a figure/animation, point them to the included demos.

## Conventions (keep consistent)

- **Geometric product**: `ab`  (mixes inner + outer)
- **Outer/wedge**: `a ∧ b`  (joins / spans)
- **Inner**: `a · b`  (projection)
- **Reversal**: `~A`  (reverse / "tilde")
- **Sandwich action** (rotors/motors): `X' = R X ~R`
- **Geometric derivative**: `∇F = Σᵢ eⁱ (∂F/∂xᵢ)` (split operator form)
- **Commutator product**: `A × B = ½(AB − BA)` (bivector rotation action)

For detailed formula tables (products, duality, rotor↔quaternion, CGA primitives), see:
- `references/index.md` — what to open (fast).
- `references/formulas.md`
- `references/algebras.md`
- `references/bivector_classification.md` — **B² determines transformation type**
- `references/cga_versors_reference.md` — **Hybrid CGA, complete versor table**
- `references/geometric_calculus_discrete.md` — **Split operator, mesh ∇, PDEs, FTGC**
- `references/leapfrog_cotan_reference.md` — **Leapfrog + cotan + FTGC synergy**
- `references/riemannian_geometry_reference.md` — **Coordinate-free Riemannian geometry** (NEW)

## What's bundled (navigation)

### Core code
- `scripts/ga_vanilla.js` — **no-dependency** GA toolkit (PGA/CGA basics, rotors/motors utilities).
- `scripts/ga_discrete_curvature.js` — DDG on triangle meshes (curvature, Steiner coefficients, flows).
- `scripts/ga_geodesics.js` — discrete geodesics (shortest/straightest, exp/log maps, transport, cut locus, Karcher mean).
- `scripts/ga_mesh.js` — PGA-style mesh measurements (length/area/volume/centroid, etc.).
- `scripts/ga_number_systems.js` — Complex/Dual/Hyperbolic number operations, bivector classification.
- `scripts/ga_utils.py` — small Python helpers.
- `scripts/discrete_geometric_calculus.py` — **Discrete ∇ operator, split form, leapfrog** (cubic grids)
- `scripts/leapfrog_cotan_mesh.py` — **Leapfrog + cotan Laplacian** (DEC formulation)
- `scripts/leapfrog_gc_mesh.py` — **FTGC on triangle meshes** (unified ∇, MeshGeometricDerivative)
- `scripts/covariant_derivative.py` — **Christoffel symbols, covariant derivative, geodesics on continuous manifolds**
- `scripts/riemannian_ga.py` — **Coordinate-free Riemannian geometry** (connection bivector, curvature 2-form)
- `src/riemannian-discrete.js` — **Discrete Riemannian on meshes** (dihedral connection, angle defect curvature, Bianchi) (NEW)

### Interactive demos (offline-friendly)
- `scripts/mesh_demo.html`
- `scripts/geodesics_demo.html`
- `scripts/geodesics_topology_demo.html` (geodesics ↔ Gauss–Bonnet ↔ Poincaré–Hopf)
- `scripts/animation_standalone.html`
- `scripts/cga_curvature_test.html`
- `scripts/cga_exponentials_demo.html` — **Point pair, circle pair, sphere pair exponentials**
- `scripts/cga_mesh_curvature.html` — **CGA circle fitting for mesh curvature**
- `scripts/number_systems_demo.html`

### Deep references (load only when needed)
- `references/ddg_reference.md` — integral + variational DDG.
- `references/geodesics_reference.md` — shortest vs straightest, exp/log, transport.
- `references/mesh_reference.md` — mesh formulas & implementation notes.
- `references/unified_topology_reference.md` — unified picture (curvature/topology/holonomy).
- `references/physics.md` — GA in physics (EM, spacetime).
- `references/cga_versors_reference.md` — **Hybrid CGA, Todd's versor classification**.
- `references/bivector_classification.md` — **B² classification, number systems, discriminant**.
- `references/geometric_calculus_discrete.md` — **Split differential operator, mesh ∇, PDEs, FTGC**
- `references/leapfrog_cotan_reference.md` — **Leapfrog + cotan + FTGC synergy**
- `references/riemannian_geometry_reference.md` — **Coordinate-free Riemannian geometry, connection bivector** (NEW)

## Decision guide by task

### A) "Compute / simplify a GA expression"
1) Confirm the algebra/signature. 2) State basis + grade selection. 3) Compute. 4) Provide a check.

### B) "Rigid motion / geometry (points, lines, planes)" (PGA)
Use `Cl(3,0,1)`; represent primitives, build a rotor/translator/motor, apply via sandwich.

### C) "Circles / spheres / conformal transforms" (CGA)
Use `Cl(4,1,0)`; embed points, build rounds/flats via outer products; apply conformal versors.

### D) "What transformation does this bivector generate?" (Hybrid CGA)
**Check B² first!**
- B² < 0 → Elliptic (rotation), invariant is circle/sphere (IPNS)
- B² = 0 → Parabolic (translation), null object
- B² > 0 → Hyperbolic (boost), invariant is **point pair** (OPNS)

Load `references/bivector_classification.md` or `references/cga_versors_reference.md`.

### E) "Exponential of point pair / circle pair / sphere pair"
- **Point pair** (PP² > 0): exp(PP) → Booster (hyperbolic flow between points)
- **Circle pair**: Form B = C₁C₂⁻¹, discriminant determines elliptic/parabolic/hyperbolic
- **Sphere pair**: Form B = S₁S₂⁻¹, gives dilation/rotation/transversion/boost

Use demos: `scripts/cga_exponentials_demo.html`

### F) "Curvature / smoothing / energies on a triangle mesh" (DDG)
Use `scripts/ga_discrete_curvature.js`. Prefer integral quantities first (angle defect / cotan Laplacian), then variational (Steiner/flows). Load `references/ddg_reference.md` for formulas.

### G) "CGA for mesh curvature estimation"
Fit circles through vertex neighbors: C = P₁ ∧ P₂ ∧ P₃, κ = 1/radius(C).
The **discriminant** of fitted circles determines surface type (elliptic/parabolic/hyperbolic = convex/flat/saddle).
Use demo: `scripts/cga_mesh_curvature.html`

### H) "Geodesics / exp-log / parallel transport"
**Discrete (meshes):** Use `scripts/ga_geodesics.js`. Clarify **shortest** (boundary value) or **straightest** (initial value). Load `references/geodesics_reference.md`.

**Continuous (analytic manifolds):** Use `scripts/covariant_derivative.py` (Christoffel-based) or `scripts/riemannian_ga.py` (coordinate-free). Provides `GeodesicSolver` for sphere, torus, hyperbolic plane.

### I) "Covariant derivative / connection / parallel transport on continuous manifolds"

**Two approaches available:**

**Coordinate-based** (traditional): Use `scripts/covariant_derivative.py`. Key classes:
- `ChristoffelSymbols(manifold)` — compute Γᵏᵢⱼ from metric
- `CovariantDerivative(manifold)` — ∇ᵢVʲ, divergence, Laplacian
- `ParallelTransport(manifold)` — transport vectors along curves
- `GeodesicSolver(manifold)` — solve d²xᵏ/dt² + Γᵏᵢⱼ vⁱ vʲ = 0

**Coordinate-free** (pure GA): Use `scripts/riemannian_ga.py`. Key classes:
- `ConnectionBivector(manifold)` — compute ωᵢ (replaces Christoffel)
- `Curvature2Form(manifold)` — compute Ω = dω + ω∧ω
- `GACovariantDerivative(manifold)` — ∇ᵤA = ∂ᵤA + ω(u) × A
- `GAGeodesicSolver(manifold)` — solve ∇ᵥv = 0
- `GAParallelTransport(manifold)` — holonomy computation

### J) "What type of transformation? / Classify a bivector"
Use `scripts/ga_number_systems.js`. Check B² to determine elliptic (−1), parabolic (0), or hyperbolic (+1). Load `references/bivector_classification.md` for theory.

### K) "Complex/Dual/Hyperbolic number operations"
Use `scripts/ga_number_systems.js`. Each class provides arithmetic, exp, and geometric utilities. Dual numbers support automatic differentiation.

### L) "IPNS vs OPNS — which should I use?"
**Use Hybrid CGA**: switch based on B².
- Negative-square objects → IPNS (circles, spheres)
- Positive-square objects → OPNS (point pairs)
- Null objects → both representations coincide

Load `references/cga_versors_reference.md` for the complete classification.

### M) "Discrete geometric derivative / PDEs / Wave equations"

Use the **split differential operator** concept:
```
∇F = Σᵢ eⁱ (∂F/∂xᵢ)
```

**Implementation choice by domain:**

| Domain | Approach | Key Feature |
|--------|----------|-------------|
| Flat cubic grid | Finite differences | Simple, fast, O(h²) |
| Curved (regular param) | Position-dependent metric | g_ij(x) field |
| Unstructured mesh (DEC) | `leapfrog_cotan_mesh.py` | DEC formulation |
| **Unstructured mesh (FTGC)** | `leapfrog_gc_mesh.py` | **Unified ∇ operator** |
| Symbolic | galgebra | Exact expressions |

**For spacetime equations** (Dirac, Maxwell, wave):
Use leapfrog time-stepping:
```
ψ(t+Δt) = ψ(t-Δt) - 2Δt·∇_spatial·ψ(t)
```

Load `references/geometric_calculus_discrete.md` for formulas and code.

### N) "FTGC on triangle meshes / Refactor DEC to geometric calculus"

Use `scripts/leapfrog_gc_mesh.py`. Key classes:

- **MeshGeometricDerivative** — unified ∇ with `.wedge()`, `.inner()`, `.grad()`, `.div()`, `.curl()`, `.laplacian()`
- **MeshMultivectorField** — staggered storage (grade 0 on vertices, grade 1 on edges, grade 2 on faces)
- **LeapfrogGC** — wave/heat/Maxwell solvers using FTGC
- **FTGCVerification** — verify ∇∧∇∧ = 0, ∇²(const) = 0

**Key insight:** Cotan weights = reciprocal basis vectors eⁱ (mesh metric).

Load `references/leapfrog_cotan_reference.md` for DEC↔FTGC translation.

### O) "Riemannian geometry / Christoffel symbols / curvature tensor" (NEW)

**Use coordinate-free GA formulation.** Christoffel symbols are coordinate artifacts — the connection bivector ω is the true geometric object.

**Master Equations:**
```
Metric:        g_ij = e_i · e_j,    e^i · e_j = δ^i_j
Connection:    ω_i = ½ e^j ∧ (∂_i e_j)
Frame:         ∂_i e_j = ω_i × e_j
Covariant:     ∇_u A = ∂_u A + ω(u) × A
Curvature:     Ω = dω + ω ∧ ω
Riemann:       [∇_u, ∇_v]w = Ω(u,v) × w
Geodesic:      ∇_v v = 0
Bianchi:       ∇ ∧ Ω = 0
```

**Translation table:**

| Classical | GA Replacement | Notes |
|-----------|---------------|-------|
| Γᵏᵢⱼ (Christoffel) | ωᵢ (connection bivector) | ω encodes frame rotation |
| Rˡᵢⱼₖ (Riemann) | Ω (curvature 2-form) | Ω = dω + ω∧ω |
| gᵢⱼ (metric) | eᵢ · eⱼ (frame inner products) | Metric implicit in frame |
| ∇ᵢVʲ | ∇ᵤv = ∂ᵤv + ω(u) × v | Commutator product rotates |

Use `scripts/riemannian_ga.py` for implementation. Load `references/riemannian_geometry_reference.md` for complete formulas.

### P) "Gaussian curvature from metric / Theorema Egregium" (NEW)

Gaussian curvature K depends only on the metric (Gauss's Theorema Egregium).

**GA approach:** K measures holonomy of tangent pseudoscalar B:
```
K · B = (∂₁n) ∧ (∂₂n)
```

For shape operator S(v) = −∇ᵥn:
```
K = det(S) = κ₁κ₂
```

**Special cases:**
- Sphere (radius R): K = 1/R² (constant positive)
- Torus: K = cos(θ) / (r(R + r cos θ)) (varies)
- Hyperbolic plane: K = −1 (constant negative)
- Isothermal coords: K = −e⁻²ᵠ∇²φ

### Q) "Shape operator / principal curvatures / pseudoscalar" (NEW)

The tangent pseudoscalar B = e₁ ∧ e₂ and normal n are dual:
```
n = −BI⁻¹     ⟺     B = nI
```

Shape operator from pseudoscalar variation:
```
S(v) = (v · ∇B) · I⁻¹
```

**Higher dimensions:** For k-manifold with tangent k-blade Bₖ:
```
∇_v B_k = ω(v) × B_k + h(v) ∧ B_k
```
- ω(v) × Bₖ: intrinsic rotation (parallel transport)
- h(v) ∧ Bₖ: extrinsic bending (second fundamental form)

## Code snippets (copy-paste ready)

### PGA 3D (Cl(3,0,1) via ga_vanilla.js)
```javascript
const GA = window.GA_Vanilla;
const pga3d = GA.PGA3D;

// Create point and plane
const p = pga3d.point(1, 2, 3);
const plane = pga3d.plane(0, 0, 1, -2);

// Reflect point through plane
const reflected = pga3d.sandwich(plane, p);

// Distance from point to plane
const dist = pga3d.distancePointPlane(p, plane);
```

### CGA 3D + Hybrid exp (Cl(4,1,0) via ga_vanilla.js)
```javascript
const cga3d = GA.CGA3D;

// Create circle through three points
const P1 = cga3d.point(1, 0, 0);
const P2 = cga3d.point(0, 1, 0);
const P3 = cga3d.point(-1, 0, 0);
const circle = cga3d.wedge3(P1, P2, P3);

// Create point pair and check sign
const PP = cga3d.wedge(P1, P2);
const B2 = cga3d.normSquared(PP);

// Hybrid exponential (works for any sign of B²)
function hybrid_exp(B, t) {
    const B2 = cga3d.normSquared(B);
    const absB = Math.sqrt(Math.abs(B2));
    if (Math.abs(B2) < 1e-10) {
        return cga3d.add(cga3d.scalar(1), cga3d.scale(B, t));  // Parabolic
    } else if (B2 < 0) {
        return cga3d.add(cga3d.scalar(Math.cos(absB*t)), 
                         cga3d.scale(B, Math.sin(absB*t)/absB));  // Elliptic
    } else {
        return cga3d.add(cga3d.scalar(Math.cosh(absB*t)), 
                         cga3d.scale(B, Math.sinh(absB*t)/absB));  // Hyperbolic
    }
}
```

### Number Systems (Bivector Classification)
```javascript
const NS = GA.NumberSystems;

// Classify a bivector
const type = NS.classifyBivectorSquare(-1);  // 'elliptic'

// Auto-differentiation with dual numbers
const f = x => x.mul(x);  // f(x) = x²
const result = f(new NS.Dual(3, 1));
// result.a = 9 (function value)
// result.b = 6 (derivative at x=3)

// Circle pair classification
const info = NS.classifyCircles2D(0, 0, 1, 1.5, 0, 0.8);
// info.type = 'elliptic' (intersecting)
```

### Discrete Geometric Derivative (Cubic Grid)
```python
# Requires: discrete_geometric_calculus.py
import numpy as np
from discrete_geometric_calculus import (
    CliffordBasis, SplitDifferentialOperator,
    EUCLIDEAN_2D, scalar_field
)

# Setup 2D Euclidean GA on 64x64 grid
basis = CliffordBasis(EUCLIDEAN_2D)
nabla = SplitDifferentialOperator(basis, spacing=0.1)

# Create scalar field f(x,y) = sin(x)cos(y)
x = np.linspace(0, 2*np.pi, 64)
y = np.linspace(0, 2*np.pi, 64)
X, Y = np.meshgrid(x, y, indexing='ij')
f = np.sin(X) * np.cos(Y)
field = scalar_field(basis, (64, 64), f)

# Compute geometric derivative (gives vector field)
grad_f = nabla.grad(field)

# Compute Laplacian
lap_f = nabla.laplacian(field)
```

### FTGC on Triangle Meshes
```python
# Requires: leapfrog_gc_mesh.py
from leapfrog_gc_mesh import (
    TriangleMesh, MeshGeometricDerivative, LeapfrogGC,
    create_grid_mesh, create_disk_mesh
)
import numpy as np

# Create mesh
mesh = create_grid_mesh(nx=32, ny=32, Lx=1.0, Ly=1.0)

# The unified geometric derivative ∇
nabla = MeshGeometricDerivative(mesh)

# Gradient of scalar field
f = mesh.vertices[:, 0]  # f = x
grad_f = nabla.grad(f)   # Vector field on edges

# Laplacian via FTGC: ∇² = ∇·∇
lap_f = nabla.laplacian(f)

# Wave equation solver
solver = LeapfrogGC(mesh)
u0 = np.exp(-50 * ((mesh.vertices[:, 0] - 0.5)**2 + 
                    (mesh.vertices[:, 1] - 0.5)**2))
v0 = np.zeros_like(u0)
dt = solver.estimate_cfl(c=1.0)
u_final = solver.wave_simulate(u0, v0, dt, n_steps=100, c=1.0)
```

### Riemannian Geometry (Coordinate-Free) (NEW)
```python
# Requires: riemannian_ga.py
from riemannian_ga import (
    Sphere2D, Torus2D, HyperbolicPlane,
    ConnectionBivector, Curvature2Form,
    GAGeodesicSolver, GAParallelTransport
)
import numpy as np

# Create sphere manifold
sphere = Sphere2D(R=1.0)

# Connection bivector (replaces Christoffel symbols)
conn = ConnectionBivector(sphere)
coords = np.array([np.pi/4, 0.0])  # θ = 45°
omega = conn.compute_at(coords)    # List of bivectors [ω_θ, ω_φ]

# Curvature via Cartan structure equation: Ω = dω + ω∧ω
curvature = Curvature2Form(sphere)
K = curvature.gaussian_curvature(coords)
print(f"K = {K:.4f} (expected: 1.0)")

# Geodesic solver: ∇_v v = 0
solver = GAGeodesicSolver(sphere)
x0 = np.array([np.pi/4, 0.0])
v0 = np.array([1.0, 0.0])  # Along meridian
t_vals, x_vals, v_vals = solver.solve(x0, v0, t_final=np.pi/2)

# Parallel transport and holonomy
transport = GAParallelTransport(sphere)
def latitude_loop(t):
    return np.array([np.pi/4, t])  # Circle at θ = π/4

w0 = np.array([1.0, 0.0])
holonomy = transport.holonomy_angle(latitude_loop, w0)
expected = 2 * np.pi * (1 - np.cos(np.pi/4))
print(f"Holonomy: {np.degrees(holonomy):.1f}° (expected: {np.degrees(expected):.1f}°)")
```

### Maxwell's Equations via FTGC
```python
# E on edges (grade 1), B on faces (grade 2)
from leapfrog_gc_mesh import LeapfrogGC, create_grid_mesh

mesh = create_grid_mesh(nx=20, ny=20)
solver = LeapfrogGC(mesh)
nabla = solver.nabla

# Initialize fields
E = np.zeros(mesh.n_edges)
B = np.zeros(mesh.n_triangles)

# Yee-style update via FTGC
dt, c = 0.01, 1.0
for step in range(100):
    # Faraday: ∂B/∂t = -∇∧E
    B = B - dt * nabla.curl(E)
    
    # Ampère: ∂E/∂t = c²∇·B  
    E = E + dt * c**2 * nabla._metric_edges_inv @ (
        nabla._wedge_1to2.T @ (nabla._metric_faces @ B)
    )
```

### galgebra Symbolic Gradient
```python
# Symbolic geometric derivative with galgebra
from sympy import symbols, sin, cos
from galgebra.ga import Ga

coords = (x, y) = symbols('x y', real=True)
ga = Ga('e', g=[1, 1], coords=coords)
grad, rgrad = ga.grads()  # Split operator: grad = Σᵢ eⁱ ∂ᵢ

f = sin(x) * cos(y)
F = ga.mv(f)
grad_F = grad * F  # Symbolic gradient
# Result: cos(x)cos(y) e_x - sin(x)sin(y) e_y
```

## Quality bar (avoid common failures)

- Always specify **(p,q,r)** when talking about signatures; don't assume.
- Distinguish **meet** vs **join** (PGA) and **round** vs **flat** (CGA).
- **Check B² before exponentiating** — it determines the formula!
- **`algebra.exp()` is elliptic-only** (cos/sin). For hyperbolic (B²>0) use `hybrid_exp()` above with cosh/sinh.
- For CGA, remember: **Hybrid CGA** = IPNS for B² < 0, OPNS for B² > 0.
- When demonstrating code, prefer the bundled scripts so examples work offline.
- For meshes: state whether results are **approximate** (Dijkstra) or geometric (straightest tracing).
- For discrete calculus: **state grid type** (cubic vs curved vs mesh) and **accuracy order** (O(h²) etc.).
- For leapfrog: **check CFL condition** (c·Δt/h < 1 for stability).
- For FTGC on meshes: verify **∇∧∇∧ = 0** and **∇²(const) = 0** as sanity checks.
- For Riemannian geometry: **use connection bivector ω, not Christoffel symbols**.
- For Riemannian: verify **Bianchi identity** ∇∧Ω = 0, check **geodesic arc-length preservation**.
- For parallel transport: verify **metric norm preservation** (transported vector length constant).

## Key Insight Summary

**The discriminant Δ = (X₁·X₂)² − X₁²X₂² determines transformation type:**

| Primitive Pair | Δ < 0 (Elliptic) | Δ = 0 (Parabolic) | Δ > 0 (Hyperbolic) |
|----------------|------------------|-------------------|---------------------|
| Circle pair | Rotation | Parabolic slide | Spiral |
| Sphere pair | Rotation | Transversion | Dilation/Boost |
| Point pair | — | Null (translation) | Booster |

**CGA ≅ Hyperbolic PGA (one dimension up)** — this explains the IPNS/OPNS duality geometrically.

**The Split Differential Operator**:
```
∇ = ∇_G ⊗ ∇_D = Σᵢ eⁱ ⊗ ∂/∂xᵢ
```

| Implementation | ∇_D (Differential) | ∇_G (Geometric) |
|----------------|-------------------|-----------------|
| galgebra | SymPy `diff()` | Symbolic eⁱ |
| Discrete (Cubic) | Finite differences | Array indexing |
| Mesh (Unstructured) | Cotan weights | Metric encoding |

**Grade change under ∇**: grade k → grades (k-1) + (k+1)
- Scalar → Vector (gradient)
- Vector → Scalar + Bivector (div + curl)

**FTGC on Meshes**:
```
Cotan weights = Reciprocal basis eⁱ = Mesh metric gⁱʲ
```

| DEC | FTGC | Mesh Element |
|-----|------|--------------|
| d (ext. derivative) | ∇∧ (outer) | Incidence matrix |
| δ = ⋆d⋆ (codifferential) | ∇· (inner) | M⁻¹ dᵀ M |
| ⋆₁ (edge Hodge) | Reciprocal basis | **COTAN WEIGHTS** |
| Δ = dδ + δd | ∇² = ∇·∇ | Cotan Laplacian |

**Riemannian Geometry in GA** (NEW):
```
Christoffel symbols Γᵏᵢⱼ  →  Connection bivector ω
Riemann tensor Rˡᵢⱼₖ     →  Curvature 2-form Ω
Metric tensor gᵢⱼ        →  Frame inner products eᵢ · eⱼ
```

**The Master Equations:**
```
Metric:        g_ij = e_i · e_j,    e^i · e_j = δ^i_j
Connection:    ω_i = ½ e^j ∧ (∂_i e_j)
Frame:         ∂_i e_j = ω_i × e_j
Covariant:     ∇_u A = ∂_u A + ω(u) × A
Curvature:     Ω = dω + ω ∧ ω
Riemann:       [∇_u, ∇_v]w = Ω(u,v) × w
Geodesic:      ∇_v v = 0
Bianchi:       ∇ ∧ Ω = 0
```

**Fundamental principle:**
```
Curvature = Derivative of Tangent Pseudoscalar
∇ ∧ ∇ B_k = Ω · B_k
```

**Holonomy interpretation:** Transport tangent pseudoscalar B around loop → rotation by ∫∫ K dA (Gauss-Bonnet).
