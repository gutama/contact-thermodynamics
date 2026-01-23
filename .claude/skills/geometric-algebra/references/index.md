# Reference Map

Use this folder as a **lookup table**: open the smallest file that answers the user's question, and link out instead of duplicating definitions.

## What to open (fast)

| User asks about… | Open… | Notes |
|---|---|---|
| GA identities, grades, products, signatures | [`formulas.md`](formulas.md), [`algebras.md`](algebras.md) | Core math + signature conventions |
| **B² classification**, elliptic/parabolic/hyperbolic, number systems | [`bivector_classification.md`](bivector_classification.md) | **Master principle: B² determines everything** |
| **CGA versors**, Hybrid CGA, IPNS vs OPNS, complete grade table | [`cga_versors_reference.md`](cga_versors_reference.md) | **Todd's Hybrid CGA, transformation-first view** |
| Exponentials of point pairs, circles, spheres | [`cga_versors_reference.md`](cga_versors_reference.md) | Booster, loxodrome, conformal transforms |
| Mesh primitives in **PGA** (points/lines/planes), simplex measures, volume/area formulas | [`mesh_reference.md`](mesh_reference.md) | Practical mesh formulas used by DDG/geodesics |
| Discrete curvature (angle defect, cotan Laplacian, Steiner/variational curvature, flows) | [`ddg_reference.md`](ddg_reference.md) | Local operators + variational view |
| CGA for mesh curvature (circle fitting, discriminant) | [`cga_versors_reference.md`](cga_versors_reference.md), [`ddg_reference.md`](ddg_reference.md) | Discriminant → surface type |
| Discrete geodesics (shortest vs straightest), exp/log maps, parallel transport, Karcher mean | [`geodesics_reference.md`](geodesics_reference.md) | Algorithms and implementation notes |
| **Discrete geometric calculus**, ∇ operator, PDEs, wave/heat equations | [`geometric_calculus_discrete.md`](geometric_calculus_discrete.md) | Split operator, leapfrog, mesh implementation |
| **Leapfrog + Cotan Laplacian**, wave equations on meshes, DEC | [`leapfrog_cotan_reference.md`](leapfrog_cotan_reference.md) | Cotan formula + leapfrog synergy |
| **FTGC on triangle meshes**, DEC↔GC translation, unified ∇ | [`leapfrog_cotan_reference.md`](leapfrog_cotan_reference.md), [`geometric_calculus_discrete.md`](geometric_calculus_discrete.md) | **Cotan weights = reciprocal basis** |
| Global constraints (Gauss–Bonnet, Poincaré–Hopf, holonomy/topology links) | [`unified_topology_reference.md`](unified_topology_reference.md) | "Why it must be true" checks |
| GA in physics (EM, spacetime, gauge framing) | [`physics.md`](physics.md) | Keep physics discussions out of DDG/geodesics unless needed |
| **Riemannian geometry**, connection bivector, curvature 2-form, coordinate-free | [`riemannian_geometry_reference.md`](riemannian_geometry_reference.md) | **Replaces Christoffel symbols with ω** (NEW) |
| **Christoffel symbols → connection bivector**, Riemann tensor → Ω | [`riemannian_geometry_reference.md`](riemannian_geometry_reference.md) | Master equations for curved spaces (NEW) |
| **Geodesics on continuous manifolds**, holonomy, parallel transport | [`riemannian_geometry_reference.md`](riemannian_geometry_reference.md), [`geodesics_reference.md`](geodesics_reference.md) | ∇ᵥv = 0 formulation |

## House rules for these references

- Prefer **links** over repetition: if a concept is defined elsewhere, point to it.
- When adding new material, place it in **one** reference file and link from the others.
- Keep "Equation Summary" sections to **just the essentials**; defer derivations to `formulas.md`.

## Key Principles (Hybrid CGA)

### The Master Classification

**B² determines transformation type:**

| B² | Type | exp(θB/2) | Representation |
|----|------|-----------|----------------|
| < 0 | Elliptic | cos + sin | IPNS (circles) |
| = 0 | Parabolic | 1 + B | Both (null) |
| > 0 | Hyperbolic | cosh + sinh | OPNS (point pairs) |

### The Discriminant

For pairs of objects (circles, spheres): **Δ = (X₁·X₂)² − X₁²X₂²**

| Δ | Configuration | Transformation |
|---|---------------|----------------|
| < 0 | Intersecting | Rotation (elliptic) |
| = 0 | Tangent | Translation (parabolic) |
| > 0 | Nested/Disjoint | Boost/Spiral (hyperbolic) |

### CGA ≅ Hyperbolic PGA

n-dimensional CGA Cl(n+1,1,0) is isomorphic to (n+1)-dimensional Hyperbolic PGA. This explains:
- Why IPNS and OPNS are both needed
- The geometric meaning of "imaginary" circles (they're real point pairs)
- PGA as a subalgebra (objects through the "eye" point)

## Key Principles (FTGC on Meshes)

### The Split Operator on Meshes

The geometric derivative ∇ = Σᵢ eⁱ ⊗ ∂/∂xᵢ where:
- **eⁱ = reciprocal basis** encoded by cotan weights
- **∂/∂xᵢ = finite differences** along edges

### DEC ↔ FTGC Translation

| DEC | Geometric Calculus | Mesh Implementation |
|-----|-------------------|---------------------|
| 0-form | Scalar (grade 0) | Vertex values |
| 1-form | Vector (grade 1) | Edge values |
| 2-form | Bivector (grade 2) | Face values |
| d (ext. derivative) | ∇∧ (outer) | Incidence matrix |
| δ = ⋆d⋆ (codifferential) | ∇· (inner) | M⁻¹ dᵀ M |
| **⋆₁ (edge Hodge)** | **Reciprocal basis** | **COTAN WEIGHTS** |
| Δ = dδ + δd | ∇² = ∇·∇ | Cotan Laplacian |

### Why Staggering Works (FTGC Explanation)

∇ changes grade: inputs and outputs have different grades → they live on different mesh elements!
- ∇(scalar on vertices) → vector on edges
- ∇(vector on edges) → scalar on vertices + bivector on faces

This is the deep reason Yee's staggered grid and DEC's primal/dual complex work.

## Key Principles (Riemannian Geometry in GA) — NEW

### The Master Equations

**Coordinate-free Riemannian geometry without Christoffel symbols:**

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

### Classical → GA Translation

| Classical | Geometric Algebra | Role |
|-----------|------------------|------|
| gᵢⱼ (metric tensor) | eᵢ · eⱼ | Implicit in frame |
| Γᵏᵢⱼ (Christoffel symbols) | ωᵢ (connection bivector) | Frame rotation |
| Rˡᵢⱼₖ (Riemann tensor) | Ω (curvature 2-form) | Holonomy per area |
| ∇ᵢVʲ (covariant derivative) | ∇ᵤA = ∂ᵤA + ω(u)×A | Multivector extension |
| Geodesic equation | ∇ᵥv = 0 | Autoparallel curve |

### Why Connection Bivector?

**Christoffel symbols are coordinate artifacts.** They exist because ∂ᵢVʲ doesn't transform tensorially.

The connection bivector ω is the **true geometric object**:
- Single bivector replaces n³ Christoffel components
- Action via commutator: ∂ᵢeⱼ = ωᵢ × eⱼ
- Curvature from Cartan: Ω = dω + ω∧ω

### Holonomy and Gauss-Bonnet

**Curvature = holonomy per unit area**

Transport tangent pseudoscalar B around infinitesimal loop:
```
δB = Ω(u,v) × B · (area)
```

**Gauss-Bonnet in GA form:**
```
∫_M Ω + ∮_∂M ω = 2πχ(M) · I
```

The integral of curvature is **topological** (Euler characteristic).

## Implementation Files

| File | Domain | Formulation |
|------|--------|-------------|
| `discrete_geometric_calculus.py` | Cubic grids | FTGC (split operator) |
| `leapfrog_cotan_mesh.py` | Triangle meshes | DEC (d, δ, ⋆) |
| `leapfrog_gc_mesh.py` | Triangle meshes | **FTGC (unified ∇)** |
| `covariant_derivative.py` | Continuous manifolds | Christoffel-based |
| `riemannian_ga.py` | Continuous manifolds | **Connection bivector** (NEW) |

## Common pitfalls (and how to avoid them)

- **Exponentiating without checking B²**
  - *Symptom:* Using cos/sin for a hyperbolic bivector, getting wrong transformation.
  - *Fix:* Always compute B² first. Use cosh/sinh for B² > 0, cos/sin for B² < 0.

- **Assuming all CGA bivectors are circles (pure IPNS thinking)**
  - *Symptom:* Point pairs appearing where circles expected; boosts misinterpreted as rotations.
  - *Fix:* Use **Hybrid CGA**: B² > 0 means the bivector is a point pair (OPNS), not a circle.

- **Meet vs join confusion (∧ vs ∨)**
  - *Symptom:* intersections look like unions; points/lines come out "at infinity".
  - *Fix:* In **PGA**, use **join (∧)** to *span* (point∧point → line, line∧point → plane) and **meet (∨)** to *intersect* (plane∨plane → line, line∨plane → point). If unsure, sanity-check the *grade* of the result.

- **Forgetting to normalize homogeneous PGA objects**
  - *Symptom:* distances/angles drift with scale; comparisons fail.
  - *Fix:* Treat PGA elements as **projective**: normalize consistently before measuring (e.g., normalize planes by Euclidean normal length; normalize points with `w=1` when appropriate).

- **Orientation/sign mismatches (winding, normals, bivector signs)**
  - *Symptom:* curvature flips sign, area becomes negative, geodesics "prefer" the wrong side.
  - *Fix:* Pick one convention (triangle winding, outward normals) and enforce it during mesh import. If you flip winding, flip any derived per-face normals/areas too.

- **Shortest vs straightest geodesics conflation**
  - *Symptom:* path looks locally straight but globally non-minimizing, or vice versa.
  - *Fix:* Use **straightest** for parallel-transport/connection-based methods; use **shortest** for metric minimization. On nonconvex surfaces, they can diverge.

- **Boundary handling omitted**
  - *Symptom:* Laplacians/flows blow up near boundaries; geodesics "leak" off the mesh.
  - *Fix:* Choose boundary conditions explicitly (Dirichlet/Neumann/mixed). For geodesics, decide whether boundaries reflect, clamp, or terminate.

- **Cotan Laplacian surprises on obtuse triangles**
  - *Symptom:* negative weights; diffusion behaves oddly; curvature flow becomes unstable.
  - *Fix:* Check mesh quality. Consider intrinsic Delaunay triangulation / edge flips or weight clamping when stability matters more than strict fidelity.

- **Degenerate / near-degenerate simplices**
  - *Symptom:* NaNs/Infs in curvature, barycentric coordinates, exp/log maps.
  - *Fix:* Detect tiny areas/edge lengths and either repair (remesh) or add eps-guards; never divide by near-zero without a fallback.

- **Unit mismatch (radians vs degrees; squared vs unsquared lengths)**
  - *Symptom:* curvature magnitudes off by ~57× or energy terms dominate unexpectedly.
  - *Fix:* Standardize: angles in **radians**, lengths in one unit system, and be explicit about whether an algorithm expects **ℓ** or **ℓ²**.

- **Assuming "flat-space" intuition for curvature/topology checks**
  - *Symptom:* expecting total Gaussian curvature ≈ 0 on a closed mesh that isn't toroidal.
  - *Fix:* Use global invariants as sanity checks: Gauss–Bonnet (total K relates to Euler characteristic) and Poincaré–Hopf for vector field index sums (see `unified_topology_reference.md`).

- **Numerical drift in rotor/exp-map composition**
  - *Symptom:* parallel transport slowly spirals; rotations gain/lose magnitude.
  - *Fix:* Re-orthonormalize rotors/frames periodically; avoid accumulating tiny floating errors in long transports.

- **Confusing DEC and FTGC sign conventions**
  - *Symptom:* Laplacian eigenvalues have wrong sign; wave equation blows up.
  - *Fix:* The cotan Laplacian is **negative semi-definite** (negative diagonal). Verify ∇²(const) = 0 and check eigenvalue signs.

- **Forgetting FTGC verification checks**
  - *Symptom:* Mysterious numerical errors in mesh calculus.
  - *Fix:* Always verify ∇∧(∇∧f) = 0 (curl of gradient is zero) and ∇²(const) = 0 as sanity checks.

- **Using Christoffel symbols when connection bivector is cleaner** (NEW)
  - *Symptom:* Index explosion, hard-to-debug tensor expressions, sign errors.
  - *Fix:* Use `riemannian_ga.py` with connection bivector ω. The single equation ∂ᵢeⱼ = ωᵢ × eⱼ replaces all Christoffel manipulations.

- **Wrong sign in connection bivector definition** (NEW)
  - *Symptom:* Geodesics curve wrong direction, curvature has wrong sign.
  - *Fix:* The correct definition is **ωᵢ = ½ eʲ ∧ (∂ᵢeⱼ)** with reciprocal frame first. This is **not** the same as ½ (∂ᵢeⱼ) ∧ eʲ (differs by sign).

- **Forgetting metric compatibility check** (NEW)
  - *Symptom:* Inner products not preserved under parallel transport.
  - *Fix:* Verify ∇(eᵢ · eⱼ) = 0, equivalently ω × (eᵢ · eⱼ) = 0.

- **Not verifying holonomy matches area integral** (NEW)
  - *Symptom:* Parallel transport around loop gives unexpected rotation.
  - *Fix:* Check that holonomy angle ≈ ∫∫ K dA (integrated Gaussian curvature) for small loops.
