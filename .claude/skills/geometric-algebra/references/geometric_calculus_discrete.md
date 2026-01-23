# Discrete Geometric Calculus Reference

## Overview

Discrete Geometric Calculus implements the **geometric derivative** (∇ operator) on discrete domains. This reference covers:
1. The split differential operator concept
2. Cubic grid implementation
3. Curved manifolds
4. Mesh implementation (with DEC translation table)
5. **FTGC on meshes** (Fundamental Theorem of Geometric Calculus)
6. Connection to galgebra

## The Split Differential Operator

### Core Concept

The geometric derivative splits into geometric and differential parts:

```
∇ = ∇_G ⊗ ∇_D
```

Where:
- **∇_G**: Reciprocal basis vectors (eⁱ)
- **∇_D**: Partial derivative operators (∂/∂xᵢ)

Full form:
```
∇F = Σᵢ eⁱ (∂F/∂xᵢ)
```

### Reciprocal Basis

For basis vectors {e₁, ..., eₙ}, the reciprocal basis {e¹, ..., eⁿ} satisfies:
```
eⁱ · eⱼ = δⁱⱼ
```

Construction:
```
eⁱ = (-1)^(i-1) (e₁ ∧ ... ∧ êᵢ ∧ ... ∧ eₙ) / Eₙ
```
where Eₙ = e₁ ∧ ... ∧ eₙ is the pseudoscalar.

For orthonormal bases: **eⁱ = g^{ii} eᵢ** (metric-weighted)

### The Unified Operator

The geometric derivative ∇ combines what in other formalisms are separate operators:

| Operation | Geometric Calculus | Result |
|-----------|-------------------|--------|
| Gradient | ∇f (scalar input) | Vector (grade 1) |
| Divergence | ∇·V = ⟨∇V⟩₀ | Scalar (grade 0) |
| Curl | ∇∧V = ⟨∇V⟩₂ | Bivector (grade 2) |
| Full derivative | ∇V | Mixed: divergence + curl |

**Key identity:** ∇V = ∇·V + ∇∧V (inner + outer = full geometric product)

## Cubic Grid Implementation

### Finite Difference Schemes

| Scheme | Formula | Order |
|--------|---------|-------|
| Forward | (f[i+1] - f[i]) / h | O(h) |
| Backward | (f[i] - f[i-1]) / h | O(h) |
| Central | (f[i+1] - f[i-1]) / (2h) | O(h²) |
| 4th order | (-f[i+2] + 8f[i+1] - 8f[i-1] + f[i-2]) / (12h) | O(h⁴) |

### Algorithm: Discrete Geometric Derivative

```python
def geometric_derivative(field, basis, spacing):
    result = zeros_like(field)
    for i in range(dim):
        # ∇_D part: finite difference
        partial_i = (field[..., i+1, ...] - field[..., i-1, ...]) / (2*h)
        
        # ∇_G part: multiply by reciprocal basis
        for each blade in field:
            sign, new_blade = basis_product(e_i, blade)
            result[new_blade] += metric[i,i] * sign * partial_i[blade]
    
    return result
```

### Grade Changes Under ∇

| Input Grade | Output Grades | Interpretation |
|-------------|---------------|----------------|
| 0 (scalar) | 1 (vector) | Gradient |
| 1 (vector) | 0 + 2 | Divergence + Curl |
| 2 (bivector) | 1 + 3 | Vector + trivector |
| k | (k-1) + (k+1) | General rule |

## Curved Manifolds

### Position-Dependent Metric

For curved spaces, the metric becomes a field:
```
g_ij = g_ij(x)
```

The geometric derivative becomes:
```
∇F = Σᵢⱼ g^{ij}(x) eⱼ (∂F/∂uⁱ)
```

> [!NOTE]
> **This is the coordinate-basis divergence formula, not a full covariant derivative.**
> For vector/multivector fields (not just scalars), a complete geometric covariant derivative
> requires connection terms (Christoffel symbols). The formulas below are valid for:
> - Scalars (gradient, Laplacian)
> - Divergence of vectors (using the √g formula which absorbs connection terms)

### Key Operators

**Gradient** (scalar → vector):
```
(∇f)^i = g^{ij} ∂ⱼf
```

**Divergence** (vector → scalar):
```
∇·V = (1/√g) ∂ᵢ(√g V^i)
```

**Laplacian** (scalar → scalar):
```
∇²f = ∇·(∇f) = (1/√g) ∂ᵢ(√g g^{ij} ∂ⱼf)
```

### Example: Sphere

Coordinates (θ, φ), metric:
```
g = R² [[1, 0], [0, sin²θ]]
```

Laplacian of f = cos(θ):
```
∇²f = -2cos(θ)/R²  (eigenvalue -2/R²)
```

## Mesh Implementation

For unstructured meshes (triangles, tetrahedra), the geometric derivative is implemented using staggered storage where different grades live on different mesh elements.

### Staggered Storage

| Grade | Mesh Location | Example |
|-------|--------------|---------|
| 0 (scalar) | Vertices | Scalar field values |
| 1 (vector) | Edges | Tangential components |
| 2 (bivector) | Faces | Area-weighted values |

### Discrete Gradient (∇ on scalars)

For a scalar field φ on vertices, the gradient along edge (i,j) is:
```
(∇φ)[ij] = φ[j] - φ[i]
```

> [!IMPORTANT]
> **This is a 1-form (edge cochain), not a metric gradient vector.**
> To get the actual gradient vector at a vertex, you must:
> 1. Apply the metric (cotan weights) to convert to a dual 1-form
> 2. Reconstruct the vector in the tangent plane
>
> In DEC terms: `grad_vector = ⋆₁(d₀φ)` with appropriate reconstruction.

### Discrete Divergence (∇· on vectors)

Using cotan weights for the mesh metric:
```
(∇·V)ᵢ = (1/Aᵢ) Σⱼ wᵢⱼ (Vⱼ - Vᵢ)
```

Where wᵢⱼ = (cot αᵢⱼ + cot βᵢⱼ)/2 are the cotan weights encoding the metric.

### Discrete Laplacian (∇²)

The **cotan Laplacian** is the mesh discretization of ∇²:
```
(∇²φ)ᵢ = (1/Aᵢ) Σⱼ wᵢⱼ (φⱼ - φᵢ)
```

This is the famous formula used throughout geometry processing.

## For DEC Users: Translation Table

If you're coming from Discrete Exterior Calculus literature, here's how the concepts map:

| DEC Notation | Geometric Calculus | Notes |
|--------------|-------------------|-------|
| 0-form | Scalar field | Grade-0 multivector |
| 1-form | Vector field (proximal) | Grade-1 multivector |
| 2-form | Bivector field | Grade-2 multivector |
| d (exterior derivative) | ∇∧ (outer part of ∇) | Grade-raising |
| δ = ⋆d⋆ (codifferential) | ∇· (inner part of ∇) | Grade-lowering |
| d + δ | ∇ | Full geometric derivative |
| Hodge star ⋆ | Duality: ⋆A = AI⁻¹ | I = pseudoscalar |
| Laplacian Δ = dδ + δd | ∇² | Same operator |
| ⋆₁ (edge Hodge) | Cotan weights | Metric on edges |
| ⋆₀ (vertex Hodge) | Vertex areas | Metric on vertices |

**Key insight:** The geometric derivative ∇ unifies d and δ into a single operator. The grade of the result tells you whether you computed the gradient (grade increases), divergence (grade decreases), or both.

---

## FTGC on Triangle Meshes (NEW)

The **Fundamental Theorem of Geometric Calculus** on meshes:
```
∫_M ∇F dV = ∮_∂M F dS
```

This section describes the FTGC-based implementation in `scripts/leapfrog_gc_mesh.py`.

### Key Classes

#### MeshGeometricDerivative
The unified ∇ operator on triangle meshes:

```python
class MeshGeometricDerivative:
    """
    The discrete geometric derivative ∇ on triangle meshes.
    
    Implements: ∇ = Σᵢ eⁱ ⊗ ∂/∂xᵢ
    
    The cotan weights encode the reciprocal basis eⁱ.
    """
    
    def wedge(self, F):    # ∇∧ (outer, grade-raising)
    def inner(self, F):    # ∇· (inner, grade-lowering)
    def apply(self, F):    # Full ∇ = ∇∧ + ∇·
    def grad(self, f):     # Gradient of scalar
    def div(self, V):      # Divergence of vector
    def curl(self, V):     # Curl of vector
    def laplacian(self, f): # ∇² = ∇·∇
```

#### MeshMultivectorField
Multivector fields with staggered storage:

```python
class MeshMultivectorField:
    """
    Grade 0 (scalars)   → Vertices  (n_vertices,)
    Grade 1 (vectors)   → Edges     (n_edges,)
    Grade 2 (bivectors) → Faces     (n_triangles,)
    """
    grade0: np.ndarray  # Vertex values
    grade1: np.ndarray  # Edge values
    grade2: np.ndarray  # Face values
```

### The Cotan Weights = Reciprocal Basis

The fundamental insight: **cotan weights encode the reciprocal basis vectors**.

In FTGC, the metric tensor gⁱʲ defines the reciprocal basis:
```
eⁱ = gⁱʲ eⱼ
```

On meshes, the cotan weights (cot α + cot β)/2 are exactly this metric integrated over dual cells!

### FTGC Identities (Verified Discretely)

The implementation verifies:

| Identity | Meaning | Numerical Check |
|----------|---------|-----------------|
| ∇∧(∇∧f) = 0 | Curl of gradient is zero | ✓ Machine precision |
| ∇²(const) = 0 | Laplacian of constant is zero | ✓ Machine precision |
| ⟨∇f, V⟩ + ⟨f, ∇·V⟩ = boundary | Integration by parts | ✓ Up to boundary terms |

### Why Staggering Works (FTGC Explanation)

The geometric derivative ∇ changes grade:
- ∇ on grade 0 → grade 1 (gradient: vertices → edges)
- ∇ on grade 1 → grade 0 + 2 (div + curl: edges → vertices + faces)
- ∇ on grade 2 → grade 1 (curl adjoint: faces → edges)

**Since inputs and outputs have different grades, they naturally live on different mesh elements!**

This is the deep reason why:
- Yee's staggered grid scheme works
- DEC places forms on different simplices
- The cotan Laplacian is so universal

### Example: Wave Equation via FTGC

```python
from leapfrog_gc_mesh import TriangleMesh, LeapfrogGC

mesh = create_grid_mesh(nx=32, ny=32)
solver = LeapfrogGC(mesh)

# Initial conditions
u0 = np.exp(-50 * ((x - 0.5)**2 + (y - 0.5)**2))
v0 = np.zeros_like(u0)

# Time stepping: u(t+dt) = 2u(t) - u(t-dt) + c²dt²∇²u
dt = solver.estimate_cfl(c=1.0)
u_final = solver.wave_simulate(u0, v0, dt, n_steps=100, c=1.0)
```

### Example: Maxwell's Equations via FTGC

```python
# E (grade 1 on edges), B (grade 2 on faces)
nabla = solver.nabla

# Faraday: ∂B/∂t = -∇∧E
B_new = B - dt * nabla.curl(E)

# Ampère: ∂E/∂t = c²∇·B
E_new = E + dt * c**2 * nabla.inner(B_new)
```

### Implementation Files

| File | Formulation | Use Case |
|------|-------------|----------|
| `discrete_geometric_calculus.py` | FTGC on cubic grids | Regular grids, PDEs |
| `leapfrog_cotan_mesh.py` | DEC on meshes | Original DEC style |
| `leapfrog_gc_mesh.py` | **FTGC on meshes** | Unified ∇, cleaner API |

All produce equivalent results; the FTGC versions provide conceptual unification.

---

## Leapfrog Time Stepping

For spacetime algebras Cl(1,n), solve wave/Dirac equations:

### Dirac Equation
```
∇ψ = 0  →  ∂_t ψ = -∇_s ψ
```

Where ∇_s is the spatial geometric derivative.

### Leapfrog Scheme
```
ψ(t+Δt) = ψ(t-Δt) + 2Δt · ∂_t ψ(t)
        = ψ(t-Δt) - 2Δt · ∇_s ψ(t)
```

**Properties**:
- Second-order accurate O(Δt²)
- Time-reversible (symplectic)
- CFL condition: c·Δt/h < 1

### Wave Equation on Mesh
```
∂²u/∂t² = c²∇²u
```

**Leapfrog update:**
```
u(t+dt) = 2u(t) - u(t-dt) + c²dt²(M⁻¹Lu)
```

Where L is the cotan Laplacian and M is the mass matrix.

## Connection to GAlgebra

| GAlgebra (Symbolic) | Discrete (Numeric) |
|---------------------|-------------------|
| `diff(F, x)` | `(F[i+1]-F[i-1])/(2h)` |
| `e^i` multivector | Array index operation |
| `grad * F` | `nabla.grad(field)` |
| Exact | O(h²) approximate |

### GAlgebra Code (Symbolic)
```python
from galgebra.ga import Ga
coords = symbols('x y z', real=True)
ga = Ga('e', g=[1,1,1], coords=coords)
grad, rgrad = ga.grads()
F = ga.mv('F', 'scalar', f=True)
grad_F = grad * F  # Symbolic ∇F
```

### Discrete Code (Numeric)
```python
from discrete_geometric_calculus import *
basis = CliffordBasis(EUCLIDEAN_3D)
nabla = SplitDifferentialOperator(basis, spacing=0.1)
field = scalar_field(basis, shape, values)
grad_F = nabla.grad(field)  # Numeric ∇F
```

## Domain Comparison

| Domain | Implementation | Metric | Best For |
|--------|----------------|--------|----------|
| Flat ℝⁿ | Cubic grid | Constant | Fast simulations |
| Curved (regular) | Curvilinear | g_ij(x) | Spheres, tori |
| Unstructured | Mesh (cotan) | Cotan weights | CAD meshes |
| Symbolic | GAlgebra | Any | Derivations |

## External Resources

### Libraries
- **galgebra**: `pip install galgebra` — Symbolic GA with `.cderiv()` for covariant derivatives
- **kingdon**: `pip install kingdon` — Numeric GA (Python)
- **numga**: Eelco's JAX-based GA with exp/log for all signatures (`pip install numga`)
- **pycomplex**: Eelco's DEC library for simplicial/cubical complexes (`pip install pycomplex`)

### Papers
- Hestenes & Sobczyk, "Clifford Algebra to Geometric Calculus" (1984) — Chapter 4: covariant derivative
- Crane, "Discrete Differential Geometry" (CMU course notes)

### Talks
- Eelco Hoogendoorn, "Catching Geometric Waves" (GAME2023) — Yee scheme in GA terms

---

## Extensions: When Connection Terms Are Needed

### ✅ Safe (no connection needed)

| Operation | Why |
|-----------|-----|
| Gradient of scalar `∇f` | No direction to transport |
| Divergence of vector `∇·V` | `(1/√g)∂ᵢ(√g V^i)` absorbs connection |
| Scalar Laplacian `∇²f` | Combines above |

### ⚠️ Requires connection (Christoffel symbols)

| Operation | Why |
|-----------|-----|
| Covariant derivative of vector `∇ᵢV^j` | Components rotate under transport |
| Parallel transport | Defines "same direction" across tangent spaces |
| Vector Laplacian | Involves second derivatives of vector components |
| Full `∇A` for grade > 0 | Each blade rotates |

### GAlgebra Covariant Derivative

```python
from sympy import symbols, sin
from galgebra.ga import Ga

# Sphere coordinates with connection
theta, phi = symbols('theta phi', real=True)
ga = Ga('e', g=[1, sin(theta)**2], coords=(theta, phi))
grad, rgrad = ga.grads()

# Covariant derivative of vector field (includes Christoffel terms):
V = ga.mv('V', 'vector', f=True)
div_V = grad | V  # Inner product = divergence
```

### Eelco's Library Ecosystem

| Library | Domain | Key Feature |
|---------|--------|-------------|
| **numga** | GA in JAX/NumPy | Lazy operator fusion, exp/log for dim ≤ 6 |
| **pycomplex** | DEC on meshes | Simplicial/cubical, multigrid, boundary handling |

**Key insight from GAME2023:** The Yee scheme (staggered grid EM) is naturally expressed in GA when you identify:
- E-field on edges ↔ grade-1 multivector
- B-field on faces ↔ grade-2 multivector
- ∂B/∂t = -∇∧E and ∂E/∂t = ∇·B

---

## Quick Reference

### Signatures
```
EUCLIDEAN_2D = Cl(2,0,0)
EUCLIDEAN_3D = Cl(3,0,0)
MINKOWSKI_4D = Cl(1,3,0)  # Spacetime
PGA_3D = Cl(3,0,1)
CGA_3D = Cl(4,1,0)
```

### Grade Extraction
```python
field.grade_select(0)  # Scalar part
field.grade_select(1)  # Vector part
field.grade_select(2)  # Bivector part
```

### Common Operations (Cubic Grid)
```python
nabla.grad(F)           # Geometric derivative ∇F
nabla.laplacian(F)      # ∇²F
nabla.partial(F, i)     # ∂F/∂xᵢ
nabla.leapfrog_dirac()  # Time step for ∇ψ = 0
```

### Common Operations (Triangle Mesh - FTGC)
```python
nabla = MeshGeometricDerivative(mesh)
nabla.grad(f)           # Gradient (scalar → edges)
nabla.div(V)            # Divergence (edges → vertices)
nabla.curl(V)           # Curl (edges → faces)
nabla.laplacian(f)      # ∇²f = div(grad(f))
nabla.wedge(F)          # ∇∧F (outer derivative)
nabla.inner(F)          # ∇·F (inner derivative)
nabla.apply(F)          # Full ∇F = ∇∧F + ∇·F
```
