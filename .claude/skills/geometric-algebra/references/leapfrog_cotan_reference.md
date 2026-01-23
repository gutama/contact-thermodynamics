# Leapfrog/Yee Scheme with Cotan Formula via FTGC

## The Question
Can we use the **leapfrog/Yee time-stepping scheme** with the **cotan formula** (the "Swiss Army" discretization) on triangle meshes?

## The Answer: Yes! They're Made for Each Other

The combination is not just possible — it's **natural and elegant**. The Fundamental Theorem of Geometric Calculus (FTGC) explains WHY:

```
Yee's Staggered Grid  ←→  Grade Staggering in ∇
        ↓                           ↓
E at edge centers    ←→  Grade-1 (vectors on edges)
B at face centers    ←→  Grade-2 (bivectors on faces)
        ↓                           ↓
   Leapfrog in time  ←→  ∇ mixes grades naturally
```

## Two Equivalent Formulations

### DEC (Discrete Exterior Calculus)
```
d (exterior derivative)  →  grade-raising
δ = ⋆d⋆ (codifferential) →  grade-lowering  
Δ = dδ + δd              →  Laplacian
⋆ (Hodge star)           →  metric/duality
```

### FTGC (Fundamental Theorem of Geometric Calculus)
```
∇∧ (outer derivative)    →  grade-raising (same as d)
∇· (inner derivative)    →  grade-lowering (same as δ)
∇² = ∇·∇                 →  Laplacian (unified, not sum)
⋆A = AI⁻¹               →  duality via pseudoscalar
```

**Key insight:** Both are equivalent! FTGC unifies them under a single operator ∇.

## The Split Differential Operator

The geometric derivative on meshes:
```
∇ = Σᵢ eⁱ ⊗ ∂/∂xᵢ
```

Where:
- **eⁱ** = reciprocal basis vectors (encoded by COTAN WEIGHTS!)
- **∂/∂xᵢ** = partial derivatives (finite differences along edges)

The cotan weights ARE the mesh metric expressed in geometric algebra form.

## Key Insight: The Cotan Formula IS the Reciprocal Basis

The "Swiss Army knife" cotan formula:
```
w_ij = (cot α_ij + cot β_ij) / 2
```

In DEC language: This is the discrete Hodge star ⋆₁.
In FTGC language: This encodes the reciprocal basis eⁱ·eʲ = gⁱʲ.

The full cotan Laplacian in both languages:

| DEC | FTGC |
|-----|------|
| Δ = ⋆₀⁻¹ d₀ᵀ ⋆₁ d₀ | ∇² = ∇·∇ = M₀⁻¹ (∇∧)ᵀ M₁ (∇∧) |

Same matrix, different conceptual framing!

## The Staggering Correspondence

| Grade | Mesh Element | DEC | FTGC | Physical Field |
|-------|--------------|-----|------|----------------|
| 0 | Vertices | 0-form | Scalar | Potential |
| 1 | Edges | 1-form | Vector | Electric field |
| 2 | Faces | 2-form | Bivector | Magnetic flux |

**Why staggering works:** The geometric derivative ∇ changes grade:
- ∇ on grade 0 → grade 1 (gradient)
- ∇ on grade 1 → grade 0 + grade 2 (divergence + curl)
- ∇ on grade 2 → grade 1 (curl adjoint)

So inputs and outputs naturally live on DIFFERENT mesh elements!

## Wave Equation: ∂²u/∂t² = c²∇²u

### Spatial Discretization (Cotan Laplacian via FTGC)
```python
# The geometric derivative ∇ with cotan metric
nabla = MeshGeometricDerivative(mesh)

# Laplacian: ∇² = ∇·(∇f) = div(grad(f))
L = nabla.laplacian_matrix()  # Uses cotan weights internally
```

### Time Discretization (Leapfrog)
```python
u_next = 2*u_curr - u_prev + c²*dt²*(L @ u_curr)
```

### Stability (CFL Condition)
```
dt < h_min / c
```

## Maxwell's Equations via FTGC

### Unified Form (Spacetime Algebra)
```
∇F = J
```
where F = E + IB is the electromagnetic bivector.

### Split into Space and Time
```
∂B/∂t = -∇∧E    (Faraday: outer derivative raises grade)
∂E/∂t = c²∇·B   (Ampère: inner derivative lowers grade)
```

### Yee-Style Update with FTGC
```python
# E (grade 1 on edges), B (grade 2 on faces)
# Staggered: E at t, B at t±dt/2

# Faraday: ∂B/∂t = -∇∧E
B_new = B - dt * nabla.curl(E)

# Ampère: ∂E/∂t = c² ∇·B  
E_new = E + dt * c² * nabla.inner(B_new)
```

## FTGC Identities (Verified Discretely)

The discrete operators preserve key identities:

| Identity | Meaning | Numerical Check |
|----------|---------|-----------------|
| ∇∧(∇∧f) = 0 | Curl of gradient is zero | Max entry ≈ 0 |
| ∇²(const) = 0 | Laplacian of constant is zero | Max entry ≈ 0 |
| ∫_M ∇F = ∮_∂M F | FTGC (Stokes generalized) | Boundary terms match |

## Why This Works So Well

### 1. Geometric Compatibility (FTGC Perspective)
- Cotan weights encode the reciprocal basis eⁱ
- The split operator ∇ = Σᵢ eⁱ ⊗ ∂ᵢ preserves grade structure
- Different grades naturally live on different mesh elements

### 2. Conservation Properties
- Leapfrog is **symplectic** → energy conservation
- FTGC preserves **topological invariants** (∇∧∇∧ = 0)
- Together: stable long-time behavior

### 3. Conceptual Unification
- One operator ∇ instead of separate d, δ, ⋆
- Grade changes explain staggering
- Cotan weights = metric = reciprocal basis

## Code Structure: DEC vs FTGC

### DEC Style (Original)
```python
class LeapfrogMesh:
    def __init__(self, mesh):
        self.dec = DECOperators(mesh)  # d₀, d₁, ⋆₀, ⋆₁, ⋆₂
    
    def maxwell_step(self, E, B, dt, c):
        B_new = B - dt * (self.dec.d1 @ E)
        E_new = E + dt*c² * (self.dec.star1_inv @ self.dec.d1.T @ self.dec.star2 @ B_new)
        return E_new, B_new
```

### FTGC Style (Refactored)
```python
class LeapfrogGC:
    def __init__(self, mesh):
        self.nabla = MeshGeometricDerivative(mesh)  # Unified ∇
    
    def maxwell_step(self, E, B, dt, c):
        B_new = B - dt * self.nabla.curl(E)      # ∇∧E
        E_new = E + dt*c² * self.nabla.inner(B_new)  # ∇·B
        return E_new, B_new
```

Same computation, cleaner conceptual model!

## Translation Table

| DEC | FTGC | Mesh Implementation |
|-----|------|---------------------|
| 0-form | Scalar (grade 0) | Vertex values |
| 1-form | Vector (grade 1) | Edge values |
| 2-form | Bivector (grade 2) | Face values |
| d (ext. derivative) | ∇∧ (outer) | Incidence matrix |
| δ = ⋆d⋆ (codifferential) | ∇· (inner) | M⁻¹ dᵀ M |
| d + δ | ∇ (full) | ∇∧ + ∇· |
| Hodge star ⋆ | Duality: ⋆A = AI⁻¹ | Metric matrices |
| **⋆₁ (edge Hodge)** | **Reciprocal basis** | **COTAN WEIGHTS** |
| Δ = dδ + δd | ∇² = ∇·∇ | Cotan Laplacian |

## Implementation Files

| File | Description |
|------|-------------|
| `scripts/leapfrog_cotan_mesh.py` | DEC-based implementation |
| `scripts/leapfrog_gc_mesh.py` | **FTGC-based implementation** (unified ∇) |

Both produce identical numerical results; the FTGC version has cleaner conceptual structure.

## The Deep Connection (FTGC Insight)

The cotan formula emerges from the geometric calculus requirement that:
```
eⁱ · eⱼ = δⁱⱼ (reciprocal basis orthonormality)
```

On a mesh, this inner product requires integrating over dual cells, which naturally produces the cotan weights! The formula is universal because it's the UNIQUE discretization that:
1. Preserves the metric structure
2. Satisfies FTGC identities
3. Is coordinate-free

## Summary

| Component | Formula | FTGC Role |
|-----------|---------|-----------|
| Cotan weights | (cot α + cot β)/2 | Reciprocal basis eⁱ |
| Mass matrix | Σ triangle_area/3 | Vertex metric M₀ |
| Gradient | ∇∧f = f[j] - f[i] | Outer derivative |
| Divergence | ∇·V = M₀⁻¹(∇∧)ᵀM₁V | Inner derivative |
| Laplacian | ∇² = ∇·∇ | div(grad) |
| Leapfrog | u(t+dt) = 2u(t) - u(t-dt) + dt²∇²u | Symplectic time |

**The cotan formula encodes the reciprocal basis; leapfrog respects the symplectic structure. Together they discretize FTGC faithfully.**

## References

1. Hestenes & Sobczyk, "Clifford Algebra to Geometric Calculus" (1984) — FTGC theory
2. Desbrun et al., "Discrete Exterior Calculus" (2005) — DEC foundation
3. Crane, "Discrete Differential Geometry" (CMU course) — Cotan formula
4. Yee, "Numerical solution of initial boundary value problems..." (1966) — Leapfrog
5. Hoogendoorn, "Catching Geometric Waves" (GAME2023) — Yee scheme in GA terms
