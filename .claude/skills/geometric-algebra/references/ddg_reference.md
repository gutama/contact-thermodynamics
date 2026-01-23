# Discrete Differential Geometry Reference

Based on Keenan Crane's "Discrete Differential Geometry: An Applied Introduction"  
Adapted for Geometric Algebra by Claude  
**Version 2.0** — Includes integral and variational approaches

## Overview

Discrete Differential Geometry (DDG) provides computational methods for geometry that preserve essential structures from smooth differential geometry. The key insight is that **the wedge product naturally computes determinants**, which are fundamental to all discrete curvature calculations.

## See also

- Mesh primitives & simplex magnitudes: [`mesh_reference.md`](mesh_reference.md)
- Global constraints (Gauss–Bonnet / topology sanity checks): [`unified_topology_reference.md`](unified_topology_reference.md)
- Core identities & operators: [`formulas.md`](formulas.md)

## The Unified Picture

The lectures present two complementary approaches to discrete curvature:

| Approach | View | Quantities | Invariant? |
|----------|------|------------|------------|
| **Integral** | Local vertex | K (angle defect), H (cotan-Laplace) | K is topological |
| **Variational** | Global derivative | Steiner coefficients, flows | Total K is topological |

**Connection:** The gradient hierarchy links scalar and vector quantities:
- ∇Volume → Area (normal vector)
- ∇Area → Mean Curvature Normal
- ∇Mean Curvature → Gaussian Curvature Normal
- ∇Total Gaussian → **0** (topological invariant!)

---

## The Wedge-Determinant Connection

The wedge product of n vectors in n-dimensional space equals the determinant:

```
a ∧ b = det([a, b]) × e₁₂                    # 2D
a ∧ b ∧ c = det([a, b, c]) × e₁₂₃            # 3D
v₀ ∨ v₁ ∨ ... ∨ vₖ = (k! × simplex measure) × I  # PGA
```

---

## Part I: Integral Approach (Local Quantities)

### Gaussian Curvature (Angle Defect)

```
Interior vertex:  K_i = 2π - Σ θⱼ
Boundary vertex:  K_i = π - Σ θⱼ
```

**Properties:** +K = sphere-like, -K = saddle-like, 0 = flat

### Mean Curvature (Cotan-Laplace)

```
H⃗_i = (1/2A_i) Σ (cot α_ij + cot β_ij)(v_j - v_i)
```

### Principal Curvatures

```
κ₁, κ₂ = H ± √(H² - K)
```

### Gauss-Bonnet Theorem

```
Σ K_i = 2π × χ
```

Examples: Sphere χ=2→4π, Torus χ=0→0, Double torus χ=-2→-4π

---

## Part II: Variational Approach (Global Derivatives)

### Steiner's Formula (Mollification)

```
V(P ⊕ B(t)) = V + A·t + H·t² + (4π/3)·t³
```

| Derivative | Result |
|------------|--------|
| V(0) | Volume |
| dV/dt\|₀ | Surface Area |
| d²V/dt²\|₀ | 2 × Total Mean Curvature |
| d³V/dt³\|₀ | 4π (constant!) |

**Key Insight:** d³V/dt³ is constant because total Gaussian curvature is topological!

### Total Mean Curvature (Edge-Based)

```
H_total = (1/2) Σ_edges l_e × θ_e
```

where θ_e = π - (dihedral angle) is the exterior dihedral angle.

**Note:** Unlike total K, total H **depends on triangulation**!

### Schläfli Formula

```
Σ l_e × dθ_e = 0
```

Constraint on polyhedral deformations.

---

## Part III: Curvature Flow

### Mean Curvature Flow
```
df/dt = -H⃗
```
Minimizes area. Smooths surfaces.

### Gaussian Curvature Flow
```
df/dt = -K·n
```

### Willmore Flow
```
Energy: W = Σ H²·A
```
Minimizes bending. Sphere minimum: W = 4π.

---

## Mixed Area (Voronoi)

```javascript
if (vertex_angle > π/2) area += tri_area / 2;
else if (other_angle > π/2) area += tri_area / 4;
else area += (cot(α) × len²₁ + cot(β) × len²₂) / 8;
```

---

## PGA Simplex Formulas

Moved to the mesh reference to avoid duplication:
- See **k-magnitudes** and **simplex carriers** in [`mesh_reference.md`](mesh_reference.md).

---
## CGA Curvature

```
Circle: C = P₀ ∧ P₁ ∧ P₂, κ = 1/radius(C)
Sphere: S = P₀ ∧ P₁ ∧ P₂ ∧ P₃
```

---

## Test Cases

| Shape | χ | Total K | Notes |
|-------|---|---------|-------|
| Sphere | 2 | 4π | K=1 everywhere |
| Torus | 0 | 0 | K>0 outer, K<0 inner |
| Cube | 2 | 4π | K=π/2 at corners |
| Tetrahedron | 2 | 4π | K≈π at vertices |

---

## Equation Summary

```
# Integral
K_i = 2π - Σ θⱼ
H⃗_i = (1/2A) Σ (cotα + cotβ)(vⱼ - vᵢ)
κ₁, κ₂ = H ± √(H² - K)
Σ K_i = 2πχ

# Variational  
V(t) = V + At + Ht² + (4π/3)t³
H_total = (1/2) Σ lₑ θₑ
W = Σ H²A
Σ lₑ dθₑ = 0

# Flow
df/dt = -H⃗  (mean curvature)
df/dt = -K·n (Gauss curvature)
```

---

## References

1. Crane, K. "Discrete Differential Geometry: An Applied Introduction" (2025)
2. Crane, K. "Discrete Curvature" lecture notes Parts I & II
3. Meyer et al. "Discrete Differential-Geometry Operators"
4. De Keninck et al. "Clean up your Mesh!" (2025)
