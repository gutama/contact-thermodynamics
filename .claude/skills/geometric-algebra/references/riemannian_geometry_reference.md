# Riemannian Geometry in Geometric Calculus

Complete coordinate-free formulation of Riemannian geometry using Geometric Algebra.

**Version 1.0** — Eliminates Christoffel symbols entirely

## Overview

This reference presents Riemannian geometry using Geometric Calculus, where:
- **Christoffel symbols** → absorbed into the **connection bivector** ω
- **Riemann tensor** → absorbed into the **curvature 2-form** Ω
- **Metric tensor** → implicit in **frame inner products**
- **Index gymnastics** → eliminated by **reciprocal frames**

## See Also

- [`ddg_reference.md`](ddg_reference.md) — Discrete curvature on meshes
- [`formulas.md`](formulas.md) — Core GA identities
- [`covariant_derivative.py`](../scripts/covariant_derivative.py) — Implementation
- [`geometric_calculus_discrete.md`](geometric_calculus_discrete.md) — Discrete ∇ operator

---

## I. Foundational Objects

### The Frame Bundle

| Object | Grade | Definition | Role |
|--------|-------|------------|------|
| Frame {eᵢ} | 1 | Basis for TₚM | Spans tangent space |
| Reciprocal {eⁱ} | 1 | eⁱ · eⱼ = δⁱⱼ | Dual basis |
| Pseudoscalar I | n | I = e₁ ∧ ··· ∧ eₙ | Volume element |
| Connection ω | 2 | ωᵢ = ½(∂ᵢeⱼ) ∧ eʲ | Replaces Γᵏᵢⱼ |
| Curvature Ω | 2 | Ω = dω + ω ∧ ω | Replaces Rˡᵢⱼₖ |

### Metric Structure

The metric is **implicit** in the frame:

```
gᵢⱼ = eᵢ · eⱼ           (covariant metric)
gⁱʲ = eⁱ · eʲ           (contravariant metric)
g = det(gᵢⱼ) = |I|²     (determinant)
```

**Index raising/lowering** is automatic:
```
eⁱ = gⁱʲ eⱼ
eᵢ = gᵢⱼ eʲ
vᵢ = v · eᵢ  (covariant)
vⁱ = v · eⁱ  (contravariant)
```

---

## II. The Connection Bivector

### Definition

The **connection bivector** encodes how the frame varies:

```
ωᵢ = ½ eʲ ∧ (∂ᵢeⱼ) = -½(∂ᵢeⱼ) ∧ eʲ
```

### Action on Frame

```
∂ᵢeⱼ = ωᵢ × eⱼ          (frame rotation)
∂ᵢeʲ = -ωᵢ × eʲ         (reciprocal rotates oppositely)
```

where × is the **commutator product**: A × B = ½(AB − BA)

### Properties

- **Metric compatibility**: ωᵢ × (eⱼ · eₖ) = 0 (preserves inner products)
- **Torsion-free**: ωᵢ × eⱼ = ωⱼ × eᵢ (symmetric in coordinate frame)
- **Bivector-valued**: Lives in ∧²TₚM

### Relation to Christoffel Symbols

```
ωᵢ = Γᵏᵢⱼ (eₖ ∧ eʲ)

Γᵏᵢⱼ = ½ gᵏˡ(∂ᵢgⱼˡ + ∂ⱼgᵢˡ − ∂ˡgᵢⱼ)
```

**Key insight**: Christoffel symbols are just **components** of ω — never needed explicitly.

---

## III. The Covariant Derivative

### The Vector Derivative

```
∇ = eⁱ ∂ᵢ = eⁱ ⊗ ∂/∂xⁱ
```

A **vector-valued differential operator**.

### Action on Scalars

```
∇f = eⁱ ∂ᵢf = (∂ᵢf) eⁱ = grad f
```

### Action on Vectors

```
∇ᵤv = ∂ᵤv + ω(u) × v
```

where:
- `∂ᵤv = uⁱ(∂ᵢvʲ)eⱼ` is the component derivative
- `ω(u) = uⁱωᵢ` is the connection along u

### Action on Multivectors

For any k-vector A:

```
∇ᵤA = ∂ᵤA + ω(u) × A
```

The commutator product **rotates** multivectors of any grade.

### Derived Operators

```
∇ · v = eⁱ · (∇_{eᵢ}v) = (1/√g) ∂ᵢ(√g vⁱ)     (divergence)
∇ ∧ v = eⁱ ∧ (∇_{eᵢ}v)                         (curl, bivector)
∇²f = ∇ · (∇f) = (1/√g) ∂ᵢ(√g gⁱʲ ∂ⱼf)        (Laplacian)
```

---

## IV. Curvature

### The Curvature 2-Form (Cartan Structure Equation)

```
Ω = dω + ω ∧ ω
```

Expanded:

```
Ωᵢⱼ = ∂ᵢωⱼ − ∂ⱼωᵢ + ωᵢ × ωⱼ
```

### Curvature as Commutator

```
[∇ᵤ, ∇ᵥ] − ∇_{[u,v]} = ℛ(u ∧ v)
```

For coordinate frame ([eᵢ, eⱼ] = 0):

```
[∇ᵢ, ∇ⱼ] = ℛ(eᵢ ∧ eⱼ)
```

### Action on Vectors

```
ℛ(u ∧ v) · w = Ω(u, v) × w
```

Curvature rotates vectors via the commutator product.

### Curvature Operator

The Riemann tensor is a **linear map on bivectors**:

```
ℛ: ∧²TₚM → ∧²TₚM
ℛ(u ∧ v) = Ω(u, v)
```

### Symmetries

**First Bianchi identity**:
```
ℛ(u∧v)·w + ℛ(v∧w)·u + ℛ(w∧u)·v = 0
```

**Pair symmetry** (self-adjoint):
```
⟨ℛ(u∧v), w∧z⟩ = ⟨ℛ(w∧z), u∧v⟩
```

**Second Bianchi identity**:
```
∇ ∧ Ω = 0
```

---

## V. Curvature Contractions

### Sectional Curvature

For 2-plane Π = u ∧ v:

```
K(Π) = Ω(u,v) · (u ∧ v) / |u ∧ v|²
```

Gaussian curvature of geodesic surface in Π-direction.

### Ricci Curvature

**Ricci operator** Ric: TₚM → TₚM:

```
Ric(v) = eⁱ · Ω(eᵢ, v)
```

**Ricci tensor**:
```
Ric(u, v) = u · Ric(v)
```

### Scalar Curvature

```
R = eⁱ · Ric(eᵢ) = tr(Ric)
```

### Einstein Tensor

```
G(v) = Ric(v) − ½Rv
```

**Einstein's equation**:
```
G(v) = 8πG T(v)
```

### Weyl Tensor

Trace-free part of Riemann (conformal curvature):

```
W = ℛ − (1/(n-2))(Ric ∧ g − R/(2(n-1)) g ∧ g)
```

---

## VI. Geodesics and Parallel Transport

### Parallel Transport

Vector w is **parallel** along curve with velocity v if:

```
∇ᵥw = 0
⟺ ẇ + ω(v) × w = 0
```

### Geodesic Equation

Velocity parallel-transports itself:

```
∇ᵥv = 0
⟺ v̇ + ω(v) × v = 0
```

### Holonomy

Transport around closed loop → rotation by curvature:

```
δR ≈ 1 + ½Ω(u, v) · (area)
```

**Curvature = holonomy per unit area**

### Geodesic Deviation (Jacobi Equation)

Separation ξ between nearby geodesics:

```
∇ᵥ∇ᵥξ = −ℛ(v ∧ ξ) · v
```

- K > 0: geodesics converge (sphere)
- K < 0: geodesics diverge (hyperbolic)
- K = 0: geodesics stay parallel (flat)

---

## VII. The Shape Operator and Pseudoscalar

### For Surfaces in ℝ³

**Tangent pseudoscalar** (bivector):
```
B = e₁ ∧ e₂
```

**Normal-bivector duality**:
```
n = −BI⁻¹     ⟺     B = nI
```

where I = e₁e₂e₃ is the 3D pseudoscalar.

### Shape Operator

```
S(v) = −∇ᵥn = (v · ∇B) · I⁻¹
```

### Gaussian Curvature

```
K · B = (∂₁n) ∧ (∂₂n)
K = det(S) = κ₁κ₂
```

### Higher Dimensions

For k-manifold in n-space with tangent k-blade Bₖ:

**Generalized Gauss-Weingarten**:
```
∇ᵥBₖ = ω(v) × Bₖ + h(v) ∧ Bₖ
```

- ω(v) × Bₖ: intrinsic (rotates within tangent space)
- h(v) ∧ Bₖ: extrinsic (tilts into normal space)

---

## VIII. Key Theorems

### Gauss-Bonnet (2D)

```
∫_M K dA + ∫_∂M κg ds = 2πχ(M)
```

GA form:
```
∫_M Ω + ∮_∂M ω = 2πχ(M) · I₂
```

### Generalized Gauss-Bonnet (2m dimensions)

```
∫_M Pf(Ω) = (2π)ᵐ χ(M)
```

where Pf(Ω) = (1/2ᵐm!) Ω ∧ ··· ∧ Ω (m times).

### Bochner Formula

```
∇²α = ∇(∇·α) + ∇·(∇∧α) + Ric(α)
```

### Cartan Structure Equations

First (torsion):
```
∇ ∧ eⁱ = ω × eⁱ
```

Second (curvature):
```
Ω = dω + ω ∧ ω
```

---

## IX. Translation Tables

### Classical → Geometric Calculus

| Classical | GA | Notes |
|-----------|-----|-------|
| gᵢⱼ | eᵢ · eⱼ | Metric implicit |
| gⁱʲ | eⁱ · eʲ | Reciprocal frame |
| Γᵏᵢⱼ | ωᵢ | Connection bivector |
| Rˡᵢⱼₖ | Ω | Curvature 2-form |
| ∇ᵢVʲ | ∇ᵤv | Covariant derivative |
| εᵢⱼₖ... | I | Pseudoscalar |

### Differential Forms → GA

| Forms | GA | Relation |
|-------|-----|----------|
| ωⁱⱼ | ω | ω = ½ωⁱⱼ eᵢ ∧ eʲ |
| Ωⁱⱼ | Ω | Ω = ½Ωⁱⱼ eᵢ ∧ eʲ |
| d | ∇∧ | Exterior derivative |
| δ = ⋆d⋆ | ∇· | Codifferential |
| Δ | ∇² | Laplacian |



---

## X. The Master Equations

All of Riemannian geometry in coordinate-free form:

```
METRIC:           gᵢⱼ = eᵢ · eⱼ
CONNECTION:       ∂ᵢeⱼ = ωᵢ × eⱼ
COVARIANT DERIV:  ∇ᵤv = ∂ᵤv + ω(u) × v
CURVATURE:        Ω = dω + ω ∧ ω
RIEMANN ACTION:   [∇ᵤ, ∇ᵥ]w = Ω(u,v) × w
GEODESIC:         v̇ + ω(v) × v = 0
PARALLEL TRANS:   ẇ + ω(v) × w = 0
BIANCHI:          ∇ ∧ Ω = 0
```

---

## XI. The Fundamental Principle

**Curvature = Derivative of Tangent Pseudoscalar**

```
∇ ∧ ∇ Bₖ = Ω · Bₖ
```

The curvature tensor is exactly the **obstruction to ∇∧∇ = 0** acting on the tangent pseudoscalar.

- 2D surface: K = rotation angle of B₂ under transport around loop
- k-manifold: ℛ = full rotation operator on Bₖ under transport

**Topology constrains its integral** (Gauss-Bonnet).

---

## XII. Implementation Notes


### Discrete Connection Bivector

On triangle meshes, the connection bivector becomes a **rotation between adjacent frames**:

| Continuous | Discrete | Mesh Element |
|------------|----------|--------------|
| ωᵢ = ½ eʲ ∧ (∂ᵢeⱼ) | Dihedral rotation | Edge (shared between triangles) |
| ∂ᵢeⱼ = ωᵢ × eⱼ | Frame alignment via unfolding | Triangle pair |

**Implementation**: When parallel transporting across edge (i,j):
1. Compute dihedral angle θ between adjacent face normals
2. Rotation axis = shared edge direction
3. Apply rotation: v' = Rθ(v) around edge axis

This is exactly the discrete analog of ω × v — the connection rotates vectors 
as they cross from one tangent plane to the next.

**Code reference**: `ga_geodesics.js` → `parallelTransport()` method

**JavaScript implementation** (NEW):
- `src/riemannian-discrete.js` → Complete discrete Riemannian module:
  - `MeshConnectionBivector` — dihedral rotation at each edge
  - `MeshCurvature2Form` — angle defect curvature K = 2π - Σθ
  - `MeshParallelTransport` — transport vectors across edges
  - `BianchiIdentityVerifier` — verify Σ Kᵥ = 2πχ(M)

### Discrete Manifolds

| Continuous | Discrete | Element |
|------------|----------|---------|
| Reciprocal eⁱ | Cotan weights | Edge |
| K (Gaussian) | 2π − Σθⱼ | Vertex |
| H (mean) | Cotan Laplacian | Vertex |
| Holonomy | Angle defect | Loop |
| Connection ω | Dihedral rotation | Edge |

### Cotan Weights = Mesh Metric

```
wᵢⱼ = ½(cot αᵢⱼ + cot βᵢⱼ)
```

This is the discrete reciprocal basis — encodes all metric information.

### Bianchi Identity on Meshes

The discrete Bianchi identity ∇ ∧ Ω = 0 becomes the **Gauss-Bonnet theorem**:

```
Σᵥ Kᵥ = 2π χ(M) = 2π (V - E + F)
```

Verified with error < 10⁻⁸ on icosahedron (χ = 2) and flat grids (χ = 1).

### Computational Advantages

- Bivector algebra parallelizes (each component independent)
- No index management
- Same algorithms at all dimensions
- Natural sparse structure

---

## References

1. Hestenes, D. & Sobczyk, G. "Clifford Algebra to Geometric Calculus" (1984)
2. Doran, C. & Lasenby, A. "Geometric Algebra for Physicists" (2003)
3. Crane, K. "Discrete Differential Geometry" (2025)
4. Cartan, É. "The Theory of Spinors" (1966)
