# CGA Versors & Hybrid CGA Reference

Based on Todd's "Comprehensive Visualization of Versors of CGA" (2024) and the Flemish school principles.

## Core Insight: Hybrid CGA

**Neither pure IPNS (sphere-based) nor pure OPNS (point-based) CGA is complete.**

The sign of B² determines which representation to use:

| B² | Representation | Object | Transformation |
|----|----------------|--------|----------------|
| **< 0** | IPNS (sphere-based) | Circle/Sphere | Rotation (elliptic) |
| **= 0** | Both superimposed | ∞ or 0-radius | Translation (parabolic) |
| **> 0** | OPNS (point-based) | **Point pair** | Boost (hyperbolic) |

This resolves the longstanding IPNS vs OPNS debate: use both, switching based on the discriminant.

---

## Flemish School Principles (Transformation-First)

1. **k-reflection A** → visualize as k hyperplanes it decomposes into
2. **Geometric product** = transform composition
3. **Geometric objects are invariants of transformations** (not vice versa)
4. **log(transformation)** generalizes "axis of rotation"

**Key implication:** Understanding transformations is primary; objects emerge as special cases.

---

## The Fundamental Isomorphism

> **n-dimensional CGA Cl(n+1,1,0) ≅ (n+1)-dimensional Hyperbolic PGA**

| CGA | Signature | = | Hyperbolic PGA |
|-----|-----------|---|----------------|
| 1D CGA | Cl(2,1,0) | = | 2D Hyperbolic PGA |
| 2D CGA | Cl(3,1,0) | = | 3D Hyperbolic PGA |
| 3D CGA | Cl(4,1,0) | = | 4D Hyperbolic PGA |

The **horosphere** (null cone) in CGA = **Klein ball boundary** in hyperbolic geometry.

### Intermediate Space Construction

1. Take the Klein ball in (n+1)D hyperbolic space
2. Place "eye" at north pole (where n∞ is tangent to null manifold)
3. Stereographic projection → CGA modelling space

This explains:
- **PGA as subalgebra:** Objects passing through the eye point
- **"Imaginary" objects:** Real objects in intermediate space with no 3D projection
- **e₋ transformation:** Hyperplane in intermediate space, "point-pair field" in modelling space

---

## Complete Versor Classification (3D CGA)

### Grade 1: 1-vectors (Reflections)

| B² | Example | Transformation | Invariant | Notes |
|----|---------|----------------|-----------|-------|
| +1 | e₁, e₊, e₃+e₀ | Sphere inversion | Sphere/plane | Extrinsic, has sidedness |
| 0 | e₀ = e₊+e₋ | Annihilation | ∞-radius sphere | Null |
| 0 | e₊−e₋ | Annihilation | 0-radius sphere | Null |
| −1 | e₋, e₋+0.5e₁ | Hyperbolic reflection | Hyperplane in intermediate space | "Point-pair field" in 3D |

### Grade 2: Bivectors (Rotors)

| B² | Example | Transformation | Invariant | Pencil Type |
|----|---------|----------------|-----------|-------------|
| −1 | e₁₂, e₂₃, e₁₂+e₁₀ | Rotation | Circle/Line | **Elliptic pencil** |
| 0 | e₀₁, e₀₂ | Translation | ∞-radius circle + ∞-radius PP | **Parabolic pencil** |
| 0 | e₁₊−e₁₋ | Parabolic motion | 0-radius circle + 0-radius PP | **Parabolic pencil** |
| +1 | e₊₋, e₁₋, e₂₋ | Boost / Scaling | **Point pair** | **Hyperbolic pencil** (dipole) |

### Grade 3: Trivectors

| B² | Example | Transformation | Invariant | Bundle Type |
|----|---------|----------------|-----------|-------------|
| −1 | e₁₂₃, e₁₂₊, e₂₃₊ | Rotoreflection | Point pair / point | Elliptic bundle |
| 0 | e₀₁₂, e₀₂₃ | Transflection | ∞-radius objects | Plane bundle |
| 0 | e₁₂₊−e₁₂₋ | Parabolic transflection | 0-radius objects | Parabolic bundle |
| +1 | e₁₂₋, e₂₃₋+e₁₃₋ | Hyperbolic transflection | **Circle/Line** | Hyperbolic bundle (spear) |

### Grade 4: Quadvectors

| B² | Example | Transformation | Invariant | Notes |
|----|---------|----------------|-----------|-------|
| −1 | e₁₂₃₊, e₁₂₊₋ | Hyperbolic screw | Sphere/plane | Poincaré ball boundary |
| 0 | e₁₂₃₀ | Euclidean screw | ∞-radius sphere | Motor (PGA) |
| 0 | e₁₂₃₊+e₁₂₃₋ | Parabolic screw | 0-radius sphere | |
| +1 | e₁₂₃₊, e₁₂₃₊+0.1e₁₂₊₋ | Isoclinic rotation | Point in intermediate space | "Circle-pair field" |

### Grade 5: Pseudoscalar

| B² | Example | Transformation | Notes |
|----|---------|----------------|-------|
| −1 | e₁₂₃₊₋ | Rotoreflection + hyperbolic translation | "Screwflection" |

---

## Exponentials of Geometric Primitives

### Point Pair Exponential → Booster

```
PP = P₁ ∧ P₂  (bivector, PP² > 0)
exp(t·PP/2) = cosh(t/2) + sinh(t/2)·P̂P
```

- Generates **hyperbolic flow** between the two points
- Points move along hyperbolic trajectories
- Special case: one point at ∞ → uniform scaling

### Circle Pair Exponential → Loxodrome

Form bivector: **B = C₁C₂⁻¹**

| Configuration | Discriminant | Type | Transformation |
|---------------|--------------|------|----------------|
| Intersecting | Δ < 0 | Elliptic | Rotation around intersection points |
| Tangent | Δ = 0 | Parabolic | Shear along tangent |
| Nested/Disjoint | Δ > 0 | Hyperbolic | Logarithmic spiral |

**Discriminant:** Δ = (C₁·C₂)² − C₁²C₂²

### Sphere Pair Exponential → Conformal Transform

Form bivector: **B = S₁S₂⁻¹**

| Configuration | Transformation | Generator |
|---------------|----------------|-----------|
| **Concentric** | Dilation | e₊₋ (hyperbolic) |
| **Intersecting** | Rotation | e_ij (elliptic) |
| **Tangent** | Transversion | e_i₀ (parabolic) |
| **Disjoint** | Boost-like | Mixed |

---

## Pencils and Bundles

### Pencils (Grade 2 — 1-parameter families of spheres)

| Type | Generator B² | Description |
|------|--------------|-------------|
| **Elliptic** | B² < 0 | Spheres through a common circle |
| **Parabolic** | B² = 0 | Spheres tangent at a point |
| **Hyperbolic** | B² > 0 | Spheres "centered" on two points |

### Bundles (Grade 3 — 2-parameter families)

| Type | Generator B² | Description |
|------|--------------|-------------|
| **Elliptic** | B² < 0 | Spheres through two points |
| **Parabolic** | B² = 0 | Planes through points at infinity |
| **Hyperbolic** | B² > 0 | Spheres orthogonal to a circle |

---

## Applications to Curved Meshes

### Curvature from Circle Fitting

```javascript
// Fit circle through 3 consecutive points
C = P₁ ∧ P₂ ∧ P₃
κ = 1 / radius(C)
```

The **discriminant of fitted circles** determines local geometry:

| Circle Config | Surface Type | Gaussian K |
|---------------|--------------|------------|
| Intersecting | Convex/Concave | K > 0 |
| Tangent | Developable | K = 0 |
| Non-intersecting | Saddle | K < 0 |

### Conformal Mesh Deformation

Sphere pair exponentials preserve **angles** (conformal), making them ideal for:
- Texture mapping without angular distortion
- Mesh parameterization
- Shape interpolation (SLERP on sphere pairs)

### Adaptive Operations

Use discriminant to guide:
- **Remeshing:** More triangles where |Δ| is large (high curvature)
- **Smoothing:** Different flows for elliptic vs hyperbolic regions
- **Feature detection:** Sign changes in Δ indicate ridges/valleys

---

## Connection to Physics

### Spacetime (Conformal Spacetime Algebra)

The intermediate-space construction extends to **Cl(2,2)** for 1+1 Minkowski space:
- Null manifold becomes **hyperboloid of one sheet** (not sphere)
- "Spheres" become hyperbolas
- Lightcones are null "circles"

Connects to: Twistor theory, de Sitter space, Penrose diagrams

### STAP Cl(3,1,1)

Space-Time-Algebra Projective: subalgebra of CSTA for electromagnetism and special relativity.

---

## Implementation Notes

### Classifying Bivector Type

```javascript
function classifyBivector(B) {
    const Bsq = B.mul(B).scalar();  // B²
    const eps = 1e-10;
    
    if (Bsq < -eps) return { type: 'elliptic', repr: 'IPNS', object: 'circle' };
    if (Bsq > eps)  return { type: 'hyperbolic', repr: 'OPNS', object: 'point_pair' };
    return { type: 'parabolic', repr: 'both', object: 'null' };
}
```

### Exponentiating by Type

```javascript
// NOTE: algebra.exp() is elliptic-only. Use hybrid_exp below:
function hybrid_exp(B, theta) {
    const Bsq = B.mul(B).scalar();
    const half = theta / 2;
    const Bhat = B.normalized();  // Use .normalized() not .normalize()
    
    if (Bsq < -1e-10) {
        // Elliptic: exp(θB/2) = cos(θ/2) + sin(θ/2)B̂
        return cga.scalar(Math.cos(half)).add(Bhat.scale(Math.sin(half)));
    } else if (Bsq > 1e-10) {
        // Hyperbolic: exp(θB/2) = cosh(θ/2) + sinh(θ/2)B̂
        return cga.scalar(Math.cosh(half)).add(Bhat.scale(Math.sinh(half)));
    } else {
        // Parabolic: exp(θB/2) = 1 + (θ/2)B
        return cga.scalar(1).add(B.scale(half));
    }
}
```

### Discriminant Calculation

```javascript
function discriminant(X1, X2) {
    const dotProd = X1.dot(X2);  // Use .dot() not .inner()
    const n1sq = X1.mul(X1).scalar();
    const n2sq = X2.mul(X2).scalar();
    return dotProd * dotProd - n1sq * n2sq;
}

function classifyPair(X1, X2) {
    const D = discriminant(X1, X2);
    if (D < -1e-10) return 'elliptic';    // Intersecting
    if (D > 1e-10)  return 'hyperbolic';  // Nested/Disjoint
    return 'parabolic';                    // Tangent
}
```

---

## References

1. Todd, H. (2024). "A Comprehensive Visualization of Versors of CGA"
2. Gunn, C. (2011). "Geometry, Kinematics, and Rigid Body Mechanics in Cayley-Klein Geometries"
3. Dorst, L. (2023). "Projective duality encodes complementary orientations"
4. Bobenko, Pinkall, Springborn (2015). "Discrete conformal maps and ideal hyperbolic polyhedra"
5. De Keninck, Roelfs et al. (2024). "From Invariant Decomposition to Spinors"
