# PGA Mesh Processing Reference

Based on "Clean up your Mesh! Part 1: Plane and Simplex"  
De Keninck, Roelfs, Dorst, Eelbode (2025)

## See also

- Map of all references: [`index.md`](index.md)
- Discrete curvature (DDG): [`ddg_reference.md`](ddg_reference.md)
- Discrete geodesics: [`geodesics_reference.md`](geodesics_reference.md)
- GA identities: [`formulas.md`](formulas.md)

## Overview

PGA (Plane-based Geometric Algebra) provides a unified, coordinate-free framework for mesh processing. This reference covers the key formulas and implementations.

## Notation

| Symbol | Meaning |
|--------|---------|
| `v` | Normalized point: `v = e₁₂₃ - xe₀₂₃ + ye₀₁₃ - ze₀₁₂` |
| `o` | Origin point: `o = e₁₂₃` |
| `∨` | Join (regressive) product |
| `∧` | Meet (outer/wedge) product |
| `‖·‖` | Euclidean norm |
| `‖·‖∞` | Ideal norm |
| `Sₖ` | k-carrier (join of k+1 points) |
| `σₖ` | k-simplex |
| `Kₖ` | k-complex |

## The Euclidean Split

Every PGA element decomposes into Euclidean and Ideal parts:

```
A = Aₑ + e₀Aᵢ
```

This split is fundamental for computing magnitudes:
- **Euclidean norm**: `‖A‖ = ‖Aₑ‖`
- **Ideal norm**: `‖A‖∞ = ‖Aᵢ‖ = ‖A ∨ o‖`

## k-Simplices and k-Carriers

### Definition

A k-simplex `σₖ = [v₀, ..., vₖ]` is the convex hull of k+1 vertices.  
Its k-carrier is the join:

```
S(σₖ) = v₀ ∨ v₁ ∨ ... ∨ vₖ
```

### Point (0-simplex)

```javascript
// PGA representation
v = e₁₂₃ - x·e₀₂₃ + y·e₀₁₃ - z·e₀₁₂

// Code
point(x, y, z) {
    mv[E123] = 1;
    mv[E023] = -x;
    mv[E013] = y;
    mv[E012] = -z;
}
```

### Edge (1-simplex)

```
v₀ ∨ v₁ = (→v₁ - →v₀)I₃ + e₀(→v₀ ∧ →v₁)I₃
```

Components:
- **Direction**: `→v₁ - →v₀` (Euclidean part)
- **Moment**: `→v₀ × →v₁` (Ideal part, Plücker coordinates)

### Triangle (2-simplex)

```
v₀ ∨ v₁ ∨ v₂ = ((→v₂ - →v₀) ∧ (→v₁ - →v₀))I₃ + e₀(→v₀ ∧ →v₁ ∧ →v₂)I₃
```

Components:
- **Normal**: `(→v₁ - →v₀) × (→v₂ - →v₀)` (Euclidean part)
- **Volume factor**: `→v₀ · (→v₁ × →v₂)` (Ideal part)

### Tetrahedron (3-simplex)

```
v₀ ∨ v₁ ∨ v₂ ∨ v₃ = det([→v₁-→v₀, →v₂-→v₀, →v₃-→v₀]) · I₃
```

This is a scalar: 6× the signed volume.

## k-Magnitudes

### Table of Magnitudes

| k-simplex | Formula | Result |
|-----------|---------|--------|
| Point | `‖v‖` | 1 (normalized) |
| Edge | `‖v₀ ∨ v₁‖` | Length |
| Triangle | `½‖v₀ ∨ v₁ ∨ v₂‖` | Area |
| Tetrahedron | `⅙‖v₀ ∨ v₁ ∨ v₂ ∨ v₃‖` | Volume |

### General Formula

```
mag(σₖ) = (1/k!) ‖Sₖ‖ = (1/k!) ‖v₀ ∨ ... ∨ vₖ‖
```

## The "Mind the Gap" Technique

Key insight: The magnitude of a k-simplex can be computed from its (k-1)-boundary:

```
‖Sₖ‖ = ‖∑ Sₖ₋₁‖∞
       ∂σₖ
```

### Example: Triangle Area from Edges

```
‖v₀ ∨ v₁ ∨ v₂‖ = ‖v₀ ∨ v₁ + v₁ ∨ v₂ + v₂ ∨ v₀‖∞
```

### Mesh Volume from Faces

For a closed mesh with faces `Fᵢ`:

```
Volume = (1/6) ‖∑ Fᵢ‖∞ = (1/6) ∑(Fᵢ ∨ o)
```

**Critical**: Sum first, then take ideal norm. The signed contributions cancel for overlapping regions.

```javascript
// Implementation
meshVolume(faces) {
    let sumD = 0;
    for (const F of facePlanes) {
        sumD += F[E0];  // d coefficient = F ∨ o
    }
    return sumD / 6;
}
```

## Mesh Gap Detection

For a mesh with non-closed boundary:

```
mag(∂̸Kₖ) = (1/(k-1)!) ‖∑ Sₖ₋₁‖
                        ∂Kₖ
```

The Euclidean norm (not ideal) gives the magnitude of missing boundary elements.

## Center of Mass

For a mesh with uniform density:

```
Ccom = (1/24V) ∑ (v₀ + v₁ + v₂ + o)(F ∨ o)
              ∂K
```

Where `V` is total volume and sum is over boundary faces.

```javascript
// Implementation (returns homogeneous point scaled by volume)
meshCenterOfMass(faces) {
    let com = zero();
    for (const face of faces) {
        const [x0,y0,z0] = pointCoords(face.v0);
        const [x1,y1,z1] = pointCoords(face.v1);
        const [x2,y2,z2] = pointCoords(face.v2);
        
        const cx = x0 + x1 + x2;
        const cy = y0 + y1 + y2;
        const cz = z0 + z1 + z2;
        
        const F = joinPPP(face.v0, face.v1, face.v2);
        const signedVol = F[E0];
        
        com[E023] += -cx * signedVol;
        com[E013] += cy * signedVol;
        com[E012] += -cz * signedVol;
        com[E123] += 4 * signedVol;
    }
    return com;  // Normalize by pointCoords() to extract position
}
```

## Plane Operations

### Plane Side Test

```
p ∨ v = ax + by + cz + d
```

Where `p = ae₁ + be₂ + ce₃ + de₀` is a plane.

- Result > 0: point above plane
- Result < 0: point below plane
- Result = 0: point on plane

### Line-Plane Intersection

```
E ∧ p = intersection point
```

Where `E` is a line bivector and `p` is a plane vector.

```javascript
// Returns point of intersection or null
meetLinePlane(E, p) {
    const dx = E[E01], dy = E[E02], dz = E[E03];
    const lx = E[E23], ly = -E[E13], lz = E[E12];
    const a = p[E1], b = p[E2], c = p[E3], d = p[E0];
    
    return point(
        (b*dz - c*dy - d*lx) / (a*lx + b*ly + c*lz),
        (c*dx - a*dz - d*ly) / (a*lx + b*ly + c*lz),
        (a*dy - b*dx - d*lz) / (a*lx + b*ly + c*lz)
    );
}
```

## Mesh Slicing (Fuel Tank Problem)

To compute volume below a cutting plane:

1. **Classify triangles**: Use plane side test on each vertex
2. **Full triangles**: Add/skip based on all vertices above/below
3. **Split triangles**: 
   - Find edge-plane intersections
   - Construct sub-triangles
   - Add/subtract contributions

```javascript
// Optimized: only accumulate d-coefficient
for (const face of faces) {
    const d0 = planeSide(plane, face.v0);
    const d1 = planeSide(plane, face.v1);
    const d2 = planeSide(plane, face.v2);
    
    if (all below) {
        sumD += facePlane[E0];
    } else if (intersecting) {
        // Split and accumulate
    }
}
return sumD / 6;
```

## Inertia Tensor

The inertia tensor is accumulated as a frame of 3 vectors:

```
I₁ = (Iy + Iz)e₁ + Ixy·e₂ + Ixz·e₃
I₂ = Ixy·e₁ + (Iz + Ix)e₂ + Iyz·e₃
I₃ = Ixz·e₁ + Iyz·e₂ + (Ix + Iy)e₃
```

Accumulated over mesh faces with signed volume weighting.

### Diagonalization

Use Jacobi eigenvalue algorithm with Givens rotations to find:
- Principal moments `[λ₁, λ₂, λ₃]`
- Principal axes (as rotor)

## Norm Summary Table

| Element | ‖·‖ (Euclidean) | ‖·‖∞ (Ideal) |
|---------|-----------------|--------------|
| Point v | 1 | Distance to origin |
| Edge v₀∨v₁ | Edge length | 2× triangle area with origin |
| Triangle v₀∨v₁∨v₂ | 2× area | 6× tetra volume with origin |
| Tetrahedron | 6× volume | 0 |

## Implementation Files

- `ga_mesh.js` - Full PGA3D mesh processing implementation
- `mesh_demo.html` - Interactive visualization demo

## References

1. De Keninck, Roelfs, Dorst, Eelbode. "Clean up your Mesh! Part 1: Plane and Simplex" (2025)
2. Dorst, De Keninck. "A Guided Tour to the Plane-Based Geometric Algebra PGA" (2022)
3. Dorst, De Keninck. "Physical Geometry by Plane-Based Geometric Algebra" (2024)
4. Roelfs, De Keninck. "Graded Symmetry Groups: Plane and Simple" (2023)
