# Discrete Geodesics Reference

Based on Keenan Crane's "Discrete Differential Geometry: An Applied Introduction"  
Geodesics I & II Lectures  
**Version 1.0**

## Overview

Geodesics generalize straight lines to curved spaces. In Euclidean space, a straight line is both the **shortest** path (minimizing length) and the **straightest** path (no turning). On curved surfaces, these two properties define geodesics locally, but their discrete realizations diverge significantly.

## See also

- Mesh primitives used by these algorithms: [`mesh_reference.md`](mesh_reference.md)
- Discrete curvature operators: [`ddg_reference.md`](ddg_reference.md)
- Global topology/holonomy checks: [`unified_topology_reference.md`](unified_topology_reference.md)

## The Two Perspectives

| Perspective | Definition | Use Case | Discrete Issues |
|-------------|------------|----------|-----------------|
| **Shortest** | Critical points of length | Boundary value (A→B) | Avoids cone vertices |
| **Straightest** | Zero geodesic curvature | Initial value (start+dir) | Exp map not surjective |

---

## Part I: Shortest Geodesics (Variational)

### Variational Characterization

A geodesic is a critical point of the length functional:

```
L(γ) = ∫|γ'(t)| dt
```

Equivalently, minimizing **Dirichlet energy** (smoothness):

```
E(γ) = ∫|γ'(t)|² dt
```

produces geodesics with constant-speed parameterization.

### Discrete Shortest Paths

**Graph Distance (Dijkstra):** Poor approximation!
- Even on refined meshes, graph distance can differ 10-20% from true geodesic
- Only follows mesh edges, missing interior shortcuts

**Exact Polyhedral Geodesics (MMP Algorithm):**
- Paths can cross triangle interiors
- Track "windows" of common geodesic paths on edges
- Propagate wavefront across triangles
- O(n² log n) complexity

### Vertex Behavior (Shortest Geodesics)

| Vertex Type | Angle Sum | Behavior |
|-------------|-----------|----------|
| **Flat** (Θ = 2π) | Sum = 2π | Path goes straight through |
| **Cone** (Θ < 2π) | Sum < 2π | Shortest paths NEVER pass through |
| **Saddle** (Θ > 2π) | Sum > 2π | Many locally shortest paths may pass through |

**Key Insight:** It's always faster to go *around* a cone vertex than through it!

### Cut Locus

The **cut locus** C(p) is the set of points q where:
1. The shortest geodesic from p is not unique, OR
2. The geodesic ceases to be shortest

Examples:
- **Sphere:** Cut locus of any point is its antipode
- **Cylinder:** Cut locus is a line opposite the source
- **Polyhedra:** Often includes saddle vertices

Related: **Medial Axis** = set of points with non-unique closest boundary points

---

## Part II: Straightest Geodesics (Geometric)

### Geometric Definition

Decompose curve curvature into:

```
κ = κ_n + κ_g
```

- **κ_n** = normal curvature (forced by surface)
- **κ_g** = geodesic curvature (lateral bending)

**A geodesic has κ_g = 0** (no sideways turning).

### Discrete Straightest Condition

When a geodesic crosses an edge, the angles on either side must be equal:

```
θ_left = θ_right
```

**Equivalently:** The path unfolds to a straight line in 2D.

At vertices, the "straightest" direction is ambiguous—cannot simply unfold around a vertex with angle defect.

### Dynamic Perspective (Covariant Derivative)

A geodesic has zero tangential acceleration:

```
∇_γ' γ' = 0
```

where ∇ is the **covariant derivative** (Levi-Civita connection).

### Geodesic Equation (Local Coordinates)

```
ẍᵏ + Γᵏᵢⱼ ẋⁱ ẋʲ = 0
```

where Γ are Christoffel symbols encoding the connection.

---

## Part III: Exponential and Logarithmic Maps

### Exponential Map

**exp_p : T_pM → M**

Maps tangent vector v at p to point on manifold by:
1. Shoot geodesic from p in direction v
2. Walk for length |v|
3. Return endpoint

```
exp_p(v) = γ(1) where γ(0)=p, γ'(0)=v, |γ|=|v|
```

### Logarithmic Map

**log_p : M → T_pM**

Inverse of exponential map:

```
log_p(q) = v such that exp_p(v) = q
```

**Properties:**
- On smooth surfaces, log is well-defined except at cut locus
- On discrete surfaces, exp map is NOT surjective (can't reach all points)

### Discrete Challenges

On polyhedral surfaces:
- Straightest geodesics from a vertex may not cover entire surface
- Some points unreachable via straightest paths
- Cut locus structure more complex

---

## Part IV: Parallel Transport

### Definition

A vector V is **parallel transported** along curve γ if:

```
∇_γ' V = 0
```

(covariant derivative along curve is zero)

### Geometric Meaning

- Vector maintains its "intrinsic direction" relative to surface
- On polyhedra: rotate by dihedral angle when crossing edges
- A geodesic parallel transports its own tangent vector

### Holonomy

Parallel transport around a closed loop may not return the original vector:

```
Holonomy = rotation after full loop
```

For a vertex: **Holonomy = Angle Defect = Gaussian Curvature**

This connects local curvature to global transport!

---

## Part V: Karcher Mean

### Problem

Find "center of mass" for points on curved surface.

Euclidean mean doesn't work—average might not be on surface!

### Definition

```
x̄ = argmin_x Σ d²(x, xᵢ)
```

Minimizes sum of squared geodesic distances.

### Algorithm

1. Start with initial guess x
2. Compute log_x(xᵢ) for all points
3. Average in tangent space: v̄ = (1/n) Σ log_x(xᵢ)
4. Update: x ← exp_x(v̄)
5. Repeat until convergence

---

## Equations Summary

### Curvature Decomposition
```
κ = κ_n + κ_g
```

### Discrete Straightest Condition
```
θ_left = θ_right
```

### Exponential/Logarithmic Maps
```
exp_p : T_pM → M
log_p : M → T_pM
```

### Karcher Mean
```
x̄ = argmin_x Σ d²(x, xᵢ)
```

Update step:
```
x ← exp_x((1/n) Σ log_x(xᵢ))
```

### Geodesic Equation
```
∇_γ' γ' = 0
```

In coordinates:
```
ẍᵏ + Γᵏᵢⱼ ẋⁱ ẋʲ = 0
```

### Curvature & topology (pointer)

This file focuses on *geodesic algorithms*. For curvature and global constraints, use:
- [`ddg_reference.md`](ddg_reference.md) — angle defect / cotan Laplacian / curvature operators  
- [`unified_topology_reference.md`](unified_topology_reference.md) — Gauss–Bonnet, Euler characteristic, holonomy

---

## Comparison: Shortest vs Straightest

| Property | Shortest | Straightest |
|----------|----------|-------------|
| **Problem Type** | Boundary value | Initial value |
| **Cone Vertices** | Never passes through | May pass through |
| **Saddle Vertices** | May pass through | May pass through |
| **On Convex Polyhedra** | Is straightest | Not necessarily shortest |
| **Exp Map** | Surjective | Not surjective |
| **Computation** | Dijkstra/MMP | Edge unfolding |

---

## Implementation Notes

### JavaScript Module: `ga_geodesics.js`

```javascript
const geo = new Geodesics.DiscreteGeodesics();
geo.loadMesh(vertices, faces);

// Shortest path (graph approximation)
const { distances, predecessors } = geo.dijkstraGraph(source);
const path = geo.reconstructPath(predecessors, target);

// Vertex classification
geo.isConeVertex(vi);    // K > 0, shortest paths avoid
geo.isSaddleVertex(vi);  // K < 0, multiple paths may cross

// Exponential map
const result = geo.exponentialMap(vertex, tangentVector);

// Logarithmic map
const { tangent, distance } = geo.logarithmicMap(source, target);

// Parallel transport
const transported = geo.parallelTransport(vector, path);

// Cut locus approximation
const { cutVertices } = geo.approximateCutLocus(source);

// Karcher mean
const mean = geo.karcherMean(vertexIndices);

// Geodesic curvature of discrete curve
const curvatures = geo.geodesicCurvature(curveVertices);
```

---

## References

1. Crane, K. "Discrete Differential Geometry: An Applied Introduction" (2025)
2. Crane, K. "Geodesics I & II" lecture notes
3. Mitchell, Mount, Papadimitriou, "The Discrete Geodesic Problem" (1987)
4. Surazhsky et al., "Fast Exact and Approximate Geodesics on Meshes" (2005)
5. Polthier & Schmies, "Straightest Geodesics on Polyhedral Surfaces"
