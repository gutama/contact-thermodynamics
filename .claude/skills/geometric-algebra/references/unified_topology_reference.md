# Unified Framework: Geodesics, Curvature, and Topology

This note is a **global “sanity layer”**: it connects *local motion* (geodesics), *local curvature* (DDG operators), and *global invariants* (Euler characteristic, index sums, holonomy).

## See also

- Discrete curvature operators and derivations: [`ddg_reference.md`](ddg_reference.md)  
- Discrete geodesic algorithms (shortest/straightest, exp/log, transport): [`geodesics_reference.md`](geodesics_reference.md)  
- Mesh primitives and simplex measures (PGA): [`mesh_reference.md`](mesh_reference.md)  
- Core GA identities: [`formulas.md`](formulas.md)  

---

## The “one picture” view

A useful mental chain (read left→right):

1) **Geodesics (local motion)**  
A curve is geodesic when its tangent is parallel-transported along itself.

2) **Connection / curvature (local geometry)**  
Curvature is the obstruction to “turning coordinates into straight lines everywhere”.

3) **Topology (global constraints)**  
Some integrals/sums of curvature are **fixed by topology**, regardless of how you bend the surface.

This is why you can use global invariants as debugging checks for meshes and operators.

---

## Local motion: the geodesic condition

Continuous form (on a smooth surface):
- “No tangential acceleration” / “tangent transported along itself”:

```
∇_v v = 0
```

Equivalent coordinate form:

```
d²xᵏ/ds² + Γᵏᵢⱼ (dxᶦ/ds)(dxʲ/ds) = 0
```

On triangle meshes you typically use **either**:
- *Shortest* geodesics (variational / distance-minimizing), or
- *Straightest* geodesics (angle-preserving under unfolding).

Details and discrete algorithms are in [`geodesics_reference.md`](geodesics_reference.md).

---

## Curvature as “failure to close”

A very practical interpretation (especially for debugging):

- Walk around a small loop.
- Transport a direction (or a frame).
- If you come back rotated, that rotation is **curvature integrated over the enclosed region**.

This shows up as **holonomy**.

### Holonomy (discrete intuition)

For a small loop around a vertex on a triangulated surface, the curvature is captured by **angle defect**:

```
Θ_v = 2π - Σ_f θ_f(v)
```

- Θ_v > 0: locally sphere-like (elliptic)
- Θ_v < 0: locally saddle-like (hyperbolic)

Angle defect and its relation to DDG operators is detailed in [`ddg_reference.md`](ddg_reference.md).

---

## Global invariant 1: Gauss–Bonnet (Euler characteristic)

For a closed surface (no boundary), discrete Gauss–Bonnet becomes:

```
Σ_v Θ_v = 2π χ
```

Where χ is the Euler characteristic:
- Sphere: χ = 2
- Torus: χ = 0
- Genus g surface: χ = 2 - 2g

### Why this matters in practice
If your curvature implementation is correct, the sum of vertex angle defects should be close to `2πχ` (up to mesh resolution / numerical error). If it’s far off, something is broken: normals, face angles, vertex neighborhoods, or inconsistent orientation.

---

## Global invariant 2: Poincaré–Hopf (index sum)

For a vector field with isolated singularities on a closed surface:

```
Σ_p index(p) = χ
```

Implications:
- On a sphere (χ=2), any continuous tangent vector field must have singularities (hairy ball).
- On a torus (χ=0), you can have a singularity-free tangent vector field.

This is most useful as a conceptual guide; for numeric checks, Gauss–Bonnet is usually simpler.

---

## Boundary case: surfaces with boundary

If the surface has boundary, Gauss–Bonnet includes boundary turning / geodesic curvature. In practice:
- If you run checks on open meshes, prefer either:
  - cap/close the mesh before computing χ checks, or
  - include boundary terms explicitly (advanced; out of scope here).

---

## Quick debugging checklist (mesh)

When topology checks fail, the usual culprits are:

- **Non-manifold edges / vertices** (broken neighborhood rings)  
- **Inconsistent orientation** (face winding mixed)  
- **Degenerate triangles** (near-zero area)  
- **Boundary treated as closed** (missing boundary terms)

For mesh hygiene and PGA representations, see [`mesh_reference.md`](mesh_reference.md).

---

## Implementation hooks in this Skill

You can validate the global picture with the bundled demos/scripts:

- `scripts/geodesics_topology_demo.html` — combined geodesics + topology visualization  
- `scripts/ga_discrete_curvature.js` — angle defect / cotan-Laplace curvature operators  
- `scripts/mesh_demo.html` — mesh primitives and numeric checks

---

## Minimal equation sheet (for copy/paste)

```
# Local
∇_v v = 0

# Discrete curvature at vertex v (triangle mesh)
Θ_v = 2π - Σ_f θ_f(v)

# Global (closed surface)
Σ_v Θ_v = 2π χ
```

---

## References (starting points)

- Crane, Keenan — Discrete Differential Geometry (applied introductions / lecture notes)
- Classical topology references on Gauss–Bonnet and Poincaré–Hopf (any standard differential geometry/topology text)
