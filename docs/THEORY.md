# Mathematical Theory

A comprehensive guide to the mathematical foundations of Contact Thermodynamics.

## Table of Contents

1. [Contact Geometry Fundamentals](#1-contact-geometry-fundamentals)
2. [1-Jet Bundles](#2-1-jet-bundles)
3. [The Canonical Contact Form](#3-the-canonical-contact-form)
4. [Contact Hamiltonian Systems](#4-contact-hamiltonian-systems)
5. [Legendrian Submanifolds](#5-legendrian-submanifolds)
6. [Thermodynamic Interpretation](#6-thermodynamic-interpretation)
7. [Gravitational Extension](#7-gravitational-extension)
8. [Riemannian Geometry via Geometric Algebra](#8-riemannian-geometry-via-geometric-algebra)
9. [Information Geometry](#9-information-geometry)

---

## 1. Contact Geometry Fundamentals

### Definition

A **contact manifold** is a pair (M, ξ) where M is a (2n+1)-dimensional manifold and ξ is a maximally non-integrable hyperplane distribution.

Locally, ξ = ker(α) for some 1-form α satisfying the **contact condition**:

$$\alpha \wedge (d\alpha)^n \neq 0$$

This form α is called a **contact form**.

### Comparison with Symplectic Geometry

| Symplectic | Contact |
|------------|---------|
| Even dimension 2n | Odd dimension 2n+1 |
| Closed 2-form ω with ω^n ≠ 0 | 1-form α with α ∧ (dα)^n ≠ 0 |
| Phase space T*Q | Extended phase space T*Q × ℝ or J¹(Q) |
| Conserves volume | Allows dissipation |

### Why Contact for Thermodynamics?

Thermodynamic systems often involve:
- An odd number of variables (intensive + extensive + potential)
- Dissipation and entropy production
- Legendre transformations between potentials

Contact geometry naturally accommodates all of these.

---

## 2. 1-Jet Bundles

### Definition

Given a manifold Q, the **1-jet bundle** J¹(Q) consists of equivalence classes of functions at points of Q, where two functions are equivalent if they have the same value and first derivatives.

Concretely, if Q has local coordinates (x¹, ..., xⁿ), then J¹(Q) has coordinates:

$$(x^a, u, p_a)$$

where:
- x^a are the base coordinates
- u represents the function value
- p_a represent the partial derivatives ∂u/∂x^a

### Dimension Formula

$$\dim J^1(Q) = n + 1 + n = 2n + 1$$

This is exactly the dimension needed for a contact manifold!

### Canonical Contact Structure

J¹(Q) carries a canonical contact form:

$$\alpha = du - p_a \, dx^a$$

This satisfies the contact condition automatically:
$$d\alpha = -dp_a \wedge dx^a$$
$$\alpha \wedge (d\alpha)^n = n! \, du \wedge dx^1 \wedge dp_1 \wedge \cdots \wedge dx^n \wedge dp_n \neq 0$$

---

## 3. The Canonical Contact Form

### Structure

$$\boxed{\alpha = du - p_a \, dx^a}$$

### Properties

1. **Kernel**: The contact distribution ξ = ker(α) is spanned by vectors satisfying:
   $$du = p_a \, dx^a$$

2. **Reeb Vector Field**: The unique vector R satisfying:
   $$\alpha(R) = 1, \qquad \iota_R \, d\alpha = 0$$
   For canonical form: R = ∂/∂u

3. **Darboux Theorem**: Any contact form can be locally transformed to the canonical form via contactomorphism.

### Physical Interpretation

| Symbol | Thermodynamics | Wave Mechanics |
|--------|----------------|----------------|
| x^a | Extensive variables | Configuration |
| u | Potential (F, G, H, ...) | Action S |
| p_a | Intensive variables | Canonical momenta |

---

## 4. Contact Hamiltonian Systems

### The Contact Hamiltonian Vector Field

Given a smooth function H: M → ℝ (the Hamiltonian), there exists a unique vector field X_H satisfying:

$$\iota_{X_H} \alpha = -H$$
$$\iota_{X_H} d\alpha = dH - (RH) \alpha$$

where R is the Reeb field and RH = α(∇H) = ∂H/∂u.

### Explicit Form

For canonical coordinates:

$$\dot{x}^a = \frac{\partial H}{\partial p_a}$$

$$\dot{p}_a = -\frac{\partial H}{\partial x^a} - p_a \cdot \frac{\partial H}{\partial u}$$

$$\dot{u} = p_a \frac{\partial H}{\partial p_a} - H$$

### Comparison with Symplectic Case

| Symplectic | Contact |
|------------|---------|
| ẋ = ∂H/∂p | ẋ = ∂H/∂p |
| ṗ = -∂H/∂x | ṗ = -∂H/∂x - p · ∂H/∂u |
| — | u̇ = p·∂H/∂p - H |

The extra term -p · ∂H/∂u in the momentum equation allows for **dissipation**.

### Energy Evolution

In symplectic mechanics, H is conserved along trajectories. In contact mechanics:

$$\frac{dH}{dt} = \{H, H\}_{contact} + (\text{Reeb term})$$

The Hamiltonian generally **changes** along the flow unless special conditions are met.

---

## 5. Legendrian Submanifolds

### Definition

A **Legendrian submanifold** L ⊂ M is an n-dimensional submanifold such that:

$$\alpha|_L = 0$$

This is the maximal dimension for which this can hold.

### Generating Functions

The most important construction uses a **generating function** A(x):

$$u = A(x), \qquad p_a = \frac{\partial A}{\partial x^a}$$

This automatically satisfies the Legendrian condition:

$$\alpha|_L = dA - \frac{\partial A}{\partial x^a} dx^a = 0 \quad ✓$$

### Hamilton-Jacobi Theory

A function A(x) generates a Legendrian that solves the Hamilton-Jacobi equation if:

$$H\left(x, A(x), \frac{\partial A}{\partial x}\right) = 0$$

This is the **Hamilton-Jacobi PDE**, fundamental to classical mechanics and optics.

### Wavefronts

In geometric optics, Legendrian submanifolds represent **wavefronts**. The generating function A is the **optical path length** (eikonal).

---

## 6. Thermodynamic Interpretation

### The Grand Model M₁₃

The full thermodynamic phase space includes:

**Base Variables Q₆:**
- Spatial position: q¹, q², q³
- Time: t
- Scale: ℓ = log(λ)
- Entropy: S

**Conjugate Momenta:**
- Wavenumber: k₁, k₂, k₃ (momentum ↔ position)
- Frequency: ω (energy ↔ time)
- Dilatation: Δ (scaling dimension ↔ scale)
- Temperature: T (heat ↔ entropy)

**Contact Form:**
$$\alpha = dA - k_i \, dq^i - \omega \, dt - \Delta \, d\ell - T \, dS$$

### Thermodynamic Potentials

Different Legendrian submanifolds correspond to different thermodynamic potentials:

| Generating Function | Physical Potential |
|--------------------|-------------------|
| A(S, V, N) | Internal Energy U |
| A(T, V, N) | Helmholtz Free Energy F |
| A(T, P, N) | Gibbs Free Energy G |
| A(S, P, N) | Enthalpy H |

### Legendre Transforms

The transition between potentials via Legendre transform corresponds to **contact transformations** that preserve the contact structure.

### The Holographic Model M₇

When spatial degrees of freedom are "integrated out" or treated as dependent, we get the reduced model:

**Base Q₃:** (t, ℓ, S)

**Emergent Space:** q^i = q^i(t, ℓ, S)

This is analogous to:
- Holographic principle in gravity
- Effective field theory after coarse-graining
- Thermodynamic limit

---

## 7. Gravitational Extension

### The Principle

The contact structure α provides **kinematics** (locally Darboux-flat).
Spacetime **curvature** enters through the Hamiltonian constraint.

### Relativistic Hamiltonian

For a particle of mass m in spacetime with metric g_μν:

$$H = \frac{1}{2} g^{\mu\nu}(x) p_\mu p_\nu - \frac{1}{2} m^2 = 0$$

This is the **mass-shell constraint**.

### With Electromagnetic Field

Adding gauge field A_μ:

$$H = \frac{1}{2} g^{\mu\nu}(x)(p_\mu - qA_\mu)(p_\nu - qA_\nu) - \frac{1}{2} m^2 = 0$$

### Hamilton-Jacobi in Curved Spacetime

The relativistic Hamilton-Jacobi equation:

$$\frac{1}{2} g^{\mu\nu}(x) \frac{\partial S}{\partial x^\mu} \frac{\partial S}{\partial x^\nu} - \frac{1}{2} m^2 = 0$$

This determines geodesics via the method of characteristics.

### Key Insight

The separation:

$$\boxed{\text{Kinematics from } \alpha \quad | \quad \text{Curvature from } g_{\mu\nu}(x) \text{ inside } H}$$

This is analogous to how:
- Contact structure → canonical phase space structure
- Metric → specific dynamics (geodesics, forces)

---

## Summary: The Three Models

| Model | Base Dim | Total Dim | Coordinates | Fields |
|-------|----------|-----------|-------------|--------|
| Grand | 6 | 13 | q^i, t, ℓ, S, A | None |
| Holographic | 3 | 7 | t, ℓ, S, A | q^i(t,ℓ,S) |
| Gauge-Extended | 7 | 15 | + (φ, I) | None |

### Physics Caveats

1. **Entropy Closure**: Using S as base coordinate assumes local thermodynamic equilibrium (LTE) or adiabatic regime.

2. **Gauge Extension**: Adding (φ, I) increases dimension by 2 (one canonical pair).

---

## Discrete Geometric Calculus on Meshes

The library includes a full implementation of the **Fundamental Theorem of Geometric Calculus (FTGC)** on triangle meshes, providing discrete differential operators.

### The Fundamental Theorem

The FTGC states:

$$\boxed{\int_M \nabla F = \oint_{\partial M} F}$$

This unifies Stokes/Gauss/Green theorems under one operator ∇.

### Discrete Exterior Calculus (DEC)

On a triangle mesh, differential forms are discretized:

| Continuous | Discrete Location | Description |
|------------|------------------|-------------|
| 0-forms (scalars) | Vertices | One value per vertex |
| 1-forms (vectors) | Edges | Signed circulation |
| 2-forms (bivectors) | Faces | Area-weighted flux |

This is **staggered storage** — the natural home for each grade.

### The Cotan Laplacian

The mesh metric is encoded by **cotan weights**:

$$w_{ij} = \frac{1}{2}(\cot \alpha_{ij} + \cot \beta_{ij})$$

where α and β are angles opposite edge (i,j) in adjacent triangles.

The discrete Laplacian:

$$(\Delta f)_i = \frac{1}{A_i} \sum_{j \sim i} w_{ij}(f_j - f_i)$$

where $A_i$ is the dual area at vertex i.

### Mixed Voronoi Dual Areas

For **non-obtuse** triangles, use circumcentric (Voronoi) areas:

$$A_i^{(k)} = \frac{1}{8}(|e_{ij}|^2 \cot\gamma_k + |e_{ik}|^2 \cot\beta_j)$$

For **obtuse** triangles, fall back to barycentric (⅓ of face area).

This ensures ∑ᵢ Aᵢ = total mesh area.

### Discrete Differential Operators

The geometric derivative splits into:

| Operator | Input | Output | Formula |
|----------|-------|--------|---------|
| grad | 0-form (vertices) | 1-form (edges) | ∇∧f |
| curl | 1-form (edges) | 2-form (faces) | ∇∧V |
| div | 1-form (edges) | 0-form (vertices) | ∇·V |
| laplacian | 0-form | 0-form | ∇·∇f |

**Key identities hold discretely:**
- curl(grad f) = 0
- laplacian(constant) = 0

### Leapfrog Time Integration

For wave/heat/Maxwell equations, we use **leapfrog** (Verlet) integration:

**Wave equation** (uₜₜ = c²∇²u):

$$u^{n+1} = 2u^n - u^{n-1} + c^2 \Delta t^2 \nabla^2 u^n$$

CFL stability: $\Delta t < h_{\min} / c$

**Heat equation** (uₜ = α∇²u):

- Explicit: $u^{n+1} = u^n + \alpha \Delta t \nabla^2 u^n$ (CFL-limited)
- Implicit: $(I - \alpha \Delta t \nabla^2) u^{n+1} = u^n$ (unconditionally stable)

**Maxwell equations** (Yee-style staggered):

$$E^{n+1} = E^n + c \Delta t \nabla \times B^n$$
$$B^{n+1} = B^n - c \Delta t \nabla \times E^{n+1}$$

### Dirichlet Boundary Conditions

Fixed boundary values are enforced via masking:

```javascript
const mask = boundaryDirichletMask(mesh);  // boundary = 1
applyDirichlet(u, mask, boundaryValues);    // fix boundary vertices
```

The Laplacian solver respects these constraints, setting ∇²u = 0 at Dirichlet vertices.

---

## 8. Riemannian Geometry via Geometric Algebra

The library implements **coordinate-free Riemannian geometry** using Geometric Algebra (GA), eliminating Christoffel symbols in favor of bivector-valued forms.

### The Connection Bivector

Instead of Christoffel symbols Γᵏᵢⱼ, we use the **connection bivector**:

$$\boxed{\omega_i = \frac{1}{2} e^j \wedge \partial_i e_j}$$

This encodes how the tangent frame rotates as you move along the manifold.

### Parallel Transport

For a vector v transported along direction u:

$$\nabla_u v = u \cdot \partial v + \omega(u) \times v$$

The commutator ω × v performs infinitesimal rotation in the plane of ω.

### Curvature 2-Form

The curvature emerges as:

$$\boxed{\Omega = d\omega + \omega \wedge \omega}$$

This is the **Cartan structure equation**. For 2D surfaces, this reduces to Gaussian curvature K.

### Key Results

| Manifold | Curvature |
|----------|-----------|
| Sphere S²(R) | K = 1/R² (constant positive) |
| Torus T² | K = cos(θ) / (r(R + r cos θ)) (variable) |
| Hyperbolic Plane H² | K = −1 (constant negative) |

### Gauss-Bonnet Theorem

The total curvature is a topological invariant:

$$\int_M K \, dA = 2\pi \chi(M)$$

where χ = V − E + F is the Euler characteristic. Verified numerically in tests.

### Holonomy

Parallel transporting around a closed loop at colatitude θ on the sphere yields rotation:

$$\text{Holonomy angle} = 2\pi(1 - \cos\theta) = \iint K \, dA$$

---

## 9. Information Geometry

Probability distributions form a **contact manifold** with coordinates:

| Variable | Description | Role |
|----------|-------------|------|
| qⁱ | Probabilities | Position (base) |
| pᵢ = −ln(qⁱ) − 1 | Surprisal | Momentum (conjugate) |
| S = −Σ qⁱ ln(qⁱ) | Shannon Entropy | Action |

### Contact Form

$$\alpha = dS - p_i \, dq^i$$

This is the canonical contact form on the probability simplex.

### Legendrian Condition

Probability distributions satisfy α|_L = 0 automatically:

$$dS = p_i \, dq^i \implies S = -\sum_i q^i \ln q^i$$

Each distribution is a point on a **Legendrian submanifold**.

### Fisher Metric

The natural Riemannian metric on probability space is the **Fisher information**:

$$g_{ij} = \mathbb{E}\left[\frac{\partial \ln p}{\partial q^i} \frac{\partial \ln p}{\partial q^j}\right] = \frac{\delta_{ij}}{q^i}$$

This is related to the Hessian of entropy.

### References (Baez)

Based on John Baez's "Information Geometry" series:
- Part 18: Legendre transform and contact geometry
- Part 19: Surprisal as conjugate variable

---

## References

1. Arnold, V.I. — *Mathematical Methods of Classical Mechanics*
2. Libermann & Marle — *Symplectic Geometry and Analytical Mechanics*
3. Geiges, H. — *An Introduction to Contact Topology*
4. Bravetti, A. — "Contact Hamiltonian Dynamics" (Entropy, 2017)
5. Mrugała, R. — "Geometric formulation of thermodynamics"
6. Crane, K. — *Discrete Differential Geometry: An Applied Introduction* (CMU)
7. Desbrun, M. et al. — "Discrete Exterior Calculus" (2005)
8. Hestenes, D. & Sobczyk, G. — *Clifford Algebra to Geometric Calculus*
9. Doran, C. & Lasenby, A. — *Geometric Algebra for Physicists*
10. Baez, J. — "Information Geometry" blog series, The n-Category Café
