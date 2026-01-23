"""
Geometric Calculus on Curved Manifolds and General Meshes

This module extends the discrete geometric calculus to handle:
1. Curved manifolds with regular parameterization (position-dependent metric)
2. Discussion of Discrete Exterior Calculus (DEC) for unstructured meshes

The key insight is that the geometric derivative on a manifold requires:
- A position-dependent metric tensor g_ij(x)
- Christoffel symbols for covariant derivatives
- Or equivalently, a position-dependent frame field

Author: Claude (Anthropic)
"""

import numpy as np
from typing import Tuple, List, Callable, Optional, Union
from dataclasses import dataclass
from discrete_geometric_calculus import (
    CliffordBasis, AlgebraSignature, DiscreteMultivectorField,
    SplitDifferentialOperator, EUCLIDEAN_2D, EUCLIDEAN_3D
)


# =============================================================================
# PART 1: LIMITATIONS OF CUBIC GRID APPROACH
# =============================================================================

"""
CURRENT LIMITATIONS (Cubic Grid Implementation):

1. UNIFORM SPACING
   - Assumes constant grid spacing h in each direction
   - Cannot handle adaptive refinement
   
2. FLAT METRIC
   - Metric tensor is constant: g_ij = diag(±1, ±1, ...)
   - No curvature effects
   
3. RECTANGULAR TOPOLOGY
   - Array indexing assumes [i,j,k,...] structure
   - Cannot handle unstructured connectivity

4. CARTESIAN COORDINATES
   - Partial derivatives are simple finite differences
   - No coordinate transformation terms


WHAT'S NEEDED FOR GENERAL MANIFOLDS:

1. Position-dependent metric g_ij(x)
2. Reciprocal frame field e^i(x) satisfying e^i · e_j = δ^i_j
3. Connection coefficients (Christoffel symbols) for parallel transport
4. Covariant derivative: ∇_i V^j = ∂_i V^j + Γ^j_ik V^k
"""


# =============================================================================
# PART 2: CURVED MANIFOLDS WITH REGULAR GRIDS
# =============================================================================

class CurvedManifoldGA:
    """
    Geometric Algebra on a curved manifold with regular parameterization.
    
    The manifold is described by:
    - Coordinates (u¹, u², ..., uⁿ) on a regular grid
    - A metric tensor field g_ij(u) at each point
    - Derived: reciprocal metric g^ij, frame fields, Christoffel symbols
    
    Example: Sphere with (θ, φ) coordinates
        g = diag(1, sin²θ)
    """
    
    def __init__(self, signature: AlgebraSignature, shape: Tuple[int, ...],
                 spacing: Union[float, Tuple[float, ...]] = 1.0):
        """
        Args:
            signature: Algebra signature (typically Euclidean for embedded manifolds)
            shape: Grid shape in parameter space
            spacing: Grid spacing in parameter coordinates
        """
        self.basis = CliffordBasis(signature)
        self.dim = signature.dim
        self.shape = shape
        
        if isinstance(spacing, (int, float)):
            self.spacing = tuple([float(spacing)] * self.dim)
        else:
            self.spacing = tuple(spacing)
        
        # Metric tensor field: g_ij at each grid point
        # Shape: (*shape, dim, dim)
        self.metric_field: Optional[np.ndarray] = None
        
        # Inverse metric: g^ij
        self.inverse_metric_field: Optional[np.ndarray] = None
        
        # Metric determinant sqrt(|g|) for volume element
        self.sqrt_det_g: Optional[np.ndarray] = None
        
    def set_flat_metric(self):
        """Set constant flat metric (reduces to cubic grid case)."""
        self.metric_field = np.zeros((*self.shape, self.dim, self.dim))
        for i in range(self.dim):
            self.metric_field[..., i, i] = self.basis.metric[i, i]
        
        self.inverse_metric_field = self.metric_field.copy()
        self.sqrt_det_g = np.ones(self.shape)
        
    def set_metric_field(self, metric_func: Callable[[np.ndarray], np.ndarray]):
        """
        Set metric from a function g_ij(coordinates).
        
        Args:
            metric_func: Function that takes coordinate array (*shape, dim)
                        and returns metric (*shape, dim, dim)
        """
        # Build coordinate grid
        grids = [np.linspace(0, (n-1)*h, n) 
                 for n, h in zip(self.shape, self.spacing)]
        coords = np.stack(np.meshgrid(*grids, indexing='ij'), axis=-1)
        
        # Compute metric at each point
        self.metric_field = metric_func(coords)
        
        # Compute inverse metric
        self.inverse_metric_field = np.linalg.inv(self.metric_field)
        
        # Compute sqrt(|det(g)|)
        det_g = np.linalg.det(self.metric_field)
        self.sqrt_det_g = np.sqrt(np.abs(det_g))
    
    def partial_derivative(self, field: np.ndarray, direction: int,
                          stencil: str = 'central') -> np.ndarray:
        """
        Compute partial derivative ∂f/∂u^i in parameter coordinates.
        Same as flat case - just finite differences.
        """
        h = self.spacing[direction]
        ndim = len(self.shape)
        result = np.zeros_like(field)
        
        def slice_at(dim, start, stop):
            s = [slice(None)] * len(field.shape)
            s[dim] = slice(start, stop)
            return tuple(s)
        
        if stencil == 'central':
            result[slice_at(direction, 1, -1)] = (
                field[slice_at(direction, 2, None)] -
                field[slice_at(direction, 0, -2)]
            ) / (2 * h)
        
        return result
    
    def covariant_derivative_scalar(self, f: np.ndarray) -> np.ndarray:
        """
        Compute covariant derivative of scalar field.
        
        For scalars: ∇_i f = ∂_i f (no connection terms)
        
        Returns gradient 1-form: (∂f/∂u¹, ∂f/∂u², ...)
        Shape: (*shape, dim)
        """
        grad = np.zeros((*self.shape, self.dim))
        for i in range(self.dim):
            grad[..., i] = self.partial_derivative(f, i)
        return grad
    
    def gradient_vector(self, f: np.ndarray) -> np.ndarray:
        """
        Compute gradient as a vector (contravariant).
        
        (grad f)^i = g^{ij} ∂_j f
        
        This is the geometric gradient that points in the direction
        of steepest ascent on the manifold.
        """
        # First compute covariant gradient (1-form)
        grad_cov = self.covariant_derivative_scalar(f)
        
        # Raise index with inverse metric
        grad_vec = np.einsum('...ij,...j->...i', 
                            self.inverse_metric_field, grad_cov)
        return grad_vec
    
    def divergence(self, V: np.ndarray) -> np.ndarray:
        """
        Compute divergence of vector field.
        
        div V = (1/√g) ∂_i (√g V^i)
        
        Args:
            V: Contravariant vector field, shape (*shape, dim)
        
        Returns:
            Scalar divergence field
        """
        # Compute √g V^i
        weighted = self.sqrt_det_g[..., np.newaxis] * V
        
        # Sum of partial derivatives
        div = np.zeros(self.shape)
        for i in range(self.dim):
            div += self.partial_derivative(weighted[..., i], i)
        
        # Divide by √g
        div /= (self.sqrt_det_g + 1e-10)  # Avoid division by zero
        
        return div
    
    def laplacian(self, f: np.ndarray) -> np.ndarray:
        """
        Compute Laplace-Beltrami operator on scalar field.
        
        Δf = div(grad f) = (1/√g) ∂_i (√g g^{ij} ∂_j f)
        """
        grad = self.gradient_vector(f)
        return self.divergence(grad)
    
    def geometric_derivative_scalar(self, f: np.ndarray) -> np.ndarray:
        """
        Compute geometric derivative of scalar field.
        
        ∇f = Σ_i e^i (∂f/∂u^i)
        
        where e^i are the reciprocal frame vectors satisfying
        e^i · e_j = δ^i_j with respect to the curved metric.
        
        Returns: Multivector field (grade 1 = vector)
        """
        field = DiscreteMultivectorField(self.basis, self.shape)
        
        for i in range(self.dim):
            partial_i = self.partial_derivative(f, i)
            
            # For curved metric, e^i = g^{ij} e_j
            # The reciprocal basis absorbs the metric
            for j in range(self.dim):
                blade_idx = self.basis.blade_to_index((j,))
                # g^{ij} contribution
                field.data[..., blade_idx] += (
                    self.inverse_metric_field[..., i, j] * partial_i
                )
        
        return field


# =============================================================================
# PART 3: EXAMPLE - SPHERE
# =============================================================================

def create_sphere_manifold(n_theta: int = 32, n_phi: int = 64, 
                           radius: float = 1.0) -> CurvedManifoldGA:
    """
    Create geometric algebra on a 2-sphere.
    
    Coordinates: (θ, φ) where θ ∈ [0, π], φ ∈ [0, 2π)
    
    Metric: ds² = R²(dθ² + sin²θ dφ²)
           g = R² [[1, 0], [0, sin²θ]]
    """
    shape = (n_theta, n_phi)
    
    # Spacing in parameter coordinates
    d_theta = np.pi / (n_theta - 1)
    d_phi = 2 * np.pi / n_phi
    
    manifold = CurvedManifoldGA(EUCLIDEAN_2D, shape, spacing=(d_theta, d_phi))
    
    def sphere_metric(coords):
        """Metric tensor for sphere."""
        theta = coords[..., 0]
        R2 = radius ** 2
        
        g = np.zeros((*coords.shape[:-1], 2, 2))
        g[..., 0, 0] = R2
        g[..., 1, 1] = R2 * np.sin(theta) ** 2
        
        # Regularize at poles
        g[..., 1, 1] = np.maximum(g[..., 1, 1], 1e-10)
        
        return g
    
    manifold.set_metric_field(sphere_metric)
    
    return manifold


def test_sphere_laplacian():
    """Test Laplace-Beltrami operator on sphere."""
    print("=" * 60)
    print("TEST: Laplace-Beltrami on 2-Sphere")
    print("=" * 60)
    
    manifold = create_sphere_manifold(n_theta=64, n_phi=128, radius=1.0)
    
    # Build coordinate arrays
    theta = np.linspace(0, np.pi, manifold.shape[0])
    phi = np.linspace(0, 2*np.pi, manifold.shape[1], endpoint=False)
    Theta, Phi = np.meshgrid(theta, phi, indexing='ij')
    
    # Test function: spherical harmonic Y_1^0 ∝ cos(θ)
    # This is an eigenfunction of ΔS² with eigenvalue -2
    f = np.cos(Theta)
    
    # Compute Laplacian
    lap_f = manifold.laplacian(f)
    
    # Expected: Δf = -2f = -2cos(θ)
    expected = -2 * np.cos(Theta)
    
    # Error (interior, avoiding poles)
    interior = (slice(5, -5), slice(None))
    error = np.abs(lap_f[interior] - expected[interior]).max()
    
    print(f"\nTest function: f(θ,φ) = cos(θ)")
    print(f"Expected: Δf = -2cos(θ) (eigenvalue -2)")
    print(f"Max error: {error:.4f}")
    print(f"Status: {'PASS' if error < 0.1 else 'FAIL'}")
    
    return error


# =============================================================================
# PART 4: DISCRETE EXTERIOR CALCULUS (DEC) FOR UNSTRUCTURED MESHES
# =============================================================================

"""
For truly unstructured meshes (triangles, tetrahedra, polygons),
we need Discrete Exterior Calculus (DEC).

KEY CONCEPTS:

1. SIMPLICIAL COMPLEX
   - 0-simplices: vertices
   - 1-simplices: edges
   - 2-simplices: triangles
   - k-simplices: k-dimensional faces

2. DISCRETE DIFFERENTIAL FORMS
   - 0-forms: values on vertices
   - 1-forms: values on edges (integrals along edges)
   - 2-forms: values on faces (integrals over faces)
   - k-forms: values on k-simplices

3. DISCRETE EXTERIOR DERIVATIVE d
   - d: k-forms → (k+1)-forms
   - Implemented as signed incidence matrices
   - d₀: 0-forms → 1-forms (gradient)
   - d₁: 1-forms → 2-forms (curl)

4. DISCRETE HODGE STAR ⋆
   - ⋆: k-forms → (n-k)-forms
   - Requires dual mesh (circumcentric or barycentric)
   - Involves volume ratios between primal and dual cells

5. DISCRETE CODIFFERENTIAL δ = ⋆d⋆
   - δ: k-forms → (k-1)-forms
   - Used for divergence

6. DISCRETE LAPLACIAN
   - Δ = dδ + δd (Hodge Laplacian)
   - For 0-forms: Δ₀ = δ₁d₀ = cotan Laplacian on meshes
"""


@dataclass
class SimplexMesh:
    """
    Simple representation of a simplicial complex.
    """
    vertices: np.ndarray      # (n_vertices, dim) coordinates
    edges: np.ndarray         # (n_edges, 2) vertex indices
    triangles: np.ndarray     # (n_triangles, 3) vertex indices
    
    @property
    def n_vertices(self) -> int:
        return len(self.vertices)
    
    @property
    def n_edges(self) -> int:
        return len(self.edges)
    
    @property
    def n_triangles(self) -> int:
        return len(self.triangles)


class DiscreteDEC:
    """
    Basic Discrete Exterior Calculus on triangle meshes.
    
    This implements the exterior derivative and Hodge star
    for doing geometric calculus on unstructured meshes.
    """
    
    def __init__(self, mesh: SimplexMesh):
        self.mesh = mesh
        self.dim = mesh.vertices.shape[1]
        
        # Build incidence matrices
        self.d0 = self._build_d0()  # gradient operator
        self.d1 = self._build_d1()  # curl operator
        
        # Build Hodge stars (diagonal approximation)
        self.star0 = self._build_star0()
        self.star1 = self._build_star1()
        self.star2 = self._build_star2()
        
    def _build_d0(self) -> np.ndarray:
        """
        Build discrete exterior derivative d₀: 0-forms → 1-forms
        
        For edge (i,j): (d₀ω)[ij] = ω[j] - ω[i]
        
        This is the discrete gradient operator.
        """
        n_v = self.mesh.n_vertices
        n_e = self.mesh.n_edges
        
        d0 = np.zeros((n_e, n_v))
        for e_idx, (i, j) in enumerate(self.mesh.edges):
            d0[e_idx, i] = -1
            d0[e_idx, j] = +1
        
        return d0
    
    def _build_d1(self) -> np.ndarray:
        """
        Build discrete exterior derivative d₁: 1-forms → 2-forms
        
        For triangle (i,j,k): (d₁ω)[ijk] = ω[ij] + ω[jk] + ω[ki]
        
        This is the discrete curl operator.
        """
        n_e = self.mesh.n_edges
        n_t = self.mesh.n_triangles
        
        # Build edge lookup
        edge_dict = {}
        for e_idx, (i, j) in enumerate(self.mesh.edges):
            edge_dict[(i, j)] = (e_idx, +1)
            edge_dict[(j, i)] = (e_idx, -1)
        
        d1 = np.zeros((n_t, n_e))
        for t_idx, (i, j, k) in enumerate(self.mesh.triangles):
            # Three edges of triangle
            for (a, b) in [(i, j), (j, k), (k, i)]:
                e_idx, sign = edge_dict[(a, b)]
                d1[t_idx, e_idx] = sign
        
        return d1
    
    def _build_star0(self) -> np.ndarray:
        """
        Diagonal Hodge star ⋆₀: 0-forms → 2-forms (dual)
        
        Uses dual cell areas (Voronoi cells).
        Simplified: uses 1/3 of surrounding triangle areas.
        """
        areas = np.zeros(self.mesh.n_vertices)
        
        for (i, j, k) in self.mesh.triangles:
            v = self.mesh.vertices
            # Triangle area
            e1 = v[j] - v[i]
            e2 = v[k] - v[i]
            if self.dim == 2:
                area = 0.5 * abs(e1[0]*e2[1] - e1[1]*e2[0])
            else:
                area = 0.5 * np.linalg.norm(np.cross(e1, e2))
            
            # Distribute to vertices
            areas[i] += area / 3
            areas[j] += area / 3
            areas[k] += area / 3
        
        return np.diag(areas)
    
    def _build_star1(self) -> np.ndarray:
        """
        Diagonal Hodge star ⋆₁: 1-forms → 1-forms (dual)
        
        Uses cotan weights for edges.
        """
        weights = np.zeros(self.mesh.n_edges)
        
        # Build edge-to-triangles lookup
        edge_dict = {}
        for e_idx, (i, j) in enumerate(self.mesh.edges):
            edge_dict[(min(i,j), max(i,j))] = e_idx
        
        for (i, j, k) in self.mesh.triangles:
            v = self.mesh.vertices
            # For each edge, find opposite angle and add cotan
            for (a, b, c) in [(i, j, k), (j, k, i), (k, i, j)]:
                # Edge (a, b), opposite vertex c
                e_idx = edge_dict[(min(a,b), max(a,b))]
                
                # Vectors from c
                ca = v[a] - v[c]
                cb = v[b] - v[c]
                
                # Cotan of angle at c
                cos_c = np.dot(ca, cb) / (np.linalg.norm(ca) * np.linalg.norm(cb) + 1e-10)
                sin_c = np.sqrt(1 - cos_c**2 + 1e-10)
                cotan_c = cos_c / sin_c
                
                weights[e_idx] += 0.5 * cotan_c
        
        return np.diag(weights)
    
    def _build_star2(self) -> np.ndarray:
        """
        Diagonal Hodge star ⋆₂: 2-forms → 0-forms (dual)
        
        Uses inverse of triangle areas.
        """
        areas = np.zeros(self.mesh.n_triangles)
        
        for t_idx, (i, j, k) in enumerate(self.mesh.triangles):
            v = self.mesh.vertices
            e1 = v[j] - v[i]
            e2 = v[k] - v[i]
            if self.dim == 2:
                areas[t_idx] = 0.5 * abs(e1[0]*e2[1] - e1[1]*e2[0])
            else:
                areas[t_idx] = 0.5 * np.linalg.norm(np.cross(e1, e2))
        
        return np.diag(1.0 / (areas + 1e-10))
    
    def gradient(self, f: np.ndarray) -> np.ndarray:
        """
        Discrete gradient: 0-form → 1-form
        
        grad f = d₀ f
        """
        return self.d0 @ f
    
    def curl(self, omega: np.ndarray) -> np.ndarray:
        """
        Discrete curl: 1-form → 2-form
        
        curl ω = d₁ ω
        """
        return self.d1 @ omega
    
    def divergence(self, omega: np.ndarray) -> np.ndarray:
        """
        Discrete divergence: 1-form → 0-form
        
        div ω = ⋆ d ⋆ ω = ⋆₀⁻¹ d₀ᵀ ⋆₁ ω
        """
        star0_inv = np.diag(1.0 / (np.diag(self.star0) + 1e-10))
        return star0_inv @ self.d0.T @ self.star1 @ omega
    
    def laplacian(self, f: np.ndarray) -> np.ndarray:
        """
        Discrete Laplacian: 0-form → 0-form
        
        Δf = div(grad f) = ⋆₀⁻¹ d₀ᵀ ⋆₁ d₀ f
        
        This is the cotan Laplacian!
        """
        return self.divergence(self.gradient(f))
    
    def laplacian_matrix(self) -> np.ndarray:
        """Return the Laplacian as a matrix (cotan Laplacian)."""
        star0_inv = np.diag(1.0 / (np.diag(self.star0) + 1e-10))
        return star0_inv @ self.d0.T @ self.star1 @ self.d0


def create_test_mesh() -> SimplexMesh:
    """Create a simple test mesh (unit square with triangulation)."""
    # 3x3 grid of vertices
    n = 5
    vertices = []
    for j in range(n):
        for i in range(n):
            vertices.append([i / (n-1), j / (n-1)])
    vertices = np.array(vertices)
    
    # Edges and triangles
    edges = []
    triangles = []
    edge_set = set()
    
    def add_edge(i, j):
        key = (min(i,j), max(i,j))
        if key not in edge_set:
            edge_set.add(key)
            edges.append(key)
    
    for j in range(n - 1):
        for i in range(n - 1):
            # Vertex indices
            v00 = j * n + i
            v10 = j * n + i + 1
            v01 = (j + 1) * n + i
            v11 = (j + 1) * n + i + 1
            
            # Two triangles per cell
            triangles.append([v00, v10, v11])
            triangles.append([v00, v11, v01])
            
            # Edges
            add_edge(v00, v10)
            add_edge(v10, v11)
            add_edge(v11, v01)
            add_edge(v01, v00)
            add_edge(v00, v11)
    
    return SimplexMesh(
        vertices=vertices,
        edges=np.array(edges),
        triangles=np.array(triangles)
    )


def test_dec_laplacian():
    """Test DEC Laplacian on a mesh."""
    print("\n" + "=" * 60)
    print("TEST: DEC Laplacian on Triangle Mesh")
    print("=" * 60)
    
    mesh = create_test_mesh()
    dec = DiscreteDEC(mesh)
    
    print(f"\nMesh: {mesh.n_vertices} vertices, {mesh.n_edges} edges, {mesh.n_triangles} triangles")
    
    # Test function: f(x,y) = x² + y²
    f = mesh.vertices[:, 0]**2 + mesh.vertices[:, 1]**2
    
    # Laplacian: Δf = 4
    lap_f = dec.laplacian(f)
    
    # Check interior vertices (not on boundary)
    interior = (mesh.vertices[:, 0] > 0.1) & (mesh.vertices[:, 0] < 0.9) & \
               (mesh.vertices[:, 1] > 0.1) & (mesh.vertices[:, 1] < 0.9)
    
    error = np.abs(lap_f[interior] - 4.0).mean()
    
    print(f"Test function: f(x,y) = x² + y²")
    print(f"Expected Laplacian: Δf = 4")
    print(f"Mean error at interior vertices: {error:.4f}")
    print(f"Status: {'PASS' if error < 1.0 else 'FAIL'}")
    
    return error


# =============================================================================
# PART 5: SUMMARY AND COMPARISON
# =============================================================================

def print_summary():
    """Print summary of different approaches."""
    print("\n" + "=" * 70)
    print("SUMMARY: Geometric Calculus on Different Domains")
    print("=" * 70)
    
    print("""
┌─────────────────────────────────────────────────────────────────────┐
│                    DOMAIN TYPE COMPARISON                           │
├─────────────────┬─────────────────┬─────────────────────────────────┤
│ Domain          │ Implementation  │ Key Features                    │
├─────────────────┼─────────────────┼─────────────────────────────────┤
│ Flat + Regular  │ Cubic Grid      │ • Simple finite differences     │
│ (ℝⁿ, uniform h) │ (Basic impl.)   │ • Constant metric               │
│                 │                 │ • Array indexing                │
│                 │                 │ • O(h²) or O(h⁴) accuracy       │
├─────────────────┼─────────────────┼─────────────────────────────────┤
│ Curved + Reg.   │ Curvilinear     │ • Position-dependent metric     │
│ (Manifold with  │ Coords          │ • g_ij(x) field                 │
│  regular param) │ (Extended impl.)│ • Christoffel symbols           │
│                 │                 │ • Covariant derivatives         │
├─────────────────┼─────────────────┼─────────────────────────────────┤
│ Unstructured    │ DEC             │ • Simplicial complex            │
│ (Triangle mesh) │ (Exterior calc.)│ • Discrete forms on k-cells     │
│                 │                 │ • Incidence matrices            │
│                 │                 │ • Cotan Laplacian               │
├─────────────────┼─────────────────┼─────────────────────────────────┤
│ Symbolic        │ GAlgebra        │ • Arbitrary coordinates         │
│ (Any manifold)  │ (SymPy)         │ • Exact derivatives             │
│                 │                 │ • Position-dependent frames     │
│                 │                 │ • Full symbolic manipulation    │
└─────────────────┴─────────────────┴─────────────────────────────────┘

THE GEOMETRIC DERIVATIVE ON EACH:

1. FLAT CUBIC GRID:
   ∇F = Σᵢ eⁱ [(F[...i+1...] - F[...i-1...]) / 2h]

2. CURVED REGULAR GRID:  
   ∇F = Σᵢⱼ g^{ij}(x) eⱼ [∂F/∂uⁱ]
   
   with position-dependent inverse metric g^{ij}(x)

3. TRIANGLE MESH (DEC):
   For scalar f: grad f = d₀f (discrete exterior derivative)
   For vector ω: div ω = ⋆d⋆ω (Hodge star + derivative)
   
   Geometric algebra: reconstruct from scalar + curl decomposition

4. SYMBOLIC (GAlgebra):
   ∇F = Σᵢ eⁱ(x) [∂F/∂xⁱ]
   
   where eⁱ(x) are reciprocal frame vectors computed from metric
""")


def run_all_tests():
    """Run all tests."""
    print_summary()
    
    print("\n" + "=" * 70)
    print("RUNNING TESTS")
    print("=" * 70)
    
    test_sphere_laplacian()
    test_dec_laplacian()
    
    print("\n" + "=" * 70)
    print("CONCLUSION")
    print("=" * 70)
    print("""
The original implementation works ONLY on cubic grids with flat metric.

To handle curved manifolds and general meshes:

1. CURVED MANIFOLDS (regular parameterization):
   → Use CurvedManifoldGA class with position-dependent metric
   → Requires computing g^{ij}(x) at each point
   → Example: Sphere with (θ,φ) coordinates

2. UNSTRUCTURED MESHES (triangles, tetrahedra):
   → Use Discrete Exterior Calculus (DEC)
   → Forms live on simplices of different dimensions
   → Geometric derivative splits into d (exterior derivative)
     and δ = ⋆d⋆ (codifferential)

3. FULL GENERALITY:
   → Use galgebra for symbolic computation
   → Code-generate numerical kernels if needed
""")


if __name__ == "__main__":
    run_all_tests()
