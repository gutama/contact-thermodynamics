"""
Coordinate-Free Riemannian Geometry via Geometric Calculus

This module extends covariant_derivative.py with pure GA formulations that
eliminate Christoffel symbols entirely. All curvature computations use:
- Connection bivector ω (replaces Γᵏᵢⱼ)
- Curvature 2-form Ω (replaces Rˡᵢⱼₖ)
- Commutator product × for bivector action

Key insight: Christoffel symbols are coordinate-dependent components of the
connection bivector. By working with ω directly, we gain:
- Geometric clarity (bivector = rotation generator)
- Computational efficiency (bivector algebra parallelizes)
- Automatic covariance (no index transformation rules)

Master Equations:
    ∂ᵢeⱼ = ωᵢ × eⱼ           (frame rotation)
    Ω = dω + ω ∧ ω           (Cartan structure equation)
    ∇ᵤv = ∂ᵤv + ω(u) × v     (covariant derivative)
    v̇ + ω(v) × v = 0         (geodesic equation)

Author: Claude (Anthropic) — Geometric Algebra Skill Extension
Version: 1.0
"""

import numpy as np
from typing import Tuple, Callable, Optional, List, Union
from dataclasses import dataclass, field
from enum import Enum


# =============================================================================
# BIVECTOR ALGEBRA
# =============================================================================

class Bivector2D:
    """
    Bivector in 2D (single component: e₁₂).
    
    In 2D, bivectors are pseudoscalars — just one component representing
    oriented area in the e₁∧e₂ plane.
    """
    def __init__(self, e12: float = 0.0):
        self.e12 = e12
    
    def __repr__(self):
        return f"Bivector2D({self.e12:.6f} e₁₂)"
    
    def __add__(self, other: 'Bivector2D') -> 'Bivector2D':
        return Bivector2D(self.e12 + other.e12)
    
    def __sub__(self, other: 'Bivector2D') -> 'Bivector2D':
        return Bivector2D(self.e12 - other.e12)
    
    def __mul__(self, scalar: float) -> 'Bivector2D':
        return Bivector2D(self.e12 * scalar)
    
    def __rmul__(self, scalar: float) -> 'Bivector2D':
        return self.__mul__(scalar)
    
    def __neg__(self) -> 'Bivector2D':
        return Bivector2D(-self.e12)
    
    @property
    def norm(self) -> float:
        return abs(self.e12)
    
    def commutator_with_vector(self, v: np.ndarray) -> np.ndarray:
        """
        Compute ω × v = ½(ωv - vω) for 2D.
        
        In 2D, the bivector e₁₂ acts on vectors by 90° rotation:
        e₁₂ × e₁ = e₂,  e₁₂ × e₂ = -e₁
        
        So for ω = ω₁₂ e₁₂ and v = (v¹, v²):
        ω × v = ω₁₂ (v² e₁ - v¹ e₂) = ω₁₂ (v², -v¹)
        
        Wait, let me recalculate. The commutator product:
        B × v = ½(Bv - vB)
        
        For B = b e₁₂ and v = v¹e₁ + v²e₂:
        Bv = b e₁₂ (v¹e₁ + v²e₂) = b(v¹ e₁₂e₁ + v² e₁₂e₂)
           = b(v¹ (-e₂) + v² e₁) = b(-v¹e₂ + v²e₁)
        vB = (v¹e₁ + v²e₂) b e₁₂ = b(v¹ e₁e₁₂ + v² e₂e₁₂)
           = b(v¹ e₂ + v² (-e₁)) = b(v¹e₂ - v²e₁)
        
        B × v = ½(Bv - vB) = ½ b [(-v¹e₂ + v²e₁) - (v¹e₂ - v²e₁)]
              = ½ b [-2v¹e₂ + 2v²e₁] = b(v²e₁ - v¹e₂)
              = b(v², -v¹)
        
        Actually in 2D: e₁₂ × e₁ = ½(e₁₂e₁ - e₁e₁₂) = ½(-e₂ - e₂) = -e₂
        Hmm, let me be more careful.
        
        e₁e₂ = e₁₂, so e₁₂e₁ = e₁e₂e₁ = -e₁e₁e₂ = -e₂
        and e₁e₁₂ = e₁e₁e₂ = e₂
        
        So e₁₂ × e₁ = ½(e₁₂e₁ - e₁e₁₂) = ½(-e₂ - e₂) = -e₂
        
        Similarly e₁₂e₂ = e₁e₂e₂ = e₁
        e₂e₁₂ = e₂e₁e₂ = -e₁e₂e₂ = -e₁
        e₁₂ × e₂ = ½(e₁ - (-e₁)) = e₁
        
        So for B = b e₁₂:
        B × (v¹e₁ + v²e₂) = v¹(b × e₁) + v²(b × e₂)
                          = v¹ b (-e₂) + v² b (e₁)
                          = b(v²e₁ - v¹e₂)
                          = b(v², -v¹)
        """
        return self.e12 * np.array([v[1], -v[0]])


class Bivector3D:
    """
    Bivector in 3D (three components: e₂₃, e₃₁, e₁₂).
    
    Represents oriented planes. The components correspond to:
    - e₂₃: YZ plane (rotation around X)
    - e₃₁: ZX plane (rotation around Y)  
    - e₁₂: XY plane (rotation around Z)
    
    Isomorphic to rotations / angular velocity / magnetic field.
    """
    def __init__(self, e23: float = 0.0, e31: float = 0.0, e12: float = 0.0):
        self.e23 = e23  # YZ plane
        self.e31 = e31  # ZX plane
        self.e12 = e12  # XY plane
    
    @classmethod
    def from_array(cls, arr: np.ndarray) -> 'Bivector3D':
        return cls(arr[0], arr[1], arr[2])
    
    @classmethod
    def from_wedge(cls, u: np.ndarray, v: np.ndarray) -> 'Bivector3D':
        """Create bivector from wedge product u ∧ v."""
        return cls(
            e23=u[1]*v[2] - u[2]*v[1],
            e31=u[2]*v[0] - u[0]*v[2],
            e12=u[0]*v[1] - u[1]*v[0]
        )
    
    def to_array(self) -> np.ndarray:
        return np.array([self.e23, self.e31, self.e12])
    
    def __repr__(self):
        return f"Bivector3D({self.e23:.4f} e₂₃ + {self.e31:.4f} e₃₁ + {self.e12:.4f} e₁₂)"
    
    def __add__(self, other: 'Bivector3D') -> 'Bivector3D':
        return Bivector3D(
            self.e23 + other.e23,
            self.e31 + other.e31,
            self.e12 + other.e12
        )
    
    def __sub__(self, other: 'Bivector3D') -> 'Bivector3D':
        return Bivector3D(
            self.e23 - other.e23,
            self.e31 - other.e31,
            self.e12 - other.e12
        )
    
    def __mul__(self, scalar: float) -> 'Bivector3D':
        return Bivector3D(
            self.e23 * scalar,
            self.e31 * scalar,
            self.e12 * scalar
        )
    
    def __rmul__(self, scalar: float) -> 'Bivector3D':
        return self.__mul__(scalar)
    
    def __neg__(self) -> 'Bivector3D':
        return Bivector3D(-self.e23, -self.e31, -self.e12)
    
    @property
    def norm_squared(self) -> float:
        return self.e23**2 + self.e31**2 + self.e12**2
    
    @property
    def norm(self) -> float:
        return np.sqrt(self.norm_squared)
    
    def commutator_with_vector(self, v: np.ndarray) -> np.ndarray:
        """
        Compute B × v = ½(Bv - vB) in 3D.
        
        The bivector acts on vectors by rotation in its plane.
        This is equivalent to cross product with the dual vector:
        
        B × v = B* × v where B* is the Hodge dual (axial vector)
        
        For B = b₂₃ e₂₃ + b₃₁ e₃₁ + b₁₂ e₁₂:
        B × v = (b₂₃, b₃₁, b₁₂) × (v¹, v², v³)  [as cross product]
        """
        b = np.array([self.e23, self.e31, self.e12])
        return np.cross(b, v)
    
    def commutator_with_bivector(self, other: 'Bivector3D') -> 'Bivector3D':
        """
        Compute B₁ × B₂ = ½(B₁B₂ - B₂B₁).
        
        For bivectors in 3D, this gives another bivector (Lie bracket).
        """
        # In 3D, bivector commutator is like cross product of dual vectors
        b1 = self.to_array()
        b2 = other.to_array()
        result = np.cross(b1, b2)
        return Bivector3D.from_array(result)
    
    def inner_with_bivector(self, other: 'Bivector3D') -> float:
        """Compute B₁ · B₂ (scalar inner product)."""
        return self.e23 * other.e23 + self.e31 * other.e31 + self.e12 * other.e12
    
    def wedge_with_bivector(self, other: 'Bivector3D') -> float:
        """
        Compute B₁ ∧ B₂ in 3D.
        
        Two bivectors wedge to a 4-vector, but in 3D this is zero
        (no grade-4 elements). Returns the pseudoscalar coefficient.
        """
        # In 3D: B ∧ B' projects to grade 4, which doesn't exist
        # But B₁B₂ has a grade-0 part (scalar) which is -B₁·B₂
        # The wedge B₁ ∧ B₂ would be grade-4, so it's 0 in 3D
        return 0.0


# =============================================================================
# TANGENT FRAME
# =============================================================================

@dataclass
class TangentFrame:
    """
    Tangent frame at a point on a manifold.
    
    Stores both the frame vectors {eᵢ} and reciprocal frame {eⁱ}
    satisfying eⁱ · eⱼ = δⁱⱼ.
    
    The metric is implicit: gᵢⱼ = eᵢ · eⱼ
    """
    vectors: np.ndarray  # Shape (dim, ambient_dim) — frame vectors eᵢ
    
    def __post_init__(self):
        self.dim = self.vectors.shape[0]
        self.ambient_dim = self.vectors.shape[1]
        self._compute_reciprocal()
        self._compute_metric()
    
    def _compute_reciprocal(self):
        """Compute reciprocal frame eⁱ such that eⁱ · eⱼ = δⁱⱼ."""
        # Reciprocal frame: eⁱ = gⁱʲ eⱼ
        # First compute metric gᵢⱼ = eᵢ · eⱼ
        g = self.vectors @ self.vectors.T
        g_inv = np.linalg.inv(g)
        self.reciprocal = g_inv @ self.vectors
    
    def _compute_metric(self):
        """Compute metric tensor and its inverse."""
        self.g = self.vectors @ self.vectors.T  # gᵢⱼ
        self.g_inv = np.linalg.inv(self.g)       # gⁱʲ
        self.sqrt_det_g = np.sqrt(np.abs(np.linalg.det(self.g)))
    
    def metric(self, i: int, j: int) -> float:
        """Return gᵢⱼ = eᵢ · eⱼ."""
        return self.g[i, j]
    
    def inverse_metric(self, i: int, j: int) -> float:
        """Return gⁱʲ = eⁱ · eʲ."""
        return self.g_inv[i, j]
    
    def pseudoscalar_2d(self) -> Bivector3D:
        """
        Compute tangent pseudoscalar B = e₁ ∧ e₂ for 2D surface in 3D.
        
        This bivector represents the oriented tangent plane.
        """
        if self.dim != 2 or self.ambient_dim != 3:
            raise ValueError("pseudoscalar_2d requires 2D surface in 3D")
        
        return Bivector3D.from_wedge(self.vectors[0], self.vectors[1])
    
    def normal_3d(self) -> np.ndarray:
        """
        Compute unit normal for 2D surface in 3D.
        
        n = B · I⁻¹ where B is tangent bivector, I is 3D pseudoscalar.
        Equivalently: n = (e₁ × e₂) / |e₁ × e₂|
        """
        if self.dim != 2 or self.ambient_dim != 3:
            raise ValueError("normal_3d requires 2D surface in 3D")
        
        cross = np.cross(self.vectors[0], self.vectors[1])
        norm = np.linalg.norm(cross)
        if norm < 1e-10:
            raise ValueError("Degenerate frame: vectors are parallel")
        return cross / norm


# =============================================================================
# CONNECTION BIVECTOR
# =============================================================================

class ConnectionBivector:
    """
    The connection bivector ω encodes frame rotation.
    
    Replaces Christoffel symbols Γᵏᵢⱼ with a geometric object:
        ωᵢ = ½(∂ᵢeⱼ) ∧ eʲ
    
    The frame derivative is then:
        ∂ᵢeⱼ = ωᵢ × eⱼ
    
    Properties:
    - Metric compatible: ωᵢ preserves inner products (generates rotations)
    - Torsion-free: ωᵢ × eⱼ = ωⱼ × eᵢ for coordinate frames
    """
    
    def __init__(self, manifold_dim: int, ambient_dim: int = 3, h: float = 1e-6):
        """
        Args:
            manifold_dim: Dimension of the manifold
            ambient_dim: Dimension of ambient space
            h: Step size for numerical differentiation
        """
        self.dim = manifold_dim
        self.ambient_dim = ambient_dim
        self.h = h
    
    def compute_at(
        self,
        coords: np.ndarray,
        frame_func: Callable[[np.ndarray], TangentFrame]
    ) -> List[Bivector3D]:
        """
        Compute connection bivectors ωᵢ at given coordinates.
        
        Args:
            coords: Coordinate values
            frame_func: Function (coords) -> TangentFrame
        
        Returns:
            List of bivectors [ω₁, ω₂, ...] (one per coordinate direction)
        """
        frame = frame_func(coords)
        omega_list = []
        
        for i in range(self.dim):
            # Compute ∂ᵢeⱼ using central differences
            coords_plus = coords.copy()
            coords_minus = coords.copy()
            coords_plus[i] += self.h
            coords_minus[i] -= self.h
            
            frame_plus = frame_func(coords_plus)
            frame_minus = frame_func(coords_minus)
            
            # ∂ᵢeⱼ for each j
            d_frame = (frame_plus.vectors - frame_minus.vectors) / (2 * self.h)
            
            # ωᵢ = ½ Σⱼ (∂ᵢeⱼ) ∧ eʲ
            omega_i = Bivector3D()
            for j in range(self.dim):
                d_ej = d_frame[j]  # ∂ᵢeⱼ
                e_j_recip = frame.reciprocal[j]  # eʲ
                
                # Add eʲ ∧ (∂ᵢeⱼ) to omega (correct sign: reciprocal frame first)
                wedge_jk = Bivector3D.from_wedge(e_j_recip, d_ej)
                omega_i = omega_i + wedge_jk * 0.5
            
            omega_list.append(omega_i)
        
        return omega_list
    
    def frame_derivative(
        self,
        omega_i: Bivector3D,
        frame: TangentFrame,
        j: int
    ) -> np.ndarray:
        """
        Compute ∂ᵢeⱼ = ωᵢ × eⱼ.
        
        This is the fundamental equation relating connection to frame variation.
        """
        ej = frame.vectors[j]
        return omega_i.commutator_with_vector(ej)


# =============================================================================
# CURVATURE 2-FORM
# =============================================================================

class Curvature2Form:
    """
    The curvature 2-form Ω via Cartan's structure equation:
        Ω = dω + ω ∧ ω
    
    This replaces the Riemann tensor Rˡᵢⱼₖ.
    
    Key properties:
    - Ωᵢⱼ = ∂ᵢωⱼ - ∂ⱼωᵢ + ωᵢ × ωⱼ
    - [∇ᵢ, ∇ⱼ]v = Ωᵢⱼ × v (Riemann action on vectors)
    - Holonomy around infinitesimal loop ≈ Ωᵢⱼ × area
    """
    
    def __init__(self, connection: ConnectionBivector, h: float = 1e-6):
        self.connection = connection
        self.dim = connection.dim
        self.h = h
    
    def compute_at(
        self,
        coords: np.ndarray,
        frame_func: Callable[[np.ndarray], TangentFrame]
    ) -> np.ndarray:
        """
        Compute curvature 2-form components Ωᵢⱼ.
        
        Returns:
            Array of shape (dim, dim) containing Bivector3D objects
            where result[i,j] = Ωᵢⱼ
        """
        # Get connection at this point and nearby points
        omega = self.connection.compute_at(coords, frame_func)
        
        # Compute partial derivatives of connection
        d_omega = []  # d_omega[i][j] = ∂ᵢωⱼ
        for i in range(self.dim):
            coords_plus = coords.copy()
            coords_minus = coords.copy()
            coords_plus[i] += self.h
            coords_minus[i] -= self.h
            
            omega_plus = self.connection.compute_at(coords_plus, frame_func)
            omega_minus = self.connection.compute_at(coords_minus, frame_func)
            
            d_omega_i = []
            for j in range(self.dim):
                d_ij = (omega_plus[j] - omega_minus[j]) * (1.0 / (2 * self.h))
                d_omega_i.append(d_ij)
            d_omega.append(d_omega_i)
        
        # Compute Ωᵢⱼ = ∂ᵢωⱼ - ∂ⱼωᵢ + ωᵢ × ωⱼ
        Omega = np.empty((self.dim, self.dim), dtype=object)
        for i in range(self.dim):
            for j in range(self.dim):
                # dω part: ∂ᵢωⱼ - ∂ⱼωᵢ
                d_part = d_omega[i][j] - d_omega[j][i]
                
                # ω ∧ ω part: ωᵢ × ωⱼ
                wedge_part = omega[i].commutator_with_bivector(omega[j])
                
                Omega[i, j] = d_part + wedge_part
        
        return Omega
    
    def sectional_curvature(
        self,
        Omega_ij: Bivector3D,
        frame: TangentFrame,
        i: int,
        j: int
    ) -> float:
        """
        Compute sectional curvature K(eᵢ ∧ eⱼ).
        
        K = Ω(eᵢ, eⱼ) · (eᵢ ∧ eⱼ) / |eᵢ ∧ eⱼ|²
        """
        # Tangent bivector for the plane
        ei, ej = frame.vectors[i], frame.vectors[j]
        plane = Bivector3D.from_wedge(ei, ej)
        
        # Inner product
        numerator = Omega_ij.inner_with_bivector(plane)
        denominator = plane.norm_squared
        
        if denominator < 1e-20:
            return 0.0
        
        return numerator / denominator


# =============================================================================
# COORDINATE-FREE COVARIANT DERIVATIVE
# =============================================================================

class GACovariantDerivative:
    """
    Covariant derivative using connection bivector:
        ∇ᵤv = ∂ᵤv + ω(u) × v
    
    This is the coordinate-free formulation — no Christoffel symbols.
    """
    
    def __init__(self, connection: ConnectionBivector):
        self.connection = connection
        self.h = connection.h
    
    def of_vector(
        self,
        v_func: Callable[[np.ndarray], np.ndarray],
        u: np.ndarray,
        coords: np.ndarray,
        frame_func: Callable[[np.ndarray], TangentFrame]
    ) -> np.ndarray:
        """
        Compute ∇ᵤv at given coordinates.
        
        Args:
            v_func: Vector field function (coords) -> v
            u: Direction vector (in coordinate basis)
            coords: Point coordinates
            frame_func: Function (coords) -> TangentFrame
        
        Returns:
            Covariant derivative ∇ᵤv (ambient space vector)
        """
        # Get connection
        omega = self.connection.compute_at(coords, frame_func)
        
        # Compute ω(u) = uⁱωᵢ
        omega_u = Bivector3D()
        for i in range(len(u)):
            omega_u = omega_u + omega[i] * u[i]
        
        # Compute partial derivative ∂ᵤv = uⁱ∂ᵢv
        partial_u_v = np.zeros(self.connection.ambient_dim)
        for i in range(len(u)):
            coords_plus = coords.copy()
            coords_minus = coords.copy()
            coords_plus[i] += self.h
            coords_minus[i] -= self.h
            
            v_plus = v_func(coords_plus)
            v_minus = v_func(coords_minus)
            
            partial_i_v = (v_plus - v_minus) / (2 * self.h)
            partial_u_v += u[i] * partial_i_v
        
        # Compute v at this point
        v = v_func(coords)
        
        # ∇ᵤv = ∂ᵤv + ω(u) × v
        return partial_u_v + omega_u.commutator_with_vector(v)
    
    def divergence(
        self,
        v_func: Callable[[np.ndarray], np.ndarray],
        coords: np.ndarray,
        frame_func: Callable[[np.ndarray], TangentFrame]
    ) -> float:
        """
        Compute divergence: div v = ∇ · v = (1/√g) ∂ᵢ(√g vⁱ)
        """
        frame = frame_func(coords)
        sqrt_g = frame.sqrt_det_g
        
        div = 0.0
        for i in range(self.connection.dim):
            coords_plus = coords.copy()
            coords_minus = coords.copy()
            coords_plus[i] += self.h
            coords_minus[i] -= self.h
            
            frame_plus = frame_func(coords_plus)
            frame_minus = frame_func(coords_minus)
            
            v_plus = v_func(coords_plus)
            v_minus = v_func(coords_minus)
            
            # vⁱ = v · eⁱ (contravariant component)
            vi_plus = np.dot(v_plus, frame_plus.reciprocal[i])
            vi_minus = np.dot(v_minus, frame_minus.reciprocal[i])
            
            sqrt_g_vi_plus = frame_plus.sqrt_det_g * vi_plus
            sqrt_g_vi_minus = frame_minus.sqrt_det_g * vi_minus
            
            div += (sqrt_g_vi_plus - sqrt_g_vi_minus) / (2 * self.h)
        
        return div / (sqrt_g + 1e-10)


# =============================================================================
# GEODESIC SOLVER (COORDINATE-FREE)
# =============================================================================

class GAGeodesicSolver:
    """
    Solve geodesic equation using connection bivector:
        v̇ + ω(v) × v = 0
    
    This is equivalent to ∇ᵥv = 0 — velocity parallel-transports itself.
    """
    
    def __init__(self, connection: ConnectionBivector):
        self.connection = connection
    
    def solve(
        self,
        x0: np.ndarray,
        v0: np.ndarray,
        frame_func: Callable[[np.ndarray], TangentFrame],
        t_final: float,
        n_steps: int = 100
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Integrate geodesic equation from initial conditions.
        
        Args:
            x0: Initial coordinates
            v0: Initial velocity (in coordinate basis)
            frame_func: Function (coords) -> TangentFrame
            t_final: Final time
            n_steps: Number of integration steps
        
        Returns:
            (t_values, x_values, v_values)
        """
        dt = t_final / n_steps
        dim = len(x0)
        
        t_values = np.linspace(0, t_final, n_steps + 1)
        x_values = np.zeros((n_steps + 1, dim))
        v_values = np.zeros((n_steps + 1, dim))
        
        x_values[0] = x0.copy()
        v_values[0] = v0.copy()
        
        x = x0.copy()
        v = v0.copy()
        
        for step in range(n_steps):
            # Get connection at current point
            omega = self.connection.compute_at(x, frame_func)
            frame = frame_func(x)
            
            # Convert coordinate velocity to ambient space
            v_ambient = np.zeros(self.connection.ambient_dim)
            for i in range(dim):
                v_ambient += v[i] * frame.vectors[i]
            
            # Compute ω(v) = vⁱωᵢ
            omega_v = Bivector3D()
            for i in range(dim):
                omega_v = omega_v + omega[i] * v[i]
            
            # Geodesic equation: v̇ = -ω(v) × v
            v_dot_ambient = -omega_v.commutator_with_vector(v_ambient)
            
            # Convert acceleration back to coordinates
            # v̇ⁱ = v̇_ambient · eⁱ
            v_dot = np.zeros(dim)
            for i in range(dim):
                v_dot[i] = np.dot(v_dot_ambient, frame.reciprocal[i])
            
            # RK4 integration would be better, but simple Euler for now
            # (The velocity update should really use the geodesic equation in coords)
            # For stability, use symplectic Euler or RK4
            
            # Simple approach: use classical geodesic equation
            # ẍᵏ + Γᵏᵢⱼ ẋⁱ ẋʲ = 0
            # But we want to stay coordinate-free...
            
            # Approximate: update position, then correct velocity to stay on geodesic
            x_new = x + v * dt
            
            # Re-project velocity to tangent space at new point
            # (This is a simplification — proper integration needs more care)
            v_new = v + v_dot * dt
            
            x = x_new
            v = v_new
            
            x_values[step + 1] = x
            v_values[step + 1] = v
        
        return t_values, x_values, v_values


# =============================================================================
# PARALLEL TRANSPORT
# =============================================================================

class GAParallelTransport:
    """
    Parallel transport using connection bivector:
        ẇ + ω(v) × w = 0
    
    Transport vector w along curve with velocity v, keeping ∇ᵥw = 0.
    """
    
    def __init__(self, connection: ConnectionBivector):
        self.connection = connection
    
    def transport(
        self,
        curve_func: Callable[[float], np.ndarray],
        w0: np.ndarray,
        frame_func: Callable[[np.ndarray], TangentFrame],
        t_start: float,
        t_end: float,
        n_steps: int = 100
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Transport vector w₀ along curve.
        
        Args:
            curve_func: Curve parameterization (t) -> coords
            w0: Initial vector (in coordinate basis)
            frame_func: Function (coords) -> TangentFrame
            t_start, t_end: Parameter range
            n_steps: Integration steps
        
        Returns:
            (t_values, w_values) — transported vector at each step
        """
        dt = (t_end - t_start) / n_steps
        dim = len(w0)
        
        t_values = np.linspace(t_start, t_end, n_steps + 1)
        w_values = np.zeros((n_steps + 1, dim))
        w_values[0] = w0.copy()
        
        w = w0.copy()
        
        for step in range(n_steps):
            t = t_start + step * dt
            coords = curve_func(t)
            
            # Curve velocity (numerical)
            coords_plus = curve_func(t + dt/10)
            coords_minus = curve_func(t - dt/10)
            v = (coords_plus - coords_minus) / (dt/5)
            
            # Get connection and frame
            omega = self.connection.compute_at(coords, frame_func)
            frame = frame_func(coords)
            
            # Convert w to ambient space
            w_ambient = np.zeros(self.connection.ambient_dim)
            for i in range(dim):
                w_ambient += w[i] * frame.vectors[i]
            
            # Compute ω(v) = vⁱωᵢ
            omega_v = Bivector3D()
            for i in range(dim):
                omega_v = omega_v + omega[i] * v[i]
            
            # Transport equation: ẇ = -ω(v) × w
            w_dot_ambient = -omega_v.commutator_with_vector(w_ambient)
            
            # Convert back to coordinates
            w_dot = np.zeros(dim)
            for i in range(dim):
                w_dot[i] = np.dot(w_dot_ambient, frame.reciprocal[i])
            
            # Update
            w = w + w_dot * dt
            w_values[step + 1] = w
        
        return t_values, w_values
    
    def holonomy(
        self,
        loop_func: Callable[[float], np.ndarray],
        w0: np.ndarray,
        frame_func: Callable[[np.ndarray], TangentFrame],
        n_steps: int = 200
    ) -> Tuple[np.ndarray, float]:
        """
        Compute holonomy around a closed loop.
        
        Returns:
            (w_final, angle) — transported vector and rotation angle
        """
        t_vals, w_vals = self.transport(
            loop_func, w0, frame_func,
            t_start=0.0, t_end=2*np.pi, n_steps=n_steps
        )
        
        w_final = w_vals[-1]
        
        # Compute rotation angle
        # cos(θ) = (w₀ · w_final) / (|w₀| |w_final|)
        dot = np.dot(w0, w_final)
        norm0 = np.linalg.norm(w0)
        norm_final = np.linalg.norm(w_final)
        
        if norm0 < 1e-10 or norm_final < 1e-10:
            return w_final, 0.0
        
        cos_angle = dot / (norm0 * norm_final)
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        angle = np.arccos(cos_angle)
        
        return w_final, angle


# =============================================================================
# SHAPE OPERATOR
# =============================================================================

class ShapeOperator:
    """
    Shape operator (Weingarten map) via pseudoscalar variation:
        S(v) = -∇ᵥn = (v · ∇B) · I⁻¹
    
    Where B is the tangent pseudoscalar and n is the unit normal.
    
    Principal curvatures are eigenvalues of S.
    Gaussian curvature K = det(S) = κ₁κ₂
    Mean curvature H = ½tr(S) = ½(κ₁ + κ₂)
    """
    
    def __init__(self, h: float = 1e-6):
        self.h = h
    
    def compute_at(
        self,
        coords: np.ndarray,
        frame_func: Callable[[np.ndarray], TangentFrame]
    ) -> np.ndarray:
        """
        Compute shape operator matrix Sⁱⱼ at given coordinates.
        
        For 2D surface in 3D:
            S(eⱼ) = Sⁱⱼ eᵢ
        
        Returns:
            (2, 2) array representing S in the frame basis
        """
        frame = frame_func(coords)
        
        if frame.dim != 2 or frame.ambient_dim != 3:
            raise ValueError("Shape operator requires 2D surface in 3D")
        
        # Get normal at this point
        n = frame.normal_3d()
        
        # Compute ∂ᵢn using central differences
        S_matrix = np.zeros((2, 2))
        
        for j in range(2):
            coords_plus = coords.copy()
            coords_minus = coords.copy()
            coords_plus[j] += self.h
            coords_minus[j] -= self.h
            
            frame_plus = frame_func(coords_plus)
            frame_minus = frame_func(coords_minus)
            
            n_plus = frame_plus.normal_3d()
            n_minus = frame_minus.normal_3d()
            
            # ∂ⱼn
            dn_j = (n_plus - n_minus) / (2 * self.h)
            
            # S(eⱼ) = -∂ⱼn projected to tangent space
            # Sⁱⱼ = -∂ⱼn · eⁱ
            for i in range(2):
                S_matrix[i, j] = -np.dot(dn_j, frame.reciprocal[i])
        
        return S_matrix
    
    def gaussian_curvature(self, S: np.ndarray) -> float:
        """K = det(S)"""
        return np.linalg.det(S)
    
    def mean_curvature(self, S: np.ndarray) -> float:
        """H = ½ tr(S)"""
        return 0.5 * np.trace(S)
    
    def principal_curvatures(self, S: np.ndarray) -> Tuple[float, float]:
        """Eigenvalues of S."""
        eigenvalues = np.linalg.eigvals(S)
        eigenvalues = np.real(eigenvalues)
        return tuple(sorted(eigenvalues))
    
    def principal_directions(self, S: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Eigenvectors of S (in coordinate basis)."""
        eigenvalues, eigenvectors = np.linalg.eig(S)
        idx = np.argsort(np.real(eigenvalues))
        return eigenvectors[:, idx[0]], eigenvectors[:, idx[1]]


# =============================================================================
# EXAMPLE MANIFOLDS
# =============================================================================

def sphere_frame(R: float = 1.0) -> Callable[[np.ndarray], TangentFrame]:
    """
    Create frame function for sphere of radius R.
    
    Coordinates: (θ, φ) where θ ∈ [0,π], φ ∈ [0,2π)
    
    Frame vectors:
        e_θ = R (cos θ cos φ, cos θ sin φ, -sin θ)
        e_φ = R sin θ (-sin φ, cos φ, 0)
    """
    def frame_func(coords: np.ndarray) -> TangentFrame:
        theta, phi = coords[0], coords[1]
        
        # Avoid singularity at poles
        theta = np.clip(theta, 0.01, np.pi - 0.01)
        
        ct, st = np.cos(theta), np.sin(theta)
        cp, sp = np.cos(phi), np.sin(phi)
        
        # Tangent vectors (derivatives of embedding)
        e_theta = R * np.array([ct * cp, ct * sp, -st])
        e_phi = R * st * np.array([-sp, cp, 0.0])
        
        vectors = np.array([e_theta, e_phi])
        return TangentFrame(vectors)
    
    return frame_func


def torus_frame(R: float = 2.0, r: float = 0.5) -> Callable[[np.ndarray], TangentFrame]:
    """
    Create frame function for torus with major radius R, minor radius r.
    
    Coordinates: (θ, φ) where θ is around the tube, φ is around the hole
    
    Embedding: 
        x = (R + r cos θ) cos φ
        y = (R + r cos θ) sin φ
        z = r sin θ
    """
    def frame_func(coords: np.ndarray) -> TangentFrame:
        theta, phi = coords[0], coords[1]
        
        ct, st = np.cos(theta), np.sin(theta)
        cp, sp = np.cos(phi), np.sin(phi)
        
        # Tangent vectors
        e_theta = r * np.array([-st * cp, -st * sp, ct])
        e_phi = (R + r * ct) * np.array([-sp, cp, 0.0])
        
        vectors = np.array([e_theta, e_phi])
        return TangentFrame(vectors)
    
    return frame_func


def paraboloid_frame(a: float = 1.0) -> Callable[[np.ndarray], TangentFrame]:
    """
    Create frame function for paraboloid z = a(x² + y²).
    
    Coordinates: (r, θ) polar
    
    Embedding:
        x = r cos θ
        y = r sin θ
        z = a r²
    """
    def frame_func(coords: np.ndarray) -> TangentFrame:
        r, theta = coords[0], coords[1]
        r = max(r, 0.01)  # Avoid singularity at origin
        
        ct, st = np.cos(theta), np.sin(theta)
        
        # Tangent vectors
        e_r = np.array([ct, st, 2 * a * r])
        e_theta = r * np.array([-st, ct, 0.0])
        
        vectors = np.array([e_r, e_theta])
        return TangentFrame(vectors)
    
    return frame_func


# =============================================================================
# TESTS
# =============================================================================

def test_sphere_curvature():
    """Test Gaussian curvature on sphere using coordinate-free method."""
    print("=" * 60)
    print("TEST: Sphere Curvature via Connection Bivector")
    print("=" * 60)
    
    R = 1.0
    frame_func = sphere_frame(R)
    
    # Create coordinate-free objects
    conn = ConnectionBivector(manifold_dim=2, ambient_dim=3)
    curv = Curvature2Form(conn)
    shape = ShapeOperator()
    
    # Test at θ = π/3 (60°)
    coords = np.array([np.pi/3, 0.0])
    
    # Method 1: Via shape operator
    S = shape.compute_at(coords, frame_func)
    K_shape = shape.gaussian_curvature(S)
    H_shape = shape.mean_curvature(S)
    
    print(f"\nMethod 1: Shape Operator")
    print(f"  S = \n{S}")
    print(f"  K = {K_shape:.6f} (expected: {1/R**2:.6f})")
    print(f"  H = {H_shape:.6f} (expected: {1/R:.6f})")
    
    # Method 2: Via curvature 2-form
    Omega = curv.compute_at(coords, frame_func)
    frame = frame_func(coords)
    K_curv = curv.sectional_curvature(Omega[0, 1], frame, 0, 1)
    
    print(f"\nMethod 2: Curvature 2-Form")
    print(f"  Ω₁₂ = {Omega[0, 1]}")
    print(f"  K = {K_curv:.6f} (expected: {1/R**2:.6f})")
    
    # Check connection bivector
    omega = conn.compute_at(coords, frame_func)
    print(f"\nConnection Bivectors:")
    print(f"  ω_θ = {omega[0]}")
    print(f"  ω_φ = {omega[1]}")
    
    error_K = abs(K_shape - 1/R**2)
    error_H = abs(H_shape - 1/R)
    
    print(f"\nErrors:")
    print(f"  |K - 1/R²| = {error_K:.6f}")
    print(f"  |H - 1/R| = {error_H:.6f}")
    
    print(f"\nStatus: {'PASS' if error_K < 0.1 and error_H < 0.1 else 'FAIL'}")


def test_torus_curvature():
    """Test varying curvature on torus."""
    print("\n" + "=" * 60)
    print("TEST: Torus Curvature (Variable K)")
    print("=" * 60)
    
    R, r = 2.0, 0.5
    frame_func = torus_frame(R, r)
    shape = ShapeOperator()
    
    # Test at different points
    test_points = [
        (0.0, 0.0, "Outer equator (K > 0)"),
        (np.pi, 0.0, "Inner equator (K < 0)"),
        (np.pi/2, 0.0, "Top (K = 0)"),
    ]
    
    for theta, phi, name in test_points:
        coords = np.array([theta, phi])
        S = shape.compute_at(coords, frame_func)
        K = shape.gaussian_curvature(S)
        H = shape.mean_curvature(S)
        
        # Expected: K = cos(θ) / (r(R + r cos θ))
        ct = np.cos(theta)
        K_expected = ct / (r * (R + r * ct))
        
        print(f"\n{name}:")
        print(f"  θ = {theta:.4f}, φ = {phi:.4f}")
        print(f"  K = {K:.6f} (expected: {K_expected:.6f})")
        print(f"  H = {H:.6f}")


def test_parallel_transport_holonomy():
    """Test holonomy on sphere via parallel transport."""
    print("\n" + "=" * 60)
    print("TEST: Parallel Transport Holonomy on Sphere")
    print("=" * 60)
    
    R = 1.0
    frame_func = sphere_frame(R)
    conn = ConnectionBivector(manifold_dim=2, ambient_dim=3)
    transport = GAParallelTransport(conn)
    
    # Transport around latitude circle at θ = π/4
    theta_fixed = np.pi / 4
    
    def latitude_loop(t):
        return np.array([theta_fixed, t])
    
    # Initial vector in θ direction
    w0 = np.array([1.0, 0.0])
    
    # Transport and get holonomy
    w_final, angle = transport.holonomy(latitude_loop, w0, frame_func, n_steps=200)
    
    # Expected holonomy: 2π(1 - cos θ)
    expected_angle = 2 * np.pi * (1 - np.cos(theta_fixed))
    
    print(f"\nLatitude circle at θ = π/4 ({np.degrees(theta_fixed):.1f}°)")
    print(f"Initial vector: {w0}")
    print(f"Final vector: {w_final}")
    print(f"Holonomy angle: {np.degrees(angle):.2f}°")
    print(f"Expected angle: {np.degrees(expected_angle):.2f}°")
    print(f"Error: {np.degrees(abs(angle - expected_angle)):.2f}°")
    
    # Note: Holonomy = ∫K dA over enclosed region
    # For spherical cap: A = 2πR²(1 - cos θ), K = 1/R²
    # So ∫K dA = 2π(1 - cos θ) ✓


def run_all_tests():
    """Run all coordinate-free Riemannian geometry tests."""
    test_sphere_curvature()
    test_torus_curvature()
    test_parallel_transport_holonomy()
    
    print("\n" + "=" * 60)
    print("ALL TESTS COMPLETED")
    print("=" * 60)


if __name__ == "__main__":
    run_all_tests()
