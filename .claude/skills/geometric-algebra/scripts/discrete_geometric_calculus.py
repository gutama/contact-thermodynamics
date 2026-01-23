"""
Discrete Geometric Calculus with Split Differential Operator

This module implements the geometric derivative on discrete grids,
bridging the concepts from:
1. GAlgebra's split differential operator (symbolic)
2. Eelco Hoogendoorn's discrete geometric calculus (numerical)

Author: Claude (Anthropic)
Date: 2025-01-18
"""

import numpy as np
from typing import Tuple, List, Dict, Optional, Callable, Union
from dataclasses import dataclass
from itertools import combinations
from functools import lru_cache


@dataclass
class AlgebraSignature:
    """Signature (p,q,r) for Clifford algebra Cl(p,q,r)."""
    p: int  # Positive dimensions
    q: int  # Negative dimensions  
    r: int  # Null dimensions
    
    @property
    def dim(self) -> int:
        return self.p + self.q + self.r
    
    @property
    def metric_diagonal(self) -> List[float]:
        return [1.0] * self.p + [-1.0] * self.q + [0.0] * self.r
    
    def __str__(self):
        return f"Cl({self.p},{self.q},{self.r})"


# =============================================================================
# COMMON ALGEBRA SIGNATURES
# =============================================================================

EUCLIDEAN_2D = AlgebraSignature(2, 0, 0)      # Cl(2,0,0) - 2D Euclidean
EUCLIDEAN_3D = AlgebraSignature(3, 0, 0)      # Cl(3,0,0) - 3D Euclidean
PGA_2D = AlgebraSignature(2, 0, 1)            # Cl(2,0,1) - 2D Projective
PGA_3D = AlgebraSignature(3, 0, 1)            # Cl(3,0,1) - 3D Projective
CGA_2D = AlgebraSignature(3, 1, 0)            # Cl(3,1,0) - 2D Conformal
CGA_3D = AlgebraSignature(4, 1, 0)            # Cl(4,1,0) - 3D Conformal
MINKOWSKI_2D = AlgebraSignature(1, 1, 0)      # Cl(1,1,0) - 1+1D Spacetime
MINKOWSKI_4D = AlgebraSignature(1, 3, 0)      # Cl(1,3,0) - 1+3D Spacetime
SPACETIME_ALGEBRA = AlgebraSignature(1, 3, 0) # Dirac algebra


class CliffordBasis:
    """
    Manages the basis blades of a Clifford algebra.
    
    Basis blades are indexed by tuples, e.g.:
    - () : scalar (grade 0)
    - (0,) : e₀ (grade 1)
    - (0,1) : e₀₁ = e₀∧e₁ (grade 2)
    - (0,1,2) : e₀₁₂ (grade 3)
    """
    
    def __init__(self, signature: AlgebraSignature):
        self.signature = signature
        self.dim = signature.dim
        self.metric = np.diag(signature.metric_diagonal)
        
        # Build all basis blades
        self.blades: List[Tuple[int, ...]] = self._build_blades()
        self.n_blades = len(self.blades)  # = 2^dim
        
        # Build lookup table
        self._blade_to_idx: Dict[Tuple[int, ...], int] = {
            b: i for i, b in enumerate(self.blades)
        }
        
        # Precompute geometric product table for basis vectors
        self._product_cache: Dict[Tuple[int, Tuple[int, ...]], Tuple[float, Tuple[int, ...]]] = {}
        
    def _build_blades(self) -> List[Tuple[int, ...]]:
        """Build all 2^n basis blades in canonical order."""
        blades = [()]
        for grade in range(1, self.dim + 1):
            blades.extend(combinations(range(self.dim), grade))
        return blades
    
    def blade_to_index(self, blade: Tuple[int, ...]) -> int:
        """Convert blade tuple to linear index."""
        return self._blade_to_idx[blade]
    
    def index_to_blade(self, idx: int) -> Tuple[int, ...]:
        """Convert linear index to blade tuple."""
        return self.blades[idx]
    
    def grade(self, blade_or_idx: Union[Tuple[int, ...], int]) -> int:
        """Get grade of a blade."""
        if isinstance(blade_or_idx, int):
            return len(self.blades[blade_or_idx])
        return len(blade_or_idx)
    
    def basis_vector_product(self, i: int, blade: Tuple[int, ...]) -> Tuple[float, Tuple[int, ...]]:
        """
        Compute eᵢ * blade (geometric product of basis vector with blade).
        
        This is the core operation for the ∇_G part of the split operator.
        
        Returns: (coefficient, result_blade)
        """
        cache_key = (i, blade)
        if cache_key in self._product_cache:
            return self._product_cache[cache_key]
        
        if i in blade:
            # Contraction case: eᵢ is already in blade
            # eᵢ * (... ∧ eᵢ ∧ ...) = metric[i,i] * (... ∧ ...)
            pos = blade.index(i)
            sign = (-1) ** pos * self.metric[i, i]
            result = tuple(b for b in blade if b != i)
        else:
            # Extension case: add eᵢ to blade
            # eᵢ * (... ∧ ...) = eᵢ ∧ (... ∧ ...)
            result = list(blade)
            pos = 0
            for j, b in enumerate(blade):
                if b > i:
                    pos = j
                    break
                pos = j + 1
            sign = (-1) ** pos
            result.insert(pos, i)
            result = tuple(result)
        
        self._product_cache[cache_key] = (sign, result)
        return (sign, result)
    
    def reciprocal_basis_product(self, i: int, blade: Tuple[int, ...]) -> Tuple[float, Tuple[int, ...]]:
        """
        Compute eⁱ * blade where eⁱ is reciprocal basis vector.
        
        For orthonormal basis: eⁱ = metric[i,i] * eᵢ
        """
        sign, result = self.basis_vector_product(i, blade)
        # For orthonormal basis, reciprocal vector is metric-weighted
        metric_factor = self.metric[i, i]
        if metric_factor == 0:
            return (0.0, ())
        return (sign * metric_factor, result)


class DiscreteMultivectorField:
    """
    A multivector-valued field on a discrete grid.
    
    Shape: (*grid_shape, n_blades)
    """
    
    def __init__(self, basis: CliffordBasis, shape: Tuple[int, ...], 
                 data: Optional[np.ndarray] = None):
        self.basis = basis
        self.shape = shape  # Grid shape
        
        if data is not None:
            assert data.shape == (*shape, basis.n_blades)
            self.data = data
        else:
            self.data = np.zeros((*shape, basis.n_blades))
    
    def __add__(self, other: 'DiscreteMultivectorField') -> 'DiscreteMultivectorField':
        return DiscreteMultivectorField(self.basis, self.shape, self.data + other.data)
    
    def __sub__(self, other: 'DiscreteMultivectorField') -> 'DiscreteMultivectorField':
        return DiscreteMultivectorField(self.basis, self.shape, self.data - other.data)
    
    def __mul__(self, scalar: float) -> 'DiscreteMultivectorField':
        return DiscreteMultivectorField(self.basis, self.shape, self.data * scalar)
    
    def __rmul__(self, scalar: float) -> 'DiscreteMultivectorField':
        return self.__mul__(scalar)
    
    def grade_select(self, k: int) -> 'DiscreteMultivectorField':
        """Extract grade-k part of the field."""
        result = DiscreteMultivectorField(self.basis, self.shape)
        for idx, blade in enumerate(self.basis.blades):
            if len(blade) == k:
                result.data[..., idx] = self.data[..., idx]
        return result
    
    def scalar_part(self) -> np.ndarray:
        """Extract scalar (grade-0) component."""
        return self.data[..., 0]
    
    def vector_part(self) -> np.ndarray:
        """Extract vector (grade-1) components as array."""
        indices = [self.basis.blade_to_index((i,)) for i in range(self.basis.dim)]
        return self.data[..., indices]
    
    def norm_squared(self) -> np.ndarray:
        """Compute |F|² = F * ~F for each point."""
        # For simple cases, approximate with sum of squares
        return np.sum(self.data ** 2, axis=-1)
    
    def copy(self) -> 'DiscreteMultivectorField':
        return DiscreteMultivectorField(self.basis, self.shape, self.data.copy())


class SplitDifferentialOperator:
    """
    The Split Differential Operator for discrete geometric calculus.
    
    Implements: ∇ = ∇_G ⊗ ∇_D
    
    Where:
    - ∇_D: Partial derivatives (finite differences)
    - ∇_G: Geometric product with reciprocal basis vectors
    
    The full geometric derivative is:
        ∇F = Σᵢ eⁱ (∂F/∂xᵢ)
    """
    
    def __init__(self, basis: CliffordBasis, spacing: Union[float, Tuple[float, ...]] = 1.0):
        self.basis = basis
        
        if isinstance(spacing, (int, float)):
            self.spacing = tuple([float(spacing)] * basis.dim)
        else:
            assert len(spacing) == basis.dim
            self.spacing = tuple(spacing)
    
    # =========================================================================
    # ∇_D: DIFFERENTIAL OPERATORS (Finite Differences)
    # =========================================================================
    
    def partial(self, field: DiscreteMultivectorField, direction: int,
                stencil: str = 'central') -> DiscreteMultivectorField:
        """
        Compute partial derivative ∂F/∂xᵢ using finite differences.
        
        Args:
            field: Input multivector field
            direction: Direction index (0, 1, ..., dim-1)
            stencil: 'forward', 'backward', or 'central'
        
        Returns:
            New field containing partial derivative
        """
        h = self.spacing[direction]
        data = field.data
        ndim = len(field.shape)
        result_data = np.zeros_like(data)
        
        def slice_at(dim, start, stop):
            """Helper to create slice tuple."""
            s = [slice(None)] * (ndim + 1)
            s[dim] = slice(start, stop)
            return tuple(s)
        
        if stencil == 'forward':
            # (f[i+1] - f[i]) / h
            result_data[slice_at(direction, 0, -1)] = (
                data[slice_at(direction, 1, None)] - 
                data[slice_at(direction, 0, -1)]
            ) / h
            
        elif stencil == 'backward':
            # (f[i] - f[i-1]) / h
            result_data[slice_at(direction, 1, None)] = (
                data[slice_at(direction, 1, None)] -
                data[slice_at(direction, 0, -1)]
            ) / h
            
        else:  # central (default, second-order accurate)
            # (f[i+1] - f[i-1]) / (2h)
            result_data[slice_at(direction, 1, -1)] = (
                data[slice_at(direction, 2, None)] -
                data[slice_at(direction, 0, -2)]
            ) / (2 * h)
        
        return DiscreteMultivectorField(self.basis, field.shape, result_data)
    
    def partial_4th_order(self, field: DiscreteMultivectorField, 
                          direction: int) -> DiscreteMultivectorField:
        """
        Fourth-order accurate central difference.
        
        (-f[i+2] + 8f[i+1] - 8f[i-1] + f[i-2]) / (12h)
        """
        h = self.spacing[direction]
        data = field.data
        ndim = len(field.shape)
        result_data = np.zeros_like(data)
        
        def slice_at(dim, start, stop):
            s = [slice(None)] * (ndim + 1)
            s[dim] = slice(start, stop)
            return tuple(s)
        
        result_data[slice_at(direction, 2, -2)] = (
            -data[slice_at(direction, 4, None)] +
            8 * data[slice_at(direction, 3, -1)] -
            8 * data[slice_at(direction, 1, -3)] +
            data[slice_at(direction, 0, -4)]
        ) / (12 * h)
        
        return DiscreteMultivectorField(self.basis, field.shape, result_data)
    
    # =========================================================================
    # ∇_G: GEOMETRIC OPERATORS (Basis Vector Products)
    # =========================================================================
    
    def left_multiply_basis(self, field: DiscreteMultivectorField, 
                           direction: int) -> DiscreteMultivectorField:
        """
        Left-multiply field by basis vector eᵢ.
        
        Computes: eᵢ * F
        """
        result_data = np.zeros_like(field.data)
        
        for idx, blade in enumerate(self.basis.blades):
            sign, new_blade = self.basis.basis_vector_product(direction, blade)
            if sign != 0:
                new_idx = self.basis.blade_to_index(new_blade)
                result_data[..., new_idx] += sign * field.data[..., idx]
        
        return DiscreteMultivectorField(self.basis, field.shape, result_data)
    
    def left_multiply_reciprocal(self, field: DiscreteMultivectorField,
                                 direction: int) -> DiscreteMultivectorField:
        """
        Left-multiply field by reciprocal basis vector eⁱ.
        
        Computes: eⁱ * F
        
        For orthonormal bases: eⁱ = metric[i,i] * eᵢ
        """
        result_data = np.zeros_like(field.data)
        
        for idx, blade in enumerate(self.basis.blades):
            sign, new_blade = self.basis.reciprocal_basis_product(direction, blade)
            if sign != 0:
                new_idx = self.basis.blade_to_index(new_blade)
                result_data[..., new_idx] += sign * field.data[..., idx]
        
        return DiscreteMultivectorField(self.basis, field.shape, result_data)
    
    # =========================================================================
    # FULL GEOMETRIC DERIVATIVE: ∇ = ∇_G ⊗ ∇_D
    # =========================================================================
    
    def grad(self, field: DiscreteMultivectorField, 
             stencil: str = 'central',
             order: int = 2) -> DiscreteMultivectorField:
        """
        Compute the geometric derivative (left gradient).
        
        ∇F = Σᵢ eⁱ (∂F/∂xᵢ)
        
        This is the split operator in action:
        1. ∇_D: Compute partial derivatives via finite differences
        2. ∇_G: Multiply by reciprocal basis vectors
        3. Sum over all directions
        
        Args:
            field: Input multivector field
            stencil: Differencing scheme ('forward', 'backward', 'central')
            order: Accuracy order (2 or 4)
        
        Returns:
            Geometric derivative of field
        """
        result = DiscreteMultivectorField(self.basis, field.shape)
        
        for i in range(self.basis.dim):
            # ∇_D part: ∂F/∂xᵢ
            if order == 4:
                partial_i = self.partial_4th_order(field, i)
            else:
                partial_i = self.partial(field, i, stencil)
            
            # ∇_G part: eⁱ * (∂F/∂xᵢ)
            ei_times_partial = self.left_multiply_reciprocal(partial_i, i)
            
            # Accumulate
            result.data += ei_times_partial.data
        
        return result
    
    def rgrad(self, field: DiscreteMultivectorField,
              stencil: str = 'central') -> DiscreteMultivectorField:
        """
        Compute the right gradient: F∇̄
        
        F∇̄ = Σᵢ (∂F/∂xᵢ) * eⁱ
        """
        # For now, implement as simple sum (correct for scalar/vector fields)
        # Full implementation would need right multiplication
        return self.grad(field, stencil)  # Simplified
    
    def laplacian(self, field: DiscreteMultivectorField,
                  stencil: str = 'central') -> DiscreteMultivectorField:
        """
        Compute the Laplacian ∇²F = ∇·(∇F).
        
        For scalar fields, this reduces to: Σᵢ ∂²F/∂xᵢ²
        """
        result = DiscreteMultivectorField(self.basis, field.shape)
        
        for i in range(self.basis.dim):
            # Second derivative: ∂²F/∂xᵢ²
            first = self.partial(field, i, stencil)
            second = self.partial(first, i, stencil)
            
            # Weight by metric
            metric = self.basis.metric[i, i]
            result.data += metric * second.data
        
        return result
    
    # =========================================================================
    # LEAPFROG TIME STEPPING (Eelco's approach)
    # =========================================================================
    
    def leapfrog_dirac(self, phi_prev: DiscreteMultivectorField,
                       phi_curr: DiscreteMultivectorField,
                       dt: float,
                       time_dim: int = 0,
                       mass: float = 0.0) -> DiscreteMultivectorField:
        """
        Leapfrog step for the Dirac equation: ∂φ + mφ = 0
        
        Splits spacetime derivative into temporal and spatial parts:
            ∂φ = ∂ₜφ + ∂ₛφ = 0
        
        Solving for time derivative:
            ∂ₜφ = -∂ₛφ - mφ
        
        Leapfrog update:
            φ(t+dt) = φ(t-dt) + 2dt·∂ₜφ
                    = φ(t-dt) - 2dt·(∂ₛφ + mφ)
        
        Args:
            phi_prev: Field at t - dt
            phi_curr: Field at t  
            dt: Time step
            time_dim: Which dimension is time (default 0)
            mass: Mass term (default 0 for massless)
        
        Returns:
            Field at t + dt
        """
        # Compute spatial geometric derivative (all dims except time)
        spatial_grad = DiscreteMultivectorField(self.basis, phi_curr.shape)
        
        for i in range(self.basis.dim):
            if i != time_dim:
                partial_i = self.partial(phi_curr, i, 'central')
                ei_times = self.left_multiply_reciprocal(partial_i, i)
                spatial_grad.data += ei_times.data
        
        # Time metric factor
        time_metric = self.basis.metric[time_dim, time_dim]
        
        # Leapfrog update: φ_next = φ_prev - 2dt * metric * (∂ₛφ + mφ)
        rhs = spatial_grad.data + mass * phi_curr.data
        phi_next_data = phi_prev.data - 2 * dt * time_metric * rhs
        
        return DiscreteMultivectorField(self.basis, phi_curr.shape, phi_next_data)
    
    def leapfrog_wave(self, u_prev: np.ndarray, u_curr: np.ndarray,
                      dt: float, c: float = 1.0) -> np.ndarray:
        """
        Leapfrog for scalar wave equation: ∂²u/∂t² = c²∇²u
        
        This is a special case that doesn't need full multivector machinery.
        
        Args:
            u_prev: Scalar field at t - dt
            u_curr: Scalar field at t
            dt: Time step
            c: Wave speed
        
        Returns:
            Scalar field at t + dt
        """
        # Compute Laplacian
        field = DiscreteMultivectorField(self.basis, u_curr.shape)
        field.data[..., 0] = u_curr  # Put in scalar component
        
        laplacian = self.laplacian(field).scalar_part()
        
        # u(t+dt) = 2u(t) - u(t-dt) + c²dt²∇²u
        return 2 * u_curr - u_prev + c**2 * dt**2 * laplacian


# =============================================================================
# FIELD CONSTRUCTORS
# =============================================================================

def scalar_field(basis: CliffordBasis, shape: Tuple[int, ...], 
                 values: np.ndarray) -> DiscreteMultivectorField:
    """Create a multivector field from scalar values."""
    field = DiscreteMultivectorField(basis, shape)
    field.data[..., 0] = values
    return field


def vector_field(basis: CliffordBasis, shape: Tuple[int, ...],
                 components: List[np.ndarray]) -> DiscreteMultivectorField:
    """Create a multivector field from vector components."""
    field = DiscreteMultivectorField(basis, shape)
    for i, comp in enumerate(components):
        idx = basis.blade_to_index((i,))
        field.data[..., idx] = comp
    return field


# =============================================================================
# EXAMPLE AND TESTING
# =============================================================================

def test_gradient_2d():
    """Test gradient computation on 2D scalar field."""
    print("=" * 60)
    print("TEST: 2D Gradient of f(x,y) = sin(x)cos(y)")
    print("=" * 60)
    
    # Setup
    basis = CliffordBasis(EUCLIDEAN_2D)
    shape = (64, 64)
    h = 0.1
    nabla = SplitDifferentialOperator(basis, spacing=h)
    
    # Create grid
    x = np.linspace(0, 2*np.pi, shape[0])
    y = np.linspace(0, 2*np.pi, shape[1])
    X, Y = np.meshgrid(x, y, indexing='ij')
    
    # Test function: f = sin(x)cos(y)
    f = np.sin(X) * np.cos(Y)
    field = scalar_field(basis, shape, f)
    
    # Compute gradient
    grad_f = nabla.grad(field, stencil='central')
    
    # Extract components
    grad_x = grad_f.data[..., basis.blade_to_index((0,))]
    grad_y = grad_f.data[..., basis.blade_to_index((1,))]
    
    # Analytical: ∇f = cos(x)cos(y)·e₁ - sin(x)sin(y)·e₂
    analytical_x = np.cos(X) * np.cos(Y)
    analytical_y = -np.sin(X) * np.sin(Y)
    
    # Error (interior points)
    err_x = np.abs(grad_x[5:-5, 5:-5] - analytical_x[5:-5, 5:-5]).max()
    err_y = np.abs(grad_y[5:-5, 5:-5] - analytical_y[5:-5, 5:-5]).max()
    
    print(f"Grid: {shape[0]}×{shape[1]}, spacing h={h}")
    print(f"Max error in ∂f/∂x: {err_x:.2e}")
    print(f"Max error in ∂f/∂y: {err_y:.2e}")
    print(f"Expected O(h²) ≈ {h**2:.2e}")
    print(f"Status: {'PASS' if max(err_x, err_y) < 10*h**2 else 'FAIL'}")
    
    return err_x, err_y


def test_laplacian_2d():
    """Test Laplacian on 2D scalar field."""
    print("\n" + "=" * 60)
    print("TEST: 2D Laplacian of f(x,y) = sin(x)cos(y)")
    print("=" * 60)
    
    basis = CliffordBasis(EUCLIDEAN_2D)
    shape = (64, 64)
    h = 0.1
    nabla = SplitDifferentialOperator(basis, spacing=h)
    
    x = np.linspace(0, 2*np.pi, shape[0])
    y = np.linspace(0, 2*np.pi, shape[1])
    X, Y = np.meshgrid(x, y, indexing='ij')
    
    f = np.sin(X) * np.cos(Y)
    field = scalar_field(basis, shape, f)
    
    # Compute Laplacian
    lap_f = nabla.laplacian(field)
    computed = lap_f.scalar_part()
    
    # Analytical: ∇²f = -2sin(x)cos(y)
    analytical = -2 * np.sin(X) * np.cos(Y)
    
    # Error
    err = np.abs(computed[5:-5, 5:-5] - analytical[5:-5, 5:-5]).max()
    
    print(f"Grid: {shape[0]}×{shape[1]}, spacing h={h}")
    print(f"Max error: {err:.2e}")
    print(f"Expected O(h²) ≈ {h**2:.2e}")
    print(f"Status: {'PASS' if err < 10*h**2 else 'FAIL'}")
    
    return err


def test_wave_equation():
    """Test leapfrog wave equation solver."""
    print("\n" + "=" * 60)
    print("TEST: 1D Wave Equation with Leapfrog")
    print("=" * 60)
    
    # 1D wave equation: ∂²u/∂t² = c²∂²u/∂x²
    basis = CliffordBasis(AlgebraSignature(1, 0, 0))  # 1D Euclidean
    
    nx = 100
    h = 0.1
    c = 1.0
    dt = 0.05  # CFL: c*dt/h < 1
    
    print(f"CFL number: {c*dt/h:.2f} (should be < 1)")
    
    nabla = SplitDifferentialOperator(basis, spacing=h)
    
    # Initial conditions: Gaussian pulse
    x = np.linspace(0, 10, nx)
    u0 = np.exp(-(x - 5)**2)
    u1 = np.exp(-((x - 5 - c*dt))**2)  # Approximate u at t=dt
    
    # Run simulation
    u_prev = u0.copy()
    u_curr = u1.copy()
    
    for step in range(50):
        u_next = nabla.leapfrog_wave(u_prev, u_curr, dt, c)
        u_prev = u_curr
        u_curr = u_next
    
    # Check that pulse has moved to the right
    peak_idx = np.argmax(u_curr)
    expected_peak = 5 + c * 50 * dt
    actual_peak = x[peak_idx]
    
    print(f"After 50 steps (t = {50*dt:.1f}):")
    print(f"Expected peak position: {expected_peak:.2f}")
    print(f"Actual peak position: {actual_peak:.2f}")
    print(f"Error: {abs(actual_peak - expected_peak):.3f}")
    print(f"Status: {'PASS' if abs(actual_peak - expected_peak) < 0.5 else 'FAIL'}")


def test_geometric_product():
    """Test basis vector geometric products."""
    print("\n" + "=" * 60)
    print("TEST: Geometric Product Rules")
    print("=" * 60)
    
    basis = CliffordBasis(EUCLIDEAN_3D)
    
    # Test: e₀ * e₁ should give (0,1) with sign +1
    sign, result = basis.basis_vector_product(0, (1,))
    print(f"e₀ * e₁ = {'+' if sign >= 0 else '-'}e{''.join(map(str, result))}")
    assert result == (0, 1), f"Expected (0,1), got {result}"
    assert sign == 1, f"Expected sign +1, got {sign}"
    
    # Test: e₁ * e₀ should give (0,1) with sign -1
    sign, result = basis.basis_vector_product(1, (0,))
    print(f"e₁ * e₀ = {'+' if sign >= 0 else '-'}e{''.join(map(str, result))}")
    assert result == (0, 1), f"Expected (0,1), got {result}"
    assert sign == -1, f"Expected sign -1, got {sign}"
    
    # Test: e₀ * e₀ should give scalar with sign +1 (Euclidean)
    sign, result = basis.basis_vector_product(0, (0,))
    print(f"e₀ * e₀ = {'+' if sign >= 0 else '-'}{result if result else '1'}")
    assert result == (), f"Expected (), got {result}"
    assert sign == 1, f"Expected sign +1, got {sign}"
    
    print("Status: PASS")


def run_all_tests():
    """Run all tests."""
    test_geometric_product()
    test_gradient_2d()
    test_laplacian_2d()
    test_wave_equation()
    
    print("\n" + "=" * 60)
    print("ALL TESTS COMPLETED")
    print("=" * 60)


if __name__ == "__main__":
    run_all_tests()
