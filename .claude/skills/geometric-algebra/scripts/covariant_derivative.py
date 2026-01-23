"""
Covariant Derivative for Geometric Calculus on Manifolds

This module provides covariant derivatives for vector/multivector fields
on curved manifolds, inspired by Cartan.jl's approach.

Key concepts:
- Christoffel symbols Γʲᵢₖ computed from metric tensor
- Covariant derivative: ∇ᵢVʲ = ∂ᵢVʲ + Γʲᵢₖ Vᵏ
- Parallel transport along curves
- Geodesic equation solver

Author: Antigravity (based on Grassmann.jl/Cartan.jl concepts)
"""

import numpy as np
from typing import Tuple, Callable, Optional, Union
from dataclasses import dataclass


@dataclass
class ManifoldMetric:
    """
    Metric tensor on a manifold with coordinate chart.
    
    Attributes:
        dim: Dimension of the manifold
        metric_func: Function (coords) -> g_ij matrix at that point
        coord_names: Names of coordinates (for display)
    """
    dim: int
    metric_func: Callable[[np.ndarray], np.ndarray]
    coord_names: Tuple[str, ...] = None
    
    def __post_init__(self):
        if self.coord_names is None:
            self.coord_names = tuple(f"x{i}" for i in range(self.dim))
    
    def metric_at(self, coords: np.ndarray) -> np.ndarray:
        """Compute metric tensor g_ij at given coordinates."""
        return self.metric_func(coords)
    
    def inverse_metric_at(self, coords: np.ndarray) -> np.ndarray:
        """Compute inverse metric g^ij at given coordinates."""
        g = self.metric_at(coords)
        return np.linalg.inv(g)
    
    def sqrt_det_g(self, coords: np.ndarray) -> float:
        """Compute √|det(g)| for volume element."""
        g = self.metric_at(coords)
        return np.sqrt(np.abs(np.linalg.det(g)))


class ChristoffelSymbols:
    """
    Compute Christoffel symbols of the second kind from a metric.
    
    Γᵏᵢⱼ = (1/2) gᵏˡ (∂ᵢgⱼˡ + ∂ⱼgᵢˡ - ∂ˡgᵢⱼ)
    """
    
    def __init__(self, manifold: ManifoldMetric, h: float = 1e-6):
        """
        Args:
            manifold: The manifold with its metric
            h: Step size for numerical differentiation
        """
        self.manifold = manifold
        self.dim = manifold.dim
        self.h = h
    
    def _partial_metric(self, coords: np.ndarray, l: int) -> np.ndarray:
        """
        Compute ∂g_ij/∂x^l using central differences.
        
        Returns: (dim, dim) array of ∂g_ij/∂x^l
        """
        coords_plus = coords.copy()
        coords_minus = coords.copy()
        coords_plus[l] += self.h
        coords_minus[l] -= self.h
        
        g_plus = self.manifold.metric_at(coords_plus)
        g_minus = self.manifold.metric_at(coords_minus)
        
        return (g_plus - g_minus) / (2 * self.h)
    
    def compute_at(self, coords: np.ndarray) -> np.ndarray:
        """
        Compute all Christoffel symbols Γᵏᵢⱼ at given coordinates.
        
        Returns: (dim, dim, dim) array where result[k,i,j] = Γᵏᵢⱼ
        """
        dim = self.dim
        g_inv = self.manifold.inverse_metric_at(coords)
        
        # Compute all partial derivatives of metric
        dg = np.zeros((dim, dim, dim))  # dg[l,i,j] = ∂g_ij/∂x^l
        for l in range(dim):
            dg[l] = self._partial_metric(coords, l)
        
        # Christoffel symbols: Γᵏᵢⱼ = (1/2) gᵏˡ (∂ᵢgⱼˡ + ∂ⱼgᵢˡ - ∂ˡgᵢⱼ)
        Gamma = np.zeros((dim, dim, dim))
        for k in range(dim):
            for i in range(dim):
                for j in range(dim):
                    for l in range(dim):
                        Gamma[k, i, j] += 0.5 * g_inv[k, l] * (
                            dg[i, j, l] + dg[j, i, l] - dg[l, i, j]
                        )
        
        return Gamma


class CovariantDerivative:
    """
    Covariant derivative operator on a manifold.
    
    For a vector field V^j:
        ∇ᵢV^j = ∂ᵢV^j + Γʲᵢₖ V^k
    
    For a covector field (1-form) ω_j:
        ∇ᵢω_j = ∂ᵢω_j - Γᵏᵢⱼ ω_k
    """
    
    def __init__(self, manifold: ManifoldMetric, h: float = 1e-6):
        self.manifold = manifold
        self.christoffel = ChristoffelSymbols(manifold, h)
        self.dim = manifold.dim
        self.h = h
    
    def _partial_vector(self, V_func: Callable, coords: np.ndarray, i: int) -> np.ndarray:
        """Compute ∂V^j/∂x^i using central differences."""
        coords_plus = coords.copy()
        coords_minus = coords.copy()
        coords_plus[i] += self.h
        coords_minus[i] -= self.h
        
        V_plus = V_func(coords_plus)
        V_minus = V_func(coords_minus)
        
        return (V_plus - V_minus) / (2 * self.h)
    
    def of_vector(self, V_func: Callable[[np.ndarray], np.ndarray], 
                  coords: np.ndarray) -> np.ndarray:
        """
        Compute covariant derivative of vector field.
        
        Args:
            V_func: Function (coords) -> V^j (contravariant vector)
            coords: Point at which to compute
        
        Returns: (dim, dim) array where result[i,j] = ∇ᵢV^j
        """
        dim = self.dim
        Gamma = self.christoffel.compute_at(coords)
        V = V_func(coords)
        
        nabla_V = np.zeros((dim, dim))
        for i in range(dim):
            # Partial derivative ∂ᵢV^j
            partial_i_V = self._partial_vector(V_func, coords, i)
            
            for j in range(dim):
                # ∇ᵢV^j = ∂ᵢV^j + Γʲᵢₖ V^k
                nabla_V[i, j] = partial_i_V[j]
                for k in range(dim):
                    nabla_V[i, j] += Gamma[j, i, k] * V[k]
        
        return nabla_V
    
    def divergence(self, V_func: Callable[[np.ndarray], np.ndarray],
                   coords: np.ndarray) -> float:
        """
        Compute divergence of vector field.
        
        div V = ∇ᵢV^i = (1/√g) ∂ᵢ(√g V^i)
        
        Note: This formula automatically includes the connection terms.
        """
        dim = self.dim
        sqrt_g = self.manifold.sqrt_det_g(coords)
        V = V_func(coords)
        
        div = 0.0
        for i in range(dim):
            # ∂ᵢ(√g V^i) using central differences
            coords_plus = coords.copy()
            coords_minus = coords.copy()
            coords_plus[i] += self.h
            coords_minus[i] -= self.h
            
            sqrt_g_V_plus = self.manifold.sqrt_det_g(coords_plus) * V_func(coords_plus)[i]
            sqrt_g_V_minus = self.manifold.sqrt_det_g(coords_minus) * V_func(coords_minus)[i]
            
            div += (sqrt_g_V_plus - sqrt_g_V_minus) / (2 * self.h)
        
        return div / (sqrt_g + 1e-10)
    
    def laplacian_scalar(self, f_func: Callable[[np.ndarray], float],
                         coords: np.ndarray) -> float:
        """
        Compute Laplace-Beltrami operator on scalar field.
        
        Δf = div(grad f) = (1/√g) ∂ᵢ(√g g^ij ∂ⱼf)
        """
        dim = self.dim
        g_inv = self.manifold.inverse_metric_at(coords)
        sqrt_g = self.manifold.sqrt_det_g(coords)
        
        # First compute gradient (contravariant)
        def grad_f(c):
            """Gradient as contravariant vector: (grad f)^i = g^ij ∂ⱼf"""
            partial_f = np.zeros(dim)
            for j in range(dim):
                c_plus = c.copy()
                c_minus = c.copy()
                c_plus[j] += self.h
                c_minus[j] -= self.h
                partial_f[j] = (f_func(c_plus) - f_func(c_minus)) / (2 * self.h)
            
            g_inv_c = self.manifold.inverse_metric_at(c)
            return g_inv_c @ partial_f
        
        # Then compute divergence of gradient
        return self.divergence(grad_f, coords)


class ParallelTransport:
    """
    Parallel transport of vectors along curves on a manifold.
    
    Solves: dV^k/dt + Γᵏᵢⱼ (dx^i/dt) V^j = 0
    """
    
    def __init__(self, manifold: ManifoldMetric, h: float = 1e-6):
        self.manifold = manifold
        self.christoffel = ChristoffelSymbols(manifold, h)
        self.dim = manifold.dim
    
    def _parallel_transport_rhs(self, t: float, V: np.ndarray, 
                                   curve: Callable, dt: float) -> np.ndarray:
        """Compute dV^k/dt = -Γᵏᵢⱼ (dx^i/dt) V^j"""
        coords = curve(t)
        
        # Compute curve velocity dx^i/dt using central difference
        coords_plus = curve(t + dt/10)
        coords_minus = curve(t - dt/10)
        velocity = (coords_plus - coords_minus) / (dt/5)
        
        # Compute Christoffel symbols
        Gamma = self.christoffel.compute_at(coords)
        
        # dV^k/dt = -Γᵏᵢⱼ (dx^i/dt) V^j
        dV = np.zeros(self.dim)
        for k in range(self.dim):
            for i in range(self.dim):
                for j in range(self.dim):
                    dV[k] -= Gamma[k, i, j] * velocity[i] * V[j]
        
        return dV
    
    def transport(self, curve: Callable[[float], np.ndarray],
                  V_initial: np.ndarray,
                  t_start: float, t_end: float, 
                  n_steps: int = 100) -> Tuple[np.ndarray, np.ndarray]:
        """
        Transport vector V along curve from t_start to t_end using RK4.
        
        Args:
            curve: Parameterized curve γ(t) -> coordinates
            V_initial: Initial vector at t_start
            t_start, t_end: Parameter range
            n_steps: Number of integration steps
        
        Returns:
            t_values: Array of t values
            V_values: Array of transported vectors at each t
        """
        dt = (t_end - t_start) / n_steps
        t_values = np.linspace(t_start, t_end, n_steps + 1)
        V_values = np.zeros((n_steps + 1, self.dim))
        V_values[0] = V_initial.copy()
        
        for step in range(n_steps):
            t = t_values[step]
            V = V_values[step]
            
            # RK4 integration
            k1 = self._parallel_transport_rhs(t, V, curve, dt)
            k2 = self._parallel_transport_rhs(t + dt/2, V + dt/2 * k1, curve, dt)
            k3 = self._parallel_transport_rhs(t + dt/2, V + dt/2 * k2, curve, dt)
            k4 = self._parallel_transport_rhs(t + dt, V + dt * k3, curve, dt)
            
            V_values[step + 1] = V + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)
        
        return t_values, V_values


class GeodesicSolver:
    """
    Solve geodesic equations on a manifold.
    
    Geodesic equation: d²x^k/dt² + Γᵏᵢⱼ (dx^i/dt)(dx^j/dt) = 0
    
    Geodesics are the "straightest possible" curves on a curved manifold,
    generalizing straight lines in flat space.
    """
    
    def __init__(self, manifold: ManifoldMetric, h: float = 1e-6):
        self.manifold = manifold
        self.christoffel = ChristoffelSymbols(manifold, h)
        self.dim = manifold.dim
    
    def _geodesic_rhs(self, state: np.ndarray) -> np.ndarray:
        """
        Compute right-hand side of geodesic ODE system.
        
        State = [x^0, x^1, ..., x^{n-1}, v^0, v^1, ..., v^{n-1}]
        where v^i = dx^i/dt
        
        Returns: [dx/dt, dv/dt] = [v, -Γᵏᵢⱼ v^i v^j]
        """
        dim = self.dim
        x = state[:dim]
        v = state[dim:]
        
        # Compute Christoffel symbols at current position
        Gamma = self.christoffel.compute_at(x)
        
        # dv^k/dt = -Γᵏᵢⱼ v^i v^j
        dv = np.zeros(dim)
        for k in range(dim):
            for i in range(dim):
                for j in range(dim):
                    dv[k] -= Gamma[k, i, j] * v[i] * v[j]
        
        # dx/dt = v, dv/dt = -Γ v v
        return np.concatenate([v, dv])
    
    def solve(self, x0: np.ndarray, v0: np.ndarray, 
              t_final: float, n_steps: int = 100) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Solve geodesic equation with initial position and velocity.
        
        Args:
            x0: Initial position (coordinates)
            v0: Initial velocity (tangent vector)
            t_final: Final time parameter
            n_steps: Number of integration steps
        
        Returns:
            t_values: Array of time values
            x_values: Array of positions at each time (n_steps+1, dim)
            v_values: Array of velocities at each time (n_steps+1, dim)
        """
        dim = self.dim
        dt = t_final / n_steps
        t_values = np.linspace(0, t_final, n_steps + 1)
        
        x_values = np.zeros((n_steps + 1, dim))
        v_values = np.zeros((n_steps + 1, dim))
        x_values[0] = x0.copy()
        v_values[0] = v0.copy()
        
        state = np.concatenate([x0, v0])
        
        for step in range(n_steps):
            # RK4 integration
            k1 = self._geodesic_rhs(state)
            k2 = self._geodesic_rhs(state + dt/2 * k1)
            k3 = self._geodesic_rhs(state + dt/2 * k2)
            k4 = self._geodesic_rhs(state + dt * k3)
            
            state = state + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)
            
            x_values[step + 1] = state[:dim]
            v_values[step + 1] = state[dim:]
        
        return t_values, x_values, v_values
    
    def arc_length(self, x_values: np.ndarray) -> float:
        """
        Compute arc length of a curve using the metric.
        
        Length = ∫ √(g_ij dx^i dx^j)
        """
        n_points = len(x_values)
        length = 0.0
        
        for i in range(n_points - 1):
            x = x_values[i]
            dx = x_values[i + 1] - x_values[i]
            g = self.manifold.metric_at(x)
            ds = np.sqrt(dx @ g @ dx)
            length += ds
        
        return length


# =============================================================================
# STANDARD MANIFOLDS
# =============================================================================

def sphere_2d(R: float = 1.0) -> ManifoldMetric:
    """
    2-sphere S² with standard metric.
    
    Coordinates: (θ, φ) where θ ∈ [0, π], φ ∈ [0, 2π)
    Metric: ds² = R²(dθ² + sin²θ dφ²)
    """
    def metric(coords):
        theta = coords[0]
        return R**2 * np.array([
            [1.0, 0.0],
            [0.0, np.sin(theta)**2 + 1e-10]  # Regularize at poles
        ])
    
    return ManifoldMetric(dim=2, metric_func=metric, coord_names=('θ', 'φ'))


def torus_2d(R: float = 2.0, r: float = 1.0) -> ManifoldMetric:
    """
    Torus with standard metric.
    
    Coordinates: (θ, φ) where θ, φ ∈ [0, 2π)
    Metric: ds² = r²dθ² + (R + r cos θ)² dφ²
    """
    def metric(coords):
        theta = coords[0]
        return np.array([
            [r**2, 0.0],
            [0.0, (R + r * np.cos(theta))**2]
        ])
    
    return ManifoldMetric(dim=2, metric_func=metric, coord_names=('θ', 'φ'))


def hyperbolic_2d() -> ManifoldMetric:
    """
    Hyperbolic plane (Poincaré half-plane model).
    
    Coordinates: (x, y) where y > 0
    Metric: ds² = (dx² + dy²) / y²
    """
    def metric(coords):
        y = coords[1]
        y_safe = max(y, 1e-10)
        return np.array([
            [1.0 / y_safe**2, 0.0],
            [0.0, 1.0 / y_safe**2]
        ])
    
    return ManifoldMetric(dim=2, metric_func=metric, coord_names=('x', 'y'))


# =============================================================================
# TESTS
# =============================================================================

def test_sphere_christoffel():
    """Test Christoffel symbols on sphere."""
    print("=" * 60)
    print("TEST: Christoffel Symbols on 2-Sphere")
    print("=" * 60)
    
    sphere = sphere_2d(R=1.0)
    christoffel = ChristoffelSymbols(sphere)
    
    # At θ = π/4 (45°)
    coords = np.array([np.pi/4, 0.0])
    Gamma = christoffel.compute_at(coords)
    
    # Known values for sphere:
    # Γ^θ_φφ = -sin(θ)cos(θ)
    # Γ^φ_θφ = Γ^φ_φθ = cot(θ)
    theta = coords[0]
    expected_Gamma_theta_phiphi = -np.sin(theta) * np.cos(theta)
    expected_Gamma_phi_thetaphi = np.cos(theta) / np.sin(theta)  # cot(θ)
    
    print(f"\nAt θ = π/4:")
    print(f"  Γ^θ_φφ computed: {Gamma[0,1,1]:.6f}")
    print(f"  Γ^θ_φφ expected: {expected_Gamma_theta_phiphi:.6f}")
    print(f"  Error: {abs(Gamma[0,1,1] - expected_Gamma_theta_phiphi):.2e}")
    
    print(f"\n  Γ^φ_θφ computed: {Gamma[1,0,1]:.6f}")
    print(f"  Γ^φ_θφ expected: {expected_Gamma_phi_thetaphi:.6f}")
    print(f"  Error: {abs(Gamma[1,0,1] - expected_Gamma_phi_thetaphi):.2e}")
    
    error1 = abs(Gamma[0,1,1] - expected_Gamma_theta_phiphi)
    error2 = abs(Gamma[1,0,1] - expected_Gamma_phi_thetaphi)
    
    print(f"\nStatus: {'PASS' if error1 < 0.01 and error2 < 0.01 else 'FAIL'}")


def test_sphere_laplacian():
    """Test Laplacian on sphere with known eigenfunction."""
    print("\n" + "=" * 60)
    print("TEST: Laplace-Beltrami on 2-Sphere")
    print("=" * 60)
    
    sphere = sphere_2d(R=1.0)
    cov_deriv = CovariantDerivative(sphere)
    
    # Test function: Y_1^0 = cos(θ), eigenvalue = -2/R² = -2
    def f(coords):
        return np.cos(coords[0])
    
    # Test at θ = π/3 (60°)
    coords = np.array([np.pi/3, 0.0])
    
    lap_f = cov_deriv.laplacian_scalar(f, coords)
    expected = -2 * np.cos(coords[0])  # Δf = -2f for Y_1^0
    
    print(f"\nTest function: f(θ,φ) = cos(θ)")
    print(f"Expected: Δf = -2cos(θ) at θ=π/3")
    print(f"  Δf computed: {lap_f:.6f}")
    print(f"  Δf expected: {expected:.6f}")
    
    relative_error = abs(lap_f - expected) / abs(expected)
    print(f"  Relative error: {relative_error:.2%}")
    
    print(f"\nStatus: {'PASS' if relative_error < 0.05 else 'FAIL'}")


def test_parallel_transport():
    """Test parallel transport on sphere."""
    print("\n" + "=" * 60)
    print("TEST: Parallel Transport on 2-Sphere")
    print("=" * 60)
    
    sphere = sphere_2d(R=1.0)
    transport = ParallelTransport(sphere)
    
    # Transport vector along latitude circle at θ = π/4
    theta_fixed = np.pi / 4
    
    def latitude_curve(t):
        """Latitude circle at fixed θ."""
        return np.array([theta_fixed, t])
    
    # Initial vector pointing in θ direction (tangent to sphere, perpendicular to curve)
    V_initial = np.array([1.0, 0.0])
    
    # Transport around full circle
    t_values, V_values = transport.transport(
        latitude_curve, V_initial, 
        t_start=0, t_end=2*np.pi, n_steps=200
    )
    
    # After full circle, vector should be rotated by holonomy angle
    # Holonomy for latitude at θ is: 2π(1 - cos(θ))
    # For θ = π/4: holonomy ≈ 2π(1 - √2/2) ≈ 1.84 radians
    
    V_final = V_values[-1]
    
    # Parallel transport preserves METRIC norm: ||V||² = g_ij V^i V^j
    # At the endpoint (same as start), compute metric norm
    coords_final = latitude_curve(2*np.pi)
    g = sphere.metric_at(coords_final)
    metric_norm_initial = np.sqrt(V_initial @ g @ V_initial)
    metric_norm_final = np.sqrt(V_final @ g @ V_final)
    
    print(f"\nCurve: Latitude circle at θ = π/4")
    print(f"Initial vector: {V_initial}")
    print(f"Final vector: {V_final}")
    print(f"Metric |V|_g initial: {metric_norm_initial:.6f}")
    print(f"Metric |V|_g final: {metric_norm_final:.6f}")
    print(f"Metric norm preserved: {abs(metric_norm_final - metric_norm_initial) < 0.05}")
    
    # Compute rotation angle in the tangent plane
    # The holonomy is the angle of rotation
    angle = np.arctan2(V_final[1] * np.sqrt(g[1,1]), V_final[0] * np.sqrt(g[0,0]))
    expected_holonomy = 2 * np.pi * (1 - np.cos(theta_fixed))
    
    print(f"\nHolonomy angle computed: {angle:.4f} rad")
    print(f"Holonomy angle expected: {expected_holonomy:.4f} rad")
    
    norm_error = abs(metric_norm_final - metric_norm_initial)
    print(f"\nStatus: {'PASS' if norm_error < 0.05 else 'FAIL'}")


def test_geodesic():
    """Test geodesic solver on sphere."""
    print("\n" + "=" * 60)
    print("TEST: Geodesics on 2-Sphere")
    print("=" * 60)
    
    sphere = sphere_2d(R=1.0)
    solver = GeodesicSolver(sphere)
    
    # Test 1: Geodesic along meridian (great circle through poles)
    # Starting at θ = π/4, φ = 0, moving in θ direction
    x0 = np.array([np.pi/4, 0.0])
    v0 = np.array([1.0, 0.0])  # Unit velocity in θ direction
    
    t_final = np.pi / 2  # Travel quarter circle
    t_values, x_values, v_values = solver.solve(x0, v0, t_final, n_steps=100)
    
    # For meridian geodesic: θ should increase linearly, φ should stay constant
    theta_final = x_values[-1, 0]
    phi_final = x_values[-1, 1]
    theta_expected = np.pi/4 + np.pi/2  # = 3π/4
    
    print(f"\nTest 1: Meridian geodesic (great circle)")
    print(f"  Start: θ = π/4 = {np.pi/4:.4f}")
    print(f"  Final θ computed: {theta_final:.4f}")
    print(f"  Final θ expected: {theta_expected:.4f}")
    print(f"  Final φ (should stay ~0): {phi_final:.6f}")
    
    theta_error = abs(theta_final - theta_expected)
    phi_drift = abs(phi_final)
    
    # Test 2: Verify arc length matches distance traveled
    arc_len = solver.arc_length(x_values)
    expected_len = t_final  # For unit speed geodesic on unit sphere
    
    print(f"\nTest 2: Arc length preservation")
    print(f"  Arc length computed: {arc_len:.4f}")
    print(f"  Arc length expected: {expected_len:.4f}")
    print(f"  Error: {abs(arc_len - expected_len):.4f}")
    
    arc_error = abs(arc_len - expected_len)
    
    # Test 3: Velocity magnitude preserved
    g_start = sphere.metric_at(x_values[0])
    g_end = sphere.metric_at(x_values[-1])
    speed_start = np.sqrt(v_values[0] @ g_start @ v_values[0])
    speed_end = np.sqrt(v_values[-1] @ g_end @ v_values[-1])
    
    print(f"\nTest 3: Speed preservation (geodesic affinely parameterized)")
    print(f"  Initial speed: {speed_start:.6f}")
    print(f"  Final speed: {speed_end:.6f}")
    
    speed_error = abs(speed_end - speed_start)
    
    passed = theta_error < 0.05 and phi_drift < 0.01 and arc_error < 0.1 and speed_error < 0.01
    print(f"\nStatus: {'PASS' if passed else 'FAIL'}")


def run_all_tests():
    """Run all tests."""
    test_sphere_christoffel()
    test_sphere_laplacian()
    test_parallel_transport()
    test_geodesic()
    
    print("\n" + "=" * 60)
    print("ALL TESTS COMPLETED")
    print("=" * 60)


if __name__ == "__main__":
    run_all_tests()



# =============================================================================
# GA INTERFACE (Coordinate-Free)
# =============================================================================

def get_ga_interface(manifold: ManifoldMetric):
    """
    Get coordinate-free GA interface for a manifold.
    
    This wraps the traditional implementation with GA semantics.
    
    Usage:
        sphere = sphere_2d(R=1.0)
        ga = get_ga_interface(sphere)
        
        omega = ga.connection_bivector(coords)  # Instead of Christoffel
        Omega = ga.curvature_2form(coords)      # Instead of Riemann tensor
    """
    from riemannian_ga import (
        ConnectionBivector, Curvature2Form, 
        GACovariantDerivative, GAGeodesicSolver, GAParallelTransport
    )
    
    class GAInterface:
        def __init__(self, manifold):
            # Create RiemannianManifold wrapper
            self.manifold = manifold
            self.dim = manifold.dim
        
        def connection_bivector(self) -> ConnectionBivector:
            """Get GA connection bivector operator."""
            return ConnectionBivector(self.manifold)
        
        def curvature_2form(self) -> Curvature2Form:
            """Get GA curvature 2-form operator."""
            return Curvature2Form(self.manifold)

        def covariant_derivative(self) -> GACovariantDerivative:
            """Get GA covariant derivative operator."""
            return GACovariantDerivative(self.manifold)
        
        def geodesic_solver(self) -> GAGeodesicSolver:
            """Get GA geodesic solver."""
            return GAGeodesicSolver(self.manifold)
        
        def parallel_transport(self) -> GAParallelTransport:
            """Get GA parallel transport operator."""
            return GAParallelTransport(self.manifold)
        
        
    return GAInterface(manifold)
