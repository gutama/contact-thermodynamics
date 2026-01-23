"""
Sanity test for connection bivector sign fix.

Tests:
1. Sphere (R=1): K = 1/R² = 1
2. Hyperbolic half-plane: K = -1

For each, we verify:
- Connection bivector ω
- Gaussian curvature K (via shape operator and curvature 2-form)
- Holonomy angle matches ∫K dA
"""

import numpy as np
from riemannian_ga import (
    TangentFrame, Bivector3D,
    ConnectionBivector, Curvature2Form,
    GAGeodesicSolver, GAParallelTransport, ShapeOperator,
    sphere_frame
)
from typing import Callable

# =============================================================================
# HYPERBOLIC HALF-PLANE (Poincaré upper half-plane)
# =============================================================================

def hyperbolic_frame(scale: float = 1.0) -> Callable[[np.ndarray], TangentFrame]:
    """
    Create frame function for hyperbolic half-plane (Poincaré model).
    
    Metric: ds² = (dx² + dy²) / y²
    
    Coordinates: (x, y) where y > 0
    
    The frame is embedded in 3D as (x, y, 0) with orthonormal tangent vectors
    scaled by 1/y to give the hyperbolic metric.
    
    Gaussian curvature: K = -1
    """
    def frame_func(coords: np.ndarray) -> TangentFrame:
        x, y = coords[0], coords[1]
        
        # Ensure y > 0 (hyperbolic plane)
        y = max(y, 0.1)
        
        # In Poincaré metric: ds² = (dx² + dy²)/y²
        # So frame vectors are: e_x = (1/y, 0, 0), e_y = (0, 1/y, 0)
        # These give: g_11 = e_x·e_x = 1/y², g_22 = 1/y², g_12 = 0
        
        # But we need to embed in 3D. Use z=0 plane.
        # The frame vectors in ambient R³:
        e_x = np.array([1.0/y, 0.0, 0.0]) * scale
        e_y = np.array([0.0, 1.0/y, 0.0]) * scale
        
        vectors = np.array([e_x, e_y])
        return TangentFrame(vectors)
    
    return frame_func


def test_sphere():
    """Test sphere curvature with sign fix."""
    print("=" * 70)
    print("TEST: Unit Sphere (R=1)")
    print("Expected: K = 1.0, H = 1.0")
    print("=" * 70)
    
    R = 1.0
    frame_func = sphere_frame(R)
    
    # Create objects
    conn = ConnectionBivector(manifold_dim=2, ambient_dim=3, h=1e-5)
    shape = ShapeOperator(h=1e-5)
    
    # Test at several points away from poles
    test_points = [
        (np.pi/4, 0.0, "θ=45°, φ=0°"),
        (np.pi/3, np.pi/4, "θ=60°, φ=45°"),
        (np.pi/2, 0.0, "θ=90° (equator)"),
    ]
    
    for theta, phi, name in test_points:
        coords = np.array([theta, phi])
        frame = frame_func(coords)
        
        # Connection bivector
        omega = conn.compute_at(coords, frame_func)
        
        # Shape operator
        S = shape.compute_at(coords, frame_func)
        K = shape.gaussian_curvature(S)
        H = shape.mean_curvature(S)
        kappa1, kappa2 = shape.principal_curvatures(S)
        
        # Check frame vectors are correct
        e_theta_norm = np.linalg.norm(frame.vectors[0])
        e_phi_norm = np.linalg.norm(frame.vectors[1])
        
        print(f"\n{name}:")
        print(f"  |e_θ| = {e_theta_norm:.4f} (expected: {R:.4f})")
        print(f"  |e_φ| = {e_phi_norm:.4f} (expected: {R*np.sin(theta):.4f})")
        print(f"  K = {K:+.6f} (expected: +1.000000)")
        print(f"  H = {H:+.6f} (expected: +1.000000 or -1.000000)")
        print(f"  κ₁ = {kappa1:+.6f}, κ₂ = {kappa2:+.6f}")
        print(f"  ω_θ = {omega[0]}")
        print(f"  ω_φ = {omega[1]}")
        
        # Verify |K - 1| < tolerance
        if abs(K - 1.0) < 0.01:
            print(f"  ✓ K is correct")
        else:
            print(f"  ✗ K error = {abs(K - 1.0):.6f}")


def test_hyperbolic():
    """Test hyperbolic half-plane curvature."""
    print("\n" + "=" * 70)
    print("TEST: Hyperbolic Half-Plane (Poincaré)")
    print("Expected: K = -1.0")
    print("=" * 70)
    
    frame_func = hyperbolic_frame()
    
    # Create objects
    conn = ConnectionBivector(manifold_dim=2, ambient_dim=3, h=1e-5)
    shape = ShapeOperator(h=1e-5)
    
    # Test at several points (y > 0 always)
    test_points = [
        (0.0, 1.0, "x=0, y=1"),
        (1.0, 1.0, "x=1, y=1"),
        (0.0, 2.0, "x=0, y=2"),
        (0.0, 0.5, "x=0, y=0.5"),
    ]
    
    for x_val, y_val, name in test_points:
        coords = np.array([x_val, y_val])
        frame = frame_func(coords)
        
        # Connection bivector
        omega = conn.compute_at(coords, frame_func)
        
        # Check metric
        g11 = np.dot(frame.vectors[0], frame.vectors[0])
        g22 = np.dot(frame.vectors[1], frame.vectors[1])
        g12 = np.dot(frame.vectors[0], frame.vectors[1])
        
        expected_g = 1.0 / (y_val ** 2)
        
        print(f"\n{name}:")
        print(f"  g₁₁ = {g11:.6f} (expected: {expected_g:.6f})")
        print(f"  g₂₂ = {g22:.6f} (expected: {expected_g:.6f})")
        print(f"  g₁₂ = {g12:.6f} (expected: 0)")
        print(f"  ω_x = {omega[0]}")
        print(f"  ω_y = {omega[1]}")
        
        # For the hyperbolic plane, the connection components should satisfy:
        # In the Poincaré half-plane with metric g = (1/y²) δ:
        # The Christoffel symbols are: Γ^x_xy = Γ^x_yx = -1/y, Γ^y_xx = 1/y, Γ^y_yy = -1/y
        # Others are zero.
        # This translates to connection bivector components.


def test_geodesic_sphere():
    """Test that geodesics on sphere are great circles."""
    print("\n" + "=" * 70)
    print("TEST: Geodesic on Sphere (should be great circle)")
    print("=" * 70)
    
    R = 1.0
    frame_func = sphere_frame(R)
    conn = ConnectionBivector(manifold_dim=2, ambient_dim=3, h=1e-5)
    solver = GAGeodesicSolver(conn)
    
    # Start at equator (θ=π/2), moving in θ direction (along meridian)
    x0 = np.array([np.pi/2, 0.0])
    v0 = np.array([1.0, 0.0])  # velocity along θ (meridian)
    
    # Solve for half a great circle
    t_vals, x_vals, v_vals = solver.solve(x0, v0, frame_func, t_final=np.pi/2, n_steps=100)
    
    print(f"\nInitial: θ = {np.degrees(x0[0]):.1f}°, φ = {np.degrees(x0[1]):.1f}°")
    print(f"Final:   θ = {np.degrees(x_vals[-1][0]):.1f}°, φ = {np.degrees(x_vals[-1][1]):.1f}°")
    
    # A meridian geodesic at φ=0 should stay at φ=0
    phi_drift = abs(x_vals[-1][1] - x_vals[0][1])
    print(f"φ drift = {np.degrees(phi_drift):.2f}° (expected: ~0°)")
    
    # θ should change by about the arc length traveled
    theta_change = abs(x_vals[-1][0] - x_vals[0][0])
    expected_theta_change = np.pi/2  # We integrated t_final = π/2
    print(f"θ change = {np.degrees(theta_change):.1f}° (traveled: {np.degrees(expected_theta_change):.1f}°)")


def test_holonomy():
    """Test holonomy matches integrated curvature."""
    print("\n" + "=" * 70)
    print("TEST: Holonomy on Sphere (should match ∫K dA)")
    print("=" * 70)
    
    R = 1.0
    frame_func = sphere_frame(R)
    conn = ConnectionBivector(manifold_dim=2, ambient_dim=3, h=1e-5)
    transport = GAParallelTransport(conn)
    
    # Transport around latitude at θ = π/4
    theta_fixed = np.pi / 4
    
    def latitude_loop(t):
        return np.array([theta_fixed, t])
    
    # Initial vector in θ direction
    w0 = np.array([1.0, 0.0])
    
    t_vals, w_vals = transport.transport(
        latitude_loop, w0, frame_func,
        t_start=0.0, t_end=2*np.pi, n_steps=400
    )
    
    w_final = w_vals[-1]
    
    # Compute angle between initial and final
    # Need to compare in same frame
    frame_start = frame_func(latitude_loop(0))
    frame_end = frame_func(latitude_loop(2*np.pi))
    
    # Convert to ambient vectors
    w0_ambient = w0[0] * frame_start.vectors[0] + w0[1] * frame_start.vectors[1]
    wf_ambient = w_final[0] * frame_end.vectors[0] + w_final[1] * frame_end.vectors[1]
    
    norm0 = np.linalg.norm(w0_ambient)
    normf = np.linalg.norm(wf_ambient)
    
    if norm0 > 1e-10 and normf > 1e-10:
        cos_angle = np.dot(w0_ambient, wf_ambient) / (norm0 * normf)
        cos_angle = np.clip(cos_angle, -1, 1)
        angle = np.arccos(cos_angle)
    else:
        angle = 0.0
    
    # Expected holonomy: 2π(1 - cos θ)
    # This is the solid angle subtended by the spherical cap
    expected_angle = 2 * np.pi * (1 - np.cos(theta_fixed))
    
    print(f"\nLatitude circle at θ = {np.degrees(theta_fixed):.1f}°")
    print(f"Initial vector (coords): {w0}")
    print(f"Final vector (coords): {w_final}")
    print(f"Holonomy angle: {np.degrees(angle):.2f}°")
    print(f"Expected (2π(1-cosθ)): {np.degrees(expected_angle):.2f}°")
    print(f"Error: {np.degrees(abs(angle - expected_angle)):.2f}°")
    
    # Also check that vector magnitude is preserved
    norm_ratio = normf / norm0 if norm0 > 1e-10 else 0
    print(f"Magnitude ratio |w_final|/|w_0| = {norm_ratio:.6f} (expected: 1.0)")


def analytical_sphere_connection():
    """
    Verify connection bivector analytically for sphere.
    
    For a unit sphere, the analytical connection should satisfy:
    ∂_θ e_θ = 0  (e_θ is unit vector, constant length R=1)
    ∂_φ e_θ = cot(θ) e_φ
    ∂_θ e_φ = 0  
    ∂_φ e_φ = -sin(θ)cos(θ) e_θ - sin²(θ) n
    
    where n is the outward normal.
    """
    print("\n" + "=" * 70)
    print("ANALYTICAL CHECK: Sphere Connection Components")
    print("=" * 70)
    
    R = 1.0
    frame_func = sphere_frame(R)
    conn = ConnectionBivector(manifold_dim=2, ambient_dim=3, h=1e-6)
    
    theta = np.pi / 3  # 60°
    phi = 0.0
    coords = np.array([theta, phi])
    
    frame = frame_func(coords)
    omega = conn.compute_at(coords, frame_func)
    
    # Get frame at this point
    e_theta = frame.vectors[0]
    e_phi = frame.vectors[1]
    n = frame.normal_3d()
    
    print(f"\nAt θ = {np.degrees(theta):.1f}°:")
    print(f"  e_θ = {e_theta}")
    print(f"  e_φ = {e_phi}")
    print(f"  n = {n}")
    
    # Check frame derivative ∂_i e_j = ω_i × e_j
    for i, name_i in [(0, 'θ'), (1, 'φ')]:
        for j, name_j in [(0, 'θ'), (1, 'φ')]:
            # Numerical derivative
            coords_plus = coords.copy()
            coords_minus = coords.copy()
            coords_plus[i] += 1e-6
            coords_minus[i] -= 1e-6
            frame_plus = frame_func(coords_plus)
            frame_minus = frame_func(coords_minus)
            
            d_ej_numerical = (frame_plus.vectors[j] - frame_minus.vectors[j]) / (2e-6)
            
            # Via connection: ω_i × e_j
            d_ej_connection = omega[i].commutator_with_vector(frame.vectors[j])
            
            error = np.linalg.norm(d_ej_numerical - d_ej_connection)
            
            print(f"\n  ∂_{name_i} e_{name_j}:")
            print(f"    Numerical:   {d_ej_numerical}")
            print(f"    ω_{name_i} × e_{name_j}: {d_ej_connection}")
            print(f"    Error: {error:.2e}")


if __name__ == "__main__":
    test_sphere()
    test_hyperbolic()
    analytical_sphere_connection()
    test_geodesic_sphere()
    test_holonomy()
    
    print("\n" + "=" * 70)
    print("ALL SANITY TESTS COMPLETED")
    print("=" * 70)
