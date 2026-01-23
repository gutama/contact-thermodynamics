"""
Alternative approach: Compute connection as rotation angle in tangent plane.

For a 2D surface, the connection 1-form ω has only one independent component
(since there's only one rotation generator in 2D).

The connection should satisfy:
  ∇_i e_1 = ω_i e_2
  ∇_i e_2 = -ω_i e_1

where {e_1, e_2} is an orthonormal frame and ω_i is a scalar (rotation rate).

Equivalently, using bivector B = e_1 ∧ e_2:
  ∇_i e_j = ω_i × e_j where ω_i = ω_i B (scalar times pseudoscalar)
"""

import numpy as np
from riemannian_ga import sphere_frame, TangentFrame

def compute_connection_2d(coords, frame_func, h=1e-6):
    """
    Compute connection 1-form for 2D surface.
    
    For orthonormal frame {e_1, e_2}:
    ω_i = (∂_i e_1) · e_2 = -(∂_i e_2) · e_1
    
    The connection bivector is then: Ω_i = ω_i (e_1 ∧ e_2)
    """
    frame = frame_func(coords)
    
    # Get orthonormal frame from coordinate frame
    g = frame.g
    sqrt_g11 = np.sqrt(g[0, 0])
    sqrt_g22 = np.sqrt(g[1, 1])
    
    # Unit tangent vectors
    hat_e1 = frame.vectors[0] / sqrt_g11
    hat_e2 = frame.vectors[1] / sqrt_g22  # Assuming orthogonal
    
    omega = np.zeros(2)
    
    for i in range(2):
        coords_plus = coords.copy()
        coords_minus = coords.copy()
        coords_plus[i] += h
        coords_minus[i] -= h
        
        frame_plus = frame_func(coords_plus)
        frame_minus = frame_func(coords_minus)
        
        # Unit vectors at perturbed points
        hat_e1_plus = frame_plus.vectors[0] / np.sqrt(frame_plus.g[0, 0])
        hat_e1_minus = frame_minus.vectors[0] / np.sqrt(frame_minus.g[0, 0])
        
        # ∂_i (hat e_1) projected onto hat e_2
        d_hat_e1 = (hat_e1_plus - hat_e1_minus) / (2 * h)
        
        # ω_i = (∂_i hat_e_1) · hat_e_2
        omega[i] = np.dot(d_hat_e1, hat_e2)
    
    return omega


def test_sphere_connection_2d():
    print("=" * 70)
    print("2D CONNECTION (Spin Connection)")
    print("=" * 70)
    
    R = 1.0
    frame_func = sphere_frame(R)
    
    theta = np.pi / 3
    phi = 0.0
    coords = np.array([theta, phi])
    
    ct, st = np.cos(theta), np.sin(theta)
    cot_theta = ct / st
    
    # Get frame
    frame = frame_func(coords)
    
    # Orthonormal basis:
    # Unit e_θ direction and unit e_φ direction
    # For sphere: e_θ has length R, e_φ has length R sin θ
    
    hat_theta = frame.vectors[0] / R  # Divide by |e_θ| = R
    hat_phi = frame.vectors[1] / (R * st)  # Divide by |e_φ| = R sin θ
    
    print(f"\nθ = {np.degrees(theta):.1f}°")
    print(f"Orthonormal frame:")
    print(f"  ê_θ = {hat_theta}")
    print(f"  ê_φ = {hat_phi}")
    
    # Analytical:
    # The spin connection for spherical coords with orthonormal frame is:
    # ω_θ = 0 (moving along meridian doesn't rotate the frame)
    # ω_φ = cos θ (moving along parallel rotates the frame)
    
    omega_analytic = np.array([0.0, ct])
    
    # Compute numerically
    omega_numeric = compute_connection_2d(coords, frame_func, h=1e-6)
    
    print(f"\nSpin connection (scalar ω_i such that ∇_i ê_1 = ω_i ê_2):")
    print(f"  ω_θ numeric   = {omega_numeric[0]:.6f}")
    print(f"  ω_θ analytic  = {omega_analytic[0]:.6f}")
    print(f"  ω_φ numeric   = {omega_numeric[1]:.6f}")  
    print(f"  ω_φ analytic  = {omega_analytic[1]:.6f} = cos θ")
    
    # The holonomy around a latitude circle:
    # Φ = ∮ ω_φ dφ = ∫_0^{2π} cos θ dφ = 2π cos θ
    # (measured in the rotating frame)
    # The actual rotation angle in inertial frame = 2π - 2π cos θ = 2π(1 - cos θ)
    
    holonomy_defect = 2 * np.pi * (1 - ct)
    print(f"\nHolonomy around latitude:")
    print(f"  2π(1 - cos θ) = {np.degrees(holonomy_defect):.2f}°")
    
    # Verify Christoffel symbol relation
    # For coordinate frame (not orthonormal):
    # The torsion-free condition with Christoffel symbols gives:
    # Γ^φ_θφ = cot θ = derivative of ln(sin θ) w.r.t θ
    
    print(f"\nChristoffel comparison:")
    print(f"  Γ^φ_θφ = cot θ = {cot_theta:.6f}")
    # The spin connection ω_φ = cos θ relates via the metric:
    # For orthonormal frame rotation rate: ω = g^{-1/2} Γ g^{1/2} type relation


if __name__ == "__main__":
    test_sphere_connection_2d()
