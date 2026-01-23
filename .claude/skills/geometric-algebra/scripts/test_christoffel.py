"""
Compare connection bivector with Christoffel symbols for sphere.

For unit sphere in spherical coords (θ, φ):
  g_θθ = 1, g_φφ = sin²θ, g_θφ = 0

Christoffel symbols (nonzero only):
  Γ^θ_φφ = -sin θ cos θ
  Γ^φ_θφ = Γ^φ_φθ = cot θ = cos θ / sin θ

So:
  ∂_θ e_θ = Γ^θ_θθ e_θ + Γ^φ_θθ e_φ + b_θθ n = 0 + 0 + b_θθ n  (purely normal)
  ∂_θ e_φ = Γ^θ_θφ e_θ + Γ^φ_θφ e_φ + b_θφ n = 0 + cot θ e_φ + 0
  ∂_φ e_θ = Γ^θ_φθ e_θ + Γ^φ_φθ e_φ + b_φθ n = 0 + cot θ e_φ + 0
  ∂_φ e_φ = Γ^θ_φφ e_θ + Γ^φ_φφ e_φ + b_φφ n = -sin θ cos θ e_θ + 0 + b_φφ n

The tangent-space connection 1-form is:
  ω^i_j = Γ^i_kj dx^k

For the bivector approach, we have:
  ω_i = Γ^k_ij (e_k ∧ e^j)

But this needs to be done carefully with metric.
"""

import numpy as np
from riemannian_ga import sphere_frame, ConnectionBivector, TangentFrame

def christoffel_comparison():
    print("=" * 70)
    print("CHRISTOFFEL vs CONNECTION BIVECTOR")
    print("=" * 70)
    
    R = 1.0
    theta = np.pi / 3
    phi = 0.0
    coords = np.array([theta, phi])
    
    ct, st = np.cos(theta), np.sin(theta)
    cot_theta = ct / st
    
    # Frame vectors (coordinate basis, NOT unit vectors)
    e_theta = R * np.array([ct * np.cos(phi), ct * np.sin(phi), -st])
    e_phi = R * st * np.array([-np.sin(phi), np.cos(phi), 0.0])
    
    # Christoffel symbols for sphere
    Gamma_theta_thetatheta = 0
    Gamma_phi_thetatheta = 0
    Gamma_theta_thetaphi = 0
    Gamma_phi_thetaphi = cot_theta
    Gamma_theta_phitheta = 0
    Gamma_phi_phitheta = cot_theta
    Gamma_theta_phiphi = -st * ct
    Gamma_phi_phiphi = 0
    
    print(f"\nθ = {np.degrees(theta):.1f}°")
    print(f"cot θ = {cot_theta:.4f}")
    print(f"sin θ cos θ = {st * ct:.4f}")
    
    print(f"\nChristoffel symbols:")
    print(f"  Γ^θ_θθ = {Gamma_theta_thetatheta}")
    print(f"  Γ^φ_θθ = {Gamma_phi_thetatheta}")
    print(f"  Γ^θ_θφ = {Gamma_theta_thetaphi}")
    print(f"  Γ^φ_θφ = {Gamma_phi_thetaphi:.4f} = cot θ")
    print(f"  Γ^θ_φφ = {Gamma_theta_phiphi:.4f} = -sin θ cos θ")
    print(f"  Γ^φ_φφ = {Gamma_phi_phiphi}")
    
    # Expected tangential derivatives from Christoffel:
    d_theta_e_theta_tang = Gamma_theta_thetatheta * e_theta + Gamma_phi_thetatheta * e_phi
    d_theta_e_phi_tang = Gamma_theta_thetaphi * e_theta + Gamma_phi_thetaphi * e_phi
    d_phi_e_theta_tang = Gamma_theta_phitheta * e_theta + Gamma_phi_phitheta * e_phi
    d_phi_e_phi_tang = Gamma_theta_phiphi * e_theta + Gamma_phi_phiphi * e_phi
    
    print(f"\nTangential parts from Christoffel (Γ^k_ij e_k):")
    print(f"  ∂_θ e_θ|_tang = {d_theta_e_theta_tang}")
    print(f"  ∂_θ e_φ|_tang = {d_theta_e_phi_tang}")
    print(f"  ∂_φ e_θ|_tang = {d_phi_e_theta_tang}")
    print(f"  ∂_φ e_φ|_tang = {d_phi_e_phi_tang}")
    
    # Get connection bivector
    frame_func = sphere_frame(R)
    conn = ConnectionBivector(manifold_dim=2, ambient_dim=3, h=1e-7)
    omega = conn.compute_at(coords, frame_func)
    
    omega_theta = omega[0].to_array()
    omega_phi = omega[1].to_array()
    
    # Compute ω × e
    print(f"\nFrom ω × e:")
    print(f"  ω_θ × e_θ = {np.cross(omega_theta, e_theta)}")
    print(f"  ω_θ × e_φ = {np.cross(omega_theta, e_phi)}")
    print(f"  ω_φ × e_θ = {np.cross(omega_phi, e_theta)}")
    print(f"  ω_φ × e_φ = {np.cross(omega_phi, e_phi)}")
    
    print(f"\nComparison:")
    for name, christoffel_val, omega_val in [
        ("∂_θ e_θ|_tang", d_theta_e_theta_tang, np.cross(omega_theta, e_theta)),
        ("∂_θ e_φ|_tang", d_theta_e_phi_tang, np.cross(omega_theta, e_phi)),
        ("∂_φ e_θ|_tang", d_phi_e_theta_tang, np.cross(omega_phi, e_theta)),
        ("∂_φ e_φ|_tang", d_phi_e_phi_tang, np.cross(omega_phi, e_phi)),
    ]:
        error = np.linalg.norm(christoffel_val - omega_val)
        match = "✓" if error < 0.01 else "✗"
        print(f"  {name}: error = {error:.6f} {match}")


if __name__ == "__main__":
    christoffel_comparison()
