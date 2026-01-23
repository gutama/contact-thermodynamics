"""
Focused test for the connection bivector sign and factor.

The key equation is: ∂ᵢeⱼ = ωᵢ × eⱼ

For the sphere at θ = π/3, φ = 0:
  e_θ = (cos θ, 0, -sin θ) = (0.5, 0, -0.866)
  e_φ = sin θ (-sin φ, cos φ, 0) = (0, 0.866, 0)
  
Analytical derivatives:
  ∂_θ e_θ = (-sin θ, 0, -cos θ) = (-0.866, 0, -0.5)
  ∂_θ e_φ = cos θ (0, 1, 0) = (0, 0.5, 0)
  ∂_φ e_θ = (0, cos θ, 0) = (0, 0.5, 0)  # Wait, let me recalculate
  ∂_φ e_φ = sin θ (-cos φ, -sin φ, 0) = (-0.866, 0, 0)
  
Actually for e_θ = (cos θ cos φ, cos θ sin φ, -sin θ):
  ∂_θ e_θ = (-sin θ cos φ, -sin θ sin φ, -cos θ)
  ∂_φ e_θ = (-cos θ sin φ, cos θ cos φ, 0) = cos θ (-sin φ, cos φ, 0) = cos θ * (e_φ / sin θ)
          = cot θ * e_φ (for normalized)
          
For e_φ = sin θ (-sin φ, cos φ, 0):
  ∂_θ e_φ = cos θ (-sin φ, cos φ, 0)
  ∂_φ e_φ = sin θ (-cos φ, -sin φ, 0)
"""

import numpy as np
from riemannian_ga import (
    TangentFrame, Bivector3D,
    ConnectionBivector, sphere_frame
)

def detailed_sphere_test():
    """Detailed test of connection bivector on sphere."""
    print("=" * 70)
    print("DETAILED SPHERE CONNECTION TEST")
    print("=" * 70)
    
    R = 1.0
    theta = np.pi / 3  # 60°
    phi = 0.0
    coords = np.array([theta, phi])
    
    ct, st = np.cos(theta), np.sin(theta)
    cp, sp = np.cos(phi), np.sin(phi)
    
    # Analytical frame vectors
    e_theta = R * np.array([ct * cp, ct * sp, -st])
    e_phi = R * st * np.array([-sp, cp, 0.0])
    
    print(f"\nθ = {np.degrees(theta):.1f}°, φ = {np.degrees(phi):.1f}°")
    print(f"cos θ = {ct:.4f}, sin θ = {st:.4f}")
    print(f"\nAnalytical frame vectors:")
    print(f"  e_θ = {e_theta}")
    print(f"  e_φ = {e_phi}")
    
    # Analytical derivatives
    # ∂_θ e_θ = R * (-sin θ cos φ, -sin θ sin φ, -cos θ)
    d_theta_e_theta = R * np.array([-st * cp, -st * sp, -ct])
    # ∂_θ e_φ = R * cos θ * (-sin φ, cos φ, 0)
    d_theta_e_phi = R * ct * np.array([-sp, cp, 0.0])
    # ∂_φ e_θ = R * (-cos θ sin φ, cos θ cos φ, 0)
    d_phi_e_theta = R * np.array([-ct * sp, ct * cp, 0.0])
    # ∂_φ e_φ = R sin θ * (-cos φ, -sin φ, 0)
    d_phi_e_phi = R * st * np.array([-cp, -sp, 0.0])
    
    print(f"\nAnalytical frame derivatives:")
    print(f"  ∂_θ e_θ = {d_theta_e_theta}")
    print(f"  ∂_θ e_φ = {d_theta_e_phi}")
    print(f"  ∂_φ e_θ = {d_phi_e_theta}")
    print(f"  ∂_φ e_φ = {d_phi_e_phi}")
    
    # Compute reciprocal frame
    g = np.array([
        [np.dot(e_theta, e_theta), np.dot(e_theta, e_phi)],
        [np.dot(e_phi, e_theta), np.dot(e_phi, e_phi)]
    ])
    g_inv = np.linalg.inv(g)
    
    e_theta_recip = g_inv[0, 0] * e_theta + g_inv[0, 1] * e_phi
    e_phi_recip = g_inv[1, 0] * e_theta + g_inv[1, 1] * e_phi
    
    print(f"\nReciprocal frame (eⁱ such that eⁱ·eⱼ = δⁱⱼ):")
    print(f"  eᶿ = {e_theta_recip}")
    print(f"  eᵠ = {e_phi_recip}")
    print(f"  Check eᶿ·e_θ = {np.dot(e_theta_recip, e_theta):.6f} (should be 1)")
    print(f"  Check eᶿ·e_φ = {np.dot(e_theta_recip, e_phi):.6f} (should be 0)")
    print(f"  Check eᵠ·e_θ = {np.dot(e_phi_recip, e_theta):.6f} (should be 0)")
    print(f"  Check eᵠ·e_φ = {np.dot(e_phi_recip, e_phi):.6f} (should be 1)")
    
    # Compute connection bivector analytically
    # ω_θ = ½ Σⱼ eʲ ∧ (∂_θ eⱼ) = ½ [eᶿ ∧ (∂_θ e_θ) + eᵠ ∧ (∂_θ e_φ)]
    
    def wedge(u, v):
        """Compute wedge product as bivector components (e23, e31, e12)."""
        return np.array([
            u[1]*v[2] - u[2]*v[1],  # e23
            u[2]*v[0] - u[0]*v[2],  # e31
            u[0]*v[1] - u[1]*v[0]   # e12
        ])
    
    omega_theta = 0.5 * (wedge(e_theta_recip, d_theta_e_theta) + 
                         wedge(e_phi_recip, d_theta_e_phi))
    omega_phi = 0.5 * (wedge(e_theta_recip, d_phi_e_theta) + 
                       wedge(e_phi_recip, d_phi_e_phi))
    
    print(f"\nAnalytical connection bivectors:")
    print(f"  ω_θ = {omega_theta} (e23, e31, e12)")
    print(f"  ω_φ = {omega_phi} (e23, e31, e12)")
    
    # Now compute using the code
    frame_func = sphere_frame(R)
    conn = ConnectionBivector(manifold_dim=2, ambient_dim=3, h=1e-7)
    
    omega = conn.compute_at(coords, frame_func)
    omega_theta_code = omega[0].to_array()
    omega_phi_code = omega[1].to_array()
    
    print(f"\nCode-computed connection bivectors:")
    print(f"  ω_θ = {omega_theta_code} (e23, e31, e12)")
    print(f"  ω_φ = {omega_phi_code} (e23, e31, e12)")
    
    print(f"\nComparison:")
    print(f"  ω_θ error = {np.linalg.norm(omega_theta_code - omega_theta):.2e}")
    print(f"  ω_φ error = {np.linalg.norm(omega_phi_code - omega_phi):.2e}")
    
    # Check ∂_i e_j = ω_i × e_j using analytical ω
    def commutator_with_vec(biv, v):
        """B × v = (b23, b31, b12) cross (v1, v2, v3)"""
        return np.cross(biv, v)
    
    print(f"\nVerify ∂ᵢeⱼ = ωᵢ × eⱼ (analytical ω):")
    print(f"  ∂_θ e_θ analytical = {d_theta_e_theta}")
    print(f"  ω_θ × e_θ          = {commutator_with_vec(omega_theta, e_theta)}")
    print(f"  Match: {np.allclose(d_theta_e_theta, commutator_with_vec(omega_theta, e_theta), atol=1e-10)}")
    
    print(f"\n  ∂_θ e_φ analytical = {d_theta_e_phi}")
    print(f"  ω_θ × e_φ          = {commutator_with_vec(omega_theta, e_phi)}")
    print(f"  Match: {np.allclose(d_theta_e_phi, commutator_with_vec(omega_theta, e_phi), atol=1e-10)}")
    
    print(f"\n  ∂_φ e_θ analytical = {d_phi_e_theta}")
    print(f"  ω_φ × e_θ          = {commutator_with_vec(omega_phi, e_theta)}")
    print(f"  Match: {np.allclose(d_phi_e_theta, commutator_with_vec(omega_phi, e_theta), atol=1e-10)}")
    
    print(f"\n  ∂_φ e_φ analytical = {d_phi_e_phi}")
    print(f"  ω_φ × e_φ          = {commutator_with_vec(omega_phi, e_phi)}")
    print(f"  Match: {np.allclose(d_phi_e_phi, commutator_with_vec(omega_phi, e_phi), atol=1e-10)}")


if __name__ == "__main__":
    detailed_sphere_test()
