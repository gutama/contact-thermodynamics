"""
Test: Does ω × e give TANGENTIAL projection of ∂e?

For embedded surfaces, ∂ᵢeⱼ has both:
1. Tangential component (in span of e_1, e_2)
2. Normal component (along n)

The connection bivector only captures the tangential rotation.
"""

import numpy as np
from riemannian_ga import sphere_frame, ConnectionBivector

def test_tangential_projection():
    print("=" * 70)
    print("TEST: Connection only captures TANGENTIAL component of ∂ᵢeⱼ")
    print("=" * 70)
    
    R = 1.0
    theta = np.pi / 3
    phi = 0.0
    coords = np.array([theta, phi])
    
    ct, st = np.cos(theta), np.sin(theta)
    cp, sp = np.cos(phi), np.sin(phi)
    
    # Frame vectors
    e_theta = R * np.array([ct * cp, ct * sp, -st])
    e_phi = R * st * np.array([-sp, cp, 0.0])
    
    # Normal (outward)
    normal = np.array([st * cp, st * sp, ct])
    
    # Analytical derivatives
    d_theta_e_theta = R * np.array([-st * cp, -st * sp, -ct])
    d_theta_e_phi = R * ct * np.array([-sp, cp, 0.0])
    d_phi_e_theta = R * np.array([-ct * sp, ct * cp, 0.0])
    d_phi_e_phi = R * st * np.array([-cp, -sp, 0.0])
    
    # Decompose each derivative into tangent + normal parts
    print(f"\nθ = {np.degrees(theta):.1f}°")
    print(f"Normal n = {normal}")
    
    derivatives = [
        ("∂_θ e_θ", d_theta_e_theta),
        ("∂_θ e_φ", d_theta_e_phi),
        ("∂_φ e_θ", d_phi_e_theta),
        ("∂_φ e_φ", d_phi_e_phi),
    ]
    
    # Get connection from code
    frame_func = sphere_frame(R)
    conn = ConnectionBivector(manifold_dim=2, ambient_dim=3, h=1e-7)
    omega = conn.compute_at(coords, frame_func)
    
    omega_theta = omega[0].to_array()
    omega_phi = omega[1].to_array()
    
    omegas = [omega_theta, omega_theta, omega_phi, omega_phi]
    e_js = [e_theta, e_phi, e_theta, e_phi]
    
    print("\n" + "-" * 70)
    
    for i, ((name, deriv), om, ej) in enumerate(zip(derivatives, omegas, e_js)):
        # Normal component
        normal_comp = np.dot(deriv, normal) * normal
        
        # Tangential component
        tangent_comp = deriv - normal_comp
        
        # What ω × e gives
        omega_cross_e = np.cross(om, ej)
        
        print(f"\n{name}:")
        print(f"  Full derivative: {deriv}")
        print(f"  Normal component: {normal_comp}")
        print(f"  Tangent component: {tangent_comp}")
        print(f"  ω × e:            {omega_cross_e}")
        print(f"  Match tangent? {np.allclose(tangent_comp, omega_cross_e, atol=1e-6)}")
        
        if not np.allclose(tangent_comp, omega_cross_e, atol=1e-6):
            print(f"  ERROR: {np.linalg.norm(tangent_comp - omega_cross_e):.6f}")


if __name__ == "__main__":
    test_tangential_projection()
