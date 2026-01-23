"""
Hyperbolic surface test: Pseudosphere (tractricoid)

The pseudosphere is a surface of revolution with constant negative Gaussian curvature K = -1.
It's the only complete surface in R³ with K = -1 (actually, it has a singularity).

Parametrization:
  x(u, v) = sech(u) cos(v)
  y(u, v) = sech(u) sin(v)  
  z(u, v) = u - tanh(u)

where u > 0 and v ∈ [0, 2π)

Gaussian curvature: K = -1 everywhere
"""

import numpy as np
from riemannian_ga import TangentFrame, ShapeOperator
from typing import Callable


def pseudosphere_frame() -> Callable[[np.ndarray], TangentFrame]:
    """
    Create frame function for pseudosphere (tractricoid).
    
    The pseudosphere has constant negative curvature K = -1.
    
    Parametrization:
        x = sech(u) cos(v)
        y = sech(u) sin(v)
        z = u - tanh(u)
    
    where sech(u) = 1/cosh(u)
    """
    def frame_func(coords: np.ndarray) -> TangentFrame:
        u, v = coords[0], coords[1]
        
        # Avoid u = 0 (singularity at the "cusp")
        u = max(u, 0.5)
        
        sechu = 1.0 / np.cosh(u)
        tanhu = np.tanh(u)
        cv, sv = np.cos(v), np.sin(v)
        
        # Derivatives of sech(u) and tanh(u)
        d_sechu = -sechu * tanhu
        d_tanhu = sechu ** 2
        
        # Tangent vectors (derivatives of embedding)
        # e_u = (∂x/∂u, ∂y/∂u, ∂z/∂u)
        e_u = np.array([
            d_sechu * cv,       # ∂x/∂u = -sech(u)tanh(u) cos(v)
            d_sechu * sv,       # ∂y/∂u = -sech(u)tanh(u) sin(v)
            1 - d_tanhu         # ∂z/∂u = 1 - sech²(u)
        ])
        
        # e_v = (∂x/∂v, ∂y/∂v, ∂z/∂v)  
        e_v = np.array([
            -sechu * sv,        # ∂x/∂v = -sech(u) sin(v)
            sechu * cv,         # ∂y/∂v = sech(u) cos(v)
            0.0                 # ∂z/∂v = 0
        ])
        
        vectors = np.array([e_u, e_v])
        return TangentFrame(vectors)
    
    return frame_func


def test_pseudosphere():
    print("=" * 70)
    print("PSEUDOSPHERE (TRACTRICOID) TEST")
    print("Expected: K = -1 (constant negative curvature)")
    print("=" * 70)
    
    frame_func = pseudosphere_frame()
    shape = ShapeOperator(h=1e-5)
    
    test_points = [
        (1.0, 0.0, "u=1, v=0"),
        (1.0, np.pi/4, "u=1, v=45°"),
        (2.0, 0.0, "u=2, v=0"),
        (0.5, np.pi/2, "u=0.5, v=90°"),
    ]
    
    all_pass = True
    
    for u_val, v_val, name in test_points:
        coords = np.array([u_val, v_val])
        
        try:
            S = shape.compute_at(coords, frame_func)
            K = shape.gaussian_curvature(S)
            H = shape.mean_curvature(S)
            
            expected_K = -1.0
            error = abs(K - expected_K)
            status = "✓" if error < 0.1 else "✗"
            all_pass = all_pass and (error < 0.1)
            
            print(f"\n{name}:")
            print(f"  K = {K:+.6f} (expected: {expected_K:.6f}) {status}")
            print(f"  H = {H:+.6f}")
            print(f"  error = {error:.6f}")
            
        except Exception as e:
            print(f"\n{name}: ERROR - {e}")
            all_pass = False
    
    print("\n" + "-" * 70)
    print(f"RESULT: {'✓ ALL PASS' if all_pass else '✗ SOME FAILED'}")
    
    return all_pass


def test_saddle_surface():
    """
    Test a simple saddle surface z = xy (hyperbolic paraboloid).
    
    At origin: K = -1 / (1 + x² + y²)² → K(0,0) = -1
    """
    print("\n" + "=" * 70)
    print("SADDLE SURFACE TEST (z = xy)")
    print("Expected: K < 0 (negative curvature)")
    print("=" * 70)
    
    def saddle_frame(coords):
        x, y = coords[0], coords[1]
        
        # Surface z = xy
        # e_x = (1, 0, ∂z/∂x) = (1, 0, y)
        # e_y = (0, 1, ∂z/∂y) = (0, 1, x)
        
        e_x = np.array([1.0, 0.0, y])
        e_y = np.array([0.0, 1.0, x])
        
        vectors = np.array([e_x, e_y])
        return TangentFrame(vectors)
    
    shape = ShapeOperator(h=1e-5)
    
    # At origin
    coords = np.array([0.0, 0.0])
    S = shape.compute_at(coords, saddle_frame)
    K = shape.gaussian_curvature(S)
    H = shape.mean_curvature(S)
    
    expected_K = -1.0  # At origin
    error = abs(K - expected_K)
    status = "✓" if error < 0.1 else "✗"
    
    print(f"\nAt origin (0, 0):")
    print(f"  K = {K:+.6f} (expected: {expected_K:.6f}) {status}")
    print(f"  H = {H:+.6f} (expected: 0)")
    print(f"  error = {error:.6f}")


if __name__ == "__main__":
    test_pseudosphere()
    test_saddle_surface()
