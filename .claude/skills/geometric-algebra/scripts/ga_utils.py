#!/usr/bin/env python3
"""
Geometric Algebra utility functions for common operations.
Works with kingdon library.

Usage:
    from ga_utils import *
    
    alg = create_pga3d()
    A, B = point3d(1,0,0), point3d(0,1,0)
    line = join(A, B)
"""

import math
from typing import Tuple, Optional

try:
    from kingdon import Algebra
except ImportError:
    raise ImportError("Please install kingdon: pip install kingdon")


# ============================================================================
# Algebra Factory Functions
# ============================================================================

def create_complex():
    """Complex numbers: Cl(0,1)"""
    return Algebra(0, 1)

def create_dual():
    """Dual numbers: Cl(0,0,1) - for automatic differentiation"""
    return Algebra(0, 0, 1)

def create_quaternions():
    """Quaternions: Cl(0,2)"""
    return Algebra(0, 2)

def create_cl3():
    """3D Euclidean Clifford algebra: Cl(3,0)"""
    return Algebra(3, 0, 0)

def create_pga2d():
    """2D Projective GA: Cl(2,0,1)"""
    return Algebra(2, 0, 1)

def create_pga3d():
    """3D Projective GA: Cl(3,0,1)"""
    return Algebra(3, 0, 1)

def create_cga2d():
    """2D Conformal GA: Cl(3,1)"""
    return Algebra(3, 1, 0)

def create_cga3d():
    """3D Conformal GA: Cl(4,1)"""
    return Algebra(4, 1, 0)

def create_sta():
    """Spacetime Algebra: Cl(1,3)"""
    return Algebra(1, 3, 0)


# ============================================================================
# PGA 2D Primitives and Operations
# ============================================================================

def pga2d_primitives(alg):
    """
    Get PGA 2D primitive constructors.
    
    Returns dict with: point, line, origin, ideal_line, EX, EY
    """
    locals().update(alg.blades)
    
    origin = 1*e12
    EX = -1*e02  # Direction of x-axis
    EY = 1*e01   # Direction of y-axis
    ideal_line = 1*e0  # Line at infinity
    
    def point(x, y):
        """Create a point at (x, y)"""
        return origin + x*EX + y*EY
    
    def line(a, b, c):
        """Create line ax + by + c = 0"""
        return a*e1 + b*e2 + c*e0
    
    return {
        'point': point,
        'line': line,
        'origin': origin,
        'ideal_line': ideal_line,
        'EX': EX,
        'EY': EY
    }


# ============================================================================
# PGA 3D Primitives and Operations
# ============================================================================

def pga3d_primitives(alg):
    """
    Get PGA 3D primitive constructors.
    
    Returns dict with: point, plane, origin, ideal_plane, EX, EY, EZ
    """
    locals().update(alg.blades)
    
    origin = 1*e123
    EX = -1*e023  # Direction of x-axis
    EY = 1*e013   # Direction of y-axis
    EZ = -1*e012  # Direction of z-axis
    ideal_plane = 1*e0  # Plane at infinity
    
    def point(x, y, z):
        """Create a point at (x, y, z)"""
        return origin + x*EX + y*EY + z*EZ
    
    def plane(a, b, c, d):
        """Create plane ax + by + cz + d = 0"""
        return a*e1 + b*e2 + c*e3 + d*e0
    
    def direction(dx, dy, dz):
        """Create ideal point (direction)"""
        return dx*EX + dy*EY + dz*EZ
    
    return {
        'point': point,
        'plane': plane,
        'direction': direction,
        'origin': origin,
        'ideal_plane': ideal_plane,
        'EX': EX,
        'EY': EY,
        'EZ': EZ
    }


# ============================================================================
# Transformation Factories
# ============================================================================

def make_rotor(angle: float, bivector):
    """
    Create a rotor for rotation by angle around bivector plane.
    
    Args:
        angle: Rotation angle in radians
        bivector: The plane of rotation (will be normalized)
    
    Returns:
        Rotor R such that R*x*~R rotates x
    """
    B = bivector.normalized()
    return math.cos(angle/2) + math.sin(angle/2) * B


def make_translator_2d(alg, dx: float, dy: float):
    """
    Create a translator for PGA 2D.
    
    Args:
        alg: The PGA2D algebra
        dx, dy: Translation amounts
    """
    locals().update(alg.blades)
    return 1 + 0.5*(dx*e01 + dy*e02)


def make_translator_3d(alg, dx: float, dy: float, dz: float):
    """
    Create a translator for PGA 3D.
    
    Args:
        alg: The PGA3D algebra  
        dx, dy, dz: Translation amounts
    """
    locals().update(alg.blades)
    return 1 + 0.5*(dx*e01 + dy*e02 + dz*e03)


# ============================================================================
# CGA Utilities
# ============================================================================

def cga3d_primitives(alg):
    """
    Get CGA 3D primitive constructors.
    
    Returns dict with embedding functions for points, spheres, circles, etc.
    """
    locals().update(alg.blades)
    
    # Null basis vectors
    ninf = e4 + e5      # Point at infinity
    no = 0.5*(e5 - e4)  # Origin
    
    def embed_point(x, y, z):
        """Embed Euclidean point into CGA"""
        p = x*e1 + y*e2 + z*e3
        return p + 0.5*(x*x + y*y + z*z)*ninf + no
    
    def sphere(cx, cy, cz, r):
        """Create sphere with center (cx,cy,cz) and radius r"""
        c = embed_point(cx, cy, cz)
        return c - 0.5*r*r*ninf
    
    def plane_from_normal(nx, ny, nz, d):
        """Create plane with normal (nx,ny,nz) and distance d from origin"""
        n = nx*e1 + ny*e2 + nz*e3
        return n + d*ninf
    
    def extract_point(P):
        """Extract Euclidean coordinates from CGA point"""
        # Normalize so coefficient of no is 1
        scale = -(P | ninf).s  # Should give -1 for normalized point
        if abs(scale) < 1e-10:
            return None  # Point at infinity
        Pn = P / scale
        return (Pn | e1).s, (Pn | e2).s, (Pn | e3).s
    
    return {
        'embed_point': embed_point,
        'sphere': sphere,
        'plane_from_normal': plane_from_normal,
        'extract_point': extract_point,
        'ninf': ninf,
        'no': no
    }


# ============================================================================
# Utility Operations
# ============================================================================

def join(A, B):
    """Join operation (vee product) - creates element spanning A and B"""
    return A & B


def meet(A, B):
    """Meet operation (wedge product) - intersection of A and B"""
    return A ^ B


def sandwich(V, X):
    """Sandwich product: V*X*~V"""
    return V >> X  # kingdon's sandwich operator


def project(A, B):
    """Project A onto B"""
    return (A << B) * B


def reject(A, B):
    """Reject A from B (orthogonal component)"""
    return (A ^ B) * B


def angle_between(A, B):
    """
    Angle between two elements (vectors or blades of same grade).
    Returns angle in radians.
    """
    An = A.normalized()
    Bn = B.normalized()
    cos_angle = (An | Bn).s
    return math.acos(max(-1, min(1, cos_angle)))


def distance_points_pga(P1, P2):
    """Distance between two PGA points (must be normalized)"""
    return (P1 & P2).norm()


# ============================================================================
# Interpolation
# ============================================================================

def slerp(R1, R2, t: float):
    """
    Spherical linear interpolation between rotors.
    
    Args:
        R1, R2: Rotors to interpolate between
        t: Parameter in [0, 1]
    
    Returns:
        Interpolated rotor at parameter t
    """
    R12 = ~R1 * R2  # R1^{-1} * R2
    return R1 * (R12 ** t)


def lerp_motors(M1, M2, t: float):
    """
    Linear interpolation of motors (approximate, for small differences).
    For large differences, use logarithmic interpolation.
    """
    return M1 * (1 - t) + M2 * t


# ============================================================================
# Example usage
# ============================================================================

if __name__ == "__main__":
    # Demo: PGA 3D basics
    print("Creating PGA 3D algebra...")
    alg = create_pga3d()
    prims = pga3d_primitives(alg)
    
    point = prims['point']
    plane = prims['plane']
    
    # Create some points
    A = point(1, 0, 0)
    B = point(0, 1, 0)
    C = point(0, 0, 1)
    
    print(f"Point A: {A}")
    print(f"Point B: {B}")
    print(f"Point C: {C}")
    
    # Line through A and B
    line_AB = join(A, B)
    print(f"Line AB: {line_AB}")
    
    # Plane through A, B, C
    plane_ABC = join(join(A, B), C)
    print(f"Plane ABC: {plane_ABC}")
    
    # Create a translator
    T = make_translator_3d(alg, 1, 2, 3)
    print(f"Translator: {T}")
    
    # Translate point A
    A_translated = sandwich(T, A)
    print(f"A translated: {A_translated}")
