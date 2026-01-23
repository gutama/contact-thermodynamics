# Geometric Algebra Reference: Algebras and Signatures

## Signature Notation

`Cl(p,q,r)` or `Algebra(p,q,r)` where:
- `p` = positive signature dimensions (e²=+1)
- `q` = negative signature dimensions (e²=-1)
- `r` = null/degenerate dimensions (e²=0)

Total dimension of algebra = 2^(p+q+r)

## Fundamental Algebras

### Real Numbers ℝ
- Signature: Cl(0,0,0)
- Basis: {1}
- Elements: 1

### Complex Numbers ℂ
- Signature: Cl(0,1,0)
- Basis: {1, e1} where e1²=-1
- Identified: e1 ↔ i
- Elements: 2

### Dual Numbers
- Signature: Cl(0,0,1)
- Basis: {1, ε} where ε²=0
- Use: Automatic differentiation
- Elements: 2
- Property: (a + bε)*(c + dε) = ac + (ad+bc)ε

### Hyperbolic Numbers (Split-Complex)
- Signature: Cl(1,0,0)
- Basis: {1, j} where j²=+1
- Elements: 2
- Contains zero divisors: (1+j)(1-j)=0

### Quaternions ℍ
- Signature: Cl(0,2,0)
- Basis: {1, i, j, k} where i²=j²=k²=ijk=-1
- Elements: 4
- Use: 3D rotations (double cover of SO(3))
- Note: Non-commutative

### Pauli Algebra / Cl(3)
- Signature: Cl(3,0,0)
- Basis: {1, e1, e2, e3, e12, e13, e23, e123}
- Elements: 8
- Even subalgebra: Quaternions
- Pseudoscalar: I = e123, I²=-1
- Use: 3D Euclidean geometry, quantum spin

## Projective Geometric Algebras (PGA)

### PGA 2D - Plane Geometry
- Signature: Cl(2,0,1) or metric [0,1,1]
- Basis: {1, e0, e1, e2, e01, e02, e12, e012}
- Elements: 8
- Vectors (grade 1): Lines (ax+by+c=0)
- Bivectors (grade 2): Points (homogeneous)
- Pseudoscalar: e012 (degenerate, e012²=0)

**Point representation:**
```
P = w*e12 + x*e02 - y*e01
  = w*e12 + x*e20 + y*e01  (alternate sign convention)
Normalized: w=1 → (x,y) Euclidean coordinates
```

**Line representation:**
```
L = a*e1 + b*e2 + c*e0
Represents line: ax + by + c = 0
```

### PGA 3D - Space Geometry
- Signature: Cl(3,0,1) or metric [0,1,1,1]
- Basis: 16 elements
- Vectors: Planes
- Bivectors: Lines (6D: 3 moment + 3 direction)
- Trivectors: Points
- Pseudoscalar: e0123 (degenerate)

**Point representation:**
```
P = w*e123 + x*(-e023) + y*(e013) + z*(-e012)
   = w*e123 + x*e032 + y*e013 + z*e021
Normalized: w=1 → (x,y,z) Euclidean coordinates
```

**Plane representation:**
```
π = a*e1 + b*e2 + c*e3 + d*e0
Normal vector: (a,b,c), distance from origin: -d/√(a²+b²+c²)
```

**Line representation (Plücker coordinates):**
```
L = d*e23 + e*e31 + f*e12 + a*e01 + b*e02 + c*e03
Direction: (d,e,f), Moment: (a,b,c)
```

## Conformal Geometric Algebras (CGA)

### CGA 2D
- Signature: Cl(3,1,0)
- Basis: 16 elements
- Extra dimensions: n∞ (infinity), no (origin)
- Use: Circles, inversions, Möbius transformations

### CGA 3D
- Signature: Cl(4,1,0)
- Basis: 32 elements
- Basis vectors: e1, e2, e3, e+, e-
- Null vectors: n∞ = e+ + e-, no = ½(e- - e+)

**Point embedding:**
```
P = x*e1 + y*e2 + z*e3 + ½(x²+y²+z²)*n∞ + no
```

**Round objects:**
```
Sphere: S = c - ½r²*n∞  (center c, radius r)
Dual sphere: S* = P1^P2^P3^P4 (through 4 points)
Circle: C* = P1^P2^P3
Line: L* = P1^P2^n∞
Plane: π* = P1^P2^P3^n∞
Point pair: PP* = P1^P2
```

## Spacetime Algebras

### Space-Time Algebra (STA)
- Signature: Cl(1,3,0)
- Basis: 16 elements
- Use: Special relativity, electromagnetism

**Basis:**
```
γ0² = +1 (timelike)
γ1² = γ2² = γ3² = -1 (spacelike)
γ0γ1γ2γ3 = I (pseudoscalar)
```

**Electromagnetic field:**
```
F = E + IB  (bivector)
E = E1*γ10 + E2*γ20 + E3*γ30
B = B1*γ23 + B2*γ31 + B3*γ12
Maxwell: ∇F = J/ε₀
```

### Algebra of Physical Space (APS)
- Signature: Cl(3,0,0)
- Relative to observer γ0
- Use: 3D physics relative to frame

## Higher-Dimensional Algebras

### Double Conformal (DCGA)
- Signature: Cl(6,2,0)
- Elements: 256
- Use: Quadric surfaces (ellipsoids, paraboloids, hyperboloids)

### Triple Conformal (TCGA)
- Signature: Cl(9,3,0)
- Elements: 4096
- Use: Cubic surfaces

### Quadric CGA (QCGA)
- Signature: Cl(9,6,0)
- Use: General quadric surfaces

## Algebra Properties

### Pseudoscalar Properties
| Algebra | Pseudoscalar | Square |
|---------|--------------|--------|
| Cl(2) | e12 | -1 |
| Cl(3) | e123 | -1 |
| PGA2D | e012 | 0 |
| PGA3D | e0123 | 0 |
| CGA3D | e12345 | +1 |
| STA | γ0123 | -1 |

### Even Subalgebras
The even subalgebra (elements of even grade) is always closed under geometric product:
- Cl(2) even → Complex numbers
- Cl(3) even → Quaternions
- PGA3D even → Dual quaternions (motors)
- CGA3D even → Versors for conformal transformations
