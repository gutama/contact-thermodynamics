# Geometric Algebra Formula Reference

## Fundamental Products

### Geometric Product
The geometric product of vectors a and b:
```
ab = a·b + a∧b
```
- Scalar part: `a·b` (inner/dot product)
- Bivector part: `a∧b` (outer/wedge product)

For unit vectors: `ab = cos(θ) + sin(θ)*B` where B is unit bivector of plane

### Inner Product (Dot/Contraction)
Left contraction `a⌋B` contracts vector a from left of B:
```
a⌋B = ½(aB - B̂a)  where B̂ = grade involution of B
```
For vector and bivector: `a⌋B = ½(aB - Ba)` (since bivectors are even)

Grade formula: `grade(A_r ⌋ B_s) = s - r` if s≥r, else 0

### Outer Product (Wedge)
```
a∧b = ½(ab - ba)        # For vectors
A∧B = <AB>_{r+s}        # General: highest grade part
```
Properties:
- Antisymmetric: `a∧b = -b∧a`
- Associative: `(a∧b)∧c = a∧(b∧c)`
- `a∧a = 0`

### Regressive Product (Vee)
```
A∨B = I⁻¹(I·A ∧ I·B) = (A* ∧ B*)*
```
Where `A* = AI⁻¹` is the dual.

In code: `A & B` or `vee(A,B)`

Use: Intersections/meets in dual representations

### Commutator Product
```
A×B = ½(AB - BA)
```
Use: Bivector acting on vector, Lie algebra structure

### Scalar Product
```
A*B = <AB>_0 = <BA>_0
```
Only the scalar (grade-0) part.

## Grade Operations

### Grade Projection
```
<M>_k = grade-k part of multivector M
```

### Grade Involution (Main Involution)
```
M̂ = Σ (-1)^k <M>_k
```
Reverses sign of odd-grade parts.

### Reversion
```
M̃ = Σ (-1)^(k(k-1)/2) <M>_k
```
Reverses order in products: `(ab)~ = ba`

Grade signs for reversion:
- Grade 0,1: +
- Grade 2,3: -
- Grade 4,5: +
- Grade 6,7: -

### Clifford Conjugate
```
M̄ = M̃̂ = M̂̃
```
Combination of reversion and grade involution.

## Inverses and Norms

### Versor Inverse
For versors (products of vectors):
```
V⁻¹ = Ṽ / (VṼ)
```

### General Inverse
Not all multivectors are invertible. When they are:
```
M⁻¹ exists iff MM̃ is a non-zero scalar
M⁻¹ = M̃ / (MM̃)
```

### Norm
```
|M|² = <MM̃>_0
|M| = √(<MM̃>_0)
```

Normalized: `M̂ = M/|M|`

## Exponentials and Logarithms

### Bivector Exponential (Rotor)
For unit bivector B (B²=-1):
```
exp(θB) = cos(θ) + sin(θ)B
```

For general bivector:
```
exp(B) = exp(θB̂) where θ=|B|, B̂=B/|B|
```

### Null Bivector Exponential (Translator)
For bivector B with B²=0:
```
exp(B) = 1 + B
```

### Motor (Combined Rotation+Translation in PGA)
```
M = exp(½(θL_rot + d*e0∧L_trans))
```

### Logarithm
For rotor R = cos(θ) + sin(θ)B:
```
log(R) = θB where θ = atan2(|<R>_2|, <R>_0)
```

## Projections and Reflections

### Orthogonal Projection
Project A onto B:
```
proj_B(A) = (A⌋B)B⁻¹ = (A⌋B)/B
```

Alternate form:
```
proj_B(A) = (A·B)B⁻¹
```

### Rejection (Orthogonal Component)
```
rej_B(A) = A - proj_B(A) = (A∧B)B⁻¹
```

### Reflection
Reflect A in hyperplane with normal n:
```
A' = -nAn⁻¹ = -nAn/n²
```

For unit normal: `A' = -nAn`

### Rotation
Rotate A in plane B by angle θ:
```
R = exp(-θB/2) = cos(θ/2) - sin(θ/2)B
A' = RAR̃
```

Double-angle: rotating by θ uses half-angle in rotor.

## PGA-Specific Formulas

### Distance Between Points (PGA)
```
d(P,Q) = |P∨Q| / (|P||Q|)
```
For normalized points: `d(P,Q) = |P∨Q|`

### Angle Between Lines (2D PGA)
```
cos(θ) = (L1·L2) / (|L1||L2|)
```

### Distance Point to Line (2D PGA)
```
d(P,L) = |P∧L| / |L|
```

### Motor Decomposition
Motor M = RT where R is rotor, T is translator:
```
Rotation part: R = <M>_even normalized
Translation part: T = MR̃
```

## CGA-Specific Formulas

### Embedding Euclidean Point
```
P = x + ½x²n∞ + no
```
where x is Euclidean position, n∞ is infinity, no is origin.

### Extracting Euclidean Point
```
x = P / (-P·n∞)
```

### Sphere (Direct Representation)
Center c, radius r:
```
S = c - ½r²n∞
```

### Distance Between Points
```
d²(P,Q) = -2(P·Q)
```
For normalized points (P·n∞ = -1).

### Circle from 3 Points
```
C = P1∧P2∧P3
```

### Intersection
```
S1 ∩ S2 = (S1∧S2)*
```
Dual of wedge gives intersection.

## Identities

### Jacobi Identity
```
[A,[B,C]] + [B,[C,A]] + [C,[A,B]] = 0
```

### Reversion of Products
```
(AB)~ = B̃Ã
```

### Contraction Identities
```
a⌋(B∧C) = (a⌋B)∧C + B̂∧(a⌋C)
(A∧B)⌋C = A⌋(B⌋C)
```

### Dual Identities
```
(A∧B)* = A*⌋B = A⌋B*
(A⌋B)* = A*∧B
A** = (-1)^(n(n-1)/2) A  (n = space dimension)
```

## Rotor-Quaternion Equivalence

The even subalgebra of Cl(3) is isomorphic to the quaternions ℍ.

### Mapping Table

| Rotor Component | Quaternion | Bivector Plane | Rotation Axis |
|-----------------|------------|----------------|---------------|
| s (scalar)      | w          | —              | —             |
| e₂₃             | i          | YZ plane       | X axis        |
| e₃₁             | j          | ZX plane       | Y axis        |
| e₁₂             | k          | XY plane       | Z axis        |

### Rotor from Axis-Angle

```
R = cos(θ/2) + sin(θ/2) · B̂

where B̂ = normalized bivector for rotation plane
      = (ax·e₂₃ + ay·e₃₁ + az·e₁₂) / |axis|
```

### Geometric Product of Rotors

```
R₁ · R₂ = (s₁ + B₁)(s₂ + B₂)
        = s₁s₂ - B₁·B₂ + s₁B₂ + s₂B₁ + B₁×B₂

Expanded:
s   = s₁s₂ - (e12₁·e12₂ + e31₁·e31₂ + e23₁·e23₂)
e12 = s₁·e12₂ + e12₁·s₂ + e31₁·e23₂ - e23₁·e31₂
e31 = s₁·e31₂ + e31₁·s₂ - e12₁·e23₂ + e23₁·e12₂
e23 = s₁·e23₂ + e23₁·s₂ + e12₁·e31₂ - e31₁·e12₂
```

### Efficient Sandwich Product

The sandwich product `R·v·R̃` can be computed efficiently using cross products:

```javascript
// Quaternion-equivalent formula for R·v·R̃
function apply(x, y, z) {
    const w = this.s;
    const i = this.e23;
    const j = -this.e31;  // Note sign!
    const k = this.e12;
    
    // First cross: q × v
    const cx = j*z - k*y;
    const cy = k*x - i*z;
    const cz = i*y - j*x;
    
    // Second cross: q × (q × v)
    const ccx = j*cz - k*cy;
    const ccy = k*cx - i*cz;
    const ccz = i*cy - j*cx;
    
    // Result: v + 2w(q×v) + 2(q×(q×v))
    return [
        x + 2*(w*cx + ccx),
        y + 2*(w*cy + ccy),
        z + 2*(w*cz + ccz)
    ];
}
```

### Sign Convention Warning

The mapping between bivector components and quaternion imaginary units
depends on handedness convention. The formulas above assume:

- Right-handed coordinate system
- Rotation direction follows right-hand rule
- e₁₂ corresponds to rotation in XY plane (around Z axis)

### SLERP for Rotors

Spherical linear interpolation:

```
slerp(R₁, R₂, t) = R₁ · (R₁⁻¹ · R₂)^t

Expanded:
θ = acos(R₁ · R₂)  // dot product of components
a = sin((1-t)θ) / sin(θ)
b = sin(tθ) / sin(θ)
R = a·R₁ + b·R₂
```

## Outermorphisms

Linear map f extends to multivectors:
```
f(a∧b∧...∧c) = f(a)∧f(b)∧...∧f(c)
```

Determinant:
```
f(I) = det(f) · I
```
