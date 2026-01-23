# Geometric Algebra in Physics

## Spacetime Algebra (STA)

### Setup
Signature: Cl(1,3) with basis {γ0, γ1, γ2, γ3}
```
γ0² = +1     (timelike)
γi² = -1     (spacelike, i=1,2,3)
γμγν = -γνγμ  for μ≠ν

Pseudoscalar: I = γ0γ1γ2γ3
I² = -1
```

### Relative Vectors (Space in γ0 frame)
```
σi = γiγ0    (Pauli matrices / relative space vectors)
σ1² = σ2² = σ3² = +1
σ1σ2σ3 = I
```

### 4-Velocity and Proper Velocity
```
v = γ(c + u)    where γ = 1/√(1-u²/c²)
Proper velocity: u = vγ0
```

### Spacetime Split
Any bivector F splits into electric + magnetic:
```
F = E + IB
E = ½(F - γ0Fγ0)   # Electric (timelike bivector part)
IB = ½(F + γ0Fγ0)  # Magnetic (spacelike bivector part)
```

## Electromagnetism

### Electromagnetic Field Bivector
```
F = E + IcB = Σ(Eiσi + cBiIσi)

Explicit:
F = E1γ10 + E2γ20 + E3γ30 + cB1γ23 + cB2γ31 + cB3γ12
```

### Maxwell's Equations (Unified)
```
∇F = J/ε0c

where ∇ = γμ∂μ = γ0(∂t/c) + γi∂i
      J = cρ + j
```

This single equation encodes all four Maxwell equations:
- ∇·E = ρ/ε0
- ∇×E = -∂B/∂t
- ∇·B = 0
- c²∇×B = j/ε0 + ∂E/∂t

### Lorentz Force
```
F = q(E + v×B)   →   dp/dτ = qF·v/c
```

### Wave Equation
```
∇²F = 0   (in vacuum)
```

### Poynting Vector
```
S = (E × B)/μ0 = c²ε0(E × B)
```
Using GA: `S = -cε0<FγoF̃>_1/2`

## Quantum Mechanics

### Pauli Equation
Spinor ψ in Cl(3):
```
σ·∇ψIσ3 = mψγ0
```

### Dirac Equation
In STA:
```
∇ψIσ3 - eAψ = mψγ0

ψ = √ρ R exp(Iβ/2)
```
where R is a rotor and ρ, β are scalar/pseudoscalar.

### Spin Bivector
```
S = (ℏ/2)ψσ3ψ̃
```

## Classical Mechanics

### Angular Momentum
```
L = r∧p = r∧mv
```
(Bivector, not pseudovector)

### Torque
```
τ = r∧F
```

### Inertia
```
L = I(ω)    # I is linear map from bivector to bivector
```

### Euler's Equations (Rigid Body)
```
İω + ω×(Iω) = τ
```
In GA with bivectors directly.

## Rotations and Angular Velocity

### Angular Velocity Bivector
```
Ω = ω1e23 + ω2e31 + ω3e12
```
(Components match traditional ω vector)

### Rotation Matrix from Rotor
For rotor R = cos(θ/2) + sin(θ/2)B:
```
v' = RvR̃
```
Each column of rotation matrix = R*ei*R̃

### Composition
```
R_total = R2 R1   (apply R1 first, then R2)
```

## Relativity

### Lorentz Transformation
Boost in direction v with rapidity φ = artanh(v/c):
```
B = exp(φn̂/2)   where n̂ = v̂∧γ0 (unit spacelike bivector)
x' = BxB̃
```

### Rotation
```
R = exp(-θB/2)   where B is spatial bivector
```

### General Lorentz Transformation
```
Λ = BR   (boost followed by rotation)
```

### Thomas Precession
Accumulated rotation from sequence of non-collinear boosts.

## Projective Mechanics

### Kinematics in PGA
Position, velocity, acceleration as points:
```
P(t) = point at time t
V = Ṗ   (velocity point/direction)
A = P̈   (acceleration)
```

### Interpolation
```
P(t) = exp(tB)P0exp(-tB)   where B encodes velocity
```

### Screw Motion
Motor M describes combined rotation + translation:
```
P(t) = M(t)P0M̃(t)
M(t) = exp(tL/2)   where L is line (bivector in PGA)
```

## Geometric Calculus

### Vector Derivative
```
∇ = ei∂i   (summed over basis)

∇F = ∇·F + ∇∧F
```

### Gradient, Divergence, Curl
```
grad(f) = ∇f           # Vector
div(F) = ∇·F           # Scalar  
curl(F) = -I(∇∧F)      # Vector (from bivector)
```

### Laplacian
```
∇²f = ∇·∇f
```

### Fundamental Theorem
```
∫_M dⁿx (∂F) = ∮_∂M dⁿ⁻¹x F
```
Generalizes Stokes, Gauss, Green's theorems.

## Code Examples

### Electromagnetic Field (kingdon)
```python
from kingdon import Algebra
import numpy as np

# Spacetime Algebra
sta = Algebra(1, 3)
locals().update(sta.blades)

# Electromagnetic field at a point
def em_field(E, B):
    """E, B are 3-vectors as lists"""
    return (E[0]*e01 + E[1]*e02 + E[2]*e03 + 
            B[0]*e23 + B[1]*e31 + B[2]*e12)

# Example: uniform E-field
F = em_field([1,0,0], [0,0,0])
```

### Rigid Body Rotation (galgebra)
```python
from galgebra.ga import Ga
from sympy import symbols, cos, sin

# 3D Euclidean
g3 = Ga('e', g=[1,1,1])
e1, e2, e3 = g3.mv()

# Angular velocity as bivector
omega = symbols('omega')
Omega = omega * e1^e2  # Rotation in xy-plane

# Rotor for finite rotation
t = symbols('t')
R = cos(omega*t/2) + sin(omega*t/2)*(e1^e2)/(e1^e2).norm()
```
