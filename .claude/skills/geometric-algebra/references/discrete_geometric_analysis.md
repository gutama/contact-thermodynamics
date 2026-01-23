# Discrete Geometric Calculus vs GAlgebra's Split Differential Operator

## Executive Summary

Both Eelco Hoogendoorn's **discrete geometric calculus** (in `geometric_calculus`) and GAlgebra's **split differential operator** address the same fundamental problem: **implementing the geometric derivative ∇** in computational settings. However, they operate in fundamentally different domains:

| Aspect | GAlgebra (Split Operator) | Eelco's Discrete Calculus |
|--------|--------------------------|---------------------------|
| Domain | Continuous (symbolic) | Discrete (numeric) |
| Backend | SymPy | NumPy/JAX/OpenCL |
| Grids | Arbitrary manifolds | Cubic/regular grids |
| Primary Use | Symbolic computation | Numerical simulation |
| Dimension | Arbitrary | Typically 2-5D |

## 1. GAlgebra's Split Differential Operator

### Concept

GAlgebra implements the geometric derivative symbolically by splitting the nabla operator:

```
∇ = ∇_G ⊗ ∇_D
```

Where:
- **∇_G** (Geometric part): Tuple of reciprocal basis vectors `(e¹, e², ..., eⁿ)`
- **∇_D** (Differential part): Tuple of partial derivative operators `(∂/∂x₁, ∂/∂x₂, ..., ∂/∂xₙ)`

The geometric derivative acting on a multivector field F is:

```
∇F = Σᵢ eⁱ (∂F/∂xᵢ)
```

### Key Features

1. **Reciprocal Frames**: For non-orthogonal bases, galgebra automatically computes:
   ```
   eⁱ = (-1)^(i-1) (e₁ ∧ ... ∧ ê_i ∧ ... ∧ eₙ) / Eₙ
   ```
   where `eⁱ · eⱼ = δⁱⱼ`

2. **Manifold Support**: Works on curved manifolds with position-dependent metrics

3. **Symbolic Manipulation**: Full SymPy integration for simplification, substitution

### Implementation (from galgebra)

```python
# The gradient is constructed as:
# grad = Σᵢ eⁱ ∂/∂xᵢ
# where eⁱ are reciprocal basis vectors

from galgebra.ga import Ga
from sympy import symbols

coords = symbols('x y z', real=True)
ga = Ga('e', g=[1,1,1], coords=coords)
grad, rgrad = ga.grads()  # left and right gradients

F = ga.mv('F', 'mv', f=True)  # multivector field
dF = grad * F  # geometric derivative
```

## 2. Eelco's Discrete Geometric Calculus

### Concept

The discrete geometric derivative operates on **cubic grids** where:
- **Fields are arrays** indexed by spatial position
- **Different grades live on different elements**:
  - Grade 0 (scalars) → vertices (0-cells)
  - Grade 1 (vectors) → edges (1-cells)  
  - Grade 2 (bivectors) → faces (2-cells)
  - Grade n (pseudoscalar) → volumes (n-cells)

The discrete derivative uses **finite differences**:

```
(∂F/∂x)ᵢ ≈ (Fᵢ₊₁ - Fᵢ) / Δx
```

### Key Innovation: Leapfrog Operators

For spacetime algebras (with mixed signature), Eelco implements **leapfrog** time-stepping:

```
φ(t+Δt) = φ(t-Δt) + 2Δt · ∂φ/∂t
```

This is second-order accurate and naturally preserves the **wave equation structure**.

### Structure (from README analysis)

```
geometric_calculus/
├── discrete/          # Discrete implementations
│   ├── numpy/         # NumPy backend
│   ├── jax/           # JAX backend (autodiff + GPU)
│   └── opencl/        # OpenCL backend (GPU)
└── continuous/        # Continuous (JAX + numga)
```

### The Geometric Equation

The central equation appears to be:

```
∂φ = 0  (Dirac equation form)
```

or for higher-order physics:

```
∂(∂(∂(∂(φ)))) = 0  (biharmonic for elasticity)
```

Where `∂` is the discrete geometric derivative.

## 3. Conceptual Relationship

### Are They the Same Concept?

**Yes, conceptually** - both implement the geometric derivative `∇ = Σᵢ eⁱ ∂ᵢ`

**No, practically** - they differ in:

| Feature | GAlgebra | Eelco's Discrete |
|---------|----------|------------------|
| ∂ᵢ implementation | Symbolic differentiation | Finite differences |
| eⁱ representation | Symbolic multivectors | Implicit in stencil structure |
| Output | Symbolic expressions | Numerical arrays |
| Accuracy | Exact (symbolic) | O(h²) or higher |

### Unified View

The geometric derivative ∇ is the foundational operator. It naturally splits into grade-raising and grade-lowering parts:

```
∇F = ∇·F + ∇∧F
     └─┬─┘   └─┬─┘
     inner   outer
```

The geometric derivative unifies:
- **Gradient**: ∇f (scalar → vector) — purely grade-raising
- **Divergence**: ∇·V = ⟨∇V⟩₀ — grade-lowering part
- **Curl**: ∇∧V = ⟨∇V⟩₂ — grade-raising part
- **Full derivative**: ∇V = divergence + curl (both grades)

## 4. Implementation: Bridging Both Approaches

Below I provide a Python implementation that unifies the concepts.

### 4.1 Discrete Split Operator (Eelco-style with galgebra concepts)

```python
import numpy as np
from typing import Tuple, List, Dict

class DiscreteGA:
    """
    Discrete Geometric Algebra on cubic grids.
    Implements the split differential operator concept from galgebra
    in a discrete numerical setting.
    """
    
    def __init__(self, signature: Tuple[int, int, int], shape: Tuple[int, ...], 
                 spacing: float = 1.0):
        """
        Args:
            signature: (p, q, r) for Cl(p,q,r) - positive, negative, null dimensions
            shape: Grid dimensions
            spacing: Grid spacing h
        """
        self.p, self.q, self.r = signature
        self.dim = self.p + self.q + self.r
        self.shape = shape
        self.h = spacing
        
        # Build basis blade indices
        self.blades = self._build_blades()
        self.n_blades = len(self.blades)
        
        # Build metric
        self.metric = self._build_metric()
        
    def _build_blades(self) -> List[Tuple[int, ...]]:
        """Build all basis blade indices."""
        from itertools import combinations
        blades = [()]  # scalar
        for grade in range(1, self.dim + 1):
            blades.extend(combinations(range(self.dim), grade))
        return blades
    
    def _build_metric(self) -> np.ndarray:
        """Build diagonal metric tensor."""
        diag = [1.0] * self.p + [-1.0] * self.q + [0.0] * self.r
        return np.diag(diag)
    
    def _blade_to_index(self, blade: Tuple[int, ...]) -> int:
        """Convert blade tuple to linear index."""
        return self.blades.index(blade)
    
    def _index_to_blade(self, idx: int) -> Tuple[int, ...]:
        """Convert linear index to blade tuple."""
        return self.blades[idx]
    
    def grade(self, idx: int) -> int:
        """Get grade of blade at index."""
        return len(self.blades[idx])
    
    def zeros(self) -> np.ndarray:
        """Create zero multivector field."""
        return np.zeros((*self.shape, self.n_blades))
    
    def scalar_field(self, values: np.ndarray) -> np.ndarray:
        """Create multivector field from scalar values."""
        field = self.zeros()
        field[..., 0] = values
        return field
    
    # ===========================================
    # SPLIT DIFFERENTIAL OPERATOR IMPLEMENTATION
    # ===========================================
    
    def partial_derivative(self, field: np.ndarray, direction: int,
                          stencil: str = 'central') -> np.ndarray:
        """
        Compute partial derivative ∂_i F.
        This is the ∇_D (differential) part of the split operator.
        
        Args:
            field: Multivector field array (..., n_blades)
            direction: Spatial direction (0, 1, ..., dim-1)
            stencil: 'forward', 'backward', or 'central'
        
        Returns:
            Partial derivative field
        """
        result = np.zeros_like(field)
        
        if stencil == 'forward':
            # (f[i+1] - f[i]) / h
            slices_plus = [slice(None)] * (len(self.shape) + 1)
            slices_here = [slice(None)] * (len(self.shape) + 1)
            slices_plus[direction] = slice(1, None)
            slices_here[direction] = slice(None, -1)
            
            result_slice = [slice(None)] * (len(self.shape) + 1)
            result_slice[direction] = slice(None, -1)
            
            result[tuple(result_slice)] = (
                field[tuple(slices_plus)] - field[tuple(slices_here)]
            ) / self.h
            
        elif stencil == 'backward':
            # (f[i] - f[i-1]) / h
            slices_here = [slice(None)] * (len(self.shape) + 1)
            slices_minus = [slice(None)] * (len(self.shape) + 1)
            slices_here[direction] = slice(1, None)
            slices_minus[direction] = slice(None, -1)
            
            result_slice = [slice(None)] * (len(self.shape) + 1)
            result_slice[direction] = slice(1, None)
            
            result[tuple(result_slice)] = (
                field[tuple(slices_here)] - field[tuple(slices_minus)]
            ) / self.h
            
        else:  # central
            # (f[i+1] - f[i-1]) / (2h)
            slices_plus = [slice(None)] * (len(self.shape) + 1)
            slices_minus = [slice(None)] * (len(self.shape) + 1)
            slices_plus[direction] = slice(2, None)
            slices_minus[direction] = slice(None, -2)
            
            result_slice = [slice(None)] * (len(self.shape) + 1)
            result_slice[direction] = slice(1, -1)
            
            result[tuple(result_slice)] = (
                field[tuple(slices_plus)] - field[tuple(slices_minus)]
            ) / (2 * self.h)
        
        return result
    
    def geometric_product_basis(self, e_i: int, blade: Tuple[int, ...]) -> Tuple[float, Tuple[int, ...]]:
        """
        Compute e_i * blade where e_i is the i-th basis vector.
        Returns (sign, result_blade).
        
        This is the ∇_G (geometric) part of the split operator.
        """
        if e_i in blade:
            # Contract: e_i is already in the blade
            pos = blade.index(e_i)
            sign = (-1) ** pos * self.metric[e_i, e_i]
            result = tuple(b for b in blade if b != e_i)
            return (sign, result)
        else:
            # Extend: add e_i to the blade
            # Find position to insert while maintaining order
            result = list(blade)
            pos = 0
            for i, b in enumerate(blade):
                if b > e_i:
                    pos = i
                    break
                pos = i + 1
            sign = (-1) ** pos
            result.insert(pos, e_i)
            return (sign, tuple(result))
    
    def apply_basis_vector(self, field: np.ndarray, direction: int) -> np.ndarray:
        """
        Apply basis vector e_direction to multivector field via geometric product.
        This computes e^i * F (left multiplication by reciprocal basis).
        """
        result = np.zeros_like(field)
        
        for idx, blade in enumerate(self.blades):
            sign, new_blade = self.geometric_product_basis(direction, blade)
            new_idx = self._blade_to_index(new_blade)
            # e^i = metric[i,i] * e_i for orthonormal basis
            metric_factor = self.metric[direction, direction]
            if metric_factor != 0:
                result[..., new_idx] += sign * metric_factor * field[..., idx]
        
        return result
    
    def geometric_derivative(self, field: np.ndarray, 
                            stencil: str = 'central') -> np.ndarray:
        """
        Compute the full geometric derivative ∇F.
        
        ∇F = Σᵢ eⁱ ∂ᵢF
        
        This combines the split operator parts:
        - ∇_D: partial derivatives (finite differences)
        - ∇_G: geometric product with basis vectors
        
        Args:
            field: Multivector field
            stencil: Differencing scheme
            
        Returns:
            Geometric derivative of field
        """
        result = np.zeros_like(field)
        
        for i in range(self.dim):
            # ∂F/∂xᵢ (the ∇_D part)
            partial_i = self.partial_derivative(field, i, stencil)
            # eⁱ * (∂F/∂xᵢ) (the ∇_G part)
            result += self.apply_basis_vector(partial_i, i)
        
        return result
    
    # ===========================================
    # LEAPFROG TIME STEPPING (Eelco's approach)
    # ===========================================
    
    def split_even_odd(self, field: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Split multivector into even and odd grade parts."""
        even = np.zeros_like(field)
        odd = np.zeros_like(field)
        
        for idx, blade in enumerate(self.blades):
            if len(blade) % 2 == 0:
                even[..., idx] = field[..., idx]
            else:
                odd[..., idx] = field[..., idx]
        
        return even, odd
    
    def leapfrog_step(self, phi_prev: np.ndarray, phi_curr: np.ndarray,
                     dt: float, time_dir: int = 0) -> np.ndarray:
        """
        Leapfrog time step for spacetime geometric equation ∂φ = 0.
        
        For Minkowski spacetime, splits into spatial and temporal parts:
        ∂φ = ∂_t φ + ∂_s φ = 0
        
        where ∂_s is the spatial geometric derivative.
        
        Args:
            phi_prev: Field at t - dt
            phi_curr: Field at t
            dt: Time step
            time_dir: Which dimension is time (default 0)
            
        Returns:
            Field at t + dt
        """
        # Compute spatial geometric derivative (excluding time direction)
        d_spatial = np.zeros_like(phi_curr)
        for i in range(self.dim):
            if i != time_dir:
                partial_i = self.partial_derivative(phi_curr, i)
                d_spatial += self.apply_basis_vector(partial_i, i)
        
        # For ∂φ = 0 → ∂_t φ = -∂_s φ
        # Leapfrog: φ(t+dt) = φ(t-dt) + 2*dt*∂_t φ
        #         = φ(t-dt) - 2*dt*∂_s φ
        
        # Apply time direction metric
        time_metric = self.metric[time_dir, time_dir]
        
        phi_next = phi_prev - 2 * dt * time_metric * d_spatial
        
        return phi_next


# ===========================================
# COMPARISON EXAMPLE
# ===========================================

def demonstrate_equivalence():
    """
    Demonstrate the conceptual equivalence between
    galgebra's split operator and discrete implementation.
    """
    print("=" * 60)
    print("Discrete Geometric Derivative vs GAlgebra Split Operator")
    print("=" * 60)
    
    # Create 2D Euclidean algebra on 32x32 grid
    ga = DiscreteGA(signature=(2, 0, 0), shape=(32, 32), spacing=0.1)
    
    # Create a test scalar field: f(x,y) = sin(x)*cos(y)
    x = np.linspace(0, 2*np.pi, 32)
    y = np.linspace(0, 2*np.pi, 32)
    X, Y = np.meshgrid(x, y, indexing='ij')
    
    scalar_values = np.sin(X) * np.cos(Y)
    field = ga.scalar_field(scalar_values)
    
    # Compute geometric derivative
    grad_field = ga.geometric_derivative(field, stencil='central')
    
    # Extract components
    # For scalar input, ∇f gives vector (grade-1)
    # grad f = e₁ ∂f/∂x + e₂ ∂f/∂y
    
    grad_x = grad_field[..., ga._blade_to_index((0,))]  # e₁ component
    grad_y = grad_field[..., ga._blade_to_index((1,))]  # e₂ component
    
    # Compare with analytical gradient
    # ∂f/∂x = cos(x)*cos(y)
    # ∂f/∂y = -sin(x)*sin(y)
    analytical_dx = np.cos(X) * np.cos(Y)
    analytical_dy = -np.sin(X) * np.sin(Y)
    
    # Compute errors (interior points only due to boundary effects)
    error_x = np.abs(grad_x[5:-5, 5:-5] - analytical_dx[5:-5, 5:-5]).max()
    error_y = np.abs(grad_y[5:-5, 5:-5] - analytical_dy[5:-5, 5:-5]).max()
    
    print(f"\nTest: Gradient of f(x,y) = sin(x)*cos(y)")
    print(f"Grid: 32x32, spacing h=0.1")
    print(f"\n∇f = e₁ ∂f/∂x + e₂ ∂f/∂y")
    print(f"\nMax error in ∂f/∂x: {error_x:.6f}")
    print(f"Max error in ∂f/∂y: {error_y:.6f}")
    print(f"\nExpected O(h²) = {0.1**2:.6f}")
    
    print("\n" + "=" * 60)
    print("The discrete implementation matches the split operator concept:")
    print("  ∇ = ∇_G ⊗ ∇_D")
    print("  ∇F = Σᵢ eⁱ (∂F/∂xᵢ)")
    print("where:")
    print("  ∇_D → finite difference partial derivatives")
    print("  ∇_G → geometric product with basis vectors")
    print("=" * 60)


if __name__ == "__main__":
    demonstrate_equivalence()
```

## 5. Key Insights

### 5.1 The Split is Fundamental

Both approaches inherently split the geometric derivative into:
1. **Scalar differentiation** (∂/∂xᵢ or finite differences)
2. **Algebraic multiplication** (eⁱ * ... or array index manipulation)

### 5.2 Staggered Grids

Eelco's approach naturally leads to **staggered grids** where different grades
live at different locations:

```
Grade 0 (vertices):     ●───●───●
                        │   │   │  
Grade 1 (edges):        ─   ─   ─
                        │   │   │
Grade 2 (faces):        ■   ■   ■
```

This is analogous to the **Whitney forms** in Finite Element Exterior Calculus.

### 5.3 Spacetime and Leapfrog

For spacetime algebras Cl(1,3) (Minkowski), the geometric equation ∂φ = 0 
naturally leads to wave equations. The leapfrog scheme:
- Is second-order accurate
- Preserves energy (symplectic)
- Matches the natural time-reversal symmetry of wave equations

### 5.4 Can We Implement Eelco's Approach in GAlgebra?

**Partially**. GAlgebra is symbolic, so we could:
1. Generate the discrete stencil symbolically
2. Code-generate the numerical kernels

But direct numerical simulation would require interfacing with NumPy/JAX.

## 6. Recommendations

### For Symbolic Analysis
Use **galgebra** when you need:
- Exact symbolic expressions
- Proofs and verifications
- LaTeX output
- Non-cubic domains

### For Numerical Simulation
Use **Eelco's approach** (or implement the above DiscreteGA) when you need:
- Fast numerical solutions
- Large-scale simulations
- GPU acceleration (JAX/OpenCL)
- Real-time applications

### For Both
Consider a **code generation** approach:
1. Define problem symbolically in galgebra
2. Generate optimized numerical kernels
3. Run with NumPy/JAX backend

## 7. Conclusion

The split differential operator concept from galgebra and Eelco's discrete 
geometric calculus are **the same mathematical idea** expressed in different 
computational paradigms:

```
         GAlgebra                    Discrete GA
            ↓                            ↓
    Symbolic SymPy expr        Numerical NumPy arrays
            ↓                            ↓
    ∇ = Σᵢ eⁱ ∂ᵢ              ∂ = Σᵢ eⁱ Δᵢ/h
            ↓                            ↓
    Exact differentiation     Finite differences
```

The implementation provided above bridges these worlds by implementing the 
galgebra-style split operator concept using discrete numerical methods.
