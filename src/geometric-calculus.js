/**
 * Discrete Geometric Calculus for Contact Thermodynamics
 * 
 * Provides the split differential operator and discrete field operations
 * for PDE solving on grids, based on the geometric algebra skill.
 * 
 * Key concept: ∇ = Σᵢ eⁱ ⊗ ∂/∂xᵢ (split form)
 * 
 * @module geometric-calculus
 * @license MIT
 */

(function (global) {
    'use strict';

    const EPSILON = 1e-10;
    const abs = Math.abs;
    const sqrt = Math.sqrt;

    // ============================================================================
    // SCALAR FIELD
    // ============================================================================

    /**
     * Discrete scalar field on a regular grid
     */
    class ScalarField {
        /**
         * @param {number[]} shape - Grid dimensions [nx, ny] or [nx, ny, nz]
         * @param {number[]} spacing - Grid spacing [dx, dy] or [dx, dy, dz]
         * @param {Float64Array|null} data - Optional initial data
         */
        constructor(shape, spacing, data = null) {
            this.shape = shape;
            this.spacing = spacing;
            this.dim = shape.length;

            const totalSize = shape.reduce((a, b) => a * b, 1);
            this.data = data || new Float64Array(totalSize);
        }

        /**
         * Get linear index from grid indices
         */
        index(...coords) {
            let idx = 0;
            let stride = 1;
            for (let d = this.dim - 1; d >= 0; d--) {
                idx += coords[d] * stride;
                stride *= this.shape[d];
            }
            return idx;
        }

        /**
         * Get grid indices from linear index
         */
        coords(idx) {
            const result = new Array(this.dim);
            for (let d = this.dim - 1; d >= 0; d--) {
                result[d] = idx % this.shape[d];
                idx = Math.floor(idx / this.shape[d]);
            }
            return result;
        }

        /**
         * Get value at grid indices
         */
        get(...coords) {
            return this.data[this.index(...coords)];
        }

        /**
         * Set value at grid indices
         */
        set(value, ...coords) {
            this.data[this.index(...coords)] = value;
        }

        /**
         * Initialize from function f(x, y) or f(x, y, z)
         */
        initFromFunction(f) {
            for (let i = 0; i < this.data.length; i++) {
                const coords = this.coords(i);
                const position = coords.map((c, d) => c * this.spacing[d]);
                this.data[i] = f(...position);
            }
        }

        clone() {
            return new ScalarField(
                [...this.shape],
                [...this.spacing],
                new Float64Array(this.data)
            );
        }

        /**
         * Compute L2 norm
         */
        norm() {
            let sum = 0;
            for (let i = 0; i < this.data.length; i++) {
                sum += this.data[i] * this.data[i];
            }
            return sqrt(sum);
        }
    }

    // ============================================================================
    // VECTOR FIELD
    // ============================================================================

    /**
     * Discrete vector field on a regular grid
     * Stores one vector per grid point
     */
    class VectorField {
        /**
         * @param {number[]} shape - Grid dimensions
         * @param {number[]} spacing - Grid spacing
         * @param {number} vectorDim - Vector dimension (usually same as grid dim)
         */
        constructor(shape, spacing, vectorDim = null) {
            this.shape = shape;
            this.spacing = spacing;
            this.dim = shape.length;
            this.vectorDim = vectorDim || shape.length;

            const totalSize = shape.reduce((a, b) => a * b, 1);
            // Store all components contiguously
            this.components = [];
            for (let i = 0; i < this.vectorDim; i++) {
                this.components.push(new Float64Array(totalSize));
            }
        }

        /**
         * Get linear index from grid indices
         */
        index(...coords) {
            let idx = 0;
            let stride = 1;
            for (let d = this.dim - 1; d >= 0; d--) {
                idx += coords[d] * stride;
                stride *= this.shape[d];
            }
            return idx;
        }

        /**
         * Get vector at grid indices
         */
        get(...coords) {
            const idx = this.index(...coords);
            return this.components.map(comp => comp[idx]);
        }

        /**
         * Set vector at grid indices
         */
        set(vector, ...coords) {
            const idx = this.index(...coords);
            for (let i = 0; i < this.vectorDim; i++) {
                this.components[i][idx] = vector[i];
            }
        }

        /**
         * Get specific component at grid indices
         */
        getComponent(comp, ...coords) {
            return this.components[comp][this.index(...coords)];
        }

        /**
         * Set specific component at grid indices
         */
        setComponent(comp, value, ...coords) {
            this.components[comp][this.index(...coords)] = value;
        }

        clone() {
            const result = new VectorField([...this.shape], [...this.spacing], this.vectorDim);
            for (let i = 0; i < this.vectorDim; i++) {
                result.components[i] = new Float64Array(this.components[i]);
            }
            return result;
        }
    }

    // ============================================================================
    // SPLIT DIFFERENTIAL OPERATOR
    // ============================================================================

    /**
     * SplitDifferentialOperator: The discrete geometric derivative
     * 
     * ∇ = Σᵢ eⁱ ⊗ ∂/∂xᵢ
     * 
     * Uses central differences for spatial derivatives.
     */
    class SplitDifferentialOperator {
        /**
         * @param {number[]} shape - Grid dimensions
         * @param {number[]} spacing - Grid spacing per dimension
         */
        constructor(shape, spacing) {
            this.shape = shape;
            this.spacing = spacing;
            this.dim = shape.length;
        }

        /**
         * Gradient of scalar field: ∇f = Σᵢ eⁱ ∂f/∂xᵢ
         * Returns VectorField
         */
        grad(scalarField) {
            const result = new VectorField(
                [...this.shape],
                [...this.spacing],
                this.dim
            );

            for (let i = 0; i < scalarField.data.length; i++) {
                const coords = scalarField.coords(i);

                for (let d = 0; d < this.dim; d++) {
                    // Central difference with periodic boundary
                    const coordsPlus = [...coords];
                    const coordsMinus = [...coords];

                    coordsPlus[d] = (coords[d] + 1) % this.shape[d];
                    coordsMinus[d] = (coords[d] - 1 + this.shape[d]) % this.shape[d];

                    const fPlus = scalarField.get(...coordsPlus);
                    const fMinus = scalarField.get(...coordsMinus);

                    result.components[d][i] = (fPlus - fMinus) / (2 * this.spacing[d]);
                }
            }

            return result;
        }

        /**
         * Divergence of vector field: ∇·V = Σᵢ ∂Vⁱ/∂xᵢ
         * Returns ScalarField
         */
        div(vectorField) {
            const result = new ScalarField([...this.shape], [...this.spacing]);

            for (let i = 0; i < result.data.length; i++) {
                const coords = result.coords(i);
                let divValue = 0;

                for (let d = 0; d < this.dim; d++) {
                    const coordsPlus = [...coords];
                    const coordsMinus = [...coords];

                    coordsPlus[d] = (coords[d] + 1) % this.shape[d];
                    coordsMinus[d] = (coords[d] - 1 + this.shape[d]) % this.shape[d];

                    const VPlus = vectorField.getComponent(d, ...coordsPlus);
                    const VMinus = vectorField.getComponent(d, ...coordsMinus);

                    divValue += (VPlus - VMinus) / (2 * this.spacing[d]);
                }

                result.data[i] = divValue;
            }

            return result;
        }

        /**
         * Laplacian of scalar field: ∇²f = ∇·∇f = Σᵢ ∂²f/∂xᵢ²
         * Uses 5-point stencil in 2D, 7-point in 3D (O(h²) accuracy)
         */
        laplacian(scalarField) {
            const result = new ScalarField([...this.shape], [...this.spacing]);

            for (let i = 0; i < result.data.length; i++) {
                const coords = result.coords(i);
                let lapValue = 0;
                const f0 = scalarField.get(...coords);

                for (let d = 0; d < this.dim; d++) {
                    const coordsPlus = [...coords];
                    const coordsMinus = [...coords];

                    coordsPlus[d] = (coords[d] + 1) % this.shape[d];
                    coordsMinus[d] = (coords[d] - 1 + this.shape[d]) % this.shape[d];

                    const fPlus = scalarField.get(...coordsPlus);
                    const fMinus = scalarField.get(...coordsMinus);
                    const h = this.spacing[d];

                    lapValue += (fPlus - 2 * f0 + fMinus) / (h * h);
                }

                result.data[i] = lapValue;
            }

            return result;
        }

        /**
         * Curl of vector field (2D): ∇×V = ∂V²/∂x¹ - ∂V¹/∂x²
         * Returns ScalarField (the out-of-plane component)
         */
        curl2D(vectorField) {
            if (this.dim !== 2) {
                throw new Error('curl2D only works for 2D fields');
            }

            const result = new ScalarField([...this.shape], [...this.spacing]);

            for (let i = 0; i < result.data.length; i++) {
                const coords = result.coords(i);

                // ∂V²/∂x¹
                const coordsXPlus = [...coords];
                const coordsXMinus = [...coords];
                coordsXPlus[0] = (coords[0] + 1) % this.shape[0];
                coordsXMinus[0] = (coords[0] - 1 + this.shape[0]) % this.shape[0];
                const dV2dx = (vectorField.getComponent(1, ...coordsXPlus) -
                    vectorField.getComponent(1, ...coordsXMinus)) / (2 * this.spacing[0]);

                // ∂V¹/∂x²
                const coordsYPlus = [...coords];
                const coordsYMinus = [...coords];
                coordsYPlus[1] = (coords[1] + 1) % this.shape[1];
                coordsYMinus[1] = (coords[1] - 1 + this.shape[1]) % this.shape[1];
                const dV1dy = (vectorField.getComponent(0, ...coordsYPlus) -
                    vectorField.getComponent(0, ...coordsYMinus)) / (2 * this.spacing[1]);

                result.data[i] = dV2dx - dV1dy;
            }

            return result;
        }

        /**
         * Curl of vector field (3D): ∇×V
         * Returns VectorField
         */
        curl3D(vectorField) {
            if (this.dim !== 3) {
                throw new Error('curl3D only works for 3D fields');
            }

            const result = new VectorField([...this.shape], [...this.spacing], 3);

            for (let i = 0; i < vectorField.components[0].length; i++) {
                const coords = new ScalarField(this.shape, this.spacing).coords(i);
                const curl = [0, 0, 0];

                // Helper to get partial derivative
                const partial = (comp, wrt) => {
                    const coordsPlus = [...coords];
                    const coordsMinus = [...coords];
                    coordsPlus[wrt] = (coords[wrt] + 1) % this.shape[wrt];
                    coordsMinus[wrt] = (coords[wrt] - 1 + this.shape[wrt]) % this.shape[wrt];
                    return (vectorField.getComponent(comp, ...coordsPlus) -
                        vectorField.getComponent(comp, ...coordsMinus)) / (2 * this.spacing[wrt]);
                };

                // (∇×V)₀ = ∂V₂/∂x₁ - ∂V₁/∂x₂
                curl[0] = partial(2, 1) - partial(1, 2);
                // (∇×V)₁ = ∂V₀/∂x₂ - ∂V₂/∂x₀
                curl[1] = partial(0, 2) - partial(2, 0);
                // (∇×V)₂ = ∂V₁/∂x₀ - ∂V₀/∂x₁
                curl[2] = partial(1, 0) - partial(0, 1);

                for (let d = 0; d < 3; d++) {
                    result.components[d][i] = curl[d];
                }
            }

            return result;
        }
    }

    // ============================================================================
    // LEAPFROG TIME INTEGRATOR
    // ============================================================================

    /**
     * LeapfrogIntegrator: Time-stepping for wave equations
     * 
     * Solves: ∂²u/∂t² = c² ∇²u
     * 
     * Uses standard leapfrog (Verlet):
     *   u(t+Δt) = 2u(t) - u(t-Δt) + (cΔt)² ∇²u(t)
     * 
     * Stability: CFL condition c·Δt/Δx < 1/√dim
     */
    class LeapfrogIntegrator {
        /**
         * @param {number[]} shape - Grid dimensions
         * @param {number[]} spacing - Grid spacing
         */
        constructor(shape, spacing) {
            this.nabla = new SplitDifferentialOperator(shape, spacing);
            this.shape = shape;
            this.spacing = spacing;
            this.dim = shape.length;
        }

        /**
         * Estimate stable time step based on CFL condition
         */
        estimateCFL(c = 1.0, safety = 0.9) {
            const minSpacing = Math.min(...this.spacing);
            return safety * minSpacing / (c * sqrt(this.dim));
        }

        /**
         * Single wave equation step
         * 
         * @param {ScalarField} u_prev - u(t-Δt)
         * @param {ScalarField} u_curr - u(t)
         * @param {number} dt - Time step
         * @param {number} c - Wave speed
         * @returns {ScalarField} u(t+Δt)
         */
        waveStep(u_prev, u_curr, dt, c = 1.0) {
            const lap_u = this.nabla.laplacian(u_curr);
            const result = new ScalarField([...this.shape], [...this.spacing]);

            const c2dt2 = c * c * dt * dt;
            for (let i = 0; i < result.data.length; i++) {
                result.data[i] = 2 * u_curr.data[i] - u_prev.data[i] + c2dt2 * lap_u.data[i];
            }

            return result;
        }

        /**
         * Run full wave simulation
         * 
         * @param {ScalarField} u0 - Initial condition u(0)
         * @param {ScalarField|null} v0 - Initial velocity ∂u/∂t(0), or null for zero
         * @param {number} dt - Time step
         * @param {number} n_steps - Number of steps
         * @param {number} c - Wave speed
         * @returns {ScalarField} Final state u(T)
         */
        waveSimulate(u0, v0, dt, n_steps, c = 1.0) {
            // Initialize u_prev = u0 - dt*v0 (first-order startup)
            let u_prev = u0.clone();
            if (v0) {
                for (let i = 0; i < u_prev.data.length; i++) {
                    u_prev.data[i] -= dt * v0.data[i];
                }
            }

            let u_curr = u0.clone();

            for (let step = 0; step < n_steps; step++) {
                const u_next = this.waveStep(u_prev, u_curr, dt, c);
                u_prev = u_curr;
                u_curr = u_next;
            }

            return u_curr;
        }

        /**
         * Heat equation step (forward Euler)
         * ∂u/∂t = α ∇²u
         * 
         * Stability: αΔt/Δx² < 1/(2·dim)
         */
        heatStep(u_curr, dt, alpha = 1.0) {
            const lap_u = this.nabla.laplacian(u_curr);
            const result = new ScalarField([...this.shape], [...this.spacing]);

            for (let i = 0; i < result.data.length; i++) {
                result.data[i] = u_curr.data[i] + alpha * dt * lap_u.data[i];
            }

            return result;
        }

        /**
         * Run heat diffusion simulation
         */
        heatSimulate(u0, dt, n_steps, alpha = 1.0) {
            let u = u0.clone();
            for (let step = 0; step < n_steps; step++) {
                u = this.heatStep(u, dt, alpha);
            }
            return u;
        }
    }

    // ============================================================================
    // EXPORTS
    // ============================================================================

    const GeometricCalculus = {
        ScalarField,
        VectorField,
        SplitDifferentialOperator,
        LeapfrogIntegrator,
        EPSILON
    };

    // Export for different module systems
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = GeometricCalculus;
    }
    if (typeof global !== 'undefined') {
        global.GeometricCalculus = GeometricCalculus;
    }
    if (typeof window !== 'undefined') {
        window.GeometricCalculus = GeometricCalculus;
    }

})(typeof globalThis !== 'undefined' ? globalThis : this);
