/**
 * Leapfrog/Yee Solvers for PDEs on Triangle Meshes via FTGC
 * 
 * Implements time-stepping schemes for:
 *   - Wave equation: ∂²u/∂t² = c²∇²u
 *   - Heat equation: ∂u/∂t = α∇²u
 *   - Maxwell equations: ∂B/∂t = -∇∧E, ∂E/∂t = c²∇·B
 * 
 * Uses the geometric derivative ∇ from MeshGeometricDerivative.
 * Supports Dirichlet boundary conditions via mask + values.
 * 
 * @module mesh-solvers
 * @license MIT
 */

(function (global) {
    'use strict';

    // Import dependencies (for CommonJS; browser assumes global)
    let MeshGeometricDerivative, TriangleMesh;
    if (typeof require !== 'undefined') {
        const ftgc = require('./mesh-ftgc.js');
        const mesh = require('./mesh.js');
        MeshGeometricDerivative = ftgc.MeshGeometricDerivative;
        TriangleMesh = mesh.TriangleMesh;
    } else if (typeof global.ContactThermo !== 'undefined') {
        MeshGeometricDerivative = global.ContactThermo.MeshGeometricDerivative;
        TriangleMesh = global.ContactThermo.TriangleMesh;
    }

    const EPSILON = 1e-10;
    const { sqrt, min, max } = Math;

    // ============================================================================
    // DIRICHLET BOUNDARY HELPERS
    // ============================================================================

    /**
     * Apply Dirichlet boundary conditions to a field.
     * @param {Float64Array} f - Field to modify (in-place)
     * @param {Uint8Array} mask - Dirichlet mask (1 = constrained)
     * @param {Float64Array} values - Values at constrained vertices
     */
    function applyDirichlet(f, mask, values) {
        if (!mask) return;
        for (let i = 0; i < f.length; i++) {
            if (mask[i]) {
                f[i] = values ? values[i] : 0;
            }
        }
    }

    /**
     * Create Dirichlet mask from mesh boundary vertices.
     * @param {TriangleMesh} mesh 
     * @returns {Uint8Array}
     */
    function boundaryDirichletMask(mesh) {
        return new Uint8Array(mesh.boundaryVertices);
    }

    // ============================================================================
    // LEAPFROG SOLVER (WAVE + HEAT + MAXWELL)
    // ============================================================================

    /**
     * Leapfrog time-stepping solver using the geometric derivative ∇.
     * 
     * Key insight: Yee's staggered grid scheme is natural in GC because
     * different grades live on different mesh elements!
     */
    class LeapfrogGCMesh {
        /**
         * @param {TriangleMesh} mesh 
         */
        constructor(mesh) {
            this.mesh = mesh;
            this.nabla = new MeshGeometricDerivative(mesh);

            // Precompute Laplacian matrices
            this.L = this.nabla.laplacianMatrix();
            this.L_weak = this.nabla.laplacianWeak();
            this.M = this.nabla._metricVertices; // Mass matrix (vertex areas)
        }

        // ========================================================================
        // STABILITY / CFL
        // ========================================================================

        /**
         * Estimate maximum stable time step (CFL condition).
         * For wave equation: c·dt/h < 1
         * @param {number} c - Wave speed
         * @returns {number} Recommended dt
         */
        estimateCFL(c = 1) {
            const edgeLengths = this.mesh.edgeLengths();
            let hMin = Infinity;
            for (let i = 0; i < edgeLengths.length; i++) {
                if (edgeLengths[i] > EPSILON) {
                    hMin = min(hMin, edgeLengths[i]);
                }
            }
            // Safety factor 0.5
            return 0.5 * hMin / c;
        }

        /**
         * Estimate maximum stable time step for explicit heat equation.
         * For heat: α·dt/h² < 0.5
         * @param {number} alpha - Diffusion coefficient
         * @returns {number} Recommended dt
         */
        estimateCFLHeat(alpha = 1) {
            const edgeLengths = this.mesh.edgeLengths();
            let hMin = Infinity;
            for (let i = 0; i < edgeLengths.length; i++) {
                if (edgeLengths[i] > EPSILON) {
                    hMin = min(hMin, edgeLengths[i]);
                }
            }
            // Safety factor 0.25
            return 0.25 * hMin * hMin / alpha;
        }

        // ========================================================================
        // WAVE EQUATION: ∂²u/∂t² = c²∇²u
        // ========================================================================

        /**
         * One leapfrog step for the wave equation.
         * u(t+dt) = 2u(t) - u(t-dt) + c²dt²∇²u(t)
         * 
         * @param {Float64Array} uPrev - u at t-dt
         * @param {Float64Array} uCurr - u at t
         * @param {number} dt - Time step
         * @param {number} c - Wave speed
         * @param {Object} opts - Options
         * @param {Uint8Array} opts.dirichletMask - Dirichlet mask
         * @param {Float64Array} opts.dirichletValues - Dirichlet values
         * @returns {Float64Array} u at t+dt
         */
        waveStep(uPrev, uCurr, dt, c = 1, opts = {}) {
            const { dirichletMask, dirichletValues } = opts;

            // Compute Laplacian using GC: ∇² = ∇·∇
            const Lu = this.L.matvec(uCurr);

            // Leapfrog update
            const uNext = new Float64Array(uCurr.length);
            const c2dt2 = c * c * dt * dt;
            for (let i = 0; i < uNext.length; i++) {
                uNext[i] = 2 * uCurr[i] - uPrev[i] + c2dt2 * Lu[i];
            }

            // Apply Dirichlet BCs
            applyDirichlet(uNext, dirichletMask, dirichletValues);

            return uNext;
        }

        /**
         * Simulate wave equation from initial conditions.
         * 
         * @param {Float64Array} u0 - Initial displacement
         * @param {Float64Array} v0 - Initial velocity
         * @param {number} dt - Time step
         * @param {number} nSteps - Number of steps
         * @param {number} c - Wave speed
         * @param {Object} opts - Options (dirichletMask, dirichletValues, callback)
         * @returns {Float64Array} Final displacement
         */
        waveSimulate(u0, v0, dt, nSteps, c = 1, opts = {}) {
            const { dirichletMask, dirichletValues, callback } = opts;

            // Initialize using Taylor expansion: u(dt) ≈ u0 + dt·v0 + 0.5·c²·dt²·∇²u0
            const Lu0 = this.L.matvec(u0);
            const uPrev = new Float64Array(u0);
            const uCurr = new Float64Array(u0.length);
            const c2dt2_half = 0.5 * c * c * dt * dt;
            for (let i = 0; i < uCurr.length; i++) {
                uCurr[i] = u0[i] + dt * v0[i] + c2dt2_half * Lu0[i];
            }
            applyDirichlet(uCurr, dirichletMask, dirichletValues);

            // Time stepping
            for (let step = 0; step < nSteps; step++) {
                const uNext = this.waveStep(uPrev, uCurr, dt, c, opts);
                
                // Shift
                for (let i = 0; i < uPrev.length; i++) {
                    uPrev[i] = uCurr[i];
                    uCurr[i] = uNext[i];
                }

                if (callback) callback(step, uCurr);
            }

            return uCurr;
        }

        // ========================================================================
        // HEAT EQUATION: ∂u/∂t = α∇²u
        // ========================================================================

        /**
         * Explicit Euler step for heat equation.
         * u(t+dt) = u(t) + α·dt·∇²u(t)
         * 
         * @param {Float64Array} u - Current field
         * @param {number} dt - Time step
         * @param {number} alpha - Diffusion coefficient
         * @param {Object} opts - Options (dirichletMask, dirichletValues)
         * @returns {Float64Array} Updated field
         */
        heatStepExplicit(u, dt, alpha = 1, opts = {}) {
            const { dirichletMask, dirichletValues } = opts;

            const Lu = this.L.matvec(u);
            const uNext = new Float64Array(u.length);
            const adt = alpha * dt;
            for (let i = 0; i < uNext.length; i++) {
                uNext[i] = u[i] + adt * Lu[i];
            }

            applyDirichlet(uNext, dirichletMask, dirichletValues);
            return uNext;
        }

        /**
         * Implicit Euler step for heat equation (unconditionally stable).
         * (I - α·dt·∇²) u(t+dt) = u(t)
         * 
         * Uses simple Jacobi iteration for the solve.
         * 
         * @param {Float64Array} u - Current field
         * @param {number} dt - Time step
         * @param {number} alpha - Diffusion coefficient
         * @param {Object} opts - Options (dirichletMask, dirichletValues, maxIter, tol)
         * @returns {Float64Array} Updated field
         */
        heatStepImplicit(u, dt, alpha = 1, opts = {}) {
            const { dirichletMask, dirichletValues, maxIter = 100, tol = 1e-8 } = opts;

            const n = u.length;
            const adt = alpha * dt;

            // Build diagonal of (I - α·dt·L)
            const diagL = this.L.diagonal();
            const diagA = new Float64Array(n);
            for (let i = 0; i < n; i++) {
                diagA[i] = 1 - adt * diagL[i];
            }

            // Jacobi iteration: x_{k+1} = D⁻¹ (b - (A - D) x_k)
            // where A = I - α·dt·L, b = u
            const uNext = new Float64Array(u);
            const rhs = new Float64Array(u);

            for (let iter = 0; iter < maxIter; iter++) {
                // Compute residual: r = b - A·x = u - (x - α·dt·L·x) = u - x + α·dt·L·x
                const Lx = this.L.matvec(uNext);
                let maxResid = 0;

                for (let i = 0; i < n; i++) {
                    if (dirichletMask && dirichletMask[i]) continue;

                    const Ax_i = uNext[i] - adt * Lx[i];
                    const resid = rhs[i] - Ax_i;
                    maxResid = max(maxResid, Math.abs(resid));
                    uNext[i] += resid / diagA[i];
                }

                applyDirichlet(uNext, dirichletMask, dirichletValues);

                if (maxResid < tol) break;
            }

            return uNext;
        }

        /**
         * Simulate heat equation from initial conditions.
         * 
         * @param {Float64Array} u0 - Initial temperature/entropy
         * @param {number} dt - Time step
         * @param {number} nSteps - Number of steps
         * @param {number} alpha - Diffusion coefficient
         * @param {Object} opts - Options (implicit, dirichletMask, dirichletValues, callback)
         * @returns {Float64Array} Final field
         */
        heatSimulate(u0, dt, nSteps, alpha = 1, opts = {}) {
            const { implicit = false, callback } = opts;

            let u = new Float64Array(u0);

            for (let step = 0; step < nSteps; step++) {
                u = implicit 
                    ? this.heatStepImplicit(u, dt, alpha, opts)
                    : this.heatStepExplicit(u, dt, alpha, opts);

                if (callback) callback(step, u);
            }

            return u;
        }

        // ========================================================================
        // MAXWELL EQUATIONS VIA GEOMETRIC CALCULUS
        // ========================================================================

        /**
         * Yee-style leapfrog step for Maxwell equations.
         * 
         * In geometric calculus, Maxwell's equations unify:
         *   ∇F = J  (F = E + IB is the electromagnetic field bivector)
         * 
         * Splitting into time and space:
         *   ∂B/∂t = -∇∧E  (Faraday: curl of E)
         *   ∂E/∂t = c²∇·B  (Ampere-like, with B on faces)
         * 
         * Where:
         *   - E is a vector field (grade 1) on edges
         *   - B is a bivector field (grade 2) on faces
         * 
         * @param {Float64Array} E - Electric field on edges
         * @param {Float64Array} B - Magnetic field on faces
         * @param {number} dt - Time step
         * @param {number} c - Speed of light
         * @returns {{E: Float64Array, B: Float64Array}} Updated fields
         */
        maxwellStep(E, B, dt, c = 1) {
            const c2 = c * c;

            // Faraday: ∂B/∂t = -∇∧E
            const curlE = this.nabla.curl(E);
            const B_new = new Float64Array(B.length);
            for (let i = 0; i < B.length; i++) {
                B_new[i] = B[i] - dt * curlE[i];
            }

            // Ampere: ∂E/∂t = c²∇·B (using inner derivative on bivector)
            // ∇·B = M₁⁻¹ d₁ᵀ M₂ B
            const tmp = new Float64Array(this.mesh.nFaces);
            for (let i = 0; i < tmp.length; i++) {
                tmp[i] = this.nabla._metricFaces[i] * B_new[i];
            }
            const tmp2 = this.nabla._wedge_1to2.matvecT(tmp);
            const divB = new Float64Array(this.mesh.nEdges);
            for (let i = 0; i < divB.length; i++) {
                divB[i] = this.nabla._metricEdgesInv[i] * tmp2[i];
            }

            const E_new = new Float64Array(E.length);
            for (let i = 0; i < E.length; i++) {
                E_new[i] = E[i] + dt * c2 * divB[i];
            }

            return { E: E_new, B: B_new };
        }

        // ========================================================================
        // UTILITY: ENERGY COMPUTATION
        // ========================================================================

        /**
         * Compute total "energy" of a scalar field (L2 norm weighted by mass).
         * E = ∫ u² dA = Σᵢ u_i² · A_i
         * 
         * @param {Float64Array} u - Scalar field on vertices
         * @returns {number} Total energy
         */
        energy(u) {
            let E = 0;
            for (let i = 0; i < u.length; i++) {
                E += u[i] * u[i] * this.M[i];
            }
            return E;
        }

        /**
         * Compute total "mass" of a scalar field.
         * M = ∫ u dA = Σᵢ u_i · A_i
         * 
         * @param {Float64Array} u - Scalar field on vertices
         * @returns {number} Total mass
         */
        mass(u) {
            let m = 0;
            for (let i = 0; i < u.length; i++) {
                m += u[i] * this.M[i];
            }
            return m;
        }

        /**
         * Compute variance of a scalar field.
         * Var = ∫ (u - mean)² dA
         * 
         * @param {Float64Array} u - Scalar field on vertices
         * @returns {number} Variance
         */
        variance(u) {
            const totalArea = this.M.reduce((a, b) => a + b, 0);
            const mean = this.mass(u) / totalArea;
            
            let variance = 0;
            for (let i = 0; i < u.length; i++) {
                const diff = u[i] - mean;
                variance += diff * diff * this.M[i];
            }
            return variance;
        }
    }

    // ============================================================================
    // EXPORTS
    // ============================================================================

    const MeshSolversModule = {
        applyDirichlet,
        boundaryDirichletMask,
        LeapfrogGCMesh
    };

    // CommonJS
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = MeshSolversModule;
    }
    // AMD
    else if (typeof define === 'function' && define.amd) {
        define(function () { return MeshSolversModule; });
    }
    // Browser global
    else {
        global.ContactThermo = global.ContactThermo || {};
        Object.assign(global.ContactThermo, MeshSolversModule);
    }

})(typeof globalThis !== 'undefined' ? globalThis : typeof window !== 'undefined' ? window : this);
