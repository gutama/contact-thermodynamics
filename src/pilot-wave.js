/**
 * Pilot-Wave Dynamics Module
 * 
 * Implements Valentini's regularized de Broglie-Bohm pilot-wave theory:
 * - WaveFunction class with multivector phase S (geometric derivative approach)
 * - Regularized de Broglie velocity (finite at nodes)
 * - H-theorem for quantum relaxation to Born rule equilibrium
 * - Integration with contact geometry (S ↔ action variable A)
 * 
 * Key insight: The phase S in ψ = R·exp(i·S) can be a multivector,
 * enabling treatment of spinor wavefunctions and curved spaces.
 * 
 * References:
 * - Valentini, "Pilot-Wave Dynamics and the Primordial Violation of the Born Rule"
 * - de Broglie, "La mécanique ondulatoire" (1927)
 * 
 * @license MIT
 */

(function (global) {
    'use strict';

    const EPSILON = 1e-12;
    const PI = Math.PI;
    const sqrt = Math.sqrt;
    const exp = Math.exp;
    const log = Math.log;
    const cos = Math.cos;
    const sin = Math.sin;

    // ============================================================================
    // SMEARING KERNEL
    // ============================================================================

    /**
     * SmearingKernel: Regularization kernel μ(x) for Valentini regularization
     * 
     * Convolves fields to ensure |ψ|² never vanishes exactly, keeping
     * de Broglie velocity finite at wavefunction nodes.
     * 
     * Kernel satisfies: ∫ μ(x) dx = 1 (normalization)
     * 
     * For time-dependent regularization, width can be a function ε(t).
     * This enables dynamically unstable equilibria where ρ ≠ (|ψ|²)_reg
     * after short-time processes.
     */
    class SmearingKernel {
        /**
         * @param {number|Function} width - Regularization scale ε or ε(t) function
         * @param {string} type - Kernel type: 'gaussian', 'uniform', 'lorentzian'
         */
        constructor(width = 1e-10, type = 'gaussian') {
            // Support both fixed and time-dependent width
            this._widthValue = typeof width === 'function' ? null : width;
            this._widthFunction = typeof width === 'function' ? width : null;
            this.type = type;
            this._time = 0;  // Current time for time-dependent kernels
        }

        /**
         * Get current width (evaluating function if time-dependent)
         * @returns {number} Current regularization scale ε
         */
        get width() {
            if (this._widthFunction) {
                return this._widthFunction(this._time);
            }
            return this._widthValue;
        }

        /**
         * Set current time for time-dependent kernels
         * @param {number} t - Time
         */
        setTime(t) {
            this._time = t;
        }

        /**
         * Check if kernel has time-dependent width
         * @returns {boolean}
         */
        get isTimeDependant() {
            return this._widthFunction !== null;
        }

        /**
         * Evaluate kernel at distance r from center
         * @param {number} r - Distance from center
         * @returns {number} Kernel value μ(r)
         */
        evaluate(r) {
            const eps = this.width;
            switch (this.type) {
                case 'gaussian':
                    // Normalized Gaussian: (1/(ε√(2π))) exp(-r²/(2ε²))
                    return exp(-r * r / (2 * eps * eps)) / (eps * sqrt(2 * PI));
                case 'uniform':
                    // Uniform ball: θ(ε - |r|) / (2ε)  (1D)
                    return Math.abs(r) <= eps ? 0.5 / eps : 0;
                case 'lorentzian':
                    // Cauchy/Lorentzian: ε / (π(r² + ε²))
                    return eps / (PI * (r * r + eps * eps));
                default:
                    return exp(-r * r / (2 * eps * eps)) / (eps * sqrt(2 * PI));
            }
        }

        /**
         * Convolve a 1D field with the kernel (discrete approximation)
         * f_reg(x) = ∫ f(y) μ(x-y) dy
         * 
         * @param {number[]} field - Field values on uniform grid
         * @param {number} dx - Grid spacing
         * @returns {number[]} Convolved field
         */
        convolve1D(field, dx) {
            const n = field.length;
            const result = new Array(n).fill(0);

            // Kernel support (3-5 widths for Gaussian)
            const support = Math.ceil(4 * this.width / dx);

            for (let i = 0; i < n; i++) {
                let sum = 0;
                let normalization = 0;

                for (let j = -support; j <= support; j++) {
                    const idx = i + j;
                    if (idx >= 0 && idx < n) {
                        const r = j * dx;
                        const weight = this.evaluate(r) * dx;
                        sum += field[idx] * weight;
                        normalization += weight;
                    }
                }

                // Normalize to ensure ∫μ = 1 even at boundaries
                result[i] = normalization > EPSILON ? sum / normalization : field[i];
            }

            return result;
        }

        /**
         * Convolve a 2D field with the kernel
         * @param {number[][]} field - 2D array of field values
         * @param {number} dx - Grid spacing in x
         * @param {number} dy - Grid spacing in y
         * @returns {number[][]} Convolved field
         */
        convolve2D(field, dx, dy = dx) {
            const nx = field.length;
            const ny = field[0].length;
            const result = [];

            const supportX = Math.ceil(4 * this.width / dx);
            const supportY = Math.ceil(4 * this.width / dy);

            for (let i = 0; i < nx; i++) {
                result[i] = [];
                for (let j = 0; j < ny; j++) {
                    let sum = 0;
                    let normalization = 0;

                    for (let di = -supportX; di <= supportX; di++) {
                        for (let dj = -supportY; dj <= supportY; dj++) {
                            const ii = i + di;
                            const jj = j + dj;
                            if (ii >= 0 && ii < nx && jj >= 0 && jj < ny) {
                                const r = sqrt((di * dx) ** 2 + (dj * dy) ** 2);
                                const weight = this.evaluate(r) * dx * dy;
                                sum += field[ii][jj] * weight;
                                normalization += weight;
                            }
                        }
                    }

                    result[i][j] = normalization > EPSILON ? sum / normalization : field[i][j];
                }
            }

            return result;
        }

        toString() {
            return `SmearingKernel(width=${this.width.toExponential(2)}, type="${this.type}")`;
        }
    }

    // ============================================================================
    // WAVEFUNCTION CLASS
    // ============================================================================

    /**
     * WaveFunction: Complex wavefunction ψ with multivector phase
     * 
     * Represents ψ = R · exp(i·S) where:
     * - R is the real amplitude (scalar field)
     * - S is the phase (can be scalar or multivector for spinors)
     * 
     * The geometric derivative ∇ψ uses the full GA machinery when S is a multivector.
     */
    class WaveFunction {
        /**
         * Create wavefunction from amplitude and phase
         * @param {number|number[]} amplitude - Amplitude R (scalar or field)
         * @param {number|number[]|Object} phase - Phase S (scalar, field, or multivector)
         * @param {Object} options - Configuration options
         */
        constructor(amplitude, phase, options = {}) {
            this.amplitude = amplitude;  // R
            this.phase = phase;          // S (scalar or multivector)
            this.mass = options.mass || 1.0;
            this.hbar = options.hbar || 1.0;  // ℏ (set to 1 for natural units)
            this.dimension = options.dimension || 1;

            // For multivector phase (spinor case)
            this.isMultivectorPhase = typeof phase === 'object' && phase.grade !== undefined;
        }

        /**
         * Create from Cartesian form ψ = ψ_r + i·ψ_i
         * @param {number|number[]} realPart - Real part ψ_r
         * @param {number|number[]} imagPart - Imaginary part ψ_i
         * @param {Object} options - Configuration
         * @returns {WaveFunction}
         */
        static fromCartesian(realPart, imagPart, options = {}) {
            if (typeof realPart === 'number') {
                const R = sqrt(realPart * realPart + imagPart * imagPart);
                const S = Math.atan2(imagPart, realPart);
                return new WaveFunction(R, S, options);
            }

            // Array case
            const R = realPart.map((re, i) => sqrt(re * re + imagPart[i] * imagPart[i]));
            const S = realPart.map((re, i) => Math.atan2(imagPart[i], re));
            return new WaveFunction(R, S, options);
        }

        /**
         * Create plane wave ψ = exp(i·k·x)
         * @param {number|number[]} k - Wave vector
         * @param {Object} options - Configuration
         * @returns {Function} Function of position x returning WaveFunction value
         */
        static planeWave(k, options = {}) {
            const hbar = options.hbar || 1.0;
            const mass = options.mass || 1.0;

            // Phase is S = k·x (scalar), velocity is v = ℏk/m
            return function (x) {
                const kx = Array.isArray(k)
                    ? k.reduce((sum, ki, i) => sum + ki * x[i], 0)
                    : k * x;
                return new WaveFunction(1.0, kx, options);
            };
        }

        /**
         * Create Gaussian wave packet
         * @param {number} x0 - Center position
         * @param {number} sigma - Width
         * @param {number} k0 - Central wave number
         * @param {Object} options - Configuration
         * @returns {Function}
         */
        static gaussianPacket(x0, sigma, k0, options = {}) {
            return function (x) {
                const dx = x - x0;
                const envelope = exp(-dx * dx / (4 * sigma * sigma));
                const phase = k0 * x;
                return new WaveFunction(envelope, phase, options);
            };
        }

        /**
         * Get probability density |ψ|² = R²
         * @returns {number|number[]}
         */
        get probabilityDensity() {
            if (typeof this.amplitude === 'number') {
                return this.amplitude * this.amplitude;
            }
            return this.amplitude.map(R => R * R);
        }

        /**
         * Get real part: Re(ψ) = R·cos(S)
         * @returns {number|number[]}
         */
        get realPart() {
            if (typeof this.amplitude === 'number') {
                return this.amplitude * cos(this.phase);
            }
            return this.amplitude.map((R, i) => R * cos(this.phase[i]));
        }

        /**
         * Get imaginary part: Im(ψ) = R·sin(S)
         * @returns {number|number[]}
         */
        get imagPart() {
            if (typeof this.amplitude === 'number') {
                return this.amplitude * sin(this.phase);
            }
            return this.amplitude.map((R, i) => R * sin(this.phase[i]));
        }

        /**
         * Compute probability current j = (ℏ/m) R² ∇S
         * For scalar phase: j^i = (ℏ/m) |ψ|² ∂_i S
         * 
         * @param {number} dx - Grid spacing for numerical gradient
         * @returns {number[]} Current vector (or array of vectors for field)
         */
        current(dx = 0.01) {
            const hbar = this.hbar;
            const m = this.mass;

            if (typeof this.phase === 'number') {
                // Single point - need external gradient info
                throw new Error('Cannot compute current at single point without gradient');
            }

            // 1D field case
            const n = this.phase.length;
            const j = new Array(n).fill(0);

            // Central difference for ∂S/∂x
            for (let i = 1; i < n - 1; i++) {
                const gradS = (this.phase[i + 1] - this.phase[i - 1]) / (2 * dx);
                const rhoSq = this.amplitude[i] * this.amplitude[i];
                j[i] = (hbar / m) * rhoSq * gradS;
            }

            // Boundary (one-sided differences)
            j[0] = (hbar / m) * this.amplitude[0] ** 2 * (this.phase[1] - this.phase[0]) / dx;
            j[n - 1] = (hbar / m) * this.amplitude[n - 1] ** 2 * (this.phase[n - 1] - this.phase[n - 2]) / dx;

            return j;
        }

        toString() {
            if (typeof this.amplitude === 'number') {
                return `ψ = ${this.amplitude.toFixed(4)} exp(i·${this.phase.toFixed(4)})`;
            }
            return `WaveFunction(n=${this.amplitude.length}, m=${this.mass}, ℏ=${this.hbar})`;
        }
    }

    // ============================================================================
    // PILOT-WAVE SYSTEM
    // ============================================================================

    /**
     * PilotWaveSystem: Implements de Broglie-Bohm guidance with Valentini regularization
     * 
     * Core equations:
     * - Guidance: v = (ℏ/m) ∇S (de Broglie, first-order)
     * - Regularized: v_reg = j_reg / ρ_reg where ρ_reg = |ψ|² ∗ μ
     * 
     * The regularization ensures v stays finite at wavefunction nodes.
     */
    class PilotWaveSystem {
        /**
         * @param {WaveFunction} wavefunction - The pilot wave ψ
         * @param {Object} options - Configuration
         */
        constructor(wavefunction, options = {}) {
            this.psi = wavefunction;
            this.kernel = options.kernel || new SmearingKernel(1e-10, 'gaussian');
            this.dx = options.dx || 0.01;  // Grid spacing

            // Cached regularized fields
            this._rhoReg = null;
            this._jReg = null;
        }

        /**
         * Compute de Broglie velocity v = (ℏ/m) ∇S
         * WARNING: Diverges at nodes where ψ = 0
         * 
         * @param {number[]} position - Position in configuration space
         * @returns {number[]} Velocity vector
         */
        deBroglieVelocity(position) {
            const hbar = this.psi.hbar;
            const m = this.psi.mass;
            const phase = this.psi.phase;

            if (typeof phase === 'number') {
                throw new Error('Need field values for velocity computation');
            }

            // For 1D: v = (ℏ/m) ∂S/∂x
            const n = phase.length;
            const v = new Array(n).fill(0);

            for (let i = 1; i < n - 1; i++) {
                v[i] = (hbar / m) * (phase[i + 1] - phase[i - 1]) / (2 * this.dx);
            }

            v[0] = (hbar / m) * (phase[1] - phase[0]) / this.dx;
            v[n - 1] = (hbar / m) * (phase[n - 1] - phase[n - 2]) / this.dx;

            return v;
        }

        /**
         * Compute regularized de Broglie velocity
         * v_reg = j_reg / ρ_reg = (j ∗ μ) / (|ψ|² ∗ μ)
         * 
         * This is always finite, even at nodes.
         * 
         * @returns {number[]} Regularized velocity field
         */
        regularizedVelocity() {
            // Compute j and ρ
            const j = this.psi.current(this.dx);
            const rho = this.psi.probabilityDensity;

            // Convolve with smearing kernel
            const jReg = this.kernel.convolve1D(j, this.dx);
            const rhoReg = this.kernel.convolve1D(rho, this.dx);

            // Cache for H-function computation
            this._jReg = jReg;
            this._rhoReg = rhoReg;

            // Regularized velocity
            const n = jReg.length;
            const vReg = new Array(n);

            for (let i = 0; i < n; i++) {
                // rhoReg is always > 0 due to smearing (this is the point!)
                vReg[i] = jReg[i] / Math.max(rhoReg[i], EPSILON);
            }

            return vReg;
        }

        /**
         * Get regularized probability density (|ψ|²)_reg
         * This is |ψ|² + (1/2) ε² ∇²|ψ|² + O(ε⁴) near nodes
         * 
         * @returns {number[]} Regularized density
         */
        regularizedDensity() {
            if (this._rhoReg === null) {
                const rho = this.psi.probabilityDensity;
                this._rhoReg = this.kernel.convolve1D(rho, this.dx);
            }
            return this._rhoReg;
        }

        /**
         * Evolve particle position under regularized guidance
         * dx/dt = v_reg(x(t), t)
         * 
         * @param {number} x0 - Initial position
         * @param {number} dt - Time step
         * @param {number} nSteps - Number of steps
         * @returns {number[]} Trajectory x(t)
         */
        evolveTrajectory(x0, dt, nSteps) {
            const trajectory = [x0];
            const vReg = this.regularizedVelocity();

            let x = x0;
            for (let step = 0; step < nSteps; step++) {
                // Interpolate velocity at current position
                const idx = Math.floor(x / this.dx);
                const n = vReg.length;

                if (idx < 0 || idx >= n - 1) {
                    // Out of bounds
                    break;
                }

                // Linear interpolation
                const t = (x / this.dx) - idx;
                const v = vReg[idx] * (1 - t) + vReg[idx + 1] * t;

                // Euler step
                x = x + v * dt;
                trajectory.push(x);
            }

            return trajectory;
        }

        toString() {
            return `PilotWaveSystem(kernel=${this.kernel}, dx=${this.dx})`;
        }
    }

    // ============================================================================
    // QUANTUM ENSEMBLE
    // ============================================================================

    /**
     * QuantumEnsemble: Tracks distribution ρ vs |ψ|² and H-function
     * 
     * Implements Valentini's subquantum H-theorem:
     * - H̄ = ∫ ρ̄ ln(ρ̄/|ψ̄|²) ≥ 0
     * - dH̄/dt ≤ 0 (relaxation)
     * - Equilibrium: ρ = |ψ|² ⟹ H̄ = 0
     */
    class QuantumEnsemble {
        /**
         * @param {number[]} distribution - The actual distribution ρ(q)
         * @param {WaveFunction} wavefunction - The guiding wave ψ
         * @param {number} dx - Grid spacing
         * @param {SmearingKernel} kernel - Optional kernel for regularized equilibrium
         */
        constructor(distribution, wavefunction, dx = 0.01, kernel = null) {
            this.rho = distribution;
            this.psi = wavefunction;
            this.dx = dx;
            this.kernel = kernel;  // Needed for regularized equilibrium

            // Normalization check
            this._normalizeDistribution();
        }

        _normalizeDistribution() {
            const total = this.rho.reduce((sum, r) => sum + r * this.dx, 0);
            if (Math.abs(total - 1) > 0.01) {
                // Renormalize
                this.rho = this.rho.map(r => r / total);
            }
        }

        /**
         * Check if in standard quantum equilibrium: ρ ≈ |ψ|²
         * @param {number} tolerance - Relative tolerance
         * @returns {boolean}
         */
        isInEquilibrium(tolerance = 0.01) {
            const rhoQM = this.psi.probabilityDensity;
            const n = this.rho.length;

            // Normalize both distributions for fair comparison
            const rhoTotal = this.rho.reduce((s, r) => s + r, 0);
            const qmTotal = rhoQM.reduce((s, r) => s + r, 0);

            // Use L2 norm relative deviation
            let sumSqDev = 0;
            let sumSqQM = 0;
            for (let i = 0; i < n; i++) {
                const rhoNorm = this.rho[i] / rhoTotal;
                const qmNorm = rhoQM[i] / qmTotal;
                const dev = rhoNorm - qmNorm;
                sumSqDev += dev * dev;
                sumSqQM += qmNorm * qmNorm;
            }

            const relativeL2 = Math.sqrt(sumSqDev / (sumSqQM + EPSILON));
            return relativeL2 < tolerance;
        }

        /**
         * Check if in REGULARIZED quantum equilibrium: ρ ≈ (|ψ|²)_reg
         * 
         * This is the correct equilibrium for regularized pilot-wave theory:
         * - (|ψ|²)_reg is nonzero at former nodes
         * - For time-dependent μ, equilibrium can be dynamically unstable
         * 
         * @param {SmearingKernel} kernel - Smearing kernel (uses this.kernel if not provided)
         * @param {number} tolerance - Relative tolerance
         * @returns {boolean}
         */
        isInRegularizedEquilibrium(kernel = null, tolerance = 0.01) {
            const k = kernel || this.kernel;
            if (!k) {
                throw new Error('Kernel required for regularized equilibrium check');
            }

            const rhoQM = this.psi.probabilityDensity;
            const rhoQMReg = k.convolve1D(rhoQM, this.dx);  // (|ψ|²)_reg
            const n = this.rho.length;

            // Normalize both distributions
            const rhoTotal = this.rho.reduce((s, r) => s + r, 0);
            const regTotal = rhoQMReg.reduce((s, r) => s + r, 0);

            // Use L2 norm relative deviation
            let sumSqDev = 0;
            let sumSqReg = 0;
            for (let i = 0; i < n; i++) {
                const rhoNorm = this.rho[i] / rhoTotal;
                const regNorm = rhoQMReg[i] / regTotal;
                const dev = rhoNorm - regNorm;
                sumSqDev += dev * dev;
                sumSqReg += regNorm * regNorm;
            }

            const relativeL2 = Math.sqrt(sumSqDev / (sumSqReg + EPSILON));
            return relativeL2 < tolerance;
        }

        /**
         * Get regularized equilibrium distribution (|ψ|²)_reg
         * @param {SmearingKernel} kernel - Smearing kernel
         * @returns {number[]} Regularized equilibrium distribution
         */
        getRegularizedEquilibrium(kernel = null) {
            const k = kernel || this.kernel;
            if (!k) {
                throw new Error('Kernel required for regularized equilibrium');
            }
            const rhoQM = this.psi.probabilityDensity;
            return k.convolve1D(rhoQM, this.dx);
        }

        /**
         * Compute coarse-grained H-function (relative entropy)
         * H̄ = ∫ ρ̄ ln(ρ̄/|ψ̄|²) dx
         * 
         * Uses coarse-graining with cell size cellSize for numerical stability.
         * 
         * @param {number} cellSize - Coarse-graining cell size (in grid units)
         * @returns {number} H-function value (≥ 0, = 0 at equilibrium)
         */
        hFunction(cellSize = 4) {
            const rhoQM = this.psi.probabilityDensity;
            const n = this.rho.length;

            // Coarse-grain by averaging over cells
            const nCells = Math.floor(n / cellSize);
            let H = 0;

            for (let c = 0; c < nCells; c++) {
                let rhoBar = 0;
                let rhoQMBar = 0;

                for (let j = 0; j < cellSize; j++) {
                    const idx = c * cellSize + j;
                    rhoBar += this.rho[idx];
                    rhoQMBar += rhoQM[idx];
                }

                rhoBar /= cellSize;
                rhoQMBar /= cellSize;

                // Avoid log(0)
                if (rhoBar > EPSILON && rhoQMBar > EPSILON) {
                    H += rhoBar * log(rhoBar / rhoQMBar) * this.dx * cellSize;
                }
            }

            return Math.max(0, H);  // H ≥ 0 by Gibbs inequality
        }

        /**
         * Fine-grained H-function (no coarse-graining, more sensitive)
         * H = ∫ ρ ln(ρ/|ψ|²) dx
         * 
         * @returns {number} H-function value
         */
        hFunctionFine() {
            const rhoQM = this.psi.probabilityDensity;
            const n = this.rho.length;
            let H = 0;

            for (let i = 0; i < n; i++) {
                if (this.rho[i] > EPSILON && rhoQM[i] > EPSILON) {
                    H += this.rho[i] * log(this.rho[i] / rhoQM[i]) * this.dx;
                }
            }

            return Math.max(0, H);
        }

        /**
         * Evolve distribution under continuity equation
         * ∂ρ/∂t + ∇·(ρv) = 0
         * 
         * @param {PilotWaveSystem} system - The pilot-wave system (provides velocity)
         * @param {number} dt - Time step
         */
        evolve(system, dt) {
            const v = system.regularizedVelocity();
            const n = this.rho.length;
            const newRho = new Array(n);

            // Upwind scheme for conservation
            for (let i = 1; i < n - 1; i++) {
                // Flux: F_i+1/2 = ρ_{upwind} * v_{i+1/2}
                const vHalfPlus = 0.5 * (v[i] + v[i + 1]);
                const vHalfMinus = 0.5 * (v[i] + v[i - 1]);

                // Upwind selection
                const fluxPlus = vHalfPlus > 0
                    ? this.rho[i] * vHalfPlus
                    : this.rho[i + 1] * vHalfPlus;
                const fluxMinus = vHalfMinus > 0
                    ? this.rho[i - 1] * vHalfMinus
                    : this.rho[i] * vHalfMinus;

                // ∂ρ/∂t = -(∂(ρv)/∂x)
                newRho[i] = this.rho[i] - dt * (fluxPlus - fluxMinus) / this.dx;
            }

            // Boundary conditions (zero flux)
            newRho[0] = this.rho[0];
            newRho[n - 1] = this.rho[n - 1];

            // Ensure non-negative
            this.rho = newRho.map(r => Math.max(0, r));
            this._normalizeDistribution();
        }

        /**
         * Simulate relaxation toward equilibrium
         * @param {PilotWaveSystem} system - The pilot-wave system
         * @param {number} dt - Time step
         * @param {number} nSteps - Number of steps
         * @param {boolean} useFine - Use fine-grained H (more sensitive)
         * @returns {number[]} H-function values at each step
         */
        relaxation(system, dt, nSteps, useFine = false) {
            const hFunc = useFine ? () => this.hFunctionFine() : () => this.hFunction();
            const hValues = [hFunc()];

            for (let step = 0; step < nSteps; step++) {
                this.evolve(system, dt);
                hValues.push(hFunc());
            }

            return hValues;
        }

        /**
         * Create ensemble in quantum equilibrium: ρ = |ψ|²
         * @param {WaveFunction} psi - Wavefunction
         * @param {number} dx - Grid spacing
         * @returns {QuantumEnsemble}
         */
        static equilibrium(psi, dx = 0.01) {
            const rhoQM = psi.probabilityDensity;
            return new QuantumEnsemble([...rhoQM], psi, dx);
        }

        /**
         * Create ensemble in nonequilibrium (uniform distribution)
         * @param {WaveFunction} psi - Wavefunction
         * @param {number} dx - Grid spacing
         * @returns {QuantumEnsemble}
         */
        static uniform(psi, n, dx = 0.01) {
            const rho = new Array(n).fill(1.0 / (n * dx));
            return new QuantumEnsemble(rho, psi, dx);
        }

        toString() {
            const H = this.hFunction();
            const eq = this.isInEquilibrium() ? 'EQUILIBRIUM' : 'NONEQUILIBRIUM';
            return `QuantumEnsemble(n=${this.rho.length}, H=${H.toFixed(6)}, ${eq})`;
        }
    }

    // ============================================================================
    // ACTION-PHASE BRIDGE (Connection to Contact Geometry)
    // ============================================================================

    /**
     * ActionPhaseBridge: Connects pilot-wave phase S to contact geometry action A
     * 
     * In the 1-jet bundle J¹(Q):
     * - Base coordinates: x^a
     * - Action/potential: A (fiber coordinate)
     * - Momenta: p_a = ∂A/∂x^a
     * 
     * The pilot wave ψ = R·exp(iS/ℏ) gives:
     * - Phase S corresponds to action A
     * - de Broglie momentum p = ∇S = ℏk
     * 
     * Legendrian condition: α|_L = 0 where α = dA - p_a dx^a
     * This is automatically satisfied when p = ∇S!
     */
    class ActionPhaseBridge {
        /**
         * @param {WaveFunction} psi - Pilot wave
         * @param {Object} contactManifold - Optional, from contact-thermodynamics index.js
         */
        constructor(psi, contactManifold = null) {
            this.psi = psi;
            this.manifold = contactManifold;
        }

        /**
         * Get action A = ℏ·S from phase
         * @param {number|number[]} phase - Phase S
         * @returns {number|number[]} Action A
         */
        phaseToAction(phase = null) {
            const S = phase || this.psi.phase;
            const hbar = this.psi.hbar;

            if (typeof S === 'number') {
                return hbar * S;
            }
            return S.map(s => hbar * s);
        }

        /**
         * Get canonical momentum p = ∇S = (1/ℏ)∇A
         * @param {number} dx - Grid spacing
         * @returns {number[]} Momentum field
         */
        momentum(dx = 0.01) {
            const S = this.psi.phase;
            if (typeof S === 'number') {
                throw new Error('Need field for momentum gradient');
            }

            const n = S.length;
            const p = new Array(n);

            for (let i = 1; i < n - 1; i++) {
                p[i] = (S[i + 1] - S[i - 1]) / (2 * dx);
            }
            p[0] = (S[1] - S[0]) / dx;
            p[n - 1] = (S[n - 1] - S[n - 2]) / dx;

            return p;
        }

        /**
         * Verify Legendrian condition: dA - p_a dx^a = 0
         * (i.e., p = ∇A = ∇S·ℏ)
         * 
         * @param {number} dx - Grid spacing
         * @returns {number} Maximum violation of Legendrian condition
         */
        legendrianViolation(dx = 0.01) {
            // For graphs of gradients, this should always be ~0
            // Just verify that p computed from S matches the "correct" value
            const p = this.momentum(dx);
            const A = this.phaseToAction();
            const hbar = this.psi.hbar;

            // ∇A/ℏ should equal ∇S
            const n = p.length;
            let maxViolation = 0;

            for (let i = 1; i < n - 1; i++) {
                const gradA = (A[i + 1] - A[i - 1]) / (2 * dx);
                const expected = gradA / hbar;
                const violation = Math.abs(p[i] - expected);
                if (violation > maxViolation) maxViolation = violation;
            }

            return maxViolation;
        }

        /**
         * Hamilton-Jacobi: H(x, ∂S/∂x) = 0 determines S
         * For free particle: (∂S/∂t) + (1/2m)|∇S|² = 0
         * 
         * @param {Function} hamiltonian - H(x, p) function
         * @param {number} dx - Grid spacing
         * @returns {number[]} H evaluated at each point (should be ~0 on shell)
         */
        hamiltonJacobi(hamiltonian, dx = 0.01) {
            const p = this.momentum(dx);
            const n = p.length;
            const H = new Array(n);

            for (let i = 0; i < n; i++) {
                const x = i * dx;
                H[i] = hamiltonian(x, p[i]);
            }

            return H;
        }

        toString() {
            return `ActionPhaseBridge(ℏ=${this.psi.hbar})`;
        }
    }

    // ============================================================================
    // CURVED SPACE PILOT-WAVE (Riemannian GA Integration)
    // ============================================================================

    /**
     * CurvedSpacePilotWave: Pilot-wave dynamics on curved Riemannian manifolds
     * 
     * Extends flat-space pilot-wave to curved geometries using:
     * - Connection bivector ω for parallel transport of spinor frames
     * - Curved-space guidance equation: v^i = (1/m) g^{ij} ∂_j S
     * - Geodesic ball regularization for finite velocity at nodes
     * 
     * Requires riemannian-ga.js module for ConnectionBivector and manifolds.
     */
    class CurvedSpacePilotWave {
        /**
         * @param {WaveFunction} psi - Wavefunction (phase on grid points)
         * @param {RiemannianManifold} manifold - The curved manifold
         * @param {SmearingKernel} kernel - Regularization kernel
         * @param {Object} options - Configuration options
         */
        constructor(psi, manifold, kernel, options = {}) {
            this.psi = psi;
            this.manifold = manifold;
            this.kernel = kernel;
            this.mass = options.mass || 1;
            this.hbar = options.hbar || 1;

            // Lazy-load connection if available
            this._connection = null;
        }

        /**
         * Get or create connection bivector
         * @returns {ConnectionBivector|null}
         */
        get connection() {
            if (this._connection === null) {
                // Try to load ConnectionBivector from riemannian-ga
                try {
                    const RGA = require('./riemannian-ga.js');
                    this._connection = new RGA.ConnectionBivector(this.manifold);
                } catch (e) {
                    // Riemannian GA not available
                    this._connection = false;
                }
            }
            return this._connection || null;
        }

        /**
         * Compute curved-space de Broglie velocity
         * 
         * v^i = (ℏ/m) g^{ij}(x) ∂_j S
         * 
         * @param {number[]} coords - Position on manifold
         * @param {number} h - Step size for numerical derivative
         * @returns {number[]} Contravariant velocity vector v^i
         */
        curvedDeBroglieVelocity(coords, h = 1e-6) {
            const dim = coords.length;
            const gInv = this.manifold.metricInverse(coords);

            // Compute ∂_j S numerically
            const gradS = this._phaseGradient(coords, h);

            // v^i = (ℏ/m) g^{ij} ∂_j S
            const prefactor = this.hbar / this.mass;
            const v = new Array(dim).fill(0);

            for (let i = 0; i < dim; i++) {
                for (let j = 0; j < dim; j++) {
                    v[i] += gInv[i][j] * gradS[j];
                }
                v[i] *= prefactor;
            }

            return v;
        }

        /**
         * Compute phase gradient at coordinates
         * @private
         */
        _phaseGradient(coords, h) {
            const dim = coords.length;
            const gradS = new Array(dim);

            // Get phase at coords (need to interpolate from grid)
            const S0 = this._interpolatePhase(coords);

            for (let i = 0; i < dim; i++) {
                const coordsPlus = [...coords];
                const coordsMinus = [...coords];
                coordsPlus[i] += h;
                coordsMinus[i] -= h;

                const Sp = this._interpolatePhase(coordsPlus);
                const Sm = this._interpolatePhase(coordsMinus);
                gradS[i] = (Sp - Sm) / (2 * h);
            }

            return gradS;
        }

        /**
         * Interpolate phase from grid to arbitrary coordinates
         * @private
         */
        _interpolatePhase(coords) {
            // For now, assume phase is defined at coords directly (grid-based)
            // In full implementation, use bilinear/trilinear interpolation
            const phase = this.psi.phase;
            if (typeof phase === 'number') {
                return phase;
            }
            // Simple: return phase at nearest grid point (placeholder)
            const idx = Math.floor(coords[0] / 0.1 + phase.length / 2);
            return phase[Math.max(0, Math.min(idx, phase.length - 1))];
        }

        /**
         * Parallel transport of spinor phase (bivector) along a path
         * 
         * Solves: dS/dλ + ω(v) × S = 0
         * 
         * Uses the connection bivector to rotate the phase.
         * 
         * @param {Object} spinorPhase - Bivector representing spinor phase
         * @param {number[]} startCoords - Starting point
         * @param {number[]} tangentDir - Direction of transport
         * @param {number} distance - Geodesic distance to transport
         * @param {number} nSteps - Number of discretization steps
         * @returns {Object} Transported spinor phase (Bivector)
         */
        parallelTransport(spinorPhase, startCoords, tangentDir, distance = 1, nSteps = 100) {
            if (!this.connection) {
                throw new Error('ConnectionBivector required for parallel transport');
            }

            const dt = distance / nSteps;
            let coords = [...startCoords];
            let S = spinorPhase.clone ? spinorPhase.clone() : { ...spinorPhase };

            for (let step = 0; step < nSteps; step++) {
                // Get connection along tangent direction
                const omega = this.connection.along(coords, tangentDir);

                // dS = -ω × S dt (commutator product rotates spinor)
                if (omega && typeof omega.commutatorWithBivector === 'function') {
                    const dS = omega.commutatorWithBivector(S).scale(-dt);
                    S = S.add(dS);
                }

                // Move along geodesic (simplified: linear for now)
                for (let i = 0; i < coords.length; i++) {
                    coords[i] += tangentDir[i] * dt;
                }
            }

            return S;
        }

        /**
         * Regularized velocity on curved space
         * 
         * Uses geodesic ball average instead of Euclidean convolution.
         * 
         * @param {number[]} coords - Position on manifold  
         * @param {number} nSamples - Monte Carlo samples in geodesic ball
         * @returns {number[]} Regularized velocity vector
         */
        regularizedVelocityCurved(coords, nSamples = 50) {
            const epsilon = this.kernel.width;
            const dim = coords.length;

            // Monte Carlo integration over geodesic ball B_ε(x)
            let jSum = new Array(dim).fill(0);
            let rhoSum = 0;
            let totalWeight = 0;

            for (let s = 0; s < nSamples; s++) {
                // Random direction on unit sphere
                const dir = this._randomUnitVector(dim);
                // Random distance (weighted by r^(dim-1) for uniform in ball)
                const r = epsilon * Math.pow(Math.random(), 1 / dim);

                // Nearby point (geodesic ball approximated as tangent)
                const nearby = coords.map((c, i) => c + r * dir[i]);

                // Kernel weight (Gaussian in geodesic distance)
                const weight = this.kernel.evaluate(r);

                // Accumulate current and density
                const v = this.curvedDeBroglieVelocity(nearby);
                const rho = this._interpolateDensity(nearby);

                for (let i = 0; i < dim; i++) {
                    jSum[i] += rho * v[i] * weight;
                }
                rhoSum += rho * weight;
                totalWeight += weight;
            }

            // Regularized: v_reg = j_reg / rho_reg
            const vReg = new Array(dim);
            for (let i = 0; i < dim; i++) {
                vReg[i] = jSum[i] / Math.max(rhoSum, EPSILON);
            }

            return vReg;
        }

        /**
         * Random unit vector in n dimensions
         * @private
         */
        _randomUnitVector(dim) {
            const v = new Array(dim);
            let norm = 0;
            for (let i = 0; i < dim; i++) {
                // Box-Muller for Gaussian components
                const u1 = Math.random();
                const u2 = Math.random();
                v[i] = sqrt(-2 * log(u1)) * cos(2 * PI * u2);
                norm += v[i] * v[i];
            }
            norm = sqrt(norm);
            return v.map(x => x / norm);
        }

        /**
         * Interpolate density at coordinates
         * @private
         */
        _interpolateDensity(coords) {
            const rho = this.psi.probabilityDensity;
            if (typeof rho === 'number') {
                return rho;
            }
            // Simple: nearest grid point
            const idx = Math.floor(coords[0] / 0.1 + rho.length / 2);
            return rho[Math.max(0, Math.min(idx, rho.length - 1))];
        }

        /**
         * Check if wavefunction has a node near given coordinates
         * @param {number[]} coords - Position to check
         * @param {number} threshold - Density threshold for node detection
         * @returns {boolean}
         */
        hasNodeNear(coords, threshold = 1e-6) {
            const rho = this._interpolateDensity(coords);
            return rho < threshold;
        }

        /**
         * Get regularized density at coordinates
         * @param {number[]} coords
         * @param {number} nSamples
         * @returns {number}
         */
        regularizedDensityCurved(coords, nSamples = 50) {
            const epsilon = this.kernel.width;
            const dim = coords.length;

            let rhoSum = 0;
            let totalWeight = 0;

            for (let s = 0; s < nSamples; s++) {
                const dir = this._randomUnitVector(dim);
                const r = epsilon * Math.pow(Math.random(), 1 / dim);
                const nearby = coords.map((c, i) => c + r * dir[i]);
                const weight = this.kernel.evaluate(r);

                rhoSum += this._interpolateDensity(nearby) * weight;
                totalWeight += weight;
            }

            return rhoSum / Math.max(totalWeight, EPSILON);
        }

        toString() {
            return `CurvedSpacePilotWave(manifold=${this.manifold.constructor.name}, m=${this.mass})`;
        }
    }

    // ============================================================================
    // EXPORTS
    // ============================================================================

    const PilotWave = {
        SmearingKernel,
        WaveFunction,
        PilotWaveSystem,
        QuantumEnsemble,
        ActionPhaseBridge,
        CurvedSpacePilotWave,

        // Constants
        EPSILON,

        // Utility: create standard hydrogen-like node
        createNodedWavefunction: function (nodePosition, options = {}) {
            // ψ(x) ~ (x - x_node) near node
            const n = options.gridSize || 100;
            const dx = options.dx || 0.1;
            const x0 = nodePosition;

            const amplitude = [];
            const phase = [];

            for (let i = 0; i < n; i++) {
                const x = (i - n / 2) * dx;
                // Linear near node: ψ ~ (x - x0)
                const dist = x - x0;
                amplitude.push(Math.abs(dist));
                // Phase changes sign at node
                phase.push(dist >= 0 ? 0 : PI);
            }

            return new WaveFunction(amplitude, phase, options);
        }
    };

    // Export for Node.js and browser
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = PilotWave;
    }
    if (typeof global !== 'undefined') {
        global.PilotWave = PilotWave;
    }
    if (typeof window !== 'undefined') {
        window.PilotWave = PilotWave;
    }

})(typeof global !== 'undefined' ? global : (typeof window !== 'undefined' ? window : this));
