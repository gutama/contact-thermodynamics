/**
 * GA-Native Pilot-Wave Guidance Module (Phase 2)
 *
 * Bridges the geometric-algebra engine (`src/algebra/multivector.js`) and the
 * Riemannian-GA layer (`src/geometry/riemannian-ga.js`) with the pilot-wave
 * dynamics in `src/physics/pilot-wave.js`. It re-expresses the de Broglie–Bohm
 * guidance equation and Valentini/Bell node-regularization using GA's
 * rotor/bivector representation of the wavefunction phase, rather than the
 * conventional `Im(∇ψ/ψ)` / `∇S` vector calculus.
 *
 * This targets research question RQ2 of the research plan: does the GA/STA
 * multivector form of the guidance equation give a more natural or more
 * numerically stable regularization near wavefunction nodes, and does it extend
 * cleanly toward the curved-spacetime case?
 *
 * Everything here is *additive*: `pilot-wave.js` (SmearingKernel, WaveFunction,
 * PilotWaveSystem, QuantumEnsemble, CurvedSpacePilotWave) is untouched, and this
 * module reuses its `SmearingKernel` and `QuantumEnsemble` verbatim.
 *
 * ── The GA construction ─────────────────────────────────────────────────────
 * A complex wavefunction ψ = R e^{iS} is represented, following the
 * Hestenes/Doran–Lasenby GA treatment of the Pauli/Schrödinger equation, as an
 * even-grade multivector ("spinor") in the plane algebra Cl(2,0,0):
 *
 *     i  ↦  I = e₁∧e₂        (unit pseudoscalar of the plane, I² = −1)
 *     ψ  ↦  ψ_r + I ψ_i  =  R · U,   U = exp(I S)   (a unit rotor)
 *
 * The guidance velocity is then the **grade-1 part of the phase rotor's
 * logarithmic derivative, de-rotated by I⁻¹**:
 *
 *     v  =  (ħ/m) ⟨ (∇U) U⁻¹ I⁻¹ ⟩₁,        U = ψ/|ψ| = exp(I S).
 *
 * This is provably equal to the standard v = (ħ/m) ∇S in flat, non-relativistic
 * space: since ∇U = (∇S) I U (the ∇S here is the ordinary gradient *vector*
 * Σ_k e_k ∂_k S), we get (∇U)U⁻¹ = (∇S) I, and right-multiplying by I⁻¹ recovers
 * (∇S) I I⁻¹ = ∇S exactly. So analytically the GA form carries **no new
 * predictive content** in flat space — it is the same velocity field.
 *
 * The numerically interesting point (RQ2) is *how* the two forms compute that
 * field. The conventional route differences the scalar phase S = atan2(ψ_i,ψ_r),
 * which has a 2π branch cut and is singular on a whole codimension-1 locus
 * emanating from any phase vortex (e.g. ψ ∝ (x+iy)e^{−r²}). The GA route
 * differences the **unit rotor U = ψ/|ψ|**, whose components (cos S, sin S) are
 * smooth everywhere R ≠ 0 — so it never forms the wrapped-angle gradient and
 * avoids the spurious branch-cut blow-up. Both still (correctly) diverge as
 * R → 0 at the node itself, which is what the smearing regularization tames.
 *
 * Node regularization is the same Valentini/Bell construction as
 * `PilotWaveSystem.regularizedVelocity` — v_reg = (j∗μ)/(ρ∗μ) — but the current
 * j = ρ v is built from the GA rotor velocity above and the reused
 * `SmearingKernel` convolution.
 *
 * ── Alternative considered and rejected ─────────────────────────────────────
 * The "full-spinor" current j = (ħ/m)⟨ψ̃ ∇ψ I⁻¹⟩ using the *un-normalized*
 * spinor ψ = ψ_r + I ψ_i (mirroring Im(ψ*∇ψ) literally). Rejected because in the
 * 2D plane algebra both a real vector `a` and `I·b` (pseudoscalar × vector) are
 * grade-1, so the amplitude-gradient term ∇R and the phase-gradient term ρ∇S do
 * *not* separate by a clean grade projection — extracting the phase part
 * requires exactly the same reversion/normalization bookkeeping that the unit
 * rotor U = ψ/|ψ| performs up front. Normalizing to U first is simpler, removes
 * the ∇R contamination by construction, and (the whole point) is what differences
 * the smooth rotor rather than the wrapped angle. The un-normalized route offers
 * no compensating advantage here, though it is the natural object once one adds
 * the quantum-potential / amplitude dynamics (flagged for a later phase).
 *
 * @module physics/ga-pilot-wave
 * @license MIT
 */

(function (global) {
    'use strict';

    let GA, PilotWave;
    if (typeof module !== 'undefined' && module.exports) {
        GA = require('../algebra/multivector.js');
        PilotWave = require('./pilot-wave.js');
    } else {
        GA = global.GA;
        PilotWave = global.PilotWave;
    }
    const { Algebra } = GA;
    const { SmearingKernel } = PilotWave;

    const GA_EPSILON = 1e-12;
    const cos = Math.cos;
    const sin = Math.sin;
    const sqrt = Math.sqrt;

    // The plane algebra Cl(2,0,0). Its unit pseudoscalar I = e₁∧e₂ (I² = −1)
    // plays the role of the complex imaginary unit i (Hestenes identification).
    // Reused across all instances; Algebra is immutable after construction.
    const PLANE = new Algebra(2, 0, 0);
    const I2 = PLANE.e(1).mul(PLANE.e(2));   // unit phase bivector
    const I2_INV = I2.inverse();             // = −I2

    /**
     * Read the coefficient of basis vector e_i (1-indexed) from a multivector.
     * @param {Object} mv - Multivector
     * @param {number} i - 1-indexed basis vector
     * @returns {number}
     */
    function vectorComponent(mv, i) {
        return mv.coeffs[1 << (i - 1)] || 0;
    }

    /**
     * The phase rotor U = exp(I·S) = cos S + I sin S (a unit rotor in Cl(2)).
     * @param {number} S - Phase
     * @returns {Object} Multivector (even grade)
     */
    function phaseRotor(S) {
        return PLANE.scalar(cos(S)).add(I2.scale(sin(S)));
    }

    /**
     * Unit rotor built directly from Cartesian components: U = (ψ_r + I ψ_i)/R.
     * This is the smooth object the GA guidance differentiates (rather than the
     * wrapped phase angle). Returns the zero multivector at an exact node.
     * @param {number} re - ψ_r
     * @param {number} im - ψ_i
     * @returns {Object} Multivector
     */
    function unitRotorFromComponents(re, im) {
        const R = sqrt(re * re + im * im);
        if (R < GA_EPSILON) return PLANE.zero();
        return PLANE.scalar(re / R).add(I2.scale(im / R));
    }

    /**
     * Core GA guidance operation at one point, given the spatial gradient of the
     * unit rotor as an array of ∂_k U (one Multivector per spatial axis) plus the
     * rotor value U itself.
     *
     *     v_k = (ħ/m) · e_k-component of  ⟨ (Σ_k e_k ∂_k U) U⁻¹ I⁻¹ ⟩₁
     *
     * @param {Object} U - Unit rotor at the point
     * @param {Object[]} dU - [∂₁U, ∂₂U, ...] (one Multivector per axis)
     * @param {number} prefactor - ħ/m
     * @returns {number[]} Velocity components (contravariant, coordinate basis)
     */
    function guidanceFromRotor(U, dU, prefactor) {
        const dim = dU.length;
        if (dim > 2) {
            throw new Error(`Unsupported spatial dimension: ${dim}. guidanceFromRotor uses the 2-D plane algebra Cl(2,0,0), which only defines basis vectors e1, e2; a higher-dimensional Clifford algebra is required to guide in ${dim} dimensions.`);
        }
        if (U.isZero()) return new Array(dim).fill(0);
        const Uinv = U.inverse();

        // ∇U = Σ_k e_k (∂_k U)
        let gradU = PLANE.zero();
        for (let k = 0; k < dim; k++) {
            gradU = gradU.add(PLANE.e(k + 1).mul(dU[k]));
        }

        // w = (∇U) U⁻¹ I⁻¹ ; its grade-1 part holds ∂_k S in component k.
        const w = gradU.mul(Uinv).mul(I2_INV).grade(1);

        const v = new Array(dim);
        for (let k = 0; k < dim; k++) {
            v[k] = prefactor * vectorComponent(w, k + 1);
        }
        return v;
    }

    // ============================================================================
    // 1D GA PILOT-WAVE SYSTEM
    // ============================================================================

    /**
     * GAPilotWaveSystem1D: GA-native de Broglie guidance on a 1D grid.
     *
     * Drop-in compatible with the interface `QuantumEnsemble` expects
     * (`regularizedVelocity()` returns a number[]), so the existing H-theorem
     * machinery can be driven by the GA velocity field unchanged.
     */
    class GAPilotWaveSystem1D {
        /**
         * @param {number[]} amplitude - Amplitude R on the grid
         * @param {number[]} phase - Phase S on the grid
         * @param {Object} options - { dx, mass, hbar, kernel }
         */
        constructor(amplitude, phase, options = {}) {
            this.amplitude = amplitude;
            this.phase = phase;
            this.dx = options.dx || 0.01;
            this.mass = options.mass || 1.0;
            this.hbar = options.hbar || 1.0;
            this.kernel = options.kernel || new SmearingKernel(1e-10, 'gaussian');

            this._rhoReg = null;
            this._jReg = null;
        }

        /**
         * Build the unit-rotor field U_i from the (amplitude, phase) grid.
         * @private
         */
        _rotorField() {
            const n = this.phase.length;
            const U = new Array(n);
            for (let i = 0; i < n; i++) {
                U[i] = phaseRotor(this.phase[i]);
            }
            return U;
        }

        /**
         * GA de Broglie velocity v = (ħ/m) ⟨(∇U)U⁻¹I⁻¹⟩₁ on the grid.
         * Diverges at nodes (as it should) — regularize with regularizedVelocity.
         * @returns {number[]}
         */
        deBroglieVelocity() {
            const U = this._rotorField();
            const n = U.length;
            const dx = this.dx;
            const prefactor = this.hbar / this.mass;
            const v = new Array(n).fill(0);

            for (let i = 0; i < n; i++) {
                const ip = Math.min(i + 1, n - 1);
                const im = Math.max(i - 1, 0);
                const denom = (ip - im) * dx;
                const dU = U[ip].sub(U[im]).scale(1 / denom);
                v[i] = guidanceFromRotor(U[i], [dU], prefactor)[0];
            }
            return v;
        }

        /**
         * GA probability current j = ρ v (built from the rotor velocity).
         * @returns {number[]}
         */
        current() {
            const v = this.deBroglieVelocity();
            const n = v.length;
            const j = new Array(n);
            for (let i = 0; i < n; i++) {
                const R = this.amplitude[i];
                j[i] = R * R * v[i];
            }
            return j;
        }

        /**
         * Valentini-regularized velocity v_reg = (j∗μ)/(ρ∗μ), finite at nodes.
         * Uses the reused SmearingKernel convolution.
         * @returns {number[]}
         */
        regularizedVelocity() {
            const j = this.current();
            const rho = this.amplitude.map(R => R * R);

            const jReg = this.kernel.convolve1D(j, this.dx);
            const rhoReg = this.kernel.convolve1D(rho, this.dx);

            this._jReg = jReg;
            this._rhoReg = rhoReg;

            const n = jReg.length;
            const vReg = new Array(n);
            for (let i = 0; i < n; i++) {
                vReg[i] = jReg[i] / Math.max(rhoReg[i], GA_EPSILON);
            }
            return vReg;
        }

        /**
         * Regularized density (|ψ|²)_reg.
         * @returns {number[]}
         */
        regularizedDensity() {
            if (this._rhoReg === null) {
                const rho = this.amplitude.map(R => R * R);
                this._rhoReg = this.kernel.convolve1D(rho, this.dx);
            }
            return this._rhoReg;
        }

        toString() {
            return `GAPilotWaveSystem1D(n=${this.phase.length}, dx=${this.dx})`;
        }
    }

    // ============================================================================
    // 2D GA PILOT-WAVE FIELD (phase vortex / node comparison)
    // ============================================================================

    /**
     * GAPilotWaveField2D: GA guidance on a 2D grid, built to expose the
     * phase-vortex near-node behavior where the GA form differs from the
     * conventional one.
     *
     * The wavefunction is supplied by its Cartesian components ψ_r(x,y),
     * ψ_i(x,y) sampled on a uniform grid. Two velocity fields are provided:
     *
     *  - deBroglieVelocityGA()   — GA rotor form: differences U = ψ/|ψ|.
     *  - deBroglieVelocityPhase() — conventional form: differences the
     *    wrapped phase angle S = atan2(ψ_i, ψ_r). This is the 2D analogue of
     *    what `pilot-wave.js` does in 1D, included as the honest baseline.
     */
    class GAPilotWaveField2D {
        /**
         * @param {number[][]} psiRe - ψ_r[i][j]
         * @param {number[][]} psiIm - ψ_i[i][j]
         * @param {Object} options - { dx, dy, mass, hbar, kernel }
         */
        constructor(psiRe, psiIm, options = {}) {
            this.psiRe = psiRe;
            this.psiIm = psiIm;
            this.nx = psiRe.length;
            this.ny = psiRe[0].length;
            this.dx = options.dx || 0.1;
            this.dy = options.dy || this.dx;
            this.mass = options.mass || 1.0;
            this.hbar = options.hbar || 1.0;
            this.kernel = options.kernel || new SmearingKernel(0.2, 'gaussian');
        }

        density() {
            const rho = [];
            for (let i = 0; i < this.nx; i++) {
                rho[i] = [];
                for (let j = 0; j < this.ny; j++) {
                    rho[i][j] = this.psiRe[i][j] ** 2 + this.psiIm[i][j] ** 2;
                }
            }
            return rho;
        }

        /**
         * GA rotor velocity field. Returns { vx, vy } as 2D arrays.
         * @returns {{vx: number[][], vy: number[][]}}
         */
        deBroglieVelocityGA() {
            const { nx, ny, dx, dy } = this;
            const prefactor = this.hbar / this.mass;
            const U = [];
            for (let i = 0; i < nx; i++) {
                U[i] = [];
                for (let j = 0; j < ny; j++) {
                    U[i][j] = unitRotorFromComponents(this.psiRe[i][j], this.psiIm[i][j]);
                }
            }

            const vx = [], vy = [];
            for (let i = 0; i < nx; i++) {
                vx[i] = new Array(ny).fill(0);
                vy[i] = new Array(ny).fill(0);
                for (let j = 0; j < ny; j++) {
                    const ip = Math.min(i + 1, nx - 1), im = Math.max(i - 1, 0);
                    const jp = Math.min(j + 1, ny - 1), jm = Math.max(j - 1, 0);
                    const dUx = U[ip][j].sub(U[im][j]).scale(1 / ((ip - im) * dx));
                    const dUy = U[i][jp].sub(U[i][jm]).scale(1 / ((jp - jm) * dy));
                    const v = guidanceFromRotor(U[i][j], [dUx, dUy], prefactor);
                    vx[i][j] = v[0];
                    vy[i][j] = v[1];
                }
            }
            return { vx, vy };
        }

        /**
         * Conventional velocity field: central-difference the wrapped phase
         * angle S = atan2(ψ_i, ψ_r). This is the baseline that suffers the
         * branch-cut artifact at a phase vortex.
         * @returns {{vx: number[][], vy: number[][]}}
         */
        deBroglieVelocityPhase() {
            const { nx, ny, dx, dy } = this;
            const prefactor = this.hbar / this.mass;
            const S = [];
            for (let i = 0; i < nx; i++) {
                S[i] = [];
                for (let j = 0; j < ny; j++) {
                    S[i][j] = Math.atan2(this.psiIm[i][j], this.psiRe[i][j]);
                }
            }

            const vx = [], vy = [];
            for (let i = 0; i < nx; i++) {
                vx[i] = new Array(ny).fill(0);
                vy[i] = new Array(ny).fill(0);
                for (let j = 0; j < ny; j++) {
                    const ip = Math.min(i + 1, nx - 1), im = Math.max(i - 1, 0);
                    const jp = Math.min(j + 1, ny - 1), jm = Math.max(j - 1, 0);
                    vx[i][j] = prefactor * (S[ip][j] - S[im][j]) / ((ip - im) * dx);
                    vy[i][j] = prefactor * (S[i][jp] - S[i][jm]) / ((jp - jm) * dy);
                }
            }
            return { vx, vy };
        }

        /**
         * Valentini-regularized GA velocity field (componentwise 2D smearing).
         * @returns {{vx: number[][], vy: number[][]}}
         */
        regularizedVelocityGA() {
            const { vx, vy } = this.deBroglieVelocityGA();
            const rho = this.density();
            const jx = [], jy = [];
            for (let i = 0; i < this.nx; i++) {
                jx[i] = []; jy[i] = [];
                for (let j = 0; j < this.ny; j++) {
                    jx[i][j] = rho[i][j] * vx[i][j];
                    jy[i][j] = rho[i][j] * vy[i][j];
                }
            }
            const jxReg = this.kernel.convolve2D(jx, this.dx, this.dy);
            const jyReg = this.kernel.convolve2D(jy, this.dx, this.dy);
            const rhoReg = this.kernel.convolve2D(rho, this.dx, this.dy);

            const vxReg = [], vyReg = [];
            for (let i = 0; i < this.nx; i++) {
                vxReg[i] = []; vyReg[i] = [];
                for (let j = 0; j < this.ny; j++) {
                    const d = Math.max(rhoReg[i][j], GA_EPSILON);
                    vxReg[i][j] = jxReg[i][j] / d;
                    vyReg[i][j] = jyReg[i][j] / d;
                }
            }
            return { vx: vxReg, vy: vyReg };
        }

        toString() {
            return `GAPilotWaveField2D(${this.nx}x${this.ny})`;
        }
    }

    // ============================================================================
    // CURVED-SPACE GA GUIDANCE
    // ============================================================================

    /**
     * GACurvedPilotWave: GA guidance velocity on a curved Riemannian manifold.
     *
     *     v^i = (ħ/m) g^{ij}(x) ∂_j S
     *
     * where the covariant phase gradient ∂_j S is extracted by the *same* GA
     * rotor operation used in flat space (guidanceFromRotor), and the index is
     * then raised with the manifold's inverse metric g^{ij}. When g^{ij} = δ^{ij}
     * (flat) this reduces exactly to the flat GA guidance — verified in the
     * tests — which is the natural bridge toward a future relativistic/
     * curved-spacetime phase (Phase 4).
     *
     * The phase is supplied as a function S(coords) so gradients can be taken at
     * arbitrary points (not just grid nodes). Composes with the existing
     * `ConnectionBivector` from `src/geometry/riemannian-ga.js` for parallel
     * transport of the phase bivector.
     */
    class GACurvedPilotWave {
        /**
         * @param {Function} phaseFn - S(coords): scalar phase field
         * @param {Object} manifold - RiemannianManifold (needs metricInverse)
         * @param {Object} options - { mass, hbar, kernel }
         */
        constructor(phaseFn, manifold, options = {}) {
            this.phaseFn = phaseFn;
            this.manifold = manifold;
            this.mass = options.mass || 1.0;
            this.hbar = options.hbar || 1.0;
            this.kernel = options.kernel || new SmearingKernel(0.2, 'gaussian');
            this._connection = null;
        }

        /**
         * Lazily build the connection bivector (mirrors CurvedSpacePilotWave).
         * @returns {Object|null}
         */
        get connection() {
            if (this._connection === null) {
                try {
                    const RGA = require('../geometry/riemannian-ga.js');
                    this._connection = new RGA.ConnectionBivector(this.manifold);
                } catch (e) {
                    if (typeof console !== 'undefined' && console.warn) {
                        console.warn('ga-pilot-wave: connection unavailable, flat fallback:', e.message);
                    }
                    this._connection = false;
                }
            }
            return this._connection || null;
        }

        /**
         * Extract the coordinate-basis phase gradient (∂_j S) at coords using the
         * GA rotor logarithmic derivative — identical operation in flat and
         * curved cases; only the index-raising below differs.
         * @param {number[]} coords
         * @param {number} h - Step size
         * @returns {number[]} ∂_j S (covariant components)
         * @private
         */
        _phaseGradientGA(coords, h) {
            const dim = coords.length;
            const U0 = phaseRotor(this.phaseFn(coords));
            const dU = new Array(dim);
            for (let k = 0; k < dim; k++) {
                const cp = [...coords]; cp[k] += h;
                const cm = [...coords]; cm[k] -= h;
                dU[k] = phaseRotor(this.phaseFn(cp))
                    .sub(phaseRotor(this.phaseFn(cm)))
                    .scale(1 / (2 * h));
            }
            // Prefactor 1 here: we want raw ∂_j S; the ħ/m enters the raised form.
            return guidanceFromRotor(U0, dU, 1);
        }

        /**
         * Flat GA guidance velocity: v^i = (ħ/m) ∂_i S (identity metric).
         * @param {number[]} coords
         * @param {number} h
         * @returns {number[]}
         */
        flatGuidanceVelocity(coords, h = 1e-6) {
            const gradS = this._phaseGradientGA(coords, h);
            const prefactor = this.hbar / this.mass;
            return gradS.map(g => prefactor * g);
        }

        /**
         * Curved GA guidance velocity: v^i = (ħ/m) g^{ij} ∂_j S.
         * @param {number[]} coords
         * @param {number} h
         * @returns {number[]}
         */
        curvedGuidanceVelocity(coords, h = 1e-6) {
            const dim = coords.length;
            const gradS = this._phaseGradientGA(coords, h);
            const gInv = this.manifold.metricInverse(coords);
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
         * Parallel transport of a phase bivector along a tangent direction,
         * solving dB/dλ + ω(v) × B = 0 with the connection bivector — reuses the
         * same `commutatorWithBivector` machinery as CurvedSpacePilotWave. With a
         * zero connection (flat), the bivector is returned unchanged.
         * @param {Object} phaseBivector - Bivector3D/Bivector4D
         * @param {number[]} startCoords
         * @param {number[]} tangentDir
         * @param {number} distance
         * @param {number} nSteps
         * @returns {Object} Transported bivector
         */
        parallelTransportPhase(phaseBivector, startCoords, tangentDir, distance = 1, nSteps = 100) {
            if (!this.connection) {
                throw new Error('ConnectionBivector required for parallel transport');
            }
            const dt = distance / nSteps;
            const coords = [...startCoords];
            let B = phaseBivector.clone ? phaseBivector.clone()
                : (phaseBivector.scale ? phaseBivector.scale(1) : { ...phaseBivector });

            for (let step = 0; step < nSteps; step++) {
                const omega = this.connection.along(coords, tangentDir);
                if (omega && typeof omega.commutatorWithBivector === 'function') {
                    const dB = omega.commutatorWithBivector(B).scale(-dt);
                    B = B.add(dB);
                }
                for (let i = 0; i < coords.length; i++) {
                    coords[i] += tangentDir[i] * dt;
                }
            }
            return B;
        }

        toString() {
            return `GACurvedPilotWave(manifold=${this.manifold.constructor.name})`;
        }
    }

    // ============================================================================
    // EXPORTS
    // ============================================================================

    const GAPilotWave = {
        GAPilotWaveSystem1D,
        GAPilotWaveField2D,
        GACurvedPilotWave,

        // Low-level GA guidance primitives (exposed for tests / reuse)
        phaseRotor,
        unitRotorFromComponents,
        guidanceFromRotor,
        vectorComponent,
        planeAlgebra: PLANE,
        phaseBivector: I2,

        GA_EPSILON
    };

    if (typeof module !== 'undefined' && module.exports) {
        module.exports = GAPilotWave;
    }
    if (typeof global !== 'undefined') {
        global.GAPilotWave = GAPilotWave;
    }
    if (typeof window !== 'undefined') {
        window.GAPilotWave = GAPilotWave;
    }

})(typeof global !== 'undefined' ? global : (typeof window !== 'undefined' ? window : this));
