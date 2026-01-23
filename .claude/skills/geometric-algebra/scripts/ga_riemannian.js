/**
 * Riemannian Geometry in Geometric Calculus - JavaScript Implementation
 * 
 * Coordinate-Free Riemannian Geometry via Geometric Calculus
 * 
 * This module provides pure GA formulations that eliminate Christoffel symbols entirely.
 * All curvature computations use:
 * - Connection bivector ω (replaces Γᵏᵢⱼ)
 * - Curvature 2-form Ω (replaces Rˡᵢⱼₖ)
 * - Commutator product × for bivector action
 * 
 * Master Equations:
 *     ∂ᵢeⱼ = ωᵢ × eⱼ           (frame rotation)
 *     Ω = dω + ω ∧ ω           (Cartan structure equation)
 *     ∇ᵤv = ∂ᵤv + ω(u) × v     (covariant derivative)
 *     v̇ + ω(v) × v = 0         (geodesic equation)
 */

(function (global) {
    'use strict';

    // ========================================================================
    // UTILITY FUNCTIONS
    // ========================================================================

    function dot(u, v) {
        let sum = 0;
        for (let i = 0; i < u.length; i++) {
            sum += u[i] * v[i];
        }
        return sum;
    }

    function cross3(u, v) {
        return [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0]
        ];
    }

    function norm(v) {
        return Math.sqrt(dot(v, v));
    }

    function normalize(v) {
        const n = norm(v);
        if (n < 1e-10) return v.slice();
        return v.map(x => x / n);
    }

    function vecAdd(u, v) {
        return u.map((x, i) => x + v[i]);
    }

    function vecSub(u, v) {
        return u.map((x, i) => x - v[i]);
    }

    function vecScale(v, s) {
        return v.map(x => x * s);
    }

    function matMul(A, B) {
        // Matrix multiplication for 2x2 or NxN matrices
        const n = A.length;
        const m = B[0].length;
        const p = B.length;
        const C = [];
        for (let i = 0; i < n; i++) {
            C[i] = [];
            for (let j = 0; j < m; j++) {
                C[i][j] = 0;
                for (let k = 0; k < p; k++) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return C;
    }

    function matVecMul(A, v) {
        const n = A.length;
        const result = [];
        for (let i = 0; i < n; i++) {
            result[i] = 0;
            for (let j = 0; j < v.length; j++) {
                result[i] += A[i][j] * v[j];
            }
        }
        return result;
    }

    function det2x2(A) {
        return A[0][0] * A[1][1] - A[0][1] * A[1][0];
    }

    function inv2x2(A) {
        const d = det2x2(A);
        if (Math.abs(d) < 1e-12) return null;
        return [
            [A[1][1] / d, -A[0][1] / d],
            [-A[1][0] / d, A[0][0] / d]
        ];
    }

    function trace2x2(A) {
        return A[0][0] + A[1][1];
    }

    function eigenvalues2x2(A) {
        // For 2x2 matrix, eigenvalues from characteristic polynomial
        const tr = trace2x2(A);
        const dt = det2x2(A);
        const disc = tr * tr - 4 * dt;
        if (disc < 0) {
            // Complex eigenvalues
            return [tr / 2, tr / 2];
        }
        const sqrtDisc = Math.sqrt(disc);
        return [(tr - sqrtDisc) / 2, (tr + sqrtDisc) / 2];
    }

    // ========================================================================
    // BIVECTOR CLASSES
    // ========================================================================

    /**
     * Bivector in 2D (single component: e₁₂)
     * In 2D, bivectors are pseudoscalars.
     */
    class Bivector2D {
        constructor(e12 = 0.0) {
            this.e12 = e12;
        }

        toString() {
            return `Bivector2D(${this.e12.toFixed(6)} e₁₂)`;
        }

        add(other) {
            return new Bivector2D(this.e12 + other.e12);
        }

        sub(other) {
            return new Bivector2D(this.e12 - other.e12);
        }

        scale(s) {
            return new Bivector2D(this.e12 * s);
        }

        neg() {
            return new Bivector2D(-this.e12);
        }

        get norm() {
            return Math.abs(this.e12);
        }

        /**
         * Compute ω × v = ½(ωv - vω) for 2D
         * For B = b e₁₂: B × (v¹e₁ + v²e₂) = b(v², -v¹)
         */
        commutatorWithVector(v) {
            return [this.e12 * v[1], -this.e12 * v[0]];
        }
    }

    /**
     * Bivector in 3D (three components: e₂₃, e₃₁, e₁₂)
     * Represents oriented planes, isomorphic to rotations.
     */
    class Bivector3D {
        constructor(e23 = 0.0, e31 = 0.0, e12 = 0.0) {
            this.e23 = e23;  // YZ plane
            this.e31 = e31;  // ZX plane
            this.e12 = e12;  // XY plane
        }

        static fromArray(arr) {
            return new Bivector3D(arr[0], arr[1], arr[2]);
        }

        static fromWedge(u, v) {
            return new Bivector3D(
                u[1] * v[2] - u[2] * v[1],  // e23
                u[2] * v[0] - u[0] * v[2],  // e31
                u[0] * v[1] - u[1] * v[0]   // e12
            );
        }

        toArray() {
            return [this.e23, this.e31, this.e12];
        }

        toString() {
            return `Bivector3D(${this.e23.toFixed(4)} e₂₃ + ${this.e31.toFixed(4)} e₃₁ + ${this.e12.toFixed(4)} e₁₂)`;
        }

        add(other) {
            return new Bivector3D(
                this.e23 + other.e23,
                this.e31 + other.e31,
                this.e12 + other.e12
            );
        }

        sub(other) {
            return new Bivector3D(
                this.e23 - other.e23,
                this.e31 - other.e31,
                this.e12 - other.e12
            );
        }

        scale(s) {
            return new Bivector3D(
                this.e23 * s,
                this.e31 * s,
                this.e12 * s
            );
        }

        neg() {
            return new Bivector3D(-this.e23, -this.e31, -this.e12);
        }

        get normSquared() {
            return this.e23 * this.e23 + this.e31 * this.e31 + this.e12 * this.e12;
        }

        get norm() {
            return Math.sqrt(this.normSquared);
        }

        /**
         * Compute B × v = ½(Bv - vB) in 3D
         * Equivalent to cross product with the dual vector
         */
        commutatorWithVector(v) {
            const b = [this.e23, this.e31, this.e12];
            return cross3(b, v);
        }

        /**
         * Compute B₁ × B₂ = ½(B₁B₂ - B₂B₁)
         * For bivectors in 3D, this gives another bivector (Lie bracket)
         */
        commutatorWithBivector(other) {
            const b1 = this.toArray();
            const b2 = other.toArray();
            const result = cross3(b1, b2);
            return Bivector3D.fromArray(result);
        }

        /**
         * Compute B₁ · B₂ (scalar inner product)
         */
        innerWithBivector(other) {
            return this.e23 * other.e23 + this.e31 * other.e31 + this.e12 * other.e12;
        }
    }

    // Generic Bivector for arbitrary dimensions (kept for compatibility)
    class Bivector {
        constructor(dim, B = null) {
            this.dim = dim;
            this.B = B || new Float64Array(dim * dim);
        }

        static fromWedge(u, v) {
            const dim = u.length;
            const B = new Float64Array(dim * dim);
            for (let i = 0; i < dim; i++) {
                for (let j = 0; j < dim; j++) {
                    B[i * dim + j] = u[i] * v[j] - v[i] * u[j];
                }
            }
            return new Bivector(dim, B);
        }

        get(i, j) {
            return this.B[i * this.dim + j];
        }

        set(i, j, val) {
            this.B[i * this.dim + j] = val;
        }

        commutator(v) {
            const result = new Float64Array(this.dim);
            for (let i = 0; i < this.dim; i++) {
                for (let j = 0; j < this.dim; j++) {
                    result[i] += this.get(i, j) * v[j];
                }
            }
            return Array.from(result);
        }

        add(other) {
            const result = new Bivector(this.dim);
            for (let i = 0; i < this.dim * this.dim; i++) {
                result.B[i] = this.B[i] + other.B[i];
            }
            return result;
        }

        scale(s) {
            const result = new Bivector(this.dim);
            for (let i = 0; i < this.dim * this.dim; i++) {
                result.B[i] = s * this.B[i];
            }
            return result;
        }
    }

    // ========================================================================
    // TANGENT FRAME
    // ========================================================================

    /**
     * Tangent frame at a point on a manifold.
     * Stores both the frame vectors {eᵢ} and reciprocal frame {eⁱ}
     * satisfying eⁱ · eⱼ = δⁱⱼ.
     */
    class TangentFrame {
        constructor(vectors) {
            // vectors: array of arrays, shape (dim, ambientDim)
            this.vectors = vectors;
            this.dim = vectors.length;
            this.ambientDim = vectors[0].length;
            this._computeReciprocal();
            this._computeMetric();
        }

        _computeReciprocal() {
            // Compute reciprocal frame eⁱ such that eⁱ · eⱼ = δⁱⱼ
            // For 2D surface in 3D
            if (this.dim === 2) {
                const g = [
                    [dot(this.vectors[0], this.vectors[0]), dot(this.vectors[0], this.vectors[1])],
                    [dot(this.vectors[1], this.vectors[0]), dot(this.vectors[1], this.vectors[1])]
                ];
                const gInv = inv2x2(g);
                if (gInv) {
                    this.reciprocal = [
                        vecAdd(vecScale(this.vectors[0], gInv[0][0]), vecScale(this.vectors[1], gInv[0][1])),
                        vecAdd(vecScale(this.vectors[0], gInv[1][0]), vecScale(this.vectors[1], gInv[1][1]))
                    ];
                } else {
                    this.reciprocal = this.vectors.map(v => v.slice());
                }
            } else {
                this.reciprocal = this.vectors.map(v => v.slice());
            }
        }

        _computeMetric() {
            // Compute metric tensor gᵢⱼ = eᵢ · eⱼ
            this.g = [];
            for (let i = 0; i < this.dim; i++) {
                this.g[i] = [];
                for (let j = 0; j < this.dim; j++) {
                    this.g[i][j] = dot(this.vectors[i], this.vectors[j]);
                }
            }
            if (this.dim === 2) {
                this.gInv = inv2x2(this.g);
                this.sqrtDetG = Math.sqrt(Math.abs(det2x2(this.g)));
            }
        }

        metric(i, j) {
            return this.g[i][j];
        }

        inverseMetric(i, j) {
            return this.gInv ? this.gInv[i][j] : 0;
        }

        /**
         * Compute tangent pseudoscalar B = e₁ ∧ e₂ for 2D surface in 3D
         */
        pseudoscalar2D() {
            if (this.dim !== 2 || this.ambientDim !== 3) {
                throw new Error("pseudoscalar2D requires 2D surface in 3D");
            }
            return Bivector3D.fromWedge(this.vectors[0], this.vectors[1]);
        }

        /**
         * Compute unit normal for 2D surface in 3D
         */
        normal3D() {
            if (this.dim !== 2 || this.ambientDim !== 3) {
                throw new Error("normal3D requires 2D surface in 3D");
            }
            const c = cross3(this.vectors[0], this.vectors[1]);
            const n = norm(c);
            if (n < 1e-10) {
                throw new Error("Degenerate frame: vectors are parallel");
            }
            return vecScale(c, 1 / n);
        }
    }

    // ========================================================================
    // MANIFOLD CLASSES
    // ========================================================================

    class Sphere {
        constructor(R = 1.0) {
            this.R = R;
            this.dim = 2;
            this.ambientDim = 3;
        }

        embedding(theta, phi) {
            const R = this.R;
            return [
                R * Math.sin(theta) * Math.cos(phi),
                R * Math.sin(theta) * Math.sin(phi),
                R * Math.cos(theta)
            ];
        }

        frame(theta, phi) {
            const R = this.R;
            // Avoid singularity at poles
            theta = Math.max(0.01, Math.min(Math.PI - 0.01, theta));

            const ct = Math.cos(theta), st = Math.sin(theta);
            const cp = Math.cos(phi), sp = Math.sin(phi);

            const e_theta = [R * ct * cp, R * ct * sp, -R * st];
            const e_phi = [R * st * (-sp), R * st * cp, 0];

            return new TangentFrame([e_theta, e_phi]);
        }

        metric(theta, phi) {
            const R = this.R;
            return [
                [R * R, 0],
                [0, R * R * Math.sin(theta) * Math.sin(theta)]
            ];
        }

        gaussianCurvature(theta, phi) {
            return 1.0 / (this.R * this.R);
        }
    }

    class Torus {
        constructor(R = 2.0, r = 0.5) {
            this.R = R;  // Major radius
            this.r = r;  // Minor radius
            this.dim = 2;
            this.ambientDim = 3;
        }

        embedding(theta, phi) {
            const { R, r } = this;
            return [
                (R + r * Math.cos(theta)) * Math.cos(phi),
                (R + r * Math.cos(theta)) * Math.sin(phi),
                r * Math.sin(theta)
            ];
        }

        frame(theta, phi) {
            const { R, r } = this;
            const ct = Math.cos(theta), st = Math.sin(theta);
            const cp = Math.cos(phi), sp = Math.sin(phi);

            const e_theta = [r * (-st) * cp, r * (-st) * sp, r * ct];
            const e_phi = [(R + r * ct) * (-sp), (R + r * ct) * cp, 0];

            return new TangentFrame([e_theta, e_phi]);
        }

        gaussianCurvature(theta, phi) {
            const { R, r } = this;
            return Math.cos(theta) / (r * (R + r * Math.cos(theta)));
        }
    }

    class Paraboloid {
        constructor(a = 1.0) {
            this.a = a;
            this.dim = 2;
            this.ambientDim = 3;
        }

        embedding(r, theta) {
            const a = this.a;
            return [
                r * Math.cos(theta),
                r * Math.sin(theta),
                a * r * r
            ];
        }

        frame(r, theta) {
            const a = this.a;
            r = Math.max(r, 0.01);  // Avoid singularity at origin

            const ct = Math.cos(theta), st = Math.sin(theta);

            const e_r = [ct, st, 2 * a * r];
            const e_theta = [r * (-st), r * ct, 0];

            return new TangentFrame([e_r, e_theta]);
        }
    }

    // ========================================================================
    // CONNECTION BIVECTOR
    // ========================================================================

    /**
     * The connection bivector ω encodes frame rotation.
     * Replaces Christoffel symbols Γᵏᵢⱼ with a geometric object:
     *     ωᵢ = ½(∂ᵢeⱼ) ∧ eʲ
     * 
     * The frame derivative is then:
     *     ∂ᵢeⱼ = ωᵢ × eⱼ
     */
    class ConnectionBivector {
        constructor(manifoldDim = 2, ambientDim = 3, h = 1e-6) {
            this.dim = manifoldDim;
            this.ambientDim = ambientDim;
            this.h = h;
        }

        /**
         * Compute connection bivectors ωᵢ at given coordinates
         * @param {Array} coords - Coordinate values
         * @param {Function} frameFunc - Function (coords) -> TangentFrame
         * @returns {Array} List of Bivector3D [ω₁, ω₂, ...]
         */
        computeAt(coords, frameFunc) {
            const frame = frameFunc(coords);
            const omegaList = [];

            for (let i = 0; i < this.dim; i++) {
                // Compute ∂ᵢeⱼ using central differences
                const coordsPlus = coords.slice();
                const coordsMinus = coords.slice();
                coordsPlus[i] += this.h;
                coordsMinus[i] -= this.h;

                const framePlus = frameFunc(coordsPlus);
                const frameMinus = frameFunc(coordsMinus);

                // ∂ᵢeⱼ for each j
                const dFrame = [];
                for (let j = 0; j < this.dim; j++) {
                    const dej = vecScale(
                        vecSub(framePlus.vectors[j], frameMinus.vectors[j]),
                        1 / (2 * this.h)
                    );
                    dFrame.push(dej);
                }

                // ωᵢ = ½ Σⱼ eʲ ∧ (∂ᵢeⱼ) (correct sign: reciprocal frame first)
                let omegaI = new Bivector3D();
                for (let j = 0; j < this.dim; j++) {
                    const dej = dFrame[j];
                    const ejRecip = frame.reciprocal[j];
                    const wedge = Bivector3D.fromWedge(ejRecip, dej);
                    omegaI = omegaI.add(wedge.scale(0.5));
                }

                omegaList.push(omegaI);
            }

            return omegaList;
        }

        /**
         * Compute ∂ᵢeⱼ = ωᵢ × eⱼ
         */
        frameDerivative(omegaI, frame, j) {
            const ej = frame.vectors[j];
            return omegaI.commutatorWithVector(ej);
        }
    }

    // ========================================================================
    // CURVATURE 2-FORM
    // ========================================================================

    /**
     * The curvature 2-form Ω via Cartan's structure equation:
     *     Ω = dω + ω ∧ ω
     * 
     * Key properties:
     * - Ωᵢⱼ = ∂ᵢωⱼ - ∂ⱼωᵢ + ωᵢ × ωⱼ
     * - [∇ᵢ, ∇ⱼ]v = Ωᵢⱼ × v (Riemann action on vectors)
     */
    class Curvature2Form {
        constructor(connection, h = 1e-6) {
            this.connection = connection;
            this.dim = connection.dim;
            this.h = h;
        }

        /**
         * Compute curvature 2-form components Ωᵢⱼ
         * @returns {Array} 2D array of Bivector3D where result[i][j] = Ωᵢⱼ
         */
        computeAt(coords, frameFunc) {
            const omega = this.connection.computeAt(coords, frameFunc);

            // Compute partial derivatives of connection
            const dOmega = [];
            for (let i = 0; i < this.dim; i++) {
                const coordsPlus = coords.slice();
                const coordsMinus = coords.slice();
                coordsPlus[i] += this.h;
                coordsMinus[i] -= this.h;

                const omegaPlus = this.connection.computeAt(coordsPlus, frameFunc);
                const omegaMinus = this.connection.computeAt(coordsMinus, frameFunc);

                const dOmegaI = [];
                for (let j = 0; j < this.dim; j++) {
                    const dij = omegaPlus[j].sub(omegaMinus[j]).scale(1 / (2 * this.h));
                    dOmegaI.push(dij);
                }
                dOmega.push(dOmegaI);
            }

            // Compute Ωᵢⱼ = ∂ᵢωⱼ - ∂ⱼωᵢ + ωᵢ × ωⱼ
            const Omega = [];
            for (let i = 0; i < this.dim; i++) {
                Omega[i] = [];
                for (let j = 0; j < this.dim; j++) {
                    // dω part: ∂ᵢωⱼ - ∂ⱼωᵢ
                    const dPart = dOmega[i][j].sub(dOmega[j][i]);

                    // ω ∧ ω part: ωᵢ × ωⱼ
                    const wedgePart = omega[i].commutatorWithBivector(omega[j]);

                    Omega[i][j] = dPart.add(wedgePart);
                }
            }

            return Omega;
        }

        /**
         * Compute sectional curvature K(eᵢ ∧ eⱼ)
         * K = Ω(eᵢ, eⱼ) · (eᵢ ∧ eⱼ) / |eᵢ ∧ eⱼ|²
         */
        sectionalCurvature(OmegaIJ, frame, i, j) {
            const ei = frame.vectors[i];
            const ej = frame.vectors[j];
            const plane = Bivector3D.fromWedge(ei, ej);

            const numerator = OmegaIJ.innerWithBivector(plane);
            const denominator = plane.normSquared;

            if (denominator < 1e-20) return 0.0;

            return numerator / denominator;
        }
    }

    // ========================================================================
    // COVARIANT DERIVATIVE
    // ========================================================================

    /**
     * Covariant derivative using connection bivector:
     *     ∇ᵤv = ∂ᵤv + ω(u) × v
     */
    class GACovariantDerivative {
        constructor(connection) {
            this.connection = connection;
            this.h = connection.h;
        }

        /**
         * Compute ∇ᵤv at given coordinates
         */
        ofVector(vFunc, u, coords, frameFunc) {
            const omega = this.connection.computeAt(coords, frameFunc);

            // Compute ω(u) = uⁱωᵢ
            let omegaU = new Bivector3D();
            for (let i = 0; i < u.length; i++) {
                omegaU = omegaU.add(omega[i].scale(u[i]));
            }

            // Compute partial derivative ∂ᵤv = uⁱ∂ᵢv
            let partialUV = [0, 0, 0];
            for (let i = 0; i < u.length; i++) {
                const coordsPlus = coords.slice();
                const coordsMinus = coords.slice();
                coordsPlus[i] += this.h;
                coordsMinus[i] -= this.h;

                const vPlus = vFunc(coordsPlus);
                const vMinus = vFunc(coordsMinus);

                const partialIV = vecScale(vecSub(vPlus, vMinus), 1 / (2 * this.h));
                partialUV = vecAdd(partialUV, vecScale(partialIV, u[i]));
            }

            const v = vFunc(coords);

            // ∇ᵤv = ∂ᵤv + ω(u) × v
            return vecAdd(partialUV, omegaU.commutatorWithVector(v));
        }

        /**
         * Compute divergence: div v = ∇ · v
         */
        divergence(vFunc, coords, frameFunc) {
            const frame = frameFunc(coords);
            const sqrtG = frame.sqrtDetG;

            let div = 0;
            for (let i = 0; i < this.connection.dim; i++) {
                const coordsPlus = coords.slice();
                const coordsMinus = coords.slice();
                coordsPlus[i] += this.h;
                coordsMinus[i] -= this.h;

                const framePlus = frameFunc(coordsPlus);
                const frameMinus = frameFunc(coordsMinus);

                const vPlus = vFunc(coordsPlus);
                const vMinus = vFunc(coordsMinus);

                // vⁱ = v · eⁱ
                const viPlus = dot(vPlus, framePlus.reciprocal[i]);
                const viMinus = dot(vMinus, frameMinus.reciprocal[i]);

                const sqrtGViPlus = framePlus.sqrtDetG * viPlus;
                const sqrtGViMinus = frameMinus.sqrtDetG * viMinus;

                div += (sqrtGViPlus - sqrtGViMinus) / (2 * this.h);
            }

            return div / (sqrtG + 1e-10);
        }
    }

    // ========================================================================
    // GEODESIC SOLVER
    // ========================================================================

    /**
     * Solve geodesic equation using connection bivector:
     *     v̇ + ω(v) × v = 0
     * 
     * This is equivalent to ∇ᵥv = 0 — velocity parallel-transports itself.
     */
    class GAGeodesicSolver {
        constructor(connection) {
            this.connection = connection;
        }

        /**
         * Integrate geodesic equation from initial conditions
         * @param {Array} x0 - Initial coordinates
         * @param {Array} v0 - Initial velocity (in coordinate basis)
         * @param {Function} frameFunc - Function (coords) -> TangentFrame
         * @param {number} tFinal - Final time
         * @param {number} nSteps - Number of integration steps
         * @returns {Object} {t, x, v} arrays
         */
        solve(x0, v0, frameFunc, tFinal, nSteps = 100) {
            const dt = tFinal / nSteps;
            const dim = x0.length;

            const tValues = [];
            const xValues = [];
            const vValues = [];

            for (let i = 0; i <= nSteps; i++) {
                tValues.push(i * dt);
            }

            let x = x0.slice();
            let v = v0.slice();

            xValues.push(x.slice());
            vValues.push(v.slice());

            for (let step = 0; step < nSteps; step++) {
                const omega = this.connection.computeAt(x, frameFunc);
                const frame = frameFunc(x);

                // Convert coordinate velocity to ambient space
                let vAmbient = [0, 0, 0];
                for (let i = 0; i < dim; i++) {
                    vAmbient = vecAdd(vAmbient, vecScale(frame.vectors[i], v[i]));
                }

                // Compute ω(v) = vⁱωᵢ
                let omegaV = new Bivector3D();
                for (let i = 0; i < dim; i++) {
                    omegaV = omegaV.add(omega[i].scale(v[i]));
                }

                // Geodesic equation: v̇ = -ω(v) × v
                const vDotAmbient = vecScale(omegaV.commutatorWithVector(vAmbient), -1);

                // Convert acceleration back to coordinates
                const vDot = [];
                for (let i = 0; i < dim; i++) {
                    vDot[i] = dot(vDotAmbient, frame.reciprocal[i]);
                }

                // Simple Euler integration
                const xNew = vecAdd(x, vecScale(v, dt));
                const vNew = vecAdd(v, vecScale(vDot, dt));

                x = xNew;
                v = vNew;

                xValues.push(x.slice());
                vValues.push(v.slice());
            }

            return { t: tValues, x: xValues, v: vValues };
        }

        /**
         * RK4 integration for better accuracy
         */
        solveRK4(x0, v0, frameFunc, tFinal, nSteps = 100) {
            const dt = tFinal / nSteps;
            const dim = x0.length;

            const tValues = [];
            const xValues = [];
            const vValues = [];

            let x = x0.slice();
            let v = v0.slice();

            for (let step = 0; step <= nSteps; step++) {
                tValues.push(step * dt);
                xValues.push(x.slice());
                vValues.push(v.slice());

                if (step === nSteps) break;

                // RK4 step
                const [k1x, k1v] = this._geodesicRHS(x, v, frameFunc);

                const x2 = vecAdd(x, vecScale(k1x, dt / 2));
                const v2 = vecAdd(v, vecScale(k1v, dt / 2));
                const [k2x, k2v] = this._geodesicRHS(x2, v2, frameFunc);

                const x3 = vecAdd(x, vecScale(k2x, dt / 2));
                const v3 = vecAdd(v, vecScale(k2v, dt / 2));
                const [k3x, k3v] = this._geodesicRHS(x3, v3, frameFunc);

                const x4 = vecAdd(x, vecScale(k3x, dt));
                const v4 = vecAdd(v, vecScale(k3v, dt));
                const [k4x, k4v] = this._geodesicRHS(x4, v4, frameFunc);

                // Combine
                for (let i = 0; i < dim; i++) {
                    x[i] += dt / 6 * (k1x[i] + 2 * k2x[i] + 2 * k3x[i] + k4x[i]);
                    v[i] += dt / 6 * (k1v[i] + 2 * k2v[i] + 2 * k3v[i] + k4v[i]);
                }
            }

            return { t: tValues, x: xValues, v: vValues };
        }

        _geodesicRHS(x, v, frameFunc) {
            const dim = x.length;
            const omega = this.connection.computeAt(x, frameFunc);
            const frame = frameFunc(x);

            // dx/dt = v
            const xDot = v.slice();

            // Convert to ambient
            let vAmbient = [0, 0, 0];
            for (let i = 0; i < dim; i++) {
                vAmbient = vecAdd(vAmbient, vecScale(frame.vectors[i], v[i]));
            }

            // ω(v)
            let omegaV = new Bivector3D();
            for (let i = 0; i < dim; i++) {
                omegaV = omegaV.add(omega[i].scale(v[i]));
            }

            // dv/dt = -ω(v) × v
            const vDotAmbient = vecScale(omegaV.commutatorWithVector(vAmbient), -1);

            const vDot = [];
            for (let i = 0; i < dim; i++) {
                vDot[i] = dot(vDotAmbient, frame.reciprocal[i]);
            }

            return [xDot, vDot];
        }
    }

    // ========================================================================
    // PARALLEL TRANSPORT
    // ========================================================================

    /**
     * Parallel transport using connection bivector:
     *     ẇ + ω(v) × w = 0
     * 
     * Transport vector w along curve with velocity v, keeping ∇ᵥw = 0.
     */
    class GAParallelTransport {
        constructor(connection) {
            this.connection = connection;
        }

        /**
         * Transport vector w₀ along curve
         */
        transport(curveFunc, w0, frameFunc, tStart, tEnd, nSteps = 100) {
            const dt = (tEnd - tStart) / nSteps;
            const dim = w0.length;

            const tValues = [];
            const wValues = [];

            let w = w0.slice();

            for (let step = 0; step <= nSteps; step++) {
                const t = tStart + step * dt;
                tValues.push(t);
                wValues.push(w.slice());

                if (step === nSteps) break;

                const coords = curveFunc(t);

                // Curve velocity (numerical)
                const coordsPlus = curveFunc(t + dt / 10);
                const coordsMinus = curveFunc(t - dt / 10);
                const v = vecScale(vecSub(coordsPlus, coordsMinus), 5 / dt);

                const omega = this.connection.computeAt(coords, frameFunc);
                const frame = frameFunc(coords);

                // Convert w to ambient space
                let wAmbient = [0, 0, 0];
                for (let i = 0; i < dim; i++) {
                    wAmbient = vecAdd(wAmbient, vecScale(frame.vectors[i], w[i]));
                }

                // Compute ω(v) = vⁱωᵢ
                let omegaV = new Bivector3D();
                for (let i = 0; i < dim; i++) {
                    omegaV = omegaV.add(omega[i].scale(v[i]));
                }

                // Transport equation: ẇ = -ω(v) × w
                const wDotAmbient = vecScale(omegaV.commutatorWithVector(wAmbient), -1);

                // Convert back to coordinates
                const wDot = [];
                for (let i = 0; i < dim; i++) {
                    wDot[i] = dot(wDotAmbient, frame.reciprocal[i]);
                }

                // Update
                w = vecAdd(w, vecScale(wDot, dt));
            }

            return { t: tValues, w: wValues };
        }

        /**
         * Compute holonomy around a closed loop
         * @returns {Object} {wFinal, angle}
         */
        holonomy(loopFunc, w0, frameFunc, nSteps = 200) {
            const result = this.transport(
                loopFunc, w0, frameFunc,
                0.0, 2 * Math.PI, nSteps
            );

            const wFinal = result.w[result.w.length - 1];

            // Compute rotation angle
            const dotProd = dot(w0, wFinal);
            const norm0 = norm(w0);
            const normFinal = norm(wFinal);

            if (norm0 < 1e-10 || normFinal < 1e-10) {
                return { wFinal, angle: 0.0 };
            }

            let cosAngle = dotProd / (norm0 * normFinal);
            cosAngle = Math.max(-1, Math.min(1, cosAngle));
            const angle = Math.acos(cosAngle);

            return { wFinal, angle };
        }
    }

    // ========================================================================
    // SHAPE OPERATOR
    // ========================================================================

    /**
     * Shape operator (Weingarten map) via pseudoscalar variation:
     *     S(v) = -∇ᵥn = (v · ∇B) · I⁻¹
     * 
     * Principal curvatures are eigenvalues of S.
     * Gaussian curvature K = det(S) = κ₁κ₂
     * Mean curvature H = ½tr(S) = ½(κ₁ + κ₂)
     */
    class ShapeOperator {
        constructor(h = 1e-6) {
            this.h = h;
        }

        /**
         * Compute shape operator matrix Sⁱⱼ at given coordinates
         */
        computeAt(coords, frameFunc) {
            const frame = frameFunc(coords);

            if (frame.dim !== 2 || frame.ambientDim !== 3) {
                throw new Error("Shape operator requires 2D surface in 3D");
            }

            const n = frame.normal3D();

            // Compute ∂ᵢn using central differences
            const S = [[0, 0], [0, 0]];

            for (let j = 0; j < 2; j++) {
                const coordsPlus = coords.slice();
                const coordsMinus = coords.slice();
                coordsPlus[j] += this.h;
                coordsMinus[j] -= this.h;

                const framePlus = frameFunc(coordsPlus);
                const frameMinus = frameFunc(coordsMinus);

                const nPlus = framePlus.normal3D();
                const nMinus = frameMinus.normal3D();

                // ∂ⱼn
                const dnJ = vecScale(vecSub(nPlus, nMinus), 1 / (2 * this.h));

                // S(eⱼ) = -∂ⱼn projected to tangent space
                // Sⁱⱼ = -∂ⱼn · eⁱ
                for (let i = 0; i < 2; i++) {
                    S[i][j] = -dot(dnJ, frame.reciprocal[i]);
                }
            }

            return S;
        }

        gaussianCurvature(S) {
            return det2x2(S);
        }

        meanCurvature(S) {
            return 0.5 * trace2x2(S);
        }

        principalCurvatures(S) {
            return eigenvalues2x2(S).sort((a, b) => a - b);
        }
    }

    // ========================================================================
    // RIEMANNIAN GA MAIN CLASS (combines all functionality)
    // ========================================================================

    class RiemannianGA {
        constructor(manifold, h = 1e-5) {
            this.manifold = manifold;
            this.h = h;

            // Create frame function from manifold
            this.frameFunc = (coords) => manifold.frame(coords[0], coords[1]);

            // Initialize components
            this.connection = new ConnectionBivector(manifold.dim, manifold.ambientDim || 3, h);
            this.curvature = new Curvature2Form(this.connection, h);
            this.covariantDerivative = new GACovariantDerivative(this.connection);
            this.geodesicSolver = new GAGeodesicSolver(this.connection);
            this.parallelTransport = new GAParallelTransport(this.connection);
            this.shapeOperator = new ShapeOperator(h);
        }

        /**
         * Get connection bivectors at a point
         */
        getConnection(coords) {
            return this.connection.computeAt(coords, this.frameFunc);
        }

        /**
         * Get curvature 2-form at a point
         */
        getCurvature(coords) {
            return this.curvature.computeAt(coords, this.frameFunc);
        }

        /**
         * Compute Gaussian curvature via shape operator
         */
        gaussianCurvature(coords) {
            const S = this.shapeOperator.computeAt(coords, this.frameFunc);
            return this.shapeOperator.gaussianCurvature(S);
        }

        /**
         * Compute mean curvature via shape operator
         */
        meanCurvature(coords) {
            const S = this.shapeOperator.computeAt(coords, this.frameFunc);
            return this.shapeOperator.meanCurvature(S);
        }

        /**
         * Get principal curvatures
         */
        principalCurvatures(coords) {
            const S = this.shapeOperator.computeAt(coords, this.frameFunc);
            return this.shapeOperator.principalCurvatures(S);
        }

        /**
         * Compute sectional curvature via curvature 2-form
         */
        sectionalCurvature(coords) {
            const Omega = this.curvature.computeAt(coords, this.frameFunc);
            const frame = this.frameFunc(coords);
            return this.curvature.sectionalCurvature(Omega[0][1], frame, 0, 1);
        }

        /**
         * Solve geodesic equation
         */
        solveGeodesic(x0, v0, tFinal, nSteps = 100, useRK4 = true) {
            if (useRK4) {
                return this.geodesicSolver.solveRK4(x0, v0, this.frameFunc, tFinal, nSteps);
            }
            return this.geodesicSolver.solve(x0, v0, this.frameFunc, tFinal, nSteps);
        }

        /**
         * Parallel transport a vector along a curve
         */
        transportVector(curveFunc, w0, tStart, tEnd, nSteps = 100) {
            return this.parallelTransport.transport(
                curveFunc, w0, this.frameFunc, tStart, tEnd, nSteps
            );
        }

        /**
         * Compute holonomy around a loop
         */
        computeHolonomy(loopFunc, w0, nSteps = 200) {
            return this.parallelTransport.holonomy(loopFunc, w0, this.frameFunc, nSteps);
        }

        /**
         * Get embedding point in ambient space
         */
        embed(coords) {
            return this.manifold.embedding(coords[0], coords[1]);
        }

        /**
         * Get tangent frame at a point
         */
        getFrame(coords) {
            return this.frameFunc(coords);
        }
    }

    // ========================================================================
    // TEST FUNCTIONS
    // ========================================================================

    function testSphereCurvature() {
        console.log("=".repeat(60));
        console.log("TEST: Sphere Curvature via Connection Bivector");
        console.log("=".repeat(60));

        const R = 1.0;
        const sphere = new Sphere(R);
        const rga = new RiemannianGA(sphere);

        // Test at θ = π/3 (60°)
        const coords = [Math.PI / 3, 0.0];

        // Method 1: Via shape operator
        const K_shape = rga.gaussianCurvature(coords);
        const H_shape = rga.meanCurvature(coords);

        // Note: With S = -∇ᵥn and outward normal, principal curvatures are -1/R
        // So H = -1/R for this convention (Gaussian K is still 1/R²)
        const expectedH = -1 / R;

        console.log("\nMethod 1: Shape Operator");
        console.log(`  K = ${K_shape.toFixed(6)} (expected: ${(1 / R / R).toFixed(6)})`);
        console.log(`  H = ${H_shape.toFixed(6)} (expected: ${expectedH.toFixed(6)})`);
        console.log("  Note: With S = -∇ᵥn and outward normal, κ = -1/R, H = -1/R");

        // Method 2: Via curvature 2-form (intrinsic)
        const K_curv = rga.sectionalCurvature(coords);

        console.log("\nMethod 2: Curvature 2-Form (intrinsic)");
        console.log(`  K = ${K_curv.toFixed(6)} (expected: ${(1 / R / R).toFixed(6)})`);

        // Check connection bivector
        const omega = rga.getConnection(coords);
        console.log("\nConnection Bivectors:");
        console.log(`  ω_θ = ${omega[0]}`);
        console.log(`  ω_φ = ${omega[1]}`);

        const errorK = Math.abs(K_shape - 1 / R / R);
        const errorH = Math.abs(H_shape - expectedH);

        console.log("\nErrors:");
        console.log(`  |K - 1/R²| = ${errorK.toFixed(6)}`);
        console.log(`  |H - (-1/R)| = ${errorH.toFixed(6)}`);

        console.log(`\nStatus: ${errorK < 0.1 && errorH < 0.1 ? 'PASS' : 'FAIL'}`);
    }

    function testTorusCurvature() {
        console.log("\n" + "=".repeat(60));
        console.log("TEST: Torus Curvature (Variable K)");
        console.log("=".repeat(60));

        const R = 2.0, r = 0.5;
        const torus = new Torus(R, r);
        const rga = new RiemannianGA(torus);

        const testPoints = [
            [0.0, 0.0, "Outer equator (K > 0)"],
            [Math.PI, 0.0, "Inner equator (K < 0)"],
            [Math.PI / 2, 0.0, "Top (K = 0)"]
        ];

        for (const [theta, phi, name] of testPoints) {
            const coords = [theta, phi];
            const K = rga.gaussianCurvature(coords);
            const H = rga.meanCurvature(coords);

            // Expected: K = cos(θ) / (r(R + r cos θ))
            const ct = Math.cos(theta);
            const K_expected = ct / (r * (R + r * ct));

            console.log(`\n${name}:`);
            console.log(`  θ = ${theta.toFixed(4)}, φ = ${phi.toFixed(4)}`);
            console.log(`  K = ${K.toFixed(6)} (expected: ${K_expected.toFixed(6)})`);
            console.log(`  H = ${H.toFixed(6)}`);
        }
    }

    function testParallelTransportHolonomy() {
        console.log("\n" + "=".repeat(60));
        console.log("TEST: Parallel Transport Holonomy on Sphere");
        console.log("=".repeat(60));

        const R = 1.0;
        const sphere = new Sphere(R);
        const rga = new RiemannianGA(sphere);

        // Transport around latitude circle at θ = π/4
        const thetaFixed = Math.PI / 4;

        const latitudeLoop = (t) => [thetaFixed, t];

        // Initial vector in θ direction
        const w0 = [1.0, 0.0];

        const { wFinal, angle } = rga.computeHolonomy(latitudeLoop, w0, 200);

        // Expected holonomy: 2π(1 - cos θ)
        const expectedAngle = 2 * Math.PI * (1 - Math.cos(thetaFixed));

        console.log(`\nLatitude circle at θ = π/4 (${(thetaFixed * 180 / Math.PI).toFixed(1)}°)`);
        console.log(`Initial vector: [${w0.join(', ')}]`);
        console.log(`Final vector: [${wFinal.map(x => x.toFixed(4)).join(', ')}]`);
        console.log(`Holonomy angle: ${(angle * 180 / Math.PI).toFixed(2)}°`);
        console.log(`Expected angle: ${(expectedAngle * 180 / Math.PI).toFixed(2)}°`);
        console.log(`Error: ${(Math.abs(angle - expectedAngle) * 180 / Math.PI).toFixed(2)}°`);
    }

    function testGeodesic() {
        console.log("\n" + "=".repeat(60));
        console.log("TEST: Geodesic on Sphere (Great Circle)");
        console.log("=".repeat(60));

        const sphere = new Sphere(1.0);
        const rga = new RiemannianGA(sphere);

        // Start at equator, move in positive theta direction (toward south pole)
        // Note: θ=0 is north pole, θ=π is south pole
        const x0 = [Math.PI / 2, 0];
        const v0 = [1, 0];  // Initial velocity in positive θ direction

        const result = rga.solveGeodesic(x0, v0, Math.PI / 2, 50);

        console.log("\nGeodesic from equator, moving along meridian (great circle):");
        console.log(`  Start: θ=${x0[0].toFixed(4)} (equator), φ=${x0[1].toFixed(4)}`);
        console.log(`  Velocity: v0=[${v0.join(', ')}] (positive θ = toward south pole)`);

        const final = result.x[result.x.length - 1];
        console.log(`  End: θ=${final[0].toFixed(4)}, φ=${final[1].toFixed(4)}`);
        console.log(`  Expected: θ ≈ π (south pole, since v₀ points toward increasing θ)`);
        console.log(`  Status: ${final[0] > x0[0] ? 'PASS' : 'FAIL'}`);
    }

    function runAllTests() {
        testSphereCurvature();
        testTorusCurvature();
        testParallelTransportHolonomy();
        testGeodesic();

        console.log("\n" + "=".repeat(60));
        console.log("ALL TESTS COMPLETED");
        console.log("=".repeat(60));
    }

    // ========================================================================
    // EXPORTS
    // ========================================================================

    const exports = {
        // Utility functions
        dot, cross3, norm, normalize, vecAdd, vecSub, vecScale,

        // Bivector classes
        Bivector,
        Bivector2D,
        Bivector3D,

        // Frame
        TangentFrame,

        // Manifolds
        Sphere,
        Torus,
        Paraboloid,

        // Core geometric objects
        ConnectionBivector,
        Curvature2Form,
        GACovariantDerivative,
        GAGeodesicSolver,
        GAParallelTransport,
        ShapeOperator,

        // Main class
        RiemannianGA,

        // Tests
        testSphereCurvature,
        testTorusCurvature,
        testParallelTransportHolonomy,
        testGeodesic,
        runAllTests
    };

    // Export to global or module
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = exports;
    }
    global.RiemannianGA = exports;

})(typeof window !== 'undefined' ? window : global);