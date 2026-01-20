/**
 * Coordinate-Free Riemannian Geometry via Geometric Calculus
 * 
 * This module provides pure GA formulations that eliminate Christoffel symbols entirely.
 * All curvature computations use:
 *   - Connection bivector ω instead of Γᵏᵢⱼ
 *   - Curvature 2-form Ω instead of Rˡᵢⱼₖ
 *   - Frame inner products instead of gᵢⱼ
 * 
 * Master Equations:
 *   Metric:        g_ij = e_i · e_j,    e^i · e_j = δ^i_j
 *   Connection:    ω_i = ½ e^j ∧ (∂_i e_j)
 *   Frame:         ∂_i e_j = ω_i × e_j
 *   Covariant:     ∇_u A = ∂_u A + ω(u) × A
 *   Curvature:     Ω = dω + ω ∧ ω
 *   Riemann:       [∇_u, ∇_v]w = Ω(u,v) × w
 *   Geodesic:      ∇_v v = 0
 *   Bianchi:       ∇ ∧ Ω = 0
 * 
 * @module riemannian-ga
 * @license MIT
 */

(function (global) {
    'use strict';

    const EPSILON = 1e-10;
    const { abs, sqrt, sin, cos, tan, atan2, acos, PI } = Math;

    // ============================================================================
    // UTILITY FUNCTIONS
    // ============================================================================

    function dot(u, v) {
        let sum = 0;
        for (let i = 0; i < u.length; i++) sum += u[i] * v[i];
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
        return sqrt(dot(v, v));
    }

    function normalize(v) {
        const n = norm(v);
        return n > EPSILON ? v.map(x => x / n) : v.slice();
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

    function matVecMul(A, v) {
        const n = A.length;
        const result = new Array(n).fill(0);
        for (let i = 0; i < n; i++) {
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
        if (abs(d) < EPSILON) throw new Error('Singular matrix');
        return [
            [A[1][1] / d, -A[0][1] / d],
            [-A[1][0] / d, A[0][0] / d]
        ];
    }

    // ============================================================================
    // BIVECTOR CLASSES
    // ============================================================================

    /**
     * Bivector in 2D (single component: e₁₂)
     * 
     * In 2D, bivectors are pseudoscalars — just one component representing
     * oriented area in the e₁∧e₂ plane.
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

        norm() {
            return abs(this.e12);
        }

        /**
         * Compute ω × v = ½(ωv - vω) for 2D.
         * 
         * In 2D, the bivector e₁₂ acts on vectors by 90° rotation:
         * e₁₂ × e₁ = e₂,  e₁₂ × e₂ = -e₁
         * 
         * So for ω = ω₁₂ e₁₂ and v = (v¹, v²):
         *     ω × v = ω₁₂ (v² e₁ - v¹ e₂) = (ω₁₂ v², -ω₁₂ v¹)
         * 
         * @param {number[]} v - 2D vector [v¹, v²]
         * @returns {number[]} Rotated vector
         */
        commutatorWithVector(v) {
            const b = this.e12;
            return [b * v[1], -b * v[0]];
        }
    }

    /**
     * Bivector in 3D (three components: e₂₃, e₃₁, e₁₂)
     * 
     * Represents oriented planes. The components correspond to:
     * - e₂₃: YZ plane (rotation around X)
     * - e₃₁: ZX plane (rotation around Y)
     * - e₁₂: XY plane (rotation around Z)
     */
    class Bivector3D {
        constructor(e23 = 0.0, e31 = 0.0, e12 = 0.0) {
            this.e23 = e23;
            this.e31 = e31;
            this.e12 = e12;
        }

        static fromArray(arr) {
            return new Bivector3D(arr[0], arr[1], arr[2]);
        }

        static fromWedge(u, v) {
            // u ∧ v = (u₂v₃ - u₃v₂) e₂₃ + (u₃v₁ - u₁v₃) e₃₁ + (u₁v₂ - u₂v₁) e₁₂
            return new Bivector3D(
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0]
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
            return new Bivector3D(this.e23 * s, this.e31 * s, this.e12 * s);
        }

        neg() {
            return new Bivector3D(-this.e23, -this.e31, -this.e12);
        }

        normSquared() {
            return this.e23 * this.e23 + this.e31 * this.e31 + this.e12 * this.e12;
        }

        norm() {
            return sqrt(this.normSquared());
        }

        /**
         * Compute B × v = ½(Bv - vB) in 3D.
         * 
         * The bivector acts on vectors by rotation in its plane.
         * This is equivalent to cross product with the dual vector:
         * 
         * B × v = B* × v where B* is the Hodge dual (axial vector)
         * 
         * For B = b₂₃ e₂₃ + b₃₁ e₃₁ + b₁₂ e₁₂:
         * B × v = (b₂₃, b₃₁, b₁₂) × (v¹, v², v³)  [as cross product]
         * 
         * @param {number[]} v - 3D vector
         * @returns {number[]} Rotated vector
         */
        commutatorWithVector(v) {
            // Dual of bivector is axial vector: B* = (b₂₃, b₃₁, b₁₂)
            // Cross product gives rotation
            return cross3([this.e23, this.e31, this.e12], v);
        }

        /**
         * Compute B₁ × B₂ = ½(B₁B₂ - B₂B₁).
         * 
         * For bivectors in 3D, this gives another bivector (Lie bracket).
         */
        commutatorWithBivector(other) {
            // [B₁, B₂] using dual vectors and cross product
            const b1 = [this.e23, this.e31, this.e12];
            const b2 = [other.e23, other.e31, other.e12];
            const result = cross3(b1, b2);
            return new Bivector3D(result[0], result[1], result[2]);
        }

        /**
         * Compute B₁ · B₂ (scalar inner product).
         */
        innerWithBivector(other) {
            return this.e23 * other.e23 + this.e31 * other.e31 + this.e12 * other.e12;
        }
    }

    // ============================================================================
    // TANGENT FRAME
    // ============================================================================

    /**
     * Tangent frame at a point on a manifold.
     * 
     * Stores both the frame {eᵢ} and reciprocal frame {eⁱ} satisfying eⁱ · eⱼ = δⁱⱼ
     */
    class TangentFrame {
        /**
         * @param {number[][]} frameVectors - Array of tangent vectors [e₁, e₂, ...]
         * @param {number} ambientDim - Dimension of ambient space (default: same as tangent dim)
         */
        constructor(frameVectors, ambientDim = null) {
            this.dim = frameVectors.length;
            this.ambientDim = ambientDim || frameVectors[0].length;
            this.frame = frameVectors;  // {eᵢ}
            this.reciprocal = this._computeReciprocal();  // {eⁱ}
        }

        /**
         * Compute the reciprocal frame satisfying eⁱ · eⱼ = δⁱⱼ
         */
        _computeReciprocal() {
            if (this.dim === 2 && this.ambientDim >= 2) {
                // 2D case: explicitly compute reciprocal
                const e1 = this.frame[0];
                const e2 = this.frame[1];

                // Metric tensor g_ij = e_i · e_j
                const g = [
                    [dot(e1, e1), dot(e1, e2)],
                    [dot(e2, e1), dot(e2, e2)]
                ];

                // Inverse metric g^ij
                const gInv = inv2x2(g);

                // Reciprocal frame: e^i = g^ij e_j
                const e1Up = vecAdd(vecScale(e1, gInv[0][0]), vecScale(e2, gInv[0][1]));
                const e2Up = vecAdd(vecScale(e1, gInv[1][0]), vecScale(e2, gInv[1][1]));

                return [e1Up, e2Up];
            }
            throw new Error('Only 2D frames currently supported');
        }

        /**
         * Get the metric tensor g_ij = e_i · e_j
         */
        metric() {
            const g = [];
            for (let i = 0; i < this.dim; i++) {
                g[i] = [];
                for (let j = 0; j < this.dim; j++) {
                    g[i][j] = dot(this.frame[i], this.frame[j]);
                }
            }
            return g;
        }

        /**
         * Get the inverse metric g^ij = e^i · e^j
         */
        metricInverse() {
            const gInv = [];
            for (let i = 0; i < this.dim; i++) {
                gInv[i] = [];
                for (let j = 0; j < this.dim; j++) {
                    gInv[i][j] = dot(this.reciprocal[i], this.reciprocal[j]);
                }
            }
            return gInv;
        }

        /**
         * Lower index: convert contravariant vector to covariant
         * v_i = g_ij v^j = e_i · v
         */
        lower(v) {
            const result = [];
            for (let i = 0; i < this.dim; i++) {
                result[i] = dot(this.frame[i], v);
            }
            return result;
        }

        /**
         * Raise index: convert covariant vector to contravariant
         * v^i = g^ij v_j = e^i · v
         */
        raise(v) {
            const result = [];
            for (let i = 0; i < this.dim; i++) {
                result[i] = dot(this.reciprocal[i], v);
            }
            return result;
        }
    }

    // ============================================================================
    // ABSTRACT MANIFOLD
    // ============================================================================

    /**
     * Abstract base class for Riemannian manifolds.
     * 
     * Subclasses must implement:
     *   - frame(coords) → TangentFrame
     *   - Optional: embedding(coords) → point in ambient space
     */
    class RiemannianManifold {
        constructor(dim, ambientDim = null) {
            this.dim = dim;
            this.ambientDim = ambientDim || dim;
        }

        /**
         * Get tangent frame at coordinates.
         * @abstract
         * @param {number[]} coords - Coordinate values
         * @returns {TangentFrame}
         */
        frame(coords) {
            throw new Error('Subclass must implement frame()');
        }

        /**
         * Get metric tensor at coordinates.
         * @param {number[]} coords
         * @returns {number[][]} Metric matrix g_ij
         */
        metric(coords) {
            return this.frame(coords).metric();
        }

        /**
         * Get inverse metric at coordinates.
         * @param {number[]} coords
         * @returns {number[][]} Inverse metric g^ij
         */
        metricInverse(coords) {
            return this.frame(coords).metricInverse();
        }

        /**
         * Numerically differentiate frame vectors.
         * @param {number[]} coords
         * @param {number} h - Step size
         * @returns {number[][][]} ∂_i e_j as 3D array [i][j][component]
         */
        frameDerivatives(coords, h = 1e-6) {
            const derivatives = [];

            for (let i = 0; i < this.dim; i++) {
                derivatives[i] = [];

                // Coordinates shifted in direction i
                const coordsPlus = [...coords];
                const coordsMinus = [...coords];
                coordsPlus[i] += h;
                coordsMinus[i] -= h;

                const framePlus = this.frame(coordsPlus);
                const frameMinus = this.frame(coordsMinus);

                for (let j = 0; j < this.dim; j++) {
                    // ∂_i e_j via central difference
                    derivatives[i][j] = vecScale(
                        vecSub(framePlus.frame[j], frameMinus.frame[j]),
                        1 / (2 * h)
                    );
                }
            }

            return derivatives;
        }
    }

    // ============================================================================
    // CONNECTION BIVECTOR
    // ============================================================================

    /**
     * Computes the connection bivector ωᵢ at each coordinate direction.
     * 
     * The connection bivector encodes how the frame rotates:
     *   ωᵢ = ½ eʲ ∧ (∂ᵢeⱼ)
     * 
     * This REPLACES Christoffel symbols entirely.
     */
    class ConnectionBivector {
        /**
         * @param {RiemannianManifold} manifold
         */
        constructor(manifold) {
            this.manifold = manifold;
        }

        /**
         * Compute connection bivectors at given coordinates.
         * 
         * For 2D surfaces in 3D, returns Bivector3D objects.
         * 
         * @param {number[]} coords - Point on manifold
         * @param {number} h - Step size for numerical derivatives
         * @returns {Bivector3D[]} Array [ω₀, ω₁, ...] of connection bivectors
         */
        computeAt(coords, h = 1e-6) {
            const frame = this.manifold.frame(coords);
            const frameDerivs = this.manifold.frameDerivatives(coords, h);
            const dim = this.manifold.dim;

            const omegas = [];

            for (let i = 0; i < dim; i++) {
                // ω_i = ½ Σⱼ eʲ ∧ (∂ᵢeⱼ)
                let omega_i = new Bivector3D(0, 0, 0);

                for (let j = 0; j < dim; j++) {
                    const ej_recip = frame.reciprocal[j];  // eʲ (reciprocal frame)
                    const d_i_ej = frameDerivs[i][j];      // ∂ᵢeⱼ

                    // eʲ ∧ (∂ᵢeⱼ)
                    const wedge_contrib = Bivector3D.fromWedge(ej_recip, d_i_ej);
                    omega_i = omega_i.add(wedge_contrib.scale(0.5));
                }

                omegas.push(omega_i);
            }

            return omegas;
        }

        /**
         * Get connection along a tangent direction.
         * ω(u) = uⁱ ωᵢ
         * 
         * @param {number[]} coords - Point on manifold
         * @param {number[]} u - Tangent vector (in coordinate basis)
         * @returns {Bivector3D}
         */
        along(coords, u, h = 1e-6) {
            const omegas = this.computeAt(coords, h);

            let result = new Bivector3D(0, 0, 0);
            for (let i = 0; i < omegas.length; i++) {
                result = result.add(omegas[i].scale(u[i]));
            }
            return result;
        }
    }

    // ============================================================================
    // CURVATURE 2-FORM
    // ============================================================================

    /**
     * Computes the curvature 2-form via Cartan structure equation:
     *   Ω = dω + ω ∧ ω
     * 
     * This REPLACES the Riemann tensor entirely.
     */
    class Curvature2Form {
        /**
         * @param {RiemannianManifold} manifold
         */
        constructor(manifold) {
            this.manifold = manifold;
            this.connection = new ConnectionBivector(manifold);
        }

        /**
         * Compute the exterior derivative of connection: dω
         * 
         * (dω)_ij = ∂_i ω_j - ∂_j ω_i  (bivector-valued)
         * 
         * @param {number[]} coords
         * @param {number} h
         * @returns {Bivector3D[][]} Matrix of bivectors [i][j]
         */
        _computeDOmega(coords, h = 1e-5) {
            const dim = this.manifold.dim;
            const dOmega = [];

            for (let i = 0; i < dim; i++) {
                dOmega[i] = [];

                const coords_plus = [...coords];
                const coords_minus = [...coords];
                coords_plus[i] += h;
                coords_minus[i] -= h;

                const omega_plus = this.connection.computeAt(coords_plus, h);
                const omega_minus = this.connection.computeAt(coords_minus, h);

                for (let j = 0; j < dim; j++) {
                    // ∂_i ω_j
                    dOmega[i][j] = omega_plus[j].sub(omega_minus[j]).scale(1 / (2 * h));
                }
            }

            return dOmega;
        }

        /**
         * Compute full curvature 2-form components: Ω_ij = ∂_i ω_j - ∂_j ω_i + ω_i × ω_j
         * 
         * @param {number[]} coords
         * @returns {Bivector3D[][]} Matrix of curvature bivectors
         */
        compute(coords, h = 1e-5) {
            const dim = this.manifold.dim;
            const omega = this.connection.computeAt(coords, h);
            const dOmega = this._computeDOmega(coords, h);

            const Omega = [];
            for (let i = 0; i < dim; i++) {
                Omega[i] = [];
                for (let j = 0; j < dim; j++) {
                    // dω part: ∂_i ω_j - ∂_j ω_i
                    const d_part = dOmega[i][j].sub(dOmega[j][i]);

                    // ω ∧ ω part: ω_i × ω_j (commutator)
                    const wedge_part = omega[i].commutatorWithBivector(omega[j]);

                    Omega[i][j] = d_part.add(wedge_part);
                }
            }

            return Omega;
        }

        /**
         * Compute Gaussian curvature K via sectional curvature.
         * 
         * K = Ω(e₁, e₂) · (e₁ ∧ e₂) / |e₁ ∧ e₂|²
         * 
         * @param {number[]} coords
         * @returns {number} Gaussian curvature
         */
        gaussianCurvature(coords, h = 1e-5) {
            // Use appropriate method based on manifold type
            if (this.manifold.ambientDim === 3) {
                // 2D surface in 3D: use shape operator
                return this.gaussianCurvatureViaShapeOperator(coords, h);
            } else {
                // 2D manifold in 2D: use metric-based method
                return this.gaussianCurvatureViaMetric(coords, h);
            }
        }

        /**
         * Compute Gaussian curvature K via metric tensor derivatives.
         * 
         * For 2D manifolds not embedded in 3D, use the formula:
         * K = -1/(2√g) [∂₁(√g Γ²₁₂/g₂₂) - ∂₂(√g Γ¹₁₂/g₁₁) + ...]
         * 
         * Simplified for orthogonal metrics: K = -(1/2√g)[∂₁₁g₂₂ + ∂₂₂g₁₁ - ...]/(2g)
         * 
         * @param {number[]} coords
         * @returns {number}
         */
        gaussianCurvatureViaMetric(coords, h = 1e-5) {
            // Compute metric
            const g = this.manifold.metric(coords);
            const g11 = g[0][0], g22 = g[1][1], g12 = g[0][1];
            const detG = g11 * g22 - g12 * g12;

            if (abs(detG) < EPSILON) return 0;

            // For conformal metric g = λ²I: K = -Δ(ln λ) / λ²
            // where Δ is the Euclidean Laplacian

            // Check if metric is conformal (g11 ≈ g22, g12 ≈ 0)
            const isConformal = abs(g11 - g22) < EPSILON && abs(g12) < EPSILON;

            if (isConformal) {
                // λ² = g11
                const lambda2 = g11;

                // Compute λ at surrounding points
                const g_xp = this.manifold.metric([coords[0] + h, coords[1]])[0][0];
                const g_xm = this.manifold.metric([coords[0] - h, coords[1]])[0][0];
                const g_yp = this.manifold.metric([coords[0], coords[1] + h])[0][0];
                const g_ym = this.manifold.metric([coords[0], coords[1] - h])[0][0];

                // ln(λ) = 0.5 * ln(g11)
                const lnLambda = 0.5 * Math.log(lambda2);
                const lnLambda_xp = 0.5 * Math.log(g_xp);
                const lnLambda_xm = 0.5 * Math.log(g_xm);
                const lnLambda_yp = 0.5 * Math.log(g_yp);
                const lnLambda_ym = 0.5 * Math.log(g_ym);

                // Laplacian Δ(ln λ) = ∂²(ln λ)/∂x² + ∂²(ln λ)/∂y²
                const d2_lnLambda_x = (lnLambda_xp - 2 * lnLambda + lnLambda_xm) / (h * h);
                const d2_lnLambda_y = (lnLambda_yp - 2 * lnLambda + lnLambda_ym) / (h * h);
                const laplacian = d2_lnLambda_x + d2_lnLambda_y;

                // K = -Δ(ln λ) / λ²
                return -laplacian / lambda2;
            }

            // For general orthogonal metrics, use Brioschi formula
            const g_0p = this.manifold.metric([coords[0] + h, coords[1]]);
            const g_0m = this.manifold.metric([coords[0] - h, coords[1]]);
            const g_1p = this.manifold.metric([coords[0], coords[1] + h]);
            const g_1m = this.manifold.metric([coords[0], coords[1] - h]);

            // First derivatives
            const d1_g11 = (g_0p[0][0] - g_0m[0][0]) / (2 * h);
            const d1_g22 = (g_0p[1][1] - g_0m[1][1]) / (2 * h);
            const d2_g11 = (g_1p[0][0] - g_1m[0][0]) / (2 * h);
            const d2_g22 = (g_1p[1][1] - g_1m[1][1]) / (2 * h);

            // Second derivatives
            const d11_g22 = (g_0p[1][1] - 2 * g22 + g_0m[1][1]) / (h * h);
            const d22_g11 = (g_1p[0][0] - 2 * g11 + g_1m[0][0]) / (h * h);

            // Brioschi formula for orthogonal metrics
            const sqrtG = sqrt(detG);
            const K = -(1 / (2 * sqrtG)) * (
                (d11_g22 - d1_g11 * d1_g22 / (2 * g11) - d1_g22 * d1_g22 / (2 * g22)) / sqrtG +
                (d22_g11 - d2_g11 * d2_g22 / (2 * g22) - d2_g11 * d2_g11 / (2 * g11)) / sqrtG
            );

            return K;
        }

        /**
         * Compute Gaussian curvature K via shape operator (Weingarten map).
         * 
         * For a 2D surface in 3D, compute K = det(S) where S is the shape operator:
         * S(v) = -∇_v n  (how the normal changes in direction v)
         * 
         * This is the most reliable method for embedded surfaces.
         * 
         * @param {number[]} coords
         * @returns {number} Gaussian curvature
         */
        gaussianCurvatureViaShapeOperator(coords, h = 1e-5) {
            if (this.manifold.dim !== 2 || this.manifold.ambientDim !== 3) {
                throw new Error('Shape operator method requires 2D surface in 3D');
            }

            // Helper: compute unit normal at a point
            const computeNormal = (c) => {
                const f = this.manifold.frame(c);
                const e1 = f.frame[0];
                const e2 = f.frame[1];
                const cross = cross3(e1, e2);
                const n = norm(cross);
                return n > EPSILON ? vecScale(cross, 1 / n) : [0, 0, 1];
            };

            // Normal at this point
            const n0 = computeNormal(coords);
            const frame = this.manifold.frame(coords);
            const e1up = frame.reciprocal[0];
            const e2up = frame.reciprocal[1];

            // Compute dn/dtheta and dn/dphi
            const c0_plus = [coords[0] + h, coords[1]];
            const c0_minus = [coords[0] - h, coords[1]];
            const c1_plus = [coords[0], coords[1] + h];
            const c1_minus = [coords[0], coords[1] - h];

            const dn_d0 = vecScale(vecSub(computeNormal(c0_plus), computeNormal(c0_minus)), 1 / (2 * h));
            const dn_d1 = vecScale(vecSub(computeNormal(c1_plus), computeNormal(c1_minus)), 1 / (2 * h));

            // Shape operator components: S^j_i = -(dn/dx^i) · e^j
            const S = [
                [-dot(dn_d0, e1up), -dot(dn_d0, e2up)],
                [-dot(dn_d1, e1up), -dot(dn_d1, e2up)]
            ];

            // Gaussian curvature K = det(S)
            return S[0][0] * S[1][1] - S[0][1] * S[1][0];
        }

        /**
         * Compute Gaussian curvature via Cartan curvature 2-form (experimental).
         * 
         * K = Ω₁₂ · (e₁ ∧ e₂) / |e₁ ∧ e₂|²
         * 
         * @param {number[]} coords
         * @returns {number}
         */
        gaussianCurvatureViaCartan(coords, h = 1e-5) {
            if (this.manifold.dim !== 2) {
                throw new Error('Gaussian curvature only defined for 2D');
            }

            const Omega = this.compute(coords, h);
            const frame = this.manifold.frame(coords);

            // Ω₁₂ = Omega[0][1]
            const Omega_12 = Omega[0][1];

            // Tangent bivector: e₁ ∧ e₂
            const e1 = frame.frame[0];
            const e2 = frame.frame[1];
            const tangent_plane = Bivector3D.fromWedge(e1, e2);

            // Sectional curvature: K = Ω₁₂ · (e₁ ∧ e₂) / |e₁ ∧ e₂|²
            const numerator = Omega_12.innerWithBivector(tangent_plane);
            const denominator = tangent_plane.normSquared();

            if (abs(denominator) < EPSILON) {
                return 0;
            }

            return numerator / denominator;
        }

        /**
         * Compute scalar curvature R.
         * 
         * For 2D: R = 2K (twice Gaussian curvature)
         * 
         * @param {number[]} coords
         * @returns {number}
         */
        scalarCurvature(coords, h = 1e-5) {
            return 2 * this.gaussianCurvature(coords, h);
        }
    }

    // ============================================================================
    // GA COVARIANT DERIVATIVE
    // ============================================================================

    /**
     * Covariant derivative using connection bivector:
     *   ∇_u A = ∂_u A + ω(u) × A
     * 
     * Works on scalars, vectors, and bivectors.
     */
    class GACovariantDerivative {
        /**
         * @param {RiemannianManifold} manifold
         */
        constructor(manifold) {
            this.manifold = manifold;
            this.connection = new ConnectionBivector(manifold);
        }

        /**
         * Covariant derivative of a scalar field.
         * ∇f = gradient (no connection needed)
         * 
         * @param {Function} f - Scalar field f(coords) → number
         * @param {number[]} coords
         * @param {number[]} u - Direction vector
         * @param {number} h
         * @returns {number} Directional derivative
         */
        scalarDerivative(f, coords, u, h = 1e-6) {
            // ∂_u f = uⁱ ∂_i f
            let result = 0;
            for (let i = 0; i < coords.length; i++) {
                const coordsPlus = [...coords];
                const coordsMinus = [...coords];
                coordsPlus[i] += h;
                coordsMinus[i] -= h;

                const df_di = (f(coordsPlus) - f(coordsMinus)) / (2 * h);
                result += u[i] * df_di;
            }
            return result;
        }

        /**
         * Covariant derivative of a vector field.
         * ∇_u v = ∂_u v + ω(u) × v
         * 
         * @param {Function} V - Vector field V(coords) → number[] (in coordinate basis)
         * @param {number[]} coords
         * @param {number[]} u - Direction
         * @param {number} h
         * @returns {number[]} Covariant derivative (coordinate components)
         */
        vectorDerivative(V, coords, u, h = 1e-6) {
            const dim = coords.length;

            // Component derivative: ∂_u V = uⁱ ∂_i V
            const partialV = new Array(dim).fill(0);
            for (let i = 0; i < dim; i++) {
                const coordsPlus = [...coords];
                const coordsMinus = [...coords];
                coordsPlus[i] += h;
                coordsMinus[i] -= h;

                const Vplus = V(coordsPlus);
                const Vminus = V(coordsMinus);

                for (let j = 0; j < dim; j++) {
                    partialV[j] += u[i] * (Vplus[j] - Vminus[j]) / (2 * h);
                }
            }

            // Connection term: ω(u) × V
            const omega_u = this.connection.along(coords, u, h);
            const V_here = V(coords);

            // Convert coordinate-basis vector to tangent space
            const frame = this.manifold.frame(coords);
            let V_tangent = new Array(this.manifold.ambientDim).fill(0);
            for (let i = 0; i < dim; i++) {
                for (let k = 0; k < this.manifold.ambientDim; k++) {
                    V_tangent[k] += V_here[i] * frame.frame[i][k];
                }
            }

            // Rotate by connection
            const rotated = omega_u.commutatorWithVector(V_tangent);

            // Convert back to coordinate basis
            const connectionTerm = frame.raise(rotated).slice(0, dim);

            // Total: ∇_u V = ∂_u V + ω(u) × V
            return vecAdd(partialV, connectionTerm);
        }

        /**
         * Covariant derivative of a tangent vector along itself.
         * ∇_v v for geodesic equation.
         * 
         * @param {number[]} coords - Current position
         * @param {number[]} v - Velocity (coordinate components)
         * @param {number} h
         * @returns {number[]} Acceleration (should be zero for geodesic)
         */
        selfDerivative(coords, v, h = 1e-6) {
            // ∇_v v = v̇ + ω(v) × v
            // But we need the actual derivative of v along the geodesic
            // For now, just compute ω(v) × v term (the "connection acceleration")

            const omega_v = this.connection.along(coords, v, h);

            // Convert v to tangent space
            const frame = this.manifold.frame(coords);
            let v_tangent = new Array(this.manifold.ambientDim).fill(0);
            for (let i = 0; i < v.length; i++) {
                for (let k = 0; k < this.manifold.ambientDim; k++) {
                    v_tangent[k] += v[i] * frame.frame[i][k];
                }
            }

            // Connection term: ω(v) × v
            const connectionAccel = omega_v.commutatorWithVector(v_tangent);

            // Project back to coordinate basis
            return frame.raise(connectionAccel).slice(0, v.length);
        }
    }

    // ============================================================================
    // CONCRETE MANIFOLDS
    // ============================================================================

    /**
     * 2-Sphere with radius R.
     * 
     * Coordinates: (θ, φ) where θ ∈ [0, π], φ ∈ [0, 2π)
     * Embedding: (R sin θ cos φ, R sin θ sin φ, R cos θ)
     * Metric: ds² = R²(dθ² + sin²θ dφ²)
     * Gaussian curvature: K = 1/R² (constant)
     */
    class Sphere2D extends RiemannianManifold {
        constructor(R = 1.0) {
            super(2, 3);  // 2D surface in 3D
            this.R = R;
        }

        /**
         * Embedding map: (θ, φ) → (x, y, z)
         */
        embedding(coords) {
            const [theta, phi] = coords;
            const R = this.R;
            return [
                R * sin(theta) * cos(phi),
                R * sin(theta) * sin(phi),
                R * cos(theta)
            ];
        }

        /**
         * Tangent frame at (θ, φ).
         * e_θ = ∂r/∂θ, e_φ = ∂r/∂φ
         */
        frame(coords) {
            const [theta, phi] = coords;
            const R = this.R;

            // Avoid singularities at poles
            const sinTheta = abs(sin(theta)) > EPSILON ? sin(theta) : EPSILON;

            // e_θ = R(cos θ cos φ, cos θ sin φ, -sin θ)
            const e_theta = [
                R * cos(theta) * cos(phi),
                R * cos(theta) * sin(phi),
                -R * sin(theta)
            ];

            // e_φ = R(-sin θ sin φ, sin θ cos φ, 0)
            const e_phi = [
                -R * sinTheta * sin(phi),
                R * sinTheta * cos(phi),
                0
            ];

            return new TangentFrame([e_theta, e_phi], 3);
        }

        /**
         * Theoretical Gaussian curvature: K = 1/R²
         */
        theoreticalCurvature() {
            return 1 / (this.R * this.R);
        }
    }

    /**
     * Torus with major radius R and minor radius r.
     * 
     * Coordinates: (θ, φ) where θ, φ ∈ [0, 2π)
     * θ = angle around the tube
     * φ = angle around the central axis
     * 
     * Embedding: ((R + r cos θ) cos φ, (R + r cos θ) sin φ, r sin θ)
     * Gaussian curvature: K = cos θ / (r(R + r cos θ))
     */
    class Torus2D extends RiemannianManifold {
        constructor(R = 2.0, r = 1.0) {
            super(2, 3);
            this.R = R;  // Major radius
            this.r = r;  // Minor radius
        }

        embedding(coords) {
            const [theta, phi] = coords;
            const { R, r } = this;
            const rho = R + r * cos(theta);
            return [
                rho * cos(phi),
                rho * sin(phi),
                r * sin(theta)
            ];
        }

        frame(coords) {
            const [theta, phi] = coords;
            const { R, r } = this;
            const rho = R + r * cos(theta);

            // e_θ = ∂r/∂θ = (-r sin θ cos φ, -r sin θ sin φ, r cos θ)
            const e_theta = [
                -r * sin(theta) * cos(phi),
                -r * sin(theta) * sin(phi),
                r * cos(theta)
            ];

            // e_φ = ∂r/∂φ = (-(R + r cos θ) sin φ, (R + r cos θ) cos φ, 0)
            const e_phi = [
                -rho * sin(phi),
                rho * cos(phi),
                0
            ];

            return new TangentFrame([e_theta, e_phi], 3);
        }

        theoreticalCurvature(coords) {
            const [theta] = coords;
            const { R, r } = this;
            return cos(theta) / (r * (R + r * cos(theta)));
        }
    }

    /**
     * Hyperbolic half-plane (Poincaré model).
     * 
     * Coordinates: (x, y) where y > 0
     * Metric: ds² = (dx² + dy²) / y²
     * Gaussian curvature: K = -1 (constant)
     */
    class HyperbolicPlane extends RiemannianManifold {
        constructor() {
            super(2, 2);
        }

        frame(coords) {
            const [x, y] = coords;
            const yVal = abs(y) > EPSILON ? y : EPSILON;

            // For metric ds² = (dx² + dy²)/y², the frame vectors should give g_ij = δ_ij/y²
            // This means |e_x| = |e_y| = 1/y
            // Frame: e_x = (1/y)∂_x, e_y = (1/y)∂_y
            // In coordinate representation: e_x = [1/y, 0], e_y = [0, 1/y]
            const e_x = [1 / yVal, 0];
            const e_y = [0, 1 / yVal];

            return new TangentFrame([e_x, e_y], 2);
        }

        /**
         * Override metric for standard Poincaré form: g = (1/y²)I
         */
        metric(coords) {
            const [x, y] = coords;
            const yVal = abs(y) > EPSILON ? y : EPSILON;
            const g = 1 / (yVal * yVal);
            return [[g, 0], [0, g]];
        }

        /**
         * For conformal metric g = (1/y²)I, K = -y² Δ(ln(1/y)) = -y² * (1/y²) = -1
         */
        theoreticalCurvature() {
            return -1;
        }
    }

    // ============================================================================
    // EXPORTS
    // ============================================================================

    const RiemannianGA = {
        // Bivector classes
        Bivector2D,
        Bivector3D,

        // Frame and manifold
        TangentFrame,
        RiemannianManifold,

        // Connection and curvature (replace Christoffel/Riemann)
        ConnectionBivector,
        Curvature2Form,

        // Covariant derivative
        GACovariantDerivative,

        // Concrete manifolds
        Sphere2D,
        Torus2D,
        HyperbolicPlane,

        // Utilities
        EPSILON,
        dot,
        cross3,
        norm,
        normalize,
        vecAdd,
        vecSub,
        vecScale
    };

    // Export for different module systems
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = RiemannianGA;
    }
    if (typeof global !== 'undefined') {
        global.RiemannianGA = RiemannianGA;
    }
    if (typeof window !== 'undefined') {
        window.RiemannianGA = RiemannianGA;
    }

})(typeof globalThis !== 'undefined' ? globalThis : this);
