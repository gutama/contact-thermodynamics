/**
 * Coordinate-Free Geodesic Solver via Geometric Calculus
 * 
 * Solves geodesic equations without indices using connection bivectors:
 *   ∇_v v = 0  ⟺  v̇ + ω(v) × v = 0
 * 
 * Also provides parallel transport and holonomy computations.
 * 
 * Key Features:
 *   - GAGeodesicSolver: Solve geodesics via RK4 integration
 *   - GAParallelTransport: Transport vectors preserving metric
 *   - Holonomy computation for Gauss-Bonnet verification
 * 
 * @module geodesic-ga
 * @license MIT
 */

(function (global) {
    'use strict';

    const EPSILON = 1e-10;
    const { abs, sqrt, sin, cos, PI } = Math;

    // Import RiemannianGA if available
    let RiemannianGA;
    if (typeof require !== 'undefined') {
        try {
            RiemannianGA = require('./riemannian-ga');
        } catch (e) {
            // Will use global
        }
    }
    if (!RiemannianGA && typeof global !== 'undefined') {
        RiemannianGA = global.RiemannianGA;
    }
    if (!RiemannianGA && typeof window !== 'undefined') {
        RiemannianGA = window.RiemannianGA;
    }

    // Utility functions (duplicated for standalone use)
    function dot(u, v) {
        let sum = 0;
        for (let i = 0; i < u.length; i++) sum += u[i] * v[i];
        return sum;
    }

    function norm(v) {
        return sqrt(dot(v, v));
    }

    function vecAdd(u, v) {
        return u.map((x, i) => x + v[i]);
    }

    function vecScale(v, s) {
        return v.map(x => x * s);
    }

    // ============================================================================
    // GA GEODESIC SOLVER
    // ============================================================================

    /**
     * Geodesic solver using connection bivector formulation.
     * 
     * Solves the geodesic equation without indices:
     *   ∇_v v = 0  ⟺  v̇ + ω(v) × v = 0
     * 
     * Uses RK4 integration in coordinate space.
     */
    class GAGeodesicSolver {
        /**
         * @param {RiemannianManifold} manifold
         */
        constructor(manifold) {
            this.manifold = manifold;
            this.connection = new (RiemannianGA.ConnectionBivector)(manifold);
        }

        /**
         * Compute the geodesic acceleration (should be zero for geodesics).
         * 
         * The second-order ODE: ẍ + Γ(ẋ,ẋ) = 0
         * In GA terms: v̇ = -ω(v) × v
         * 
         * @param {number[]} x - Position (coordinates)
         * @param {number[]} v - Velocity (coordinate components)
         * @returns {number[]} Acceleration dv/dt
         */
        _acceleration(x, v, h = 1e-6) {
            // Get connection along velocity
            const omega_v = this.connection.along(x, v, h);

            // Convert v to tangent space vector
            const frame = this.manifold.frame(x);
            let v_tangent = new Array(this.manifold.ambientDim).fill(0);
            for (let i = 0; i < v.length; i++) {
                for (let k = 0; k < this.manifold.ambientDim; k++) {
                    v_tangent[k] += v[i] * frame.frame[i][k];
                }
            }

            // ω(v) × v (rotation of v by connection)
            const rotated = omega_v.commutatorWithVector(v_tangent);

            // Project back to coordinate basis
            const connectionAccel = frame.raise(rotated).slice(0, v.length);

            // Acceleration = -ω(v) × v
            return connectionAccel.map(x => -x);
        }

        /**
         * Single RK4 step for geodesic integration.
         * 
         * State: [x, v] where x is position, v is velocity
         * 
         * @param {number[]} x - Position
         * @param {number[]} v - Velocity
         * @param {number} dt - Time step
         * @returns {Object} {x: new position, v: new velocity}
         */
        _rk4Step(x, v, dt) {
            const dim = x.length;

            // k1
            const a1 = this._acceleration(x, v);
            const k1_x = v;
            const k1_v = a1;

            // k2 (midpoint)
            const x2 = vecAdd(x, vecScale(k1_x, dt / 2));
            const v2 = vecAdd(v, vecScale(k1_v, dt / 2));
            const a2 = this._acceleration(x2, v2);
            const k2_x = v2;
            const k2_v = a2;

            // k3 (midpoint)
            const x3 = vecAdd(x, vecScale(k2_x, dt / 2));
            const v3 = vecAdd(v, vecScale(k2_v, dt / 2));
            const a3 = this._acceleration(x3, v3);
            const k3_x = v3;
            const k3_v = a3;

            // k4 (endpoint)
            const x4 = vecAdd(x, vecScale(k3_x, dt));
            const v4 = vecAdd(v, vecScale(k3_v, dt));
            const a4 = this._acceleration(x4, v4);
            const k4_x = v4;
            const k4_v = a4;

            // Combine: y_new = y + (dt/6)(k1 + 2k2 + 2k3 + k4)
            const x_new = [];
            const v_new = [];
            for (let i = 0; i < dim; i++) {
                x_new[i] = x[i] + (dt / 6) * (k1_x[i] + 2 * k2_x[i] + 2 * k3_x[i] + k4_x[i]);
                v_new[i] = v[i] + (dt / 6) * (k1_v[i] + 2 * k2_v[i] + 2 * k3_v[i] + k4_v[i]);
            }

            return { x: x_new, v: v_new };
        }

        /**
         * Solve geodesic equation from initial conditions.
         * 
         * @param {number[]} x0 - Initial position (coordinates)
         * @param {number[]} v0 - Initial velocity (coordinate components)
         * @param {number} tFinal - Final time
         * @param {number} dt - Time step (default: auto from CFL-like estimate)
         * @returns {Object} {t: times, x: positions, v: velocities}
         */
        solve(x0, v0, tFinal, dt = null) {
            // Auto-select dt if not provided
            if (dt === null) {
                dt = tFinal / 100;  // 100 steps by default
            }

            const nSteps = Math.ceil(tFinal / dt);
            const actualDt = tFinal / nSteps;

            const ts = [0];
            const xs = [x0.slice()];
            const vs = [v0.slice()];

            let x = x0.slice();
            let v = v0.slice();

            for (let step = 0; step < nSteps; step++) {
                const result = this._rk4Step(x, v, actualDt);
                x = result.x;
                v = result.v;

                ts.push((step + 1) * actualDt);
                xs.push(x.slice());
                vs.push(v.slice());
            }

            return { t: ts, x: xs, v: vs };
        }

        /**
         * Compute the arc length of velocity (should be preserved for geodesics).
         * 
         * Uses the metric: |v|² = g_ij v^i v^j
         * 
         * @param {number[]} x - Position
         * @param {number[]} v - Velocity
         * @returns {number} Metric norm of velocity
         */
        velocityNorm(x, v) {
            const g = this.manifold.metric(x);
            let normSq = 0;
            for (let i = 0; i < v.length; i++) {
                for (let j = 0; j < v.length; j++) {
                    normSq += g[i][j] * v[i] * v[j];
                }
            }
            return sqrt(abs(normSq));
        }

        /**
         * Verify that velocity norm is preserved along geodesic.
         * 
         * @param {Object} solution - Output from solve()
         * @returns {Object} {norms: array, maxDeviation: number}
         */
        verifyArcLengthPreservation(solution) {
            const norms = [];
            for (let i = 0; i < solution.x.length; i++) {
                norms.push(this.velocityNorm(solution.x[i], solution.v[i]));
            }

            const initialNorm = norms[0];
            const maxDeviation = Math.max(...norms.map(n => abs(n - initialNorm) / initialNorm));

            return { norms, maxDeviation };
        }
    }

    // ============================================================================
    // GA PARALLEL TRANSPORT
    // ============================================================================

    /**
     * Parallel transport of vectors along curves.
     * 
     * Solves: ẇ + ω(v) × w = 0
     * where v is the curve velocity and w is the transported vector.
     * 
     * Preserves the metric norm of the transported vector.
     */
    class GAParallelTransport {
        /**
         * @param {RiemannianManifold} manifold
         */
        constructor(manifold) {
            this.manifold = manifold;
            this.connection = new (RiemannianGA.ConnectionBivector)(manifold);
        }

        /**
         * Compute the transport derivative: ∂w/∂t = -ω(v) × w
         * 
         * @param {number[]} x - Position on curve
         * @param {number[]} v - Curve velocity
         * @param {number[]} w - Vector being transported (coordinate components)
         * @returns {number[]} Time derivative of w
         */
        _transportDerivative(x, v, w, h = 1e-6) {
            const omega_v = this.connection.along(x, v, h);

            // Convert w to tangent space
            const frame = this.manifold.frame(x);
            let w_tangent = new Array(this.manifold.ambientDim).fill(0);
            for (let i = 0; i < w.length; i++) {
                for (let k = 0; k < this.manifold.ambientDim; k++) {
                    w_tangent[k] += w[i] * frame.frame[i][k];
                }
            }

            // ω(v) × w
            const rotated = omega_v.commutatorWithVector(w_tangent);

            // Project back to coordinate basis
            const result = frame.raise(rotated).slice(0, w.length);

            // dw/dt = -ω(v) × w
            return result.map(x => -x);
        }

        /**
         * RK4 step for parallel transport.
         * 
         * @param {number[]} x - Position
         * @param {number[]} v - Curve velocity
         * @param {number[]} w - Vector being transported
         * @param {number} dt
         * @returns {number[]} Updated w
         */
        _rk4Step(x, v, w, dt) {
            // k1
            const k1 = this._transportDerivative(x, v, w);

            // k2
            const w2 = vecAdd(w, vecScale(k1, dt / 2));
            const k2 = this._transportDerivative(x, v, w2);

            // k3
            const w3 = vecAdd(w, vecScale(k2, dt / 2));
            const k3 = this._transportDerivative(x, v, w3);

            // k4
            const w4 = vecAdd(w, vecScale(k3, dt));
            const k4 = this._transportDerivative(x, v, w4);

            // Combine
            const w_new = [];
            for (let i = 0; i < w.length; i++) {
                w_new[i] = w[i] + (dt / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
            }
            return w_new;
        }

        /**
         * Transport a vector along a parametrized curve.
         * 
         * @param {Function} curve - Curve γ(t) → coordinates
         * @param {number[]} w0 - Initial vector at γ(t0)
         * @param {number} t0 - Start parameter
         * @param {number} t1 - End parameter
         * @param {number} nSteps - Number of integration steps
         * @returns {number[]} Transported vector at γ(t1)
         */
        transport(curve, w0, t0, t1, nSteps = 100) {
            const dt = (t1 - t0) / nSteps;
            const h = abs(dt) * 0.01;  // For numerical derivatives

            let w = w0.slice();
            let t = t0;

            for (let step = 0; step < nSteps; step++) {
                const x = curve(t);

                // Compute curve velocity via finite difference
                const xPlus = curve(t + h);
                const xMinus = curve(t - h);
                const v = vecScale(vecAdd(xPlus, vecScale(xMinus, -1)), 1 / (2 * h));

                // Transport step
                w = this._rk4Step(x, v, w, dt);
                t += dt;
            }

            return w;
        }

        /**
         * Transport along a geodesic solution.
         * 
         * @param {Object} geodesic - Output from GAGeodesicSolver.solve()
         * @param {number[]} w0 - Initial vector
         * @returns {Object} {vectors: transported vectors at each time}
         */
        transportAlongGeodesic(geodesic, w0) {
            const { t, x, v } = geodesic;
            const vectors = [w0.slice()];

            let w = w0.slice();

            for (let i = 1; i < t.length; i++) {
                const dt = t[i] - t[i - 1];

                // Use position and velocity from geodesic
                w = this._rk4Step(x[i - 1], v[i - 1], w, dt);
                vectors.push(w.slice());
            }

            return { vectors };
        }

        /**
         * Compute holonomy angle for transport around a closed loop.
         * 
         * The holonomy measures how much a vector rotates when transported
         * around a closed curve. For a sphere, this equals the enclosed
         * solid angle times the curvature.
         * 
         * @param {Function} loop - Closed curve γ(t): [0, 1] → M with γ(0) = γ(1)
         * @param {number[]} w0 - Initial vector
         * @param {number} nSteps
         * @returns {number} Rotation angle in radians
         */
        holonomyAngle(loop, w0, nSteps = 200) {
            // Transport around the loop
            const w1 = this.transport(loop, w0, 0, 1, nSteps);

            // Get frame at starting point
            const x0 = loop(0);
            const frame = this.manifold.frame(x0);

            // Convert both vectors to tangent space
            let w0_tangent = new Array(this.manifold.ambientDim).fill(0);
            let w1_tangent = new Array(this.manifold.ambientDim).fill(0);

            for (let i = 0; i < w0.length; i++) {
                for (let k = 0; k < this.manifold.ambientDim; k++) {
                    w0_tangent[k] += w0[i] * frame.frame[i][k];
                    w1_tangent[k] += w1[i] * frame.frame[i][k];
                }
            }

            // Compute angle between w0 and w1 in tangent plane
            const norm0 = norm(w0_tangent);
            const norm1 = norm(w1_tangent);

            if (norm0 < EPSILON || norm1 < EPSILON) {
                return 0;
            }

            // Dot product gives cos(angle)
            const cosAngle = dot(w0_tangent, w1_tangent) / (norm0 * norm1);

            // For 2D, we can also get the sign from the cross product
            // In tangent plane, this projects to the normal
            if (this.manifold.dim === 2 && this.manifold.ambientDim === 3) {
                // Cross product to get signed angle
                const cross = [
                    w0_tangent[1] * w1_tangent[2] - w0_tangent[2] * w1_tangent[1],
                    w0_tangent[2] * w1_tangent[0] - w0_tangent[0] * w1_tangent[2],
                    w0_tangent[0] * w1_tangent[1] - w0_tangent[1] * w1_tangent[0]
                ];

                // Normal vector
                const e1 = frame.frame[0];
                const e2 = frame.frame[1];
                const normal = [
                    e1[1] * e2[2] - e1[2] * e2[1],
                    e1[2] * e2[0] - e1[0] * e2[2],
                    e1[0] * e2[1] - e1[1] * e2[0]
                ];

                // Sign from dot(cross, normal)
                const sign = dot(cross, normal) >= 0 ? 1 : -1;

                // Use atan2 for full range
                const sinAngle = sign * norm(cross) / (norm0 * norm1);
                return Math.atan2(sinAngle, Math.max(-1, Math.min(1, cosAngle)));
            }

            // For 2D manifolds in 2D, use atan2 directly
            if (this.manifold.dim === 2 && this.manifold.ambientDim === 2) {
                const cross2D = w0_tangent[0] * w1_tangent[1] - w0_tangent[1] * w1_tangent[0];
                return Math.atan2(cross2D, cosAngle * norm0 * norm1);
            }

            return Math.acos(Math.max(-1, Math.min(1, cosAngle)));
        }

        /**
         * Verify holonomy matches Gauss-Bonnet prediction.
         * 
         * For a region Σ on a 2D surface:
         *   holonomy = ∫∫_Σ K dA
         * 
         * For a latitude circle at θ on a unit sphere:
         *   holonomy = 2π(1 - cos θ)
         * 
         * @param {Function} loop - Closed curve
         * @param {number} expectedHolonomy - Theoretical value
         * @param {number[]} w0 - Initial vector
         * @returns {Object} {computed, expected, error}
         */
        verifyGaussBonnet(loop, expectedHolonomy, w0, nSteps = 200) {
            const computed = this.holonomyAngle(loop, w0, nSteps);
            const error = abs(computed - expectedHolonomy);

            return {
                computed,
                expected: expectedHolonomy,
                error,
                relativeError: abs(expectedHolonomy) > EPSILON ? error / abs(expectedHolonomy) : error
            };
        }
    }

    // ============================================================================
    // GAUSS-BONNET INTEGRATOR
    // ============================================================================

    /**
     * Numerical integration of Gaussian curvature for Gauss-Bonnet verification.
     * 
     * ∫∫_M K dA = 2π χ(M)
     * 
     * For sphere: ∫∫ K dA = 4π (χ = 2)
     * For torus: ∫∫ K dA = 0 (χ = 0)
     */
    class GaussBonnetIntegrator {
        /**
         * @param {RiemannianManifold} manifold
         */
        constructor(manifold) {
            this.manifold = manifold;
            this.curvature = new (RiemannianGA.Curvature2Form)(manifold);
        }

        /**
         * Integrate K dA over a rectangular region in parameter space.
         * 
         * @param {number[]} uRange - [u_min, u_max]
         * @param {number[]} vRange - [v_min, v_max]
         * @param {number} nU - Grid points in u
         * @param {number} nV - Grid points in v
         * @returns {number} Integral value
         */
        integrate(uRange, vRange, nU = 50, nV = 50) {
            const [uMin, uMax] = uRange;
            const [vMin, vMax] = vRange;
            const du = (uMax - uMin) / nU;
            const dv = (vMax - vMin) / nV;

            let integral = 0;

            for (let i = 0; i < nU; i++) {
                for (let j = 0; j < nV; j++) {
                    const u = uMin + (i + 0.5) * du;
                    const v = vMin + (j + 0.5) * dv;
                    const coords = [u, v];

                    // Get Gaussian curvature
                    const K = this.curvature.gaussianCurvature(coords);

                    // Get area element √(det g)
                    const g = this.manifold.metric(coords);
                    const detG = g[0][0] * g[1][1] - g[0][1] * g[1][0];
                    const dA = sqrt(abs(detG)) * du * dv;

                    integral += K * dA;
                }
            }

            return integral;
        }

        /**
         * Verify Gauss-Bonnet for a closed manifold.
         * 
         * @param {number} eulerCharacteristic
         * @param {number[]} uRange
         * @param {number[]} vRange
         * @param {number} nU
         * @param {number} nV
         * @returns {Object}
         */
        verifyGaussBonnet(eulerCharacteristic, uRange, vRange, nU = 50, nV = 50) {
            const integral = this.integrate(uRange, vRange, nU, nV);
            const expected = 2 * PI * eulerCharacteristic;

            return {
                integral,
                expected,
                error: abs(integral - expected),
                relativeError: abs(expected) > EPSILON ? abs(integral - expected) / abs(expected) : abs(integral)
            };
        }
    }

    // ============================================================================
    // EXPORTS
    // ============================================================================

    const GeodesicGA = {
        GAGeodesicSolver,
        GAParallelTransport,
        GaussBonnetIntegrator,

        // Utility
        EPSILON
    };

    // Export for different module systems
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = GeodesicGA;
    }
    if (typeof global !== 'undefined') {
        global.GeodesicGA = GeodesicGA;
    }
    if (typeof window !== 'undefined') {
        window.GeodesicGA = GeodesicGA;
    }

})(typeof globalThis !== 'undefined' ? globalThis : this);
