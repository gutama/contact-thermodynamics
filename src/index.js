/**
 * Geometric Mechanics of Extended Thermodynamics
 * Framework: Contact Geometry / 1-Jet Bundles
 * 
 * Implements:
 * - 1-jet bundle J¹(Q) with canonical contact form α = du - p_a dx^a
 * - Grand model M₁₃ = J¹(Q₆) with dim = 13
 * - Holographic model M₇ = J¹(Q₃) with dim = 7
 * - Contact Hamiltonian dynamics
 * - Reeb vector fields
 * - Legendrian submanifolds
 * - Gravitational extension (GR coupling)
 * 
 * Based on the reference sheet for extended thermodynamics
 * Using Geometric Algebra conventions where applicable
 * 
 * @license MIT
 */

(function (global) {
    'use strict';

    const EPSILON = 1e-12;
    const abs = Math.abs;
    const sqrt = Math.sqrt;
    const sin = Math.sin;
    const cos = Math.cos;
    const exp = Math.exp;
    const log = Math.log;

    function isZero(x) { return abs(x) < EPSILON; }

    // ============================================================================
    // I. CONTACT MANIFOLDS (imported from ./contact/manifold.js)
    // ============================================================================
    //
    // ContactManifold, ContactPoint, GrandContactManifold and
    // HolographicContactManifold are defined once in the contact namespace and
    // consumed here so there is a single source of truth for the manifold classes.

    let ContactManifold, ContactPoint, GrandContactManifold, HolographicContactManifold;
    // Gate the synchronous require() on `module.exports` (CommonJS/Node) rather
    // than `typeof require`, which is also truthy under AMD loaders (e.g.
    // RequireJS) whose `require` is not a synchronous CommonJS loader and would
    // throw here before the UMD `define.amd` branch runs.
    if (typeof module !== 'undefined' && module.exports) {
        const Manifold = require('./contact/manifold.js');
        ContactManifold = Manifold.ContactManifold;
        ContactPoint = Manifold.ContactPoint;
        GrandContactManifold = Manifold.GrandContactManifold;
        HolographicContactManifold = Manifold.HolographicContactManifold;
    } else if (typeof global !== 'undefined' && global.ContactManifold) {
        ({ ContactManifold, ContactPoint, GrandContactManifold, HolographicContactManifold } = global.ContactManifold);
    }

    // ============================================================================
    // IV. CONTACT HAMILTONIAN DYNAMICS (imported from ./contact/hamiltonian.js)
    // V.  LEGENDRIAN SUBMANIFOLDS   (imported from ./contact/legendrian.js)
    // ============================================================================
    //
    // ContactHamiltonian and LegendrianSubmanifold live in the contact namespace
    // and are consumed here so there is a single source of truth. See the note on
    // the manifold import above for why the require() is gated on module.exports.

    let ContactHamiltonian, LegendrianSubmanifold;
    if (typeof module !== 'undefined' && module.exports) {
        ContactHamiltonian = require('./contact/hamiltonian.js').ContactHamiltonian;
        LegendrianSubmanifold = require('./contact/legendrian.js').LegendrianSubmanifold;
    } else if (typeof global !== 'undefined' && global.ContactHamiltonian) {
        ({ ContactHamiltonian } = global.ContactHamiltonian);
        ({ LegendrianSubmanifold } = global.LegendrianSubmanifold);
    }

    // ============================================================================
    // VI. GRAVITATIONAL EXTENSION (GENERAL RELATIVITY)
    // ============================================================================

    /**
     * SpacetimeMetric: Encapsulates g_μν(x) for GR coupling
     * 
     * Kinematics from α (contact structure)
     * Curvature from g_μν inside H
     */
    class SpacetimeMetric {
        /**
         * @param {Function} metricFunc - (x) => 4x4 covariant metric g_μν
         * @param {Function} [inverseFunc] - (x) => 4x4 contravariant metric g^μν
         */
        constructor(metricFunc, inverseFunc = null) {
            this.g = metricFunc;
            this._gInv = inverseFunc;
        }

        /**
         * Covariant metric g_μν at point x
         */
        covariant(x) {
            return this.g(x);
        }

        /**
         * Contravariant metric g^μν at point x (numerical inversion if not provided)
         */
        contravariant(x) {
            if (this._gInv) {
                return this._gInv(x);
            }
            return this.invertMatrix(this.g(x));
        }

        /**
         * 4x4 matrix inversion
         */
        invertMatrix(m) {
            // Gaussian elimination for 4x4
            const n = 4;
            const aug = m.map((row, i) => [...row, ...Array(n).fill(0).map((_, j) => i === j ? 1 : 0)]);

            for (let col = 0; col < n; col++) {
                // Find pivot
                let maxRow = col;
                for (let row = col + 1; row < n; row++) {
                    if (abs(aug[row][col]) > abs(aug[maxRow][col])) maxRow = row;
                }
                [aug[col], aug[maxRow]] = [aug[maxRow], aug[col]];

                if (isZero(aug[col][col])) {
                    throw new Error('Singular metric');
                }

                // Scale
                const scale = aug[col][col];
                for (let j = 0; j < 2 * n; j++) aug[col][j] /= scale;

                // Eliminate
                for (let row = 0; row < n; row++) {
                    if (row !== col) {
                        const factor = aug[row][col];
                        for (let j = 0; j < 2 * n; j++) {
                            aug[row][j] -= factor * aug[col][j];
                        }
                    }
                }
            }

            return aug.map(row => row.slice(n));
        }

        /**
         * Minkowski metric (flat spacetime, signature +---)
         */
        static minkowski() {
            return new SpacetimeMetric(
                x => [
                    [1, 0, 0, 0],
                    [0, -1, 0, 0],
                    [0, 0, -1, 0],
                    [0, 0, 0, -1]
                ],
                x => [
                    [1, 0, 0, 0],
                    [0, -1, 0, 0],
                    [0, 0, -1, 0],
                    [0, 0, 0, -1]
                ]
            );
        }

        /**
         * Schwarzschild metric (spherical coords: t, r, θ, φ)
         * @param {number} M - Mass parameter (in geometric units, G=c=1)
         */
        static schwarzschild(M = 1) {
            return new SpacetimeMetric(x => {
                const [t, r, theta, phi] = x;
                const f = 1 - 2 * M / r;
                const sinTheta2 = sin(theta) ** 2;
                return [
                    [f, 0, 0, 0],
                    [0, -1 / f, 0, 0],
                    [0, 0, -r * r, 0],
                    [0, 0, 0, -r * r * sinTheta2]
                ];
            });
        }

        /**
         * FLRW metric (cosmological, comoving coords: t, χ, θ, φ)
         * @param {Function} a - Scale factor a(t)
         * @param {number} k - Curvature: +1 (closed), 0 (flat), -1 (open)
         */
        static flrw(a, k = 0) {
            return new SpacetimeMetric(x => {
                const [t, chi, theta, phi] = x;
                const aT = a(t);
                const a2 = aT * aT;
                let Sk2;
                if (k === 0) Sk2 = chi * chi;
                else if (k === 1) Sk2 = sin(chi) ** 2;
                else Sk2 = (exp(chi) - exp(-chi)) / 2; // sinh²

                const sinTheta2 = sin(theta) ** 2;
                return [
                    [1, 0, 0, 0],
                    [0, -a2, 0, 0],
                    [0, 0, -a2 * Sk2, 0],
                    [0, 0, 0, -a2 * Sk2 * sinTheta2]
                ];
            });
        }
    }

    // ============================================================================
    // VI-B. CHRISTOFFEL SYMBOLS
    // ============================================================================

    /**
     * ChristoffelSymbols: Compute connection coefficients from a metric
     * 
     * Γᵏᵢⱼ = ½ gᵏˡ(∂ᵢgⱼˡ + ∂ⱼgᵢˡ - ∂ˡgᵢⱼ)
     * 
     * Essential for:
     * - Geodesic equation: ẍ^μ + Γ^μ_αβ ẋ^α ẋ^β = 0
     * - Covariant derivatives
     * - Parallel transport
     */
    class ChristoffelSymbols {
        /**
         * @param {SpacetimeMetric} metric - The spacetime metric
         * @param {number} h - Step size for numerical differentiation
         */
        constructor(metric, h = 1e-7) {
            this.metric = metric;
            this.h = h;
            this.dim = 4; // Spacetime dimension
        }

        /**
         * Compute ∂gᵢⱼ/∂xˡ using central differences
         */
        _partialMetric(x, l) {
            const xPlus = [...x];
            const xMinus = [...x];
            xPlus[l] += this.h;
            xMinus[l] -= this.h;

            const gPlus = this.metric.covariant(xPlus);
            const gMinus = this.metric.covariant(xMinus);

            const result = [];
            for (let i = 0; i < this.dim; i++) {
                result[i] = [];
                for (let j = 0; j < this.dim; j++) {
                    result[i][j] = (gPlus[i][j] - gMinus[i][j]) / (2 * this.h);
                }
            }
            return result;
        }

        /**
         * Compute all Christoffel symbols Γᵏᵢⱼ at given coordinates
         * 
         * @param {number[]} x - Spacetime coordinates [x⁰, x¹, x², x³]
         * @returns {number[][][]} Γ[k][i][j] = Γᵏᵢⱼ
         */
        computeAt(x) {
            const gInv = this.metric.contravariant(x);

            // Precompute all metric partial derivatives
            const partials = [];
            for (let l = 0; l < this.dim; l++) {
                partials[l] = this._partialMetric(x, l);
            }

            // Compute Γᵏᵢⱼ = ½ gᵏˡ(∂ᵢgⱼˡ + ∂ⱼgᵢˡ - ∂ˡgᵢⱼ)
            const gamma = [];
            for (let k = 0; k < this.dim; k++) {
                gamma[k] = [];
                for (let i = 0; i < this.dim; i++) {
                    gamma[k][i] = [];
                    for (let j = 0; j < this.dim; j++) {
                        let sum = 0;
                        for (let l = 0; l < this.dim; l++) {
                            sum += gInv[k][l] * (
                                partials[i][j][l] +  // ∂ᵢgⱼˡ
                                partials[j][i][l] -  // ∂ⱼgᵢˡ
                                partials[l][i][j]    // ∂ˡgᵢⱼ
                            );
                        }
                        gamma[k][i][j] = 0.5 * sum;
                    }
                }
            }
            return gamma;
        }

        /**
         * Compute Riemann curvature tensor R^ρ_σμν (for verification)
         * 
         * R^ρ_σμν = ∂_μΓ^ρ_νσ - ∂_νΓ^ρ_μσ + Γ^ρ_μλ Γ^λ_νσ - Γ^ρ_νλ Γ^λ_μσ
         */
        riemannAt(x) {
            const h = this.h;
            const gammaAt = this.computeAt(x);

            // Compute partial derivatives of Christoffel symbols
            const dGamma = [];
            for (let mu = 0; mu < this.dim; mu++) {
                const xPlus = [...x]; xPlus[mu] += h;
                const xMinus = [...x]; xMinus[mu] -= h;
                const gammaPlus = this.computeAt(xPlus);
                const gammaMinus = this.computeAt(xMinus);

                dGamma[mu] = [];
                for (let rho = 0; rho < this.dim; rho++) {
                    dGamma[mu][rho] = [];
                    for (let alpha = 0; alpha < this.dim; alpha++) {
                        dGamma[mu][rho][alpha] = [];
                        for (let beta = 0; beta < this.dim; beta++) {
                            dGamma[mu][rho][alpha][beta] =
                                (gammaPlus[rho][alpha][beta] - gammaMinus[rho][alpha][beta]) / (2 * h);
                        }
                    }
                }
            }

            // Compute Riemann tensor
            const R = [];
            for (let rho = 0; rho < this.dim; rho++) {
                R[rho] = [];
                for (let sigma = 0; sigma < this.dim; sigma++) {
                    R[rho][sigma] = [];
                    for (let mu = 0; mu < this.dim; mu++) {
                        R[rho][sigma][mu] = [];
                        for (let nu = 0; nu < this.dim; nu++) {
                            let val = dGamma[mu][rho][nu][sigma] - dGamma[nu][rho][mu][sigma];
                            for (let lambda = 0; lambda < this.dim; lambda++) {
                                val += gammaAt[rho][mu][lambda] * gammaAt[lambda][nu][sigma];
                                val -= gammaAt[rho][nu][lambda] * gammaAt[lambda][mu][sigma];
                            }
                            R[rho][sigma][mu][nu] = val;
                        }
                    }
                }
            }
            return R;
        }

        /**
         * Compute Ricci tensor R_μν
         */
        ricciAt(x) {
            const R = this.riemannAt(x);
            const Ric = [];
            for (let mu = 0; mu < this.dim; mu++) {
                Ric[mu] = [];
                for (let nu = 0; nu < this.dim; nu++) {
                    let sum = 0;
                    for (let rho = 0; rho < this.dim; rho++) {
                        sum += R[rho][mu][rho][nu];
                    }
                    Ric[mu][nu] = sum;
                }
            }
            return Ric;
        }

        /**
         * Compute Ricci scalar R
         */
        ricciScalarAt(x) {
            const Ric = this.ricciAt(x);
            const gInv = this.metric.contravariant(x);
            let R = 0;
            for (let mu = 0; mu < this.dim; mu++) {
                for (let nu = 0; nu < this.dim; nu++) {
                    R += gInv[mu][nu] * Ric[mu][nu];
                }
            }
            return R;
        }
    }

    // ============================================================================
    // VI-C. COVARIANT DERIVATIVE
    // ============================================================================

    /**
     * CovariantDerivative: Covariant derivative operator on spacetime
     * 
     * For a vector field V^j:
     *   ∇ᵢV^j = ∂ᵢV^j + Γʲᵢₖ V^k
     * 
     * For a covector field (1-form) ω_j:
     *   ∇ᵢω_j = ∂ᵢω_j - Γᵏᵢⱼ ω_k
     */
    class CovariantDerivative {
        /**
         * @param {SpacetimeMetric} metric
         * @param {number} h - Step size for numerical differentiation
         */
        constructor(metric, h = 1e-7) {
            this.metric = metric;
            this.christoffel = new ChristoffelSymbols(metric, h);
            this.h = h;
            this.dim = 4;
        }

        /**
         * Compute covariant derivative of contravariant vector field
         * ∇ᵢV^j = ∂ᵢV^j + Γʲᵢₖ V^k
         * 
         * @param {Function} V_func - (x) => [V⁰, V¹, V², V³]
         * @param {number[]} x - Point to evaluate at
         * @returns {number[][]} ∇V[i][j] = ∇ᵢV^j
         */
        ofVector(V_func, x) {
            const gamma = this.christoffel.computeAt(x);
            const V = V_func(x);

            const nablaV = [];
            for (let i = 0; i < this.dim; i++) {
                nablaV[i] = [];

                // Compute ∂ᵢV^j
                const xPlus = [...x]; xPlus[i] += this.h;
                const xMinus = [...x]; xMinus[i] -= this.h;
                const VPlus = V_func(xPlus);
                const VMinus = V_func(xMinus);

                for (let j = 0; j < this.dim; j++) {
                    let term = (VPlus[j] - VMinus[j]) / (2 * this.h);  // ∂ᵢV^j
                    for (let k = 0; k < this.dim; k++) {
                        term += gamma[j][i][k] * V[k];  // Γʲᵢₖ V^k
                    }
                    nablaV[i][j] = term;
                }
            }
            return nablaV;
        }

        /**
         * Divergence of vector field: div V = ∇ᵢV^i
         * 
         * Uses formula: div V = (1/√g) ∂ᵢ(√g V^i)
         */
        divergence(V_func, x) {
            const g = this.metric.covariant(x);
            const detG = this._det4x4(g);
            const sqrtG = sqrt(abs(detG));

            let div = 0;
            for (let i = 0; i < this.dim; i++) {
                const xPlus = [...x]; xPlus[i] += this.h;
                const xMinus = [...x]; xMinus[i] -= this.h;

                const gPlus = this.metric.covariant(xPlus);
                const gMinus = this.metric.covariant(xMinus);
                const sqrtGPlus = sqrt(abs(this._det4x4(gPlus)));
                const sqrtGMinus = sqrt(abs(this._det4x4(gMinus)));

                const VPlus = V_func(xPlus);
                const VMinus = V_func(xMinus);

                div += (sqrtGPlus * VPlus[i] - sqrtGMinus * VMinus[i]) / (2 * this.h);
            }
            return div / sqrtG;
        }

        /**
         * Laplace-Beltrami operator on scalar field
         * Δf = div(grad f) = (1/√g) ∂ᵢ(√g g^ij ∂ⱼf)
         */
        laplacian(f_func, x) {
            // Gradient: (grad f)^i = g^ij ∂ⱼf
            const grad_f = (coords) => {
                const gInv = this.metric.contravariant(coords);
                const df = [];
                for (let j = 0; j < this.dim; j++) {
                    const xPlus = [...coords]; xPlus[j] += this.h;
                    const xMinus = [...coords]; xMinus[j] -= this.h;
                    df[j] = (f_func(xPlus) - f_func(xMinus)) / (2 * this.h);
                }
                const gradF = [];
                for (let i = 0; i < this.dim; i++) {
                    gradF[i] = 0;
                    for (let j = 0; j < this.dim; j++) {
                        gradF[i] += gInv[i][j] * df[j];
                    }
                }
                return gradF;
            };

            return this.divergence(grad_f, x);
        }

        _det4x4(m) {
            // 4x4 determinant using cofactor expansion
            return (
                m[0][0] * (m[1][1] * (m[2][2] * m[3][3] - m[2][3] * m[3][2]) - m[1][2] * (m[2][1] * m[3][3] - m[2][3] * m[3][1]) + m[1][3] * (m[2][1] * m[3][2] - m[2][2] * m[3][1])) -
                m[0][1] * (m[1][0] * (m[2][2] * m[3][3] - m[2][3] * m[3][2]) - m[1][2] * (m[2][0] * m[3][3] - m[2][3] * m[3][0]) + m[1][3] * (m[2][0] * m[3][2] - m[2][2] * m[3][0])) +
                m[0][2] * (m[1][0] * (m[2][1] * m[3][3] - m[2][3] * m[3][1]) - m[1][1] * (m[2][0] * m[3][3] - m[2][3] * m[3][0]) + m[1][3] * (m[2][0] * m[3][1] - m[2][1] * m[3][0])) -
                m[0][3] * (m[1][0] * (m[2][1] * m[3][2] - m[2][2] * m[3][1]) - m[1][1] * (m[2][0] * m[3][2] - m[2][2] * m[3][0]) + m[1][2] * (m[2][0] * m[3][1] - m[2][1] * m[3][0]))
            );
        }
    }

    // ============================================================================
    // VI-D. PARALLEL TRANSPORT
    // ============================================================================

    /**
     * ParallelTransport: Transport vectors along curves preserving the connection
     * 
     * Solves: dV^k/dt + Γᵏᵢⱼ (dx^i/dt) V^j = 0
     */
    class ParallelTransport {
        /**
         * @param {SpacetimeMetric} metric
         * @param {number} h - Step size for numerical differentiation
         */
        constructor(metric, h = 1e-7) {
            this.metric = metric;
            this.christoffel = new ChristoffelSymbols(metric, h);
            this.dim = 4;
        }

        /**
         * Transport vector along a parameterized curve
         * 
         * @param {Function} curve - γ(t) => [x⁰, x¹, x², x³]
         * @param {number[]} V_initial - Initial vector V^μ(t₀)
         * @param {number} t_start - Starting parameter
         * @param {number} t_end - Ending parameter
         * @param {number} n_steps - Number of integration steps
         * @returns {Object[]} Array of {t, x, V} at each step
         */
        transport(curve, V_initial, t_start, t_end, n_steps = 100) {
            const dt = (t_end - t_start) / n_steps;
            const h = 1e-6;  // For computing curve tangent

            let V = [...V_initial];
            const trajectory = [{ t: t_start, x: curve(t_start), V: [...V] }];

            for (let step = 0; step < n_steps; step++) {
                const t = t_start + step * dt;
                const x = curve(t);

                // Compute curve tangent dx^i/dt
                const xPlus = curve(t + h);
                const xMinus = curve(t - h);
                const dxdt = [];
                for (let i = 0; i < this.dim; i++) {
                    dxdt[i] = (xPlus[i] - xMinus[i]) / (2 * h);
                }

                // RK4 integration of parallel transport equation
                const k1 = this._transportRHS(x, V, dxdt);

                const x2 = curve(t + dt / 2);
                const dxdt2 = [];
                const xPlus2 = curve(t + dt / 2 + h);
                const xMinus2 = curve(t + dt / 2 - h);
                for (let i = 0; i < this.dim; i++) {
                    dxdt2[i] = (xPlus2[i] - xMinus2[i]) / (2 * h);
                }
                const V2 = V.map((v, i) => v + 0.5 * dt * k1[i]);
                const k2 = this._transportRHS(x2, V2, dxdt2);

                const V3 = V.map((v, i) => v + 0.5 * dt * k2[i]);
                const k3 = this._transportRHS(x2, V3, dxdt2);

                const x4 = curve(t + dt);
                const dxdt4 = [];
                const xPlus4 = curve(t + dt + h);
                const xMinus4 = curve(t + dt - h);
                for (let i = 0; i < this.dim; i++) {
                    dxdt4[i] = (xPlus4[i] - xMinus4[i]) / (2 * h);
                }
                const V4 = V.map((v, i) => v + dt * k3[i]);
                const k4 = this._transportRHS(x4, V4, dxdt4);

                // Combine RK4 steps
                V = V.map((v, i) => v + dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));

                trajectory.push({ t: t + dt, x: curve(t + dt), V: [...V] });
            }

            return trajectory;
        }

        /**
         * RHS of parallel transport equation: dV^k/dt = -Γᵏᵢⱼ (dx^i/dt) V^j
         */
        _transportRHS(x, V, dxdt) {
            const gamma = this.christoffel.computeAt(x);
            const dVdt = [];

            for (let k = 0; k < this.dim; k++) {
                let sum = 0;
                for (let i = 0; i < this.dim; i++) {
                    for (let j = 0; j < this.dim; j++) {
                        sum += gamma[k][i][j] * dxdt[i] * V[j];
                    }
                }
                dVdt[k] = -sum;
            }
            return dVdt;
        }

        /**
         * Compute holonomy around a closed loop
         * 
         * @param {Function} loop - Closed curve γ(t), γ(0) = γ(1)
         * @param {number} n_steps - Integration steps
         * @returns {number[][]} Holonomy transformation matrix
         */
        holonomy(loop, n_steps = 200) {
            // Transport each basis vector around the loop
            const holonomyMatrix = [];

            for (let mu = 0; mu < this.dim; mu++) {
                const e_mu = Array(this.dim).fill(0);
                e_mu[mu] = 1;

                const result = this.transport(loop, e_mu, 0, 1, n_steps);
                holonomyMatrix[mu] = result[result.length - 1].V;
            }

            // Transpose to get transformation matrix
            const H = [];
            for (let i = 0; i < this.dim; i++) {
                H[i] = [];
                for (let j = 0; j < this.dim; j++) {
                    H[i][j] = holonomyMatrix[j][i];
                }
            }
            return H;
        }
    }

    /**
     * RelativisticHamiltonian: Mass-shell constraint for contact dynamics
     * 
     * H = ½ g^μν(x)(p_μ - qA_μ)(p_ν - qA_ν) + ½m² = 0
     * 
     * Couples contact kinematics with spacetime curvature
     */
    class RelativisticHamiltonian {
        /**
         * @param {SpacetimeMetric} metric
         * @param {number} mass
         * @param {Function} [gaugePotential] - A_μ(x) for EM coupling
         * @param {number} [charge] - Charge q
         */
        constructor(metric, mass = 1, gaugePotential = null, charge = 0) {
            this.metric = metric;
            this.m = mass;
            this.A = gaugePotential || (x => [0, 0, 0, 0]);
            this.q = charge;
        }

        /**
         * Evaluate Hamiltonian at spacetime point with 4-momentum
         * H = ½(g^μν(p_μ - qA_μ)(p_ν - qA_ν) + m²) = 0 on shell
         * 
         * With signature (+,-,-,-), on shell: p₀² - |p|² = m²
         * So g^μν p_μ p_ν = p₀² - |p|² should equal m² (not -m²)
         * Thus H = ½(g^μν p_μ p_ν - m²) = 0 on shell
         * 
         * @param {number[]} x - Spacetime coordinates [x⁰, x¹, x², x³]
         * @param {number[]} p - 4-momentum [p₀, p₁, p₂, p₃]
         */
        evaluate(x, p) {
            const gInv = this.metric.contravariant(x);
            const A = this.A(x);

            let H = 0;
            for (let mu = 0; mu < 4; mu++) {
                for (let nu = 0; nu < 4; nu++) {
                    const pMinusA_mu = p[mu] - this.q * A[mu];
                    const pMinusA_nu = p[nu] - this.q * A[nu];
                    H += 0.5 * gInv[mu][nu] * pMinusA_mu * pMinusA_nu;
                }
            }
            // Mass shell: g^μν p_μ p_ν = m² → H = ½(g^μν p_μ p_ν - m²) = 0
            H -= 0.5 * this.m * this.m;

            return H;
        }

        /**
         * Hamilton-Jacobi PDE residual
         * ½ g^μν(∂_μA - qA_μ)(∂_νA - qA_ν) + ½m² = 0
         * @param {number[]} x - Spacetime point
         * @param {Function} actionFunc - A(x)
         * @param {Function} gradA - ∂_μA(x)
         */
        hjResidual(x, actionFunc, gradA) {
            const dA = gradA(x);
            return this.evaluate(x, dA);
        }

        /**
         * Geodesic equations from variational principle
         * Returns Hamilton's equations for (x^μ, p_μ)
         */
        geodesicEquations(x, p) {
            const h = 1e-7;
            const dHdx = [];
            const dHdp = [];

            // ∂H/∂x^μ
            for (let mu = 0; mu < 4; mu++) {
                const xPlus = [...x]; xPlus[mu] += h;
                const xMinus = [...x]; xMinus[mu] -= h;
                dHdx.push((this.evaluate(xPlus, p) - this.evaluate(xMinus, p)) / (2 * h));
            }

            // ∂H/∂p_μ
            for (let mu = 0; mu < 4; mu++) {
                const pPlus = [...p]; pPlus[mu] += h;
                const pMinus = [...p]; pMinus[mu] -= h;
                dHdp.push((this.evaluate(x, pPlus) - this.evaluate(x, pMinus)) / (2 * h));
            }

            // ẋ^μ = ∂H/∂p_μ, ṗ_μ = -∂H/∂x^μ
            return {
                xDot: dHdp,
                pDot: dHdx.map(d => -d)
            };
        }

        /**
         * Integrate geodesic
         * @param {number[]} x0 - Initial position
         * @param {number[]} p0 - Initial momentum
         * @param {number} dtau - Proper time step
         * @param {number} steps
         */
        integrateGeodesic(x0, p0, dtau, steps) {
            let x = [...x0];
            let p = [...p0];
            const trajectory = [{ x: [...x], p: [...p], tau: 0 }];

            for (let i = 0; i < steps; i++) {
                // RK4
                const k1 = this.geodesicEquations(x, p);

                const x2 = x.map((xi, mu) => xi + 0.5 * dtau * k1.xDot[mu]);
                const p2 = p.map((pi, mu) => pi + 0.5 * dtau * k1.pDot[mu]);
                const k2 = this.geodesicEquations(x2, p2);

                const x3 = x.map((xi, mu) => xi + 0.5 * dtau * k2.xDot[mu]);
                const p3 = p.map((pi, mu) => pi + 0.5 * dtau * k2.pDot[mu]);
                const k3 = this.geodesicEquations(x3, p3);

                const x4 = x.map((xi, mu) => xi + dtau * k3.xDot[mu]);
                const p4 = p.map((pi, mu) => pi + dtau * k3.pDot[mu]);
                const k4 = this.geodesicEquations(x4, p4);

                x = x.map((xi, mu) => xi + dtau / 6 * (k1.xDot[mu] + 2 * k2.xDot[mu] + 2 * k3.xDot[mu] + k4.xDot[mu]));
                p = p.map((pi, mu) => pi + dtau / 6 * (k1.pDot[mu] + 2 * k2.pDot[mu] + 2 * k3.pDot[mu] + k4.pDot[mu]));

                trajectory.push({ x: [...x], p: [...p], tau: (i + 1) * dtau });
            }

            return trajectory;
        }
    }

    // ============================================================================
    // VII. EXTENDED THERMODYNAMIC HAMILTONIAN
    // ============================================================================

    /**
     * ThermodynamicHamiltonian: Contact Hamiltonian for extended thermodynamics
     * 
     * Combines kinetic, potential, and thermodynamic terms
     */
    class ThermodynamicHamiltonian extends ContactHamiltonian {
        /**
         * @param {ContactManifold} manifold - Grand or Holographic manifold
         * @param {Object} params - Physical parameters
         */
        constructor(manifold, params = {}) {
            // Default parameters
            const {
                mass = 1,
                potential = () => 0,
                thermalCoupling = 1,
                scalingDimension = 0
            } = params;

            // Build Hamiltonian function
            const H = coords => {
                let result = 0;

                // Kinetic term: ½m|k|² = ½m(k₁² + k₂² + k₃²)
                if ('k1' in coords) {
                    const k1 = coords.k1 || 0;
                    const k2 = coords.k2 || 0;
                    const k3 = coords.k3 || 0;
                    result += 0.5 * mass * (k1 * k1 + k2 * k2 + k3 * k3);
                }

                // Frequency-time coupling: ω
                if ('omega' in coords) {
                    result -= coords.omega;  // -ω for E = ℏω
                }

                // Dilatation term: Δ·ℓ contribution
                if ('Delta' in coords && 'ell' in coords) {
                    result += scalingDimension * coords.Delta;
                }

                // Thermal term: T·S coupling
                if ('T' in coords && 'S' in coords) {
                    result -= thermalCoupling * coords.T * coords.S;
                }

                // External potential
                result += potential(coords);

                return result;
            };

            super(manifold, H);
            this.params = params;
        }

        /**
         * Create dispersion relation Hamiltonian
         * H = ω - c|k| (massless) or H = ω - √(c²|k|² + m²c⁴) (massive)
         */
        static dispersionRelation(manifold, c = 1, mass = 0) {
            const H = coords => {
                const k1 = coords.k1 || 0;
                const k2 = coords.k2 || 0;
                const k3 = coords.k3 || 0;
                const omega = coords.omega || 0;
                const kSq = k1 * k1 + k2 * k2 + k3 * k3;

                if (mass === 0) {
                    return omega - c * sqrt(kSq);
                } else {
                    return omega - sqrt(c * c * kSq + mass * mass * c * c * c * c);
                }
            };

            return new ContactHamiltonian(manifold, H);
        }

        /**
         * Create thermodynamic equation of state Hamiltonian
         * For ideal gas: PV = NkT → H involves (T, S, ℓ)
         */
        static equationOfState(manifold, type = 'ideal') {
            const H = coords => {
                const T = coords.T || 0;
                const S = coords.S || 0;
                const ell = coords.ell || 0;
                const lambda = exp(ell);

                if (type === 'ideal') {
                    // F = U - TS, with U ∝ T and V ∝ λ³
                    return 1.5 * T - T * S + T * log(lambda * lambda * lambda);
                } else if (type === 'van_der_waals') {
                    // Van der Waals with constants a, b
                    const a = 1, b = 0.1;
                    const V = lambda * lambda * lambda;
                    return T / (V - b) - a / (V * V) - T * S;
                }
                return 0;
            };

            return new ContactHamiltonian(manifold, H);
        }
    }

    // ============================================================================
    // VIII. GAUGE EXTENSION
    // ============================================================================

    /**
     * GaugeExtendedManifold: Adds gauge coordinate pair (φ, I)
     * 
     * For Grand model + gauge: dim = 15
     * New coordinates: φ (gauge phase), I (gauge flux)
     * 
     * Extended contact form: α = dA - ... - I dφ
     */
    class GaugeExtendedManifold extends GrandContactManifold {
        constructor() {
            super();

            // Add gauge coordinates
            this.baseCoords.push('phi');
            this.momentaCoords.push('I');
            this.n = this.baseCoords.length;  // Now 7
            this.dim = 2 * this.n + 1;  // Now 15

            this.physicalInterpretations.phi = 'gauge phase';
            this.physicalInterpretations.I = 'gauge flux/current';
        }

        contactFormSymbolic() {
            return 'dA - k₁dq¹ - k₂dq² - k₃dq³ - ω dt - Δ dℓ - T dS - I dφ';
        }

        verifyContactCondition(pt) {
            return 5040; // 7!
        }

        toString() {
            return 'Gauge-Extended Contact Manifold M₁₅ = J¹(Q₇): dim=15';
        }
    }

    // ============================================================================
    // IX. UTILITY CLASSES
    // ============================================================================

    /**
     * DifferentialForm: Symbolic differential form representation
     */
    class DifferentialForm {
        constructor(degree, terms = []) {
            this.degree = degree;
            this.terms = terms; // [{coeff, basis: ['dx1', 'dx2', ...]}]
        }

        static oneForm(coord) {
            return new DifferentialForm(1, [{ coeff: 1, basis: [`d${coord}`] }]);
        }

        wedge(other) {
            const newDegree = this.degree + other.degree;
            const newTerms = [];

            for (const t1 of this.terms) {
                for (const t2 of other.terms) {
                    // Check for repeated basis elements
                    const combined = [...t1.basis, ...t2.basis];
                    const unique = new Set(combined);
                    if (unique.size === combined.length) {
                        // Compute sign from permutation
                        let sign = 1;
                        // Bubble sort to count inversions
                        const arr = [...combined];
                        for (let i = 0; i < arr.length; i++) {
                            for (let j = i + 1; j < arr.length; j++) {
                                if (arr[i] > arr[j]) {
                                    [arr[i], arr[j]] = [arr[j], arr[i]];
                                    sign *= -1;
                                }
                            }
                        }
                        newTerms.push({
                            coeff: sign * t1.coeff * t2.coeff,
                            basis: arr
                        });
                    }
                }
            }

            return new DifferentialForm(newDegree, newTerms);
        }

        toString() {
            if (this.terms.length === 0) return '0';
            return this.terms.map(t => {
                const coeffStr = t.coeff === 1 ? '' : (t.coeff === -1 ? '-' : `${t.coeff}`);
                return coeffStr + t.basis.join('∧');
            }).join(' + ');
        }
    }

    /**
     * Summary table generator
     */
    function summaryTable() {
        return {
            columns: ['Theory Version', 'Base Dim', 'Total Contact Dim', 'Independent Coordinates', 'Dependent Fields'],
            rows: [
                ['Grand (Standard)', '6', '13', 'q¹,q²,q³,t,ℓ,S + A', 'None'],
                ['Holographic', '3', '7', 't,ℓ,S + A', 'q^i(t,ℓ,S)'],
                ['Gauge-Extended (Grand+φ)', '7', '15', 'add (φ,I)', 'None']
            ]
        };
    }

    // ============================================================================
    // MODULE IMPORTS
    // ============================================================================

    // Mesh / FTGC (Discrete Geometric Calculus on Meshes)
    let MeshModule, MeshFTGCModule, MeshSolversModule, EntropicGravityModule, MeshEntropicFlowModule;
    // Riemannian Geometry modules
    let RiemannianGAModule, GeodesicGAModule, GeometricCalculusModule;
    let RiemannianDiscreteModule, RiemannianSpacetimeModule;
    // Quantum / Information modules
    let PilotWaveModule, InformationGeometryModule;
    // New modular subpackages
    let AlgebraModule, GeometryModule, CalculusModule, ContactModule, PhysicsModule;

    // These are all first-party modules that must resolve. Requiring them
    // directly (rather than swallowing failures in an empty catch) means a
    // broken path surfaces immediately instead of silently producing a
    // half-populated public API.
    //
    // Gate on `module.exports` (CommonJS/Node) rather than `typeof require`,
    // which is also truthy under AMD loaders (e.g. RequireJS) whose `require`
    // is not a synchronous CommonJS loader and would throw here before the UMD
    // `define.amd` branch at the bottom of this file runs.
    if (typeof module !== 'undefined' && module.exports) {
        // Mesh modules (legacy paths for compatibility)
        MeshModule = require('./mesh.js');
        MeshFTGCModule = require('./mesh-ftgc.js');
        MeshSolversModule = require('./mesh-solvers.js');
        EntropicGravityModule = require('./entropic-gravity.js');
        MeshEntropicFlowModule = require('./mesh-entropic-flow.js');

        // Riemannian Geometry modules (legacy paths)
        RiemannianGAModule = require('./riemannian-ga.js');
        GeodesicGAModule = require('./geodesic-ga.js');
        GeometricCalculusModule = require('./geometric-calculus.js');
        RiemannianDiscreteModule = require('./riemannian-discrete.js');
        RiemannianSpacetimeModule = require('./riemannian-spacetime.js');

        // Quantum / Information modules (legacy paths)
        PilotWaveModule = require('./pilot-wave.js');
        InformationGeometryModule = require('./information-geometry.js');

        // New modular subpackages
        AlgebraModule = require('./algebra');
        GeometryModule = require('./geometry');
        CalculusModule = require('./calculus');
        ContactModule = require('./contact');
        PhysicsModule = require('./physics');
    }

    // ============================================================================
    // EXPORTS
    // ============================================================================

    const ExtendedThermodynamics = {
        // Core classes
        ContactManifold,
        ContactPoint,

        // Specific manifolds
        GrandContactManifold,
        HolographicContactManifold,
        GaugeExtendedManifold,

        // Dynamics
        ContactHamiltonian,
        LegendrianSubmanifold,
        ThermodynamicHamiltonian,

        // Gravitational extension
        SpacetimeMetric,
        RelativisticHamiltonian,

        // Differential geometry (from GA skill)
        ChristoffelSymbols,
        CovariantDerivative,
        ParallelTransport,

        // Mesh / FTGC (Discrete Geometric Calculus on Triangle Meshes)
        ...(MeshModule || {}),
        ...(MeshFTGCModule || {}),
        ...(MeshSolversModule || {}),

        // Entropic Gravity (Bianconi Framework)
        // Two-metric system (g, G), relative entropy action, modified Einstein equations
        ...(EntropicGravityModule || {}),

        // Entropic Gradient Flow on Meshes (Particles + Dynamic Spacetimes)
        ...(MeshEntropicFlowModule || {}),

        // Riemannian Geometry via Geometric Algebra (Coordinate-free)
        // Connection bivector ω, Curvature 2-form Ω, Sphere2D, Torus2D, HyperbolicPlane
        ...(RiemannianGAModule || {}),

        // Geodesic Solver (Coordinate-free)
        // GAGeodesicSolver, GAParallelTransport, HolonomyComputer
        ...(GeodesicGAModule || {}),

        // Geometric Calculus on Regular Grids
        // ScalarField, VectorField, SplitDifferentialOperator, LeapfrogIntegrator
        ...(GeometricCalculusModule || {}),

        // Discrete Riemannian Geometry on Meshes
        // Bivector3D, MeshCurvature2Form, MeshConnectionBivector
        ...(RiemannianDiscreteModule || {}),

        // Spacetime GA for General Relativity
        // SpacetimeManifoldGA with Cl(1,3) algebra
        ...(RiemannianSpacetimeModule || {}),

        // Pilot-Wave Theory with Valentini Regularization
        // SmearingKernel, WaveFunction, PilotWaveSystem, QuantumEnsemble
        ...(PilotWaveModule || {}),

        // Information Geometry
        // ProbabilityManifold based on Baez's Information Geometry
        ...(InformationGeometryModule || {}),

        // Utilities
        DifferentialForm,
        summaryTable,

        // Constants
        EPSILON,

        // ============================================================================
        // MODULAR NAMESPACED EXPORTS (New structure)
        // ============================================================================
        // These provide access to the new modular subpackages:
        //   algebra: Multivector, Algebra, NumberSystems (Complex, Dual, Hyperbolic)
        //   geometry: Riemannian GA (connection bivector, curvature), geodesics
        //   calculus: Discrete geometric calculus (grid + mesh operators, solvers)
        //   contact: Contact manifolds, Hamiltonians, Legendrians (extracted from index)
        //   physics: Pilot-wave, entropic gravity, spacetime
        algebra: AlgebraModule || {},
        geometry: GeometryModule || {},
        calculus: CalculusModule || {},
        contact: ContactModule || {},
        physics: PhysicsModule || {},

        // Factory functions
        grandManifold: () => new GrandContactManifold(),
        holographicManifold: () => new HolographicContactManifold(),
        gaugeExtended: () => new GaugeExtendedManifold()
    };

    // Export for different module systems
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = ExtendedThermodynamics;
    } else if (typeof define === 'function' && define.amd) {
        define([], () => ExtendedThermodynamics);
    } else {
        global.ExtendedThermodynamics = ExtendedThermodynamics;
    }

})(typeof window !== 'undefined' ? window : global);
