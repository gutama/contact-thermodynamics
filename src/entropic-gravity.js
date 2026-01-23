/**
 * Entropic Gravity Module (Bianconi Framework)
 *
 * Implements "Gravity from Entropy" (Bianconi, 2025) as a dynamical closure
 * for the GMET contact geometry framework.
 *
 * Key Concepts:
 * - Two-Metric System: spacetime metric g and matter-induced metric G
 * - Quantum Relative Entropy Action: S(G||g) = Tr[G(ln G - ln g)]
 * - Modified Einstein Equations with emergent cosmological constant Λ_G
 *
 * Integration with GMET:
 * - GMET provides the kinematic scaffold (13D/7D contact manifold)
 * - Bianconi provides the gravitational dynamics (theory selection rule)
 * - Together: GMET is the container, Bianconi is the preferred fluid
 *
 * @module entropic-gravity
 * @license MIT
 */

(function (global) {
    'use strict';

    // [NEW] GA Integration (Node.js support)
    let SpacetimeManifoldGA;
    if (typeof require === 'function') {
        try {
            const mod = require('./riemannian-spacetime.js');
            SpacetimeManifoldGA = mod.SpacetimeManifoldGA;
        } catch (e) {
            // Ignore in browser or if module missing
        }
    }

    const EPSILON = 1e-12;
    const { abs, sqrt, sin, cos, exp, log, PI } = Math;

    // ============================================================================
    // UTILITY FUNCTIONS
    // ============================================================================

    function isZero(x) { return abs(x) < EPSILON; }

    /**
     * 4x4 matrix operations for spacetime tensors
     */
    function mat4Add(A, B) {
        return A.map((row, i) => row.map((val, j) => val + B[i][j]));
    }

    function mat4Sub(A, B) {
        return A.map((row, i) => row.map((val, j) => val - B[i][j]));
    }

    function mat4Scale(A, s) {
        return A.map(row => row.map(val => val * s));
    }

    function mat4Mul(A, B) {
        const C = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]];
        for (let i = 0; i < 4; i++) {
            for (let j = 0; j < 4; j++) {
                for (let k = 0; k < 4; k++) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return C;
    }

    function mat4Trace(A) {
        return A[0][0] + A[1][1] + A[2][2] + A[3][3];
    }

    function mat4Det(m) {
        return (
            m[0][0] * (m[1][1] * (m[2][2] * m[3][3] - m[2][3] * m[3][2]) - m[1][2] * (m[2][1] * m[3][3] - m[2][3] * m[3][1]) + m[1][3] * (m[2][1] * m[3][2] - m[2][2] * m[3][1])) -
            m[0][1] * (m[1][0] * (m[2][2] * m[3][3] - m[2][3] * m[3][2]) - m[1][2] * (m[2][0] * m[3][3] - m[2][3] * m[3][0]) + m[1][3] * (m[2][0] * m[3][2] - m[2][2] * m[3][0])) +
            m[0][2] * (m[1][0] * (m[2][1] * m[3][3] - m[2][3] * m[3][1]) - m[1][1] * (m[2][0] * m[3][3] - m[2][3] * m[3][0]) + m[1][3] * (m[2][0] * m[3][1] - m[2][1] * m[3][0])) -
            m[0][3] * (m[1][0] * (m[2][1] * m[3][2] - m[2][2] * m[3][1]) - m[1][1] * (m[2][0] * m[3][2] - m[2][2] * m[3][0]) + m[1][2] * (m[2][0] * m[3][1] - m[2][1] * m[3][0]))
        );
    }

    function mat4Inv(m) {
        const n = 4;
        const aug = m.map((row, i) => [...row, ...Array(n).fill(0).map((_, j) => i === j ? 1 : 0)]);

        for (let col = 0; col < n; col++) {
            let maxRow = col;
            for (let row = col + 1; row < n; row++) {
                if (abs(aug[row][col]) > abs(aug[maxRow][col])) maxRow = row;
            }
            [aug[col], aug[maxRow]] = [aug[maxRow], aug[col]];

            if (isZero(aug[col][col])) {
                throw new Error('Singular metric in mat4Inv');
            }

            const scale = aug[col][col];
            for (let j = 0; j < 2 * n; j++) aug[col][j] /= scale;

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
     * Matrix logarithm approximation for positive-definite symmetric matrices
     * Uses series expansion: ln(I + A) ≈ A - A²/2 + A³/3 - ...
     * For metrics, we compute ln(g) via eigendecomposition approximation
     */
    function mat4Log(A, maxIter = 20) {
        // For symmetric positive-definite matrices near identity
        // ln(A) = ln(I + (A-I)) ≈ (A-I) - (A-I)²/2 + ...
        const I = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
        const AmI = mat4Sub(A, I);

        let result = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]];
        let power = AmI;

        for (let k = 1; k <= maxIter; k++) {
            const sign = (k % 2 === 1) ? 1 : -1;
            const term = mat4Scale(power, sign / k);
            result = mat4Add(result, term);
            power = mat4Mul(power, AmI);

            // Check convergence
            if (mat4Trace(mat4Mul(term, term)) < EPSILON * EPSILON) break;
        }

        return result;
    }

    /**
     * Matrix exponential via Padé approximation
     */
    function mat4Exp(A, order = 6) {
        const I = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
        let result = I;
        let term = I;

        for (let k = 1; k <= order; k++) {
            term = mat4Scale(mat4Mul(term, A), 1 / k);
            result = mat4Add(result, term);
        }

        return result;
    }

    /**
     * Identity matrix 4x4
     */
    function mat4Identity() {
        return [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
    }

    /**
     * Zero matrix 4x4
     */
    function mat4Zero() {
        return [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]];
    }

    // ============================================================================
    // I. MATTER-INDUCED METRIC
    // ============================================================================

    /**
     * MatterInducedMetric: The metric G_μν induced by topological matter fields
     *
     * G_μν = R_μν + ∂_μφ ∂_νφ + F_μρ F_ν^ρ + H_μρσ H_ν^ρσ
     *
     * Where:
     *   φ: scalar field (0-form)
     *   A_μ: gauge potential (1-form), F = dA
     *   B_μν: 2-form potential, H = dB
     *   R_μν: Ricci tensor from spacetime metric g
     */
    class MatterInducedMetric {
        /**
         * @param {Object} options
         * @param {Function} options.scalarField - φ(x) → number
         * @param {Function} options.vectorPotential - A(x) → [A_0, A_1, A_2, A_3]
         * @param {Function} options.twoFormPotential - B(x) → 4x4 antisymmetric matrix B_μν
         * @param {Function} options.ricciTensor - R(x) → 4x4 Ricci tensor R_μν
         */
        constructor(options = {}) {
            this.phi = options.scalarField || (x => 0);
            this.A = options.vectorPotential || (x => [0, 0, 0, 0]);
            this.B = options.twoFormPotential || (x => mat4Zero());
            this.R = options.ricciTensor || (x => mat4Zero());
            this.h = options.differentiationStep || 1e-7;
        }

        /**
         * Compute field strength F_μν = ∂_μ A_ν - ∂_ν A_μ
         */
        fieldStrength(x) {
            const F = mat4Zero();
            const h = this.h;

            for (let mu = 0; mu < 4; mu++) {
                for (let nu = 0; nu < 4; nu++) {
                    if (mu !== nu) {
                        // ∂_μ A_ν
                        const xPlus = [...x]; xPlus[mu] += h;
                        const xMinus = [...x]; xMinus[mu] -= h;
                        const dmuAnu = (this.A(xPlus)[nu] - this.A(xMinus)[nu]) / (2 * h);

                        // ∂_ν A_μ
                        const xPlus2 = [...x]; xPlus2[nu] += h;
                        const xMinus2 = [...x]; xMinus2[nu] -= h;
                        const dnuAmu = (this.A(xPlus2)[mu] - this.A(xMinus2)[mu]) / (2 * h);

                        F[mu][nu] = dmuAnu - dnuAmu;
                    }
                }
            }
            return F;
        }

        /**
         * Compute 3-form field strength H_μνρ = ∂_μ B_νρ + cyclic
         * Returns contracted form H_μρσ H_ν^ρσ as 4x4 matrix
         */
        threeFormContraction(x, gInv) {
            const h = this.h;
            const result = mat4Zero();

            // Simplified: compute H·H contraction
            // For full implementation, need 3-form components
            // Here we use approximation based on B field derivatives

            const B = this.B(x);

            for (let mu = 0; mu < 4; mu++) {
                for (let nu = 0; nu < 4; nu++) {
                    let sum = 0;
                    for (let rho = 0; rho < 4; rho++) {
                        for (let sigma = 0; sigma < 4; sigma++) {
                            // Approximate H_μρσ ≈ ∂_μ B_ρσ (simplified)
                            const xPlus = [...x]; xPlus[mu] += h;
                            const xMinus = [...x]; xMinus[mu] -= h;
                            const dB_mu = (this.B(xPlus)[rho][sigma] - this.B(xMinus)[rho][sigma]) / (2 * h);

                            const xPlus2 = [...x]; xPlus2[nu] += h;
                            const xMinus2 = [...x]; xMinus2[nu] -= h;
                            const dB_nu = (this.B(xPlus2)[rho][sigma] - this.B(xMinus2)[rho][sigma]) / (2 * h);

                            // Contract with inverse metric
                            for (let alpha = 0; alpha < 4; alpha++) {
                                for (let beta = 0; beta < 4; beta++) {
                                    sum += dB_mu * dB_nu * gInv[rho][alpha] * gInv[sigma][beta];
                                }
                            }
                        }
                    }
                    result[mu][nu] = sum;
                }
            }

            return result;
        }

        /**
         * Compute gradient ∂_μφ
         */
        scalarGradient(x) {
            const grad = [0, 0, 0, 0];
            const h = this.h;

            for (let mu = 0; mu < 4; mu++) {
                const xPlus = [...x]; xPlus[mu] += h;
                const xMinus = [...x]; xMinus[mu] -= h;
                grad[mu] = (this.phi(xPlus) - this.phi(xMinus)) / (2 * h);
            }
            return grad;
        }

        /**
         * Compute matter-induced metric G_μν at point x
         *
         * G_μν = R_μν + ∂_μφ ∂_νφ + F_μρ F_ν^ρ + H_μρσ H_ν^ρσ
         *
         * @param {number[]} x - Spacetime coordinates
         * @param {number[][]} gInv - Inverse spacetime metric g^μν for raising indices
         */
        covariant(x, gInv = null) {
            gInv = gInv || mat4Identity();

            // Ricci contribution
            const G = this.R(x).map(row => [...row]);

            // Scalar field contribution: ∂_μφ ∂_νφ
            const dphi = this.scalarGradient(x);
            for (let mu = 0; mu < 4; mu++) {
                for (let nu = 0; nu < 4; nu++) {
                    G[mu][nu] += dphi[mu] * dphi[nu];
                }
            }

            // Gauge field contribution: F_μρ F_ν^ρ = F_μρ g^ρσ F_νσ
            const F = this.fieldStrength(x);
            for (let mu = 0; mu < 4; mu++) {
                for (let nu = 0; nu < 4; nu++) {
                    let sum = 0;
                    for (let rho = 0; rho < 4; rho++) {
                        for (let sigma = 0; sigma < 4; sigma++) {
                            sum += F[mu][rho] * gInv[rho][sigma] * F[nu][sigma];
                        }
                    }
                    G[mu][nu] += sum;
                }
            }

            // 2-form contribution (simplified)
            const HH = this.threeFormContraction(x, gInv);
            for (let mu = 0; mu < 4; mu++) {
                for (let nu = 0; nu < 4; nu++) {
                    G[mu][nu] += HH[mu][nu];
                }
            }

            return G;
        }

        /**
         * Contravariant matter metric G^μν
         */
        contravariant(x, gInv = null) {
            const G = this.covariant(x, gInv);
            return mat4Inv(G);
        }
    }

    // ============================================================================
    // II. TWO-METRIC SYSTEM
    // ============================================================================

    /**
     * TwoMetricSystem: The (g, G) pair central to Bianconi's entropic gravity
     *
     * - g_μν: Spacetime metric (the "quantum operator")
     * - G_μν: Matter-induced metric (the "classical state")
     *
     * The interplay between g and G generates gravitational dynamics
     */
    class TwoMetricSystem {
        /**
         * @param {Object} spacetimeMetric - SpacetimeMetric instance (g)
         * @param {MatterInducedMetric} matterMetric - MatterInducedMetric instance (G)
         */
        constructor(spacetimeMetric, matterMetric) {
            if (!spacetimeMetric || typeof spacetimeMetric.covariant !== 'function') {
                throw new Error('spacetimeMetric must have covariant() method');
            }
            this.g = spacetimeMetric;
            this.G = matterMetric || new MatterInducedMetric();
        }

        /**
         * Get spacetime metric g_μν at point x
         */
        spacetimeMetric(x) {
            return this.g.covariant(x);
        }

        /**
         * Get matter-induced metric G_μν at point x
         */
        matterMetric(x) {
            const gInv = this.g.contravariant ? this.g.contravariant(x) : mat4Inv(this.g.covariant(x));
            return this.G.covariant(x, gInv);
        }

        /**
         * Metric difference: δg_μν = G_μν - g_μν
         */
        metricDifference(x) {
            const g = this.spacetimeMetric(x);
            const G = this.matterMetric(x);
            return mat4Sub(G, g);
        }

        /**
         * Metric ratio: G · g⁻¹ (for eigenvalue analysis)
         */
        metricRatio(x) {
            const gInv = this.g.contravariant ? this.g.contravariant(x) : mat4Inv(this.g.covariant(x));
            const G = this.matterMetric(x);
            return mat4Mul(G, gInv);
        }

        /**
         * Volume element ratio: √(det G / det g)
         */
        volumeRatio(x) {
            const g = this.spacetimeMetric(x);
            const G = this.matterMetric(x);
            const detG = mat4Det(G);
            const detg = mat4Det(g);

            if (abs(detg) < EPSILON) return 0;
            return sqrt(abs(detG / detg));
        }
    }

    // ============================================================================
    // III. RELATIVE ENTROPY ACTION
    // ============================================================================

    /**
     * RelativeEntropyAction: The quantum relative entropy functional S(G||g)
     *
     * S(G||g) = ∫ d⁴x √|g| Tr[G(ln G - ln g)]
     *
     * This action serves as the generating function for gravitational dynamics.
     * Minimizing S(G||g) yields the modified Einstein equations.
     */
    class RelativeEntropyAction {
        /**
         * @param {TwoMetricSystem} twoMetricSystem
         * @param {Object} options
         */
        constructor(twoMetricSystem, options = {}) {
            this.system = twoMetricSystem;
            this.h = options.integrationStep || 0.1;
            this.logIterations = options.logIterations || 20;
        }

        /**
         * Compute local entropy density: Tr[G(ln G - ln g)]
         *
         * Uses a regularized first-order approximation for numerical stability:
         * S ≈ Tr[(G - g) · g⁻¹] for small deviations
         *
         * @param {number[]} x - Spacetime point
         * @returns {number} Local entropy density
         */
        localDensity(x) {
            const g = this.system.spacetimeMetric(x);
            const G = this.system.matterMetric(x);

            const trG = mat4Trace(G);
            const trg = mat4Trace(g);

            // If G is negligible, return 0
            if (abs(trG) < EPSILON * 10) return 0;

            // For numerical stability, use first-order approximation:
            // S(G||g) ≈ Tr[(G - g) · g⁻¹] when metrics differ
            // This avoids computing log of potentially ill-conditioned matrices

            try {
                const gInv = mat4Inv(g);
                const diff = mat4Sub(G, g);
                const product = mat4Mul(diff, gInv);
                const firstOrder = mat4Trace(product);

                // Add small regularization for stability
                const normG = sqrt(abs(mat4Trace(mat4Mul(G, G))));
                const normg = sqrt(abs(mat4Trace(mat4Mul(g, g))));

                const regularization = (normg > EPSILON)
                    ? 0.001 * (normG / (normg + EPSILON))
                    : 0;

                // Return bounded value
                const result = firstOrder + regularization;
                return isFinite(result) ? result : 0;
            } catch (e) {
                return 0;
            }
        }

        /**
         * Compute the relative entropy integrand with volume measure
         *
         * s(x) = √|g| · Tr[G(ln G - ln g)]
         */
        integrand(x) {
            const g = this.system.spacetimeMetric(x);
            const sqrtDetG = sqrt(abs(mat4Det(g)));
            return sqrtDetG * this.localDensity(x);
        }

        /**
         * Compute total action over a spacetime region
         *
         * S = ∫ d⁴x √|g| Tr[G(ln G - ln g)]
         *
         * @param {number[][]} bounds - [[t_min,t_max], [x_min,x_max], [y_min,y_max], [z_min,z_max]]
         * @param {number} nPoints - Points per dimension
         */
        totalAction(bounds, nPoints = 10) {
            const h = [];
            for (let i = 0; i < 4; i++) {
                h[i] = (bounds[i][1] - bounds[i][0]) / nPoints;
            }

            let S = 0;
            for (let i0 = 0; i0 < nPoints; i0++) {
                const t = bounds[0][0] + (i0 + 0.5) * h[0];
                for (let i1 = 0; i1 < nPoints; i1++) {
                    const x = bounds[1][0] + (i1 + 0.5) * h[1];
                    for (let i2 = 0; i2 < nPoints; i2++) {
                        const y = bounds[2][0] + (i2 + 0.5) * h[2];
                        for (let i3 = 0; i3 < nPoints; i3++) {
                            const z = bounds[3][0] + (i3 + 0.5) * h[3];
                            S += this.integrand([t, x, y, z]);
                        }
                    }
                }
            }

            return S * h[0] * h[1] * h[2] * h[3];
        }

        /**
         * Variation of action with respect to spacetime metric: δS/δg^μν
         *
         * This gives the "entropic stress-energy tensor"
         */
        variationWrtMetric(x, h = 1e-6) {
            const variation = mat4Zero();

            for (let mu = 0; mu < 4; mu++) {
                for (let nu = mu; nu < 4; nu++) {
                    // Perturb g^μν
                    const originalG = this.system.g;

                    // Create perturbed metric
                    const perturbedMetricPlus = {
                        covariant: (pt) => {
                            const g = originalG.covariant(pt);
                            const gPerturbed = g.map(row => [...row]);
                            // Perturb the covariant metric (approximate)
                            gPerturbed[mu][nu] += h;
                            if (mu !== nu) gPerturbed[nu][mu] += h;
                            return gPerturbed;
                        },
                        contravariant: (pt) => mat4Inv(perturbedMetricPlus.covariant(pt))
                    };

                    const perturbedMetricMinus = {
                        covariant: (pt) => {
                            const g = originalG.covariant(pt);
                            const gPerturbed = g.map(row => [...row]);
                            gPerturbed[mu][nu] -= h;
                            if (mu !== nu) gPerturbed[nu][mu] -= h;
                            return gPerturbed;
                        },
                        contravariant: (pt) => mat4Inv(perturbedMetricMinus.covariant(pt))
                    };

                    // Compute numerical derivative
                    const systemPlus = new TwoMetricSystem(perturbedMetricPlus, this.system.G);
                    const systemMinus = new TwoMetricSystem(perturbedMetricMinus, this.system.G);
                    const actionPlus = new RelativeEntropyAction(systemPlus);
                    const actionMinus = new RelativeEntropyAction(systemMinus);

                    const dS = (actionPlus.localDensity(x) - actionMinus.localDensity(x)) / (2 * h);
                    variation[mu][nu] = dS;
                    if (mu !== nu) variation[nu][mu] = dS;
                }
            }

            return variation;
        }
    }

    // ============================================================================
    // IV. EMERGENT COSMOLOGICAL CONSTANT
    // ============================================================================

    /**
     * EmergentCosmologicalConstant: Computes Λ_G from the G-field
     *
     * In Bianconi's framework, Λ emerges from the auxiliary G-field
     * required to enforce the constraint between g and G.
     *
     * Λ_G = (1/4) Tr[G·g⁻¹ - 4·I]
     */
    class EmergentCosmologicalConstant {
        /**
         * @param {TwoMetricSystem} twoMetricSystem
         */
        constructor(twoMetricSystem) {
            this.system = twoMetricSystem;
        }

        /**
         * Compute local emergent Λ at point x
         *
         * Λ_G = (1/4) Tr[G·g⁻¹ - 4·I]
         */
        localValue(x) {
            const ratio = this.system.metricRatio(x);
            const trace = mat4Trace(ratio);
            return (trace - 4) / 4;
        }

        /**
         * Compute average Λ over a region
         */
        averageValue(bounds, nPoints = 10) {
            const h = [];
            for (let i = 0; i < 4; i++) {
                h[i] = (bounds[i][1] - bounds[i][0]) / nPoints;
            }

            let sum = 0;
            let count = 0;

            for (let i0 = 0; i0 < nPoints; i0++) {
                const t = bounds[0][0] + (i0 + 0.5) * h[0];
                for (let i1 = 0; i1 < nPoints; i1++) {
                    const x = bounds[1][0] + (i1 + 0.5) * h[1];
                    for (let i2 = 0; i2 < nPoints; i2++) {
                        const y = bounds[2][0] + (i2 + 0.5) * h[2];
                        for (let i3 = 0; i3 < nPoints; i3++) {
                            const z = bounds[3][0] + (i3 + 0.5) * h[3];
                            sum += this.localValue([t, x, y, z]);
                            count++;
                        }
                    }
                }
            }

            return sum / count;
        }

        /**
         * Check if Λ is positive (de Sitter) or negative (anti-de Sitter)
         */
        classify(x) {
            const lambda = this.localValue(x);
            if (lambda > EPSILON) return 'de Sitter (Λ > 0)';
            if (lambda < -EPSILON) return 'anti-de Sitter (Λ < 0)';
            return 'Minkowski (Λ ≈ 0)';
        }
    }

    // ============================================================================
    // V. MODIFIED EINSTEIN SOLVER
    // ============================================================================

    /**
     * ModifiedEinsteinSolver: Solves the modified Einstein equations
     *
     * G_μν + Λ_G g_μν = 8πG T_μν + T^(ent)_μν
     *
     * Where T^(ent)_μν is the entropic stress-energy tensor from δS/δg
     */
    class ModifiedEinsteinSolver {
        /**
         * @param {TwoMetricSystem} twoMetricSystem
         * @param {Object} options
         */
        constructor(twoMetricSystem, options = {}) {
            this.system = twoMetricSystem;
            this.action = new RelativeEntropyAction(twoMetricSystem);
            this.lambda = new EmergentCosmologicalConstant(twoMetricSystem);
            this.G_newton = options.newtonConstant || 1;  // Geometric units
            this.h = options.differentiationStep || 1e-7;
        }

        /**
         * Compute Einstein tensor G_μν = R_μν - (1/2)R g_μν
         *
         * @param {number[]} x - Spacetime point
         * @param {ChristoffelSymbols} christoffel - Christoffel symbols instance
         */
        einsteinTensor(x, christoffel) {
            const Ric = christoffel.ricciAt(x);
            const R = christoffel.ricciScalarAt(x);
            const g = this.system.spacetimeMetric(x);

            const G_tensor = mat4Zero();
            for (let mu = 0; mu < 4; mu++) {
                for (let nu = 0; nu < 4; nu++) {
                    G_tensor[mu][nu] = Ric[mu][nu] - 0.5 * R * g[mu][nu];
                }
            }
            return G_tensor;
        }

        /**
         * Compute entropic stress-energy tensor T^(ent)_μν
         */
        entropicStressEnergy(x) {
            return this.action.variationWrtMetric(x);
        }

        /**
         * Compute effective stress-energy: 8πG T_μν + T^(ent)_μν - Λ_G g_μν
         *
         * @param {number[]} x - Spacetime point
         * @param {Function} T_matter - Matter stress-energy T_μν(x)
         */
        effectiveStressEnergy(x, T_matter = null) {
            const T_mat = T_matter ? T_matter(x) : mat4Zero();
            const T_ent = this.entropicStressEnergy(x);
            const g = this.system.spacetimeMetric(x);
            const Lambda = this.lambda.localValue(x);

            const T_eff = mat4Zero();
            for (let mu = 0; mu < 4; mu++) {
                for (let nu = 0; nu < 4; nu++) {
                    T_eff[mu][nu] = 8 * PI * this.G_newton * T_mat[mu][nu]
                        + T_ent[mu][nu]
                        - Lambda * g[mu][nu];
                }
            }
            return T_eff;
        }

        /**
         * Check field equations residual: G_μν - T^eff_μν ≈ 0
         */
        residual(x, christoffel, T_matter = null) {
            const G_tensor = this.einsteinTensor(x, christoffel);
            const T_eff = this.effectiveStressEnergy(x, T_matter);
            return mat4Sub(G_tensor, T_eff);
        }

        /**
         * Compute residual norm (for convergence checking)
         */
        residualNorm(x, christoffel, T_matter = null) {
            const R = this.residual(x, christoffel, T_matter);
            let norm = 0;
            for (let mu = 0; mu < 4; mu++) {
                for (let nu = 0; nu < 4; nu++) {
                    norm += R[mu][nu] * R[mu][nu];
                }
            }
            return sqrt(norm);
        }
    }

    // ============================================================================
    // VI. ENTROPIC GRAVITY HAMILTONIAN
    // ============================================================================

    /**
     * EntropicGravityHamiltonian: Contact Hamiltonian for GMET + Bianconi
     *
     * H = H_geo + H_entropy
     *
     * Where:
     *   H_geo = ½ g^μν (p_μ - qA_μ)(p_ν - qA_ν) - ½m² (mass-shell)
     *   H_entropy = α · S(G||g) (entropic correction)
     *
     * This combines GMET kinematics with Bianconi dynamics
     */
    class EntropicGravityHamiltonian {
        /**
         * @param {Object} manifold - Contact manifold (GrandContactManifold)
         * @param {TwoMetricSystem} twoMetricSystem
         * @param {Object} options
         */
        constructor(manifold, twoMetricSystem, options = {}) {
            this.manifold = manifold;
            this.system = twoMetricSystem;
            this.action = new RelativeEntropyAction(twoMetricSystem);

            this.mass = options.mass || 1;
            this.charge = options.charge || 0;
            this.gaugePotential = options.gaugePotential || (x => [0, 0, 0, 0]);
            this.entropicCoupling = options.entropicCoupling || 0.1;
            this.h = options.differentiationStep || 1e-7;
            this.useGA = options.useGA || false;  // [NEW] Toggle for GA mode

            // [NEW] GA Manifold instance (lazy init)
            this.gaManifold = null;
            if (this.useGA) {
                // Wrap the metric function for GA
                // system.g.covariant must be a function of x
                const metricFunc = (x) => this.system.g.covariant(x);
                this.gaManifold = new SpacetimeManifoldGA(metricFunc);
            }
        }

        /**
         * Extract spacetime coordinates from contact point
         */
        _extractSpacetime(coords) {
            // Map GMET coords to spacetime: (t, q1, q2, q3) → (x^0, x^1, x^2, x^3)
            return [
                coords.t || 0,
                coords.q1 || 0,
                coords.q2 || 0,
                coords.q3 || 0
            ];
        }

        /**
         * Extract 4-momentum from contact point
         */
        _extractMomentum(coords) {
            // Map GMET momenta: (omega, k1, k2, k3) → (p_0, p_1, p_2, p_3)
            return [
                coords.omega || 0,
                coords.k1 || 0,
                coords.k2 || 0,
                coords.k3 || 0
            ];
        }

        /**
         * Geometric Hamiltonian (mass-shell constraint)
         * H_geo = ½ g^μν (p_μ - qA_μ)(p_ν - qA_ν) - ½m²
         */
        geometricPart(coords) {
            const x = this._extractSpacetime(coords);
            const p = this._extractMomentum(coords);
            const A = this.gaugePotential(x);
            const gInv = this.system.g.contravariant
                ? this.system.g.contravariant(x)
                : mat4Inv(this.system.g.covariant(x));

            let H = 0;
            for (let mu = 0; mu < 4; mu++) {
                for (let nu = 0; nu < 4; nu++) {
                    const pMinusA_mu = p[mu] - this.charge * A[mu];
                    const pMinusA_nu = p[nu] - this.charge * A[nu];
                    H += 0.5 * gInv[mu][nu] * pMinusA_mu * pMinusA_nu;
                }
            }
            H -= 0.5 * this.mass * this.mass;

            return H;
        }

        /**
         * Entropic Hamiltonian contribution
         * H_entropy = α · S_local(G||g)
         */
        entropicPart(coords) {
            const x = this._extractSpacetime(coords);
            return this.entropicCoupling * this.action.localDensity(x);
        }

        /**
         * Total Hamiltonian
         * H = H_geo + H_entropy
         */
        evaluate(coords) {
            if (this.useGA) {
                // In GA mode, we use the connection bivector field
                return this.geometricPart(coords) + this.entropicPartGA(coords);
            }
            return this.geometricPart(coords) + this.entropicPart(coords);
        }

        /**
         * [NEW] Entropic Hamiltonian Part via Geometric Algebra
         * H_entropy = α · < Ω · S >  (coupling curvature 2-form to entropy flux?)
         * 
         * Simplified for this stage:
         * Uses the magnitude of the Connection Bivector as a proxy for
         * the information density of the gravitational field.
         * S ≈ |ω|²
         */
        entropicPartGA(coords) {
            const x = this._extractSpacetime(coords);

            // Get Connection Bivectors [omega_0, omega_1, omega_2, omega_3]
            // These represent the local gravitational field strength potentials
            const omegas = this.gaManifold.connectionBivector(x);

            // Scalar density S = sum |omega_mu|^2
            let S = 0;
            // omega_mu is a Multivector. normSq() gives magnitude squared.
            // Contraction with metric required? 
            // S = g^uv <omega_u * ~omega_v>

            const gInv = this.system.g.contravariant
                ? this.system.g.contravariant(x)
                : mat4Inv(this.system.g.covariant(x));

            for (let mu = 0; mu < 4; mu++) {
                for (let nu = 0; nu < 4; nu++) {
                    const g_uv = gInv[mu][nu];
                    if (abs(g_uv) > 1e-9) {
                        // Inner product of bivectors: <A B~>_0
                        // Effectively dot product if orthogonal basis
                        const prod = omegas[mu].mul(omegas[nu].reverse()).scalar();
                        S += g_uv * prod;
                    }
                }
            }

            // Determine sign? Gravity is attractive -> Potential decreases?
            // "Entropy" usually maximizes.
            // If S is information density, particles might flow towards it (entropic force).

            return this.entropicCoupling * S;
        }

        /**
         * Numerical gradient of Hamiltonian
         */
        gradient(pt, h = null) {
            h = h || this.h;
            const grad = {};
            const allCoords = this.manifold.allCoords;

            for (const c of allCoords) {
                const coordsPlus = { ...pt.coords, [c]: pt.coords[c] + h };
                const coordsMinus = { ...pt.coords, [c]: pt.coords[c] - h };
                grad[c] = (this.evaluate(coordsPlus) - this.evaluate(coordsMinus)) / (2 * h);
            }
            return grad;
        }

        /**
         * Contact Hamiltonian vector field X_H
         *
         * For canonical α = du - p_a dx^a:
         *   ẋ^a = ∂H/∂p_a
         *   ṗ_a = -∂H/∂x^a - p_a · ∂H/∂u
         *   u̇ = p_a · ∂H/∂p_a - H
         */
        vectorField(pt) {
            const grad = this.gradient(pt);
            const RH = grad[this.manifold.fiberCoord];
            const Hval = this.evaluate(pt.coords);

            const X = {};
            const n = this.manifold.n;

            // ẋ^a = ∂H/∂p_a
            for (let i = 0; i < n; i++) {
                const pCoord = this.manifold.momentaCoords[i];
                const xCoord = this.manifold.baseCoords[i];
                X[xCoord] = grad[pCoord];
            }

            // ṗ_a = -∂H/∂x^a - p_a · RH
            for (let i = 0; i < n; i++) {
                const pCoord = this.manifold.momentaCoords[i];
                const xCoord = this.manifold.baseCoords[i];
                const p = pt.get(pCoord);
                X[pCoord] = -grad[xCoord] - p * RH;
            }

            // u̇ = p · ∂_p H - H
            let pDotDpH = 0;
            for (let i = 0; i < n; i++) {
                const pCoord = this.manifold.momentaCoords[i];
                pDotDpH += pt.get(pCoord) * grad[pCoord];
            }
            X[this.manifold.fiberCoord] = pDotDpH - Hval;

            return X;
        }

        /**
         * Flow under contact Hamiltonian dynamics (RK4)
         */
        flow(pt, dt, steps = 1) {
            let current = pt.clone();
            const trajectory = [current.clone()];

            for (let i = 0; i < steps; i++) {
                const k1 = this.vectorField(current);

                const pt2 = current.add(k1, dt / 2);
                const k2 = this.vectorField(pt2);

                const pt3 = current.add(k2, dt / 2);
                const k3 = this.vectorField(pt3);

                const pt4 = current.add(k3, dt);
                const k4 = this.vectorField(pt4);

                const combined = {};
                for (const c of this.manifold.allCoords) {
                    combined[c] = (k1[c] + 2 * k2[c] + 2 * k3[c] + k4[c]) / 6;
                }

                current = current.add(combined, dt);
                trajectory.push(current.clone());
            }

            return trajectory;
        }

        /**
         * Track Hamiltonian evolution along trajectory
         */
        hamiltonianEvolution(trajectory) {
            return trajectory.map(pt => this.evaluate(pt.coords));
        }

        /**
         * Track entropy evolution along trajectory
         */
        entropyEvolution(trajectory) {
            return trajectory.map(pt => this.entropicPart(pt.coords));
        }
    }

    // ============================================================================
    // VII. FACTORY FUNCTIONS
    // ============================================================================

    /**
     * Create a simple matter configuration for testing
     *
     * @param {string} type - 'vacuum', 'scalar', 'em', or 'full'
     */
    function createMatterConfig(type = 'vacuum') {
        switch (type) {
            case 'vacuum':
                return new MatterInducedMetric({});

            case 'scalar':
                // Scalar field: φ = exp(-r²/σ²)
                return new MatterInducedMetric({
                    scalarField: x => {
                        const r2 = x[1] * x[1] + x[2] * x[2] + x[3] * x[3];
                        return exp(-r2);
                    }
                });

            case 'em':
                // Electromagnetic: Coulomb-like potential
                return new MatterInducedMetric({
                    vectorPotential: x => {
                        const r = sqrt(x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + 0.01);
                        return [1 / r, 0, 0, 0];  // A_0 = φ = 1/r
                    }
                });

            case 'full':
                // Full matter content
                return new MatterInducedMetric({
                    scalarField: x => exp(-(x[1] * x[1] + x[2] * x[2] + x[3] * x[3])),
                    vectorPotential: x => {
                        const r = sqrt(x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + 0.01);
                        return [1 / r, 0, 0, 0];
                    },
                    twoFormPotential: x => mat4Zero()
                });

            default:
                return new MatterInducedMetric({});
        }
    }

    /**
     * Create standard spacetime metrics (wrapper for compatibility)
     */
    const StandardMetrics = {
        minkowski: () => ({
            covariant: x => [
                [1, 0, 0, 0],
                [0, -1, 0, 0],
                [0, 0, -1, 0],
                [0, 0, 0, -1]
            ],
            contravariant: x => [
                [1, 0, 0, 0],
                [0, -1, 0, 0],
                [0, 0, -1, 0],
                [0, 0, 0, -1]
            ]
        }),

        schwarzschild: (M = 1) => ({
            covariant: x => {
                const [t, r, theta, phi] = x;
                const f = 1 - 2 * M / Math.max(r, 0.01);
                const sinTheta2 = sin(theta) ** 2;
                return [
                    [f, 0, 0, 0],
                    [0, -1 / f, 0, 0],
                    [0, 0, -r * r, 0],
                    [0, 0, 0, -r * r * sinTheta2]
                ];
            },
            contravariant: x => {
                const [t, r, theta, phi] = x;
                const f = 1 - 2 * M / Math.max(r, 0.01);
                const sinTheta2 = sin(theta) ** 2;
                return [
                    [1 / f, 0, 0, 0],
                    [0, -f, 0, 0],
                    [0, 0, -1 / (r * r), 0],
                    [0, 0, 0, -1 / (r * r * sinTheta2 + EPSILON)]
                ];
            }
        }),

        flrw: (a_func, k = 0) => ({
            covariant: x => {
                const [t, chi, theta, phi] = x;
                const aT = a_func(t);
                const a2 = aT * aT;
                let Sk2;
                if (k === 0) Sk2 = chi * chi;
                else if (k === 1) Sk2 = sin(chi) ** 2;
                else Sk2 = ((exp(chi) - exp(-chi)) / 2) ** 2;

                const sinTheta2 = sin(theta) ** 2;
                return [
                    [1, 0, 0, 0],
                    [0, -a2, 0, 0],
                    [0, 0, -a2 * Sk2, 0],
                    [0, 0, 0, -a2 * Sk2 * sinTheta2]
                ];
            }
        })
    };

    // ============================================================================
    // VIII. EXPORTS
    // ============================================================================

    const EntropicGravity = {
        // Core classes
        MatterInducedMetric,
        TwoMetricSystem,
        RelativeEntropyAction,
        EmergentCosmologicalConstant,
        ModifiedEinsteinSolver,
        EntropicGravityHamiltonian,

        // Factory functions
        createMatterConfig,
        StandardMetrics,

        // Matrix utilities (exposed for testing)
        utils: {
            mat4Add,
            mat4Sub,
            mat4Scale,
            mat4Mul,
            mat4Trace,
            mat4Det,
            mat4Inv,
            mat4Log,
            mat4Exp,
            mat4Identity,
            mat4Zero
        },

        // Constants
        EPSILON
    };

    // Export for different module systems
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = EntropicGravity;
    } else if (typeof define === 'function' && define.amd) {
        define([], () => EntropicGravity);
    } else {
        global.EntropicGravity = EntropicGravity;
    }

})(typeof window !== 'undefined' ? window : global);
