/**
 * Information Geometry Module
 *
 * Implements the contact geometry of probability distributions based on
 * John Baez's "Information Geometry" series (Part 18, 19).
 *
 * The manifold of probability distributions is treated as a Contact Manifold.
 * Coordinates:
 *  - q: Probabilities (p_i in stats, q_i in geometric analogy to position)
 *  - S: Shannon Entropy (analogue to Action or negative Free Energy)
 *  - p: Surprisal (conjugate to probability, analogue to momentum)
 *
 * Contact Form: alpha = dS - sum(p_i dq^i)
 *
 * @module information-geometry
 * @license MIT
 */

(function (global) {
    'use strict';

    // ============================================================================
    // IMPORT DEPENDENCIES
    // ============================================================================

    let GA;
    if (typeof require !== 'undefined') {
        try {
            GA = require('./multivector.js');
        } catch (e) {
            // Will use global
        }
    }
    // Fallback to global (browser)
    if (!GA && typeof global !== 'undefined') {
        GA = global.GA;
    }

    const { abs, log, sqrt } = Math;
    const EPSILON = 1e-10;

    // ============================================================================
    // PROBABILITY MANIFOLD CLASS
    // ============================================================================

    class ProbabilityManifold {
        /**
         * @param {number} n - Number of microstates (dimension of probability vector)
         */
        constructor(n) {
            this.n = n;
            // Tangent space dimension is 2*n + 1 for the full contact bundle relative to R^n?
            // Actually, the probability simplex has dim n-1.
            // The extended phase space (thermodynamic phase space) has dim 2n+1 (including S).
            // Let's model the full ambient space first:
            // q^1...q^n (probabilities)
            // p_1...p_n (surprisals)
            // S (entropy)
            // Total dim = 2n + 1

            // Algebra: Cl(2n+1, 0, 0)
            // Basis mapping:
            // e_1 ... e_n : dq^1 ... dq^n
            // e_{n+1} ... e_{2n} : dp_1 ... dp_n
            // e_{2n+1} : dS
            this.algebra = new GA.Algebra(2 * n + 1, 0, 0);

            // Pre-compute basis vectors for convenience
            this.dqs = [];
            this.dps = [];
            for (let i = 0; i < n; i++) {
                this.dqs.push(this.algebra.e(i + 1));
                this.dps.push(this.algebra.e(n + i + 1));
            }
            this.dS = this.algebra.e(2 * n + 1);
        }

        /**
         * Compute Shannon Entropy S(q) = - sum q_i ln(q_i)
         * @param {number[]} q - Probability distribution
         */
        entropy(q) {
            let s = 0;
            for (let x of q) {
                if (x > EPSILON) {
                    s -= x * log(x);
                }
            }
            return s;
        }

        /**
         * Compute Surprisal p_i = - ln(q_i) - 1
         * Note: This is the conjugate variable such that dS = sum p_i dq^i
         * S = - sum q_i ln q_i
         * dS/dq_k = - (1 * ln q_k + q_k * 1/q_k) = - ln q_k - 1
         *
         * @param {number[]} q - Probability distribution
         */
        surprisal(q) {
            // Baez Part 19: p_i = -ln(q_i) - 1
            // (The -1 comes from the product rule)
            return q.map(x => (x > EPSILON) ? -log(x) - 1 : 0);
        }

        /**
         * Construct the Contact Form alpha at a point in the phase space.
         * alpha = dS - p_i dq^i
         *
         * This 1-form lives in the cotangent bundle of the thermodynamic phase space.
         * In our GA representation, 1-forms obtained by "d" are represented as vectors.
         *
         * @returns {Multivector} The contact form alpha
         */
        contactForm(q, p) {
            // alpha = dS - sum(p_i * dq^i)
            // In our basis: e_{2n+1} - sum(p[i] * e_{i+1})

            let alpha = this.dS;
            for (let i = 0; i < this.n; i++) {
                // Term: p_i dq^i
                const term = this.dqs[i].scale(p[i]);
                alpha = alpha.sub(term);
            }
            return alpha;
        }

        /**
         * Verify the Contact Condition: alpha ^ (d alpha)^n != 0
         *
         * Wait, for a contact manifold of dimension 2k+1, the condition is alpha ^ (d alpha)^k != 0.
         * Here our manifold dimension is 2n+1.
         *
         * d(alpha) calculation:
         * alpha = dS - p_i dq^i
         * d(alpha) = d(dS) - d(p_i dq^i)
         *          = 0 - dp_i ^ dq^i
         *          = - sum_{i=1}^n dp_i ^ dq^i
         *          = sum_{i=1}^n dq^i ^ dp_i (using anti-symmetry)
         *
         * In terms of our basis:
         * dq^i = e_{i+1}
         * dp_i = e_{n+i+1}
         * d(alpha) = sum e_{i+1} ^ e_{n+i+1}
         *
         * Since this is constant (independent of point), we can just compute it once.
         */
        exteriorDerivativeAlpha() {
            // dAlpha = sum(dq^i ^ dp_i)
            let dAlpha = this.algebra.zero();
            for (let i = 0; i < this.n; i++) {
                // dq^i ^ dp_i
                const term = this.dqs[i].wedge(this.dps[i]);
                dAlpha = dAlpha.add(term);
            }
            return dAlpha;
        }

        /**
         * Compute the volume form alpha ^ (d alpha)^n
         * This should be non-zero (specifically, a multiple of the volume element).
         */
        volumeForm(q, p) {
            const alpha = this.contactForm(q, p);
            const dAlpha = this.exteriorDerivativeAlpha();

            // Compute (dAlpha)^n
            // We multiply dAlpha by itself n times using the wedge product
            let dAlphaN = dAlpha;
            for (let k = 1; k < this.n; k++) {
                dAlphaN = dAlphaN.wedge(dAlpha);
            }

            return alpha.wedge(dAlphaN);
        }

        /**
         * Check if a set of probabilities lies on the constraint surface M.
         * M is defined by p_i = surprisal(q_i) and S = entropy(q).
         * Actually, the "Legendrian Submanifold" is defined by p_i = dS/dq_i.
         * On this submanifold, alpha should vanish (pull-back of alpha is 0).
         */
        checkLegendrianCondition(q) {
            // This is a bit subtle to check with just single points.
            // A submanifold is Legendrian if alpha restricts to 0 on it.
            // i.e., tangent vectors to the submanifold belong to the kernel of alpha.
            // Tangent vectors to the submanifold generated by q -> (q, p(q), S(q))
            // are linear combinations of v_k = d/dq_k + (dp_i/dq_k) d/dp_i + (dS/dq_k) d/dS
            //
            // We can explicitly construct these tangent vectors and check alpha(v_k) = 0.

            const results = [];
            const p = this.surprisal(q);
            const alpha = this.contactForm(q, p);

            // For each coordinate q_k, construct the tangent vector v_k
            for (let k = 0; k < this.n; k++) {
                // v_k corresponds to varying q_k
                // Component along dq_k is 1.
                // Component along dS is dS/dq_k = p_k
                // Component along dp_j is dp_j/dq_k = d(-ln q_j - 1)/dq_k
                // = -1/q_j * delta_{jk}

                // v_k = e_{q_k} + p_k e_S + sum_j (dp_j/dq_k) e_{p_j}

                let vk = this.dqs[k]; // d/dq_k
                vk = vk.add(this.dS.scale(p[k])); // + p_k d/dS

                // Add dp terms
                // dp_k/dq_k = -1/q_k
                if (q[k] > EPSILON) {
                    const dp_dq = -1.0 / q[k];
                    vk = vk.add(this.dps[k].scale(dp_dq));
                }

                // Now compute contraction alpha(v_k).
                // Since alpha is a 1-form (vector in GA dual), and v_k is a vector,
                // we can use the dot product (scalar part of geometric product).
                // alpha . v_k
                // Both are represented as vectors in the Euclidean Cl(2n+1).
                // Wait, forms and vectors are dual.
                // In Euclidean metric, we can just dot them if we are consistent.
                // alpha = e_S - p_i e_qi
                // v_k = e_qk + p_k e_S + ...
                // alpha . v_k = (e_S . e_qk) + p_k(e_S . e_S) + ... - p_i(e_qi . e_qk) ...
                // = 0 + p_k(1) - p_k(1) = 0.

                const contraction = alpha.dot(vk).scalar();
                results.push({ k, contraction, passed: abs(contraction) < 1e-6 });
            }
            return results;
        }

    }

    // ============================================================================
    // EXPORTS
    // ============================================================================

    const InformationGeometryModule = { ProbabilityManifold };

    // Export for different module systems
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = InformationGeometryModule;
    }
    if (typeof define === 'function' && define.amd) {
        define([], () => InformationGeometryModule);
    }
    if (typeof global !== 'undefined') {
        global.InformationGeometry = InformationGeometryModule;
    }

})(typeof globalThis !== 'undefined' ? globalThis : (typeof window !== 'undefined' ? window : global));
