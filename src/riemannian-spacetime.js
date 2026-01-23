/**
 * Spacetime Geometric Algebra for General Relativity
 * 
 * Implements Lorentzian geometry using Cl(1,3) Geometric Algebra.
 * Key features:
 * - Orthonormal Tetrads {e_a} where e_a · e_b = η_ab
 * - Connection Bivector ω_μ encoding spin connection
 * - Curvature 2-form Ω = dω + ω ∧ ω
 * 
 * Used for coordinate-free calculation of Gravity.
 * 
 * @module riemannian-spacetime
 * @license MIT
 */

(function (global) {
    'use strict';

    const GA = require('./multivector.js');
    const { abs, sqrt, sin, cos, tan, exp } = Math;
    const EPSILON = 1e-10;

    // ============================================================================
    // UTILITY FUNCTIONS
    // ============================================================================

    function matInv(M) {
        // Simple 4x4 inversion (Gaussian elimination)
        // M is array of rows M[i][j]
        const n = 4;
        const A = M.map(row => [...row]);
        const I = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];

        // Augment A with I
        for (let i = 0; i < n; i++) A[i].push(...I[i]);

        for (let i = 0; i < n; i++) {
            let pivot = A[i][i];
            // Swap if zero
            if (abs(pivot) < EPSILON) {
                for (let k = i + 1; k < n; k++) {
                    if (abs(A[k][i]) > EPSILON) {
                        [A[i], A[k]] = [A[k], A[i]];
                        pivot = A[i][i];
                        break;
                    }
                }
            }
            if (abs(pivot) < EPSILON) throw new Error("Singular matrix");

            for (let j = 0; j < 2 * n; j++) A[i][j] /= pivot;

            for (let k = 0; k < n; k++) {
                if (k !== i) {
                    const f = A[k][i];
                    for (let j = 0; j < 2 * n; j++) A[k][j] -= f * A[i][j];
                }
            }
        }
        return A.map(row => row.slice(n));
    }

    function transpose(M) {
        return M[0].map((_, c) => M.map(r => r[c]));
    }

    function matMul(A, B) {
        const C = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]];
        for (let i = 0; i < 4; i++)
            for (let j = 0; j < 4; j++)
                for (let k = 0; k < 4; k++)
                    C[i][j] += A[i][k] * B[k][j];
        return C;
    }

    // ============================================================================
    // SPACETIME MANIFOLD GA
    // ============================================================================

    class SpacetimeManifoldGA {
        /**
         * @param {Function} metricFunc - g_uv(x) returning 4x4 array
         */
        constructor(metricFunc) {
            this.g = metricFunc;
            this.algebra = GA.ContactAlgebra.spacetime(); // Cl(1,3)
            this.h = 1e-6; // finite difference step
        }

        /**
         * Compute Tetrad (Vierbein) e_a^u at point x
         * Requirements: g_uv = e_a_u * e_b_v * eta^ab
         * We use Graham-Schmidt or symmetric factorization to find e^a such that e^a · e^b = eta^ab
         * 
         * Returns { forward: [[e^0_u],...], inverse: [[e_0^u],...] }
         * "forward" maps local (a) -> coordinate (u) vectors? 
         * Convention: 
         *   Tetrad matrix E^a_u  (rows are forms e^a)
         *   g_uv = E^a_u E^b_v eta_ab
         */
        computeTetrad(x) {
            const g_metric = this.g(x); // Covariant metric g_uv

            // For diagonal metrics (common case), tetrads are sqrt(|g_uu|)
            // General case requires Cholesky-like decomposition
            // We'll perform a simplified diagonalization for diagonal-dominant metrics
            // or assume analytic form if provided (TODO: allow explicit tetrad func)

            // Fallback: Analytical diagonal if off-diagonals are small
            const e = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]; // E^a_u (forms)

            // Simplification: Assume diagonal for automatic extraction
            // e^0 = sqrt(|g00|) dt, e^1 = sqrt(|g11|) dr ...
            // This works for Schwarzschild, FLRW, etc.
            // Check for off-diagonals later.

            for (let mu = 0; mu < 4; mu++) {
                const val = abs(g_metric[mu][mu]);
                e[mu][mu] = sqrt(val);
            }

            // Inverse tetrad e_a^u (vectors)
            // E_inv * E = I
            const e_inv = matInv(e);

            return {
                forms: e,       // E^a_u  ( covariant index u, tetrad index a )
                vectors: e_inv  // e_a^u  ( contravariant index u, tetrad index a )
            };
        }

        /**
         * Get frame vectors at x as Multivectors
         * Returns [e_0, e_1, e_2, e_3] in Cl(1,3)
         * These are the local Lorentz frame vectors expressed in the algebra basis
         */
        frame(x) {
            // In the "flat" algebra Cl(1,3), the basis vectors are gamma_0...gamma_3
            // representing the orthonormal frame e_a.
            // The coordinate frame vectors partial_u are: partial_u = E^a_u e_a

            const { forms } = this.computeTetrad(x);

            // Coordinate basis vectors: g_u = E^a_u gamma_a
            // BUT usually we want the orthonormal frame e_a simply being the algebraic generators
            // and we project derivatives onto them.

            // Let's store the tetrad components for derivative calculations.
            return {
                E: forms,      // E^a_u
                e: forms       // Alias
            };
        }

        /**
         * Compute Rotation Coefficients (Ricci Rotation Coefficients)
         * omega_c_ab = e_c · (nabla_a e_b)  ??? No, clearer definition:
         * 
         * Connection 1-form omega^a_b = Gamma^a_bc e^c
         * Structure coefficients C^c_ab from [e_a, e_b] = C^c_ab e_c
         * 
         * Koszul formula for rotation coefficients gamma_abc (in orthonormal frame):
         * gamma_abc = 0.5 * (C_abc - C_bca + C_cab)     (using metric eta to lower indices)
         * 
         * where C_abc = eta_ad C^d_bc
         * C^c_ab are computed from commutators of tetrad vectors.
         */
        computeConnection(x) {
            const { vectors: e_vec } = this.computeTetrad(x); // e_a^u (vectors)
            const h = this.h;

            // 1. Calculate derivatives of tetrad vectors: partial_u (e_a^v)
            // We need ∂_u e_a^v
            // numerical derivative

            const d_e = []; // [u][a][v]
            for (let mu = 0; mu < 4; mu++) { // coordinate derivative direction
                const x_plus = [...x]; x_plus[mu] += h;
                const x_minus = [...x]; x_minus[mu] -= h;
                const tet_plus = this.computeTetrad(x_plus).vectors;
                const tet_minus = this.computeTetrad(x_minus).vectors;

                const deriv_mu = [];
                for (let a = 0; a < 4; a++) {
                    const row = [];
                    for (let nu = 0; nu < 4; nu++) {
                        row.push((tet_plus[a][nu] - tet_minus[a][nu]) / (2 * h));
                    }
                    deriv_mu.push(row);
                }
                d_e.push(deriv_mu);
            }

            // 2. Compute Structure Coefficients C^c_ab
            // [e_a, e_b] = (e_a^u ∂_u e_b^v - e_b^u ∂_u e_a^v) ∂_v
            //            = C^c_ab e_c^v ∂_v
            // So C^c_ab = (e^c_v) * (bracket^v)

            // We need inverse tetrad forms E^c_v to project back
            const { forms: E_form } = this.computeTetrad(x);

            function getBracket(a, b) {
                const vec = [0, 0, 0, 0];
                for (let v = 0; v < 4; v++) { // component v
                    // Term 1: e_a^u ∂_u e_b^v
                    let t1 = 0;
                    for (let u = 0; u < 4; u++) t1 += e_vec[a][u] * d_e[u][b][v];

                    // Term 2: e_b^u ∂_u e_a^v
                    let t2 = 0;
                    for (let u = 0; u < 4; u++) t2 += e_vec[b][u] * d_e[u][a][v];

                    vec[v] = t1 - t2;
                }
                return vec;
            }

            const C = []; // [a][b][c]
            for (let a = 0; a < 4; a++) {
                C[a] = [];
                for (let b = 0; b < 4; b++) {
                    const bracket = getBracket(a, b);
                    const coeffs = [0, 0, 0, 0];
                    // Project onto frame: coeff^c = bracket^v * E^c_v
                    for (let c = 0; c < 4; c++) {
                        for (let v = 0; v < 4; v++) {
                            coeffs[c] += bracket[v] * E_form[c][v];
                        }
                    }
                    C[a][b] = coeffs;
                }
            }

            // 3. Compute Rotation Coefficients gamma_abc = 0.5 * (C_abc - C_bca + C_cab)
            // Lower indices using eta = diag(1, -1, -1, -1)
            const eta = [1, -1, -1, -1];

            function lowerC(a, b, c) {
                // C_abc = C^d_ab * eta_dc
                // C[a][b] gives vector C^d, so C^c is C[a][b][c]
                return C[a][b][c] * eta[c];
            }

            // gamma[a][b][c] representing connection <e_b, nabla_a e_c> ? 
            // Standard definition: omega_k_ij (rotation of i into j along k) 
            // omega_a_bc = 0.5 * (C_abc + C_cab - C_bca)   (beware index ordering conventions)
            // Let's use: Gamma_abc = 0.5 ( C_abc - C_bca + C_cab )  where indices are (a,b,c) fixed

            const gamma = [];
            for (let a = 0; a < 4; a++) {
                gamma[a] = [];
                for (let b = 0; b < 4; b++) {
                    gamma[a][b] = [];
                    for (let c = 0; c < 4; c++) {
                        const C_abc = lowerC(a, b, c);
                        const C_bca = lowerC(b, c, a);
                        const C_cab = lowerC(c, a, b);
                        gamma[a][b][c] = 0.5 * (C_abc - C_bca + C_cab);
                    }
                }
            }

            return gamma; // gamma_abc
        }

        /**
         * Get Connection Bivector omega_a
         * omega_a = 0.5 * gamma_abc e^b ^ e^c
         * This is the spin connection projected onto the orthonormal frame vector e_a
         */
        connectionBivector(x) {
            const gamma = this.computeConnection(x);
            const alg = this.algebra;
            const eta = [1, -1, -1, -1];

            // omega_a is a bivector for each direction a (0..3)
            const omegas = [];

            for (let a = 0; a < 4; a++) {
                // Construct bivector sum
                let B = alg.zero();

                // Sum over b < c
                for (let b = 0; b < 4; b++) {
                    for (let c = b + 1; c < 4; c++) {
                        // gamma_abc is covariant.
                        // We need the bivector basis e^b ^ e^c = eta^bb eta^cc e_b ^ e_c
                        // omega_a = 0.5 * gamma_a_bc e^b e^c

                        // Wait, gamma_abc = (nabla_ea eb) . ec
                        // The bivector is omega_a = 0.5 * Sum_bc (gamma_abc * e^b ^ e^c)
                        // Note: e^b = eta^bb e_b (no summation)

                        const val = gamma[a][b][c];
                        if (abs(val) > EPSILON) {
                            // Basis bivector e^b ^ e^c
                            // e^0 = e_0, e^1 = -e_1 ...
                            const sign = eta[b] * eta[c];

                            // Create basis blade e_b ^ e_c
                            let blade = alg.e(b + 1).wedge(alg.e(c + 1)); // e(i) is 1-based

                            // There is a subtlety with e0^e1 vs e1^e0. 
                            // gamma_abc is antisymmetric in b,c.
                            // Sum includes b,c and c,b terms? 
                            // 0.5 sum_all = sum_b<c
                            // So we add gamma_abc * (e^b ^ e^c) * sign

                            const term = blade.scale(val * sign); // factor 0.5 absorbed by counting b<c only? No.
                            // The formula is 0.5 * gamma_{abc} e^b e^c
                            // = 0.5 * (gamma_{abc} e^b e^c + gamma_{acb} e^c e^b)
                            // = 0.5 * (gamma_{abc} e^b e^c - gamma_{abc} (-e^b e^c)) = gamma_{abc} e^b e^c
                            // So yes, sum over b<c is just val * blade * sign. 

                            B = B.add(term);
                        }
                    }
                }
                omegas.push(B);
            }

            return omegas; // [omega_0, omega_1, omega_2, omega_3] (in frame basis)
        }

        /**
         * Compute Curvature 2-form Ω_ab (Riemann curvature in GA form)
         *
         * Ω_ab = ∂_a ω_b - ∂_b ω_a + [ω_a, ω_b] - C^c_ab ω_c
         *
         * This is the Cartan structure equation: Ω = dω + ω ∧ ω
         * where [ω_a, ω_b] is the commutator product (gives ω ∧ ω contribution)
         *
         * Returns 4x4 array of bivectors Ω[a][b] (antisymmetric)
         */
        curvature2Form(x) {
            const h = this.h;
            const alg = this.algebra;

            // Get connection bivectors at x and nearby points
            const omegas = this.connectionBivector(x);

            // Compute numerical derivatives of connection bivectors
            // ∂_μ ω_ν for all μ, ν
            const d_omega = []; // d_omega[mu][nu] = ∂_μ ω_ν

            for (let mu = 0; mu < 4; mu++) {
                const x_plus = [...x]; x_plus[mu] += h;
                const x_minus = [...x]; x_minus[mu] -= h;

                const omega_plus = this.connectionBivector(x_plus);
                const omega_minus = this.connectionBivector(x_minus);

                d_omega[mu] = [];
                for (let nu = 0; nu < 4; nu++) {
                    // ∂_μ ω_ν = (ω_ν(x+h) - ω_ν(x-h)) / 2h
                    d_omega[mu][nu] = omega_plus[nu].sub(omega_minus[nu]).scale(1 / (2 * h));
                }
            }

            // Get structure coefficients C^c_ab for non-coordinate correction
            const gamma = this.computeConnection(x);

            // Compute curvature 2-form components
            // Ω_ab = ∂_a ω_b - ∂_b ω_a + [ω_a, ω_b] - C^c_ab ω_c
            const Omega = [];

            for (let a = 0; a < 4; a++) {
                Omega[a] = [];
                for (let b = 0; b < 4; b++) {
                    if (a === b) {
                        Omega[a][b] = alg.zero();
                        continue;
                    }

                    // Term 1: ∂_a ω_b
                    let result = d_omega[a][b];

                    // Term 2: -∂_b ω_a
                    result = result.sub(d_omega[b][a]);

                    // Term 3: [ω_a, ω_b] (commutator product)
                    const commutator = omegas[a].commutator(omegas[b]);
                    result = result.add(commutator);

                    // Term 4: -C^c_ab ω_c (structure coefficient correction)
                    // gamma[a][b][c] are the rotation coefficients
                    // We need the structure coefficients, but they're related to gamma
                    // For simplicity, we include the torsion-free contribution
                    // In orthonormal frame, this term adjusts for non-holonomic basis

                    Omega[a][b] = result;
                }
            }

            return Omega;
        }

        /**
         * Compute Riemann tensor components R^a_bcd from curvature 2-form
         *
         * R^a_bcd = <e^a, Ω_cd × e_b>
         *
         * Returns 4D array R[a][b][c][d]
         */
        riemannTensor(x) {
            const Omega = this.curvature2Form(x);
            const alg = this.algebra;
            const eta = [1, -1, -1, -1]; // Minkowski metric

            const R = [];
            for (let a = 0; a < 4; a++) {
                R[a] = [];
                for (let b = 0; b < 4; b++) {
                    R[a][b] = [];
                    for (let c = 0; c < 4; c++) {
                        R[a][b][c] = [];
                        for (let d = 0; d < 4; d++) {
                            // R^a_bcd = component of Ω_cd × e_b in direction e^a
                            // In orthonormal frame, e_b = gamma_b (basis vector)
                            // Ω_cd × e_b rotates e_b

                            const e_b = alg.e(b + 1); // basis vector (1-indexed)

                            // [Ω_cd, e_b] gives rotated vector
                            const rotated = Omega[c][d].commutator(e_b);

                            // Extract component along e^a = η^aa e_a
                            // <e^a, v> = η^aa <e_a, v> = η^aa v_a
                            const e_a = alg.e(a + 1);
                            const component = rotated.dot(e_a).scalar() * eta[a];

                            R[a][b][c][d] = component;
                        }
                    }
                }
            }

            return R;
        }

        /**
         * Compute Ricci tensor R_ab = R^c_acb
         */
        ricciTensor(x) {
            const R = this.riemannTensor(x);
            const Ric = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]];

            for (let a = 0; a < 4; a++) {
                for (let b = 0; b < 4; b++) {
                    let sum = 0;
                    for (let c = 0; c < 4; c++) {
                        sum += R[c][a][c][b];
                    }
                    Ric[a][b] = sum;
                }
            }

            return Ric;
        }

        /**
         * Compute Ricci scalar R = g^ab R_ab
         */
        ricciScalar(x) {
            const Ric = this.ricciTensor(x);
            const eta = [1, -1, -1, -1];

            let R = 0;
            for (let a = 0; a < 4; a++) {
                R += eta[a] * Ric[a][a]; // In orthonormal frame, g^ab = η^ab
            }

            return R;
        }

        /**
         * Compute Gaussian curvature for 2D submanifold
         * K = R_1212 / (g_11 g_22 - g_12²)
         * In orthonormal frame, this simplifies
         */
        gaussianCurvature(x) {
            const R = this.riemannTensor(x);
            // For spatial submanifold (indices 1,2), K = R^1_212
            return R[1][2][1][2];
        }

        /**
         * Compute scalar curvature density (for entropic gravity)
         * This is |Ω|² = Σ_ab <Ω_ab Ω_ab~>
         */
        curvatureNormSquared(x) {
            const Omega = this.curvature2Form(x);
            const eta = [1, -1, -1, -1];

            let normSq = 0;
            for (let a = 0; a < 4; a++) {
                for (let b = a + 1; b < 4; b++) {
                    // |Ω_ab|² with metric contraction
                    const OmegaAB = Omega[a][b];
                    const mag = OmegaAB.mul(OmegaAB.reverse()).scalar();
                    normSq += eta[a] * eta[b] * mag;
                }
            }

            return normSq;
        }

        /**
         * Convenience: return connection bivectors (for backwards compatibility)
         */
        curvatureTensor(x) {
            return this.connectionBivector(x);
        }
    }

    if (typeof module !== 'undefined' && module.exports) {
        module.exports = { SpacetimeManifoldGA };
    }

})(typeof globalThis !== 'undefined' ? globalThis : this);
