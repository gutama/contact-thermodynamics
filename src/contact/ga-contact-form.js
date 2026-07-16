/**
 * GA-Native Contact Form Module (Phase 1)
 *
 * Bridges the geometric-algebra layer (`src/algebra/multivector.js`) with the
 * contact-geometry layer (`src/contact/`). It represents the canonical contact
 * 1-form α = ds − p_a dq^a, its exterior derivative dα, the Reeb field R, and
 * the contact Hamiltonian vector field X_H as genuine `Multivector` objects,
 * rather than the display-string / scalar-evaluator representation used by
 * `ContactManifold`.
 *
 * This addresses research gaps G1 (no GA object for α/dα, stubbed
 * non-degeneracy check) and G2 (the GA and contact layers were disconnected).
 * It is *additive*: the conventional `ContactManifold` / `ContactHamiltonian`
 * classes are untouched.
 *
 * ── Construction ──────────────────────────────────────────────────────────
 * The cotangent basis 1-forms {dq^a, dp_a, ds} are modeled as an orthonormal
 * grade-1 basis of a Clifford/Grassmann algebra Cl(2n+1, 0, 0):
 *
 *     index 0 .. n-1    →  dq^a   (base coordinates)
 *     index n .. 2n-1   →  dp_a   (conjugate momenta)
 *     index 2n          →  ds     (fiber / Reeb axis)
 *
 * Only the outer product (`op`/`wedge`) is used for α, dα and the
 * non-degeneracy wedge, so those results are metric-independent (pure
 * Grassmann algebra). The natural pairing α(v) and the interior product ι_v
 * are realized by the left contraction (`ip`/`dot`): identifying the tangent
 * basis ∂_i with the cotangent basis dx^i under the Euclidean coordinate
 * identity metric makes `ip(v, ω)` reproduce ⟨ω, v⟩ and ι_v ω exactly in
 * Darboux coordinates (this is the coordinate musical isomorphism, faithful for
 * the pairing/contraction operations the contact structure needs).
 *
 * The Reeb axis ds is the *radical* (null direction) of the degenerate 2-form
 * dα: dα is the canonical symplectic form on the (q,p)-plane and annihilates
 * ds, so ds spans ker(dα). This is the "Heisenberg-type null/Reeb direction"
 * the foundations doc calls for — realized here as the kernel of dα rather than
 * as a separately-signed generator, which keeps a single clean algebra and
 * avoids a bespoke degenerate metric.
 *
 * Alternative considered and rejected: a degenerate algebra Cl(n, n, 1) with an
 * explicit null generator for ds and opposite signs for q vs p (so that dα
 * would be the metric bivector). Rejected because (a) the contact conditions
 * α(R)=1, ι_R dα=0, α(X_H)=−H, ι_{X_H} dα = dH − R(H)α are all defined via the
 * natural pairing / interior product and need no metric, so a Euclidean-
 * identity carrier is sufficient and simpler; and (b) the degenerate generator
 * makes `ip`-based contractions ambiguous (contracting into a null direction),
 * complicating the very operations we rely on. The null Reeb structure is still
 * present — as ker(dα) — without paying that cost.
 *
 * @module contact/ga-contact-form
 * @license MIT
 */

(function (global) {
    'use strict';

    let GA;
    if (typeof module !== 'undefined' && module.exports) {
        GA = require('../algebra/multivector.js');
    } else {
        GA = global.GA;
    }
    const { Algebra } = GA;

    const EPSILON = 1e-10;

    // Guard: the generic Algebra precomputes O(size²) = O(4^n) product tables.
    // Cl(2n+1) is only feasible for modest n. dim 11 (n=5, size 2048) is fine;
    // dim 13 (Grand manifold, n=6, size 8192) needs a sparse backend (Phase 2).
    const MAX_ALGEBRA_SIZE = 4096;

    /**
     * GA-native representation of a canonical contact manifold's structure.
     */
    class GAContactForm {
        /**
         * @param {ContactManifold} manifold - conventional contact manifold
         *   (supplies baseCoords, momentaCoords, fiberCoord, n).
         */
        constructor(manifold) {
            this.manifold = manifold;
            this.n = manifold.n;
            this.dim = 2 * this.n + 1;

            const size = 1 << this.dim;
            if (size > MAX_ALGEBRA_SIZE) {
                throw new Error(
                    `GAContactForm: base dim n=${this.n} → Cl(${this.dim}) has ${size} ` +
                    `blades, exceeding MAX_ALGEBRA_SIZE=${MAX_ALGEBRA_SIZE}. The generic ` +
                    `Algebra precompute is O(4^n); large manifolds (e.g. Grand M₁₃) need ` +
                    `a sparse exterior-algebra backend (deferred to Phase 2).`
                );
            }

            // Basis ordering: [dq^a] , [dp_a] , ds
            this._coordOrder = [
                ...manifold.baseCoords,
                ...manifold.momentaCoords,
                manifold.fiberCoord
            ];
            this._indexOf = {};
            const basisNames = [];
            this._coordOrder.forEach((c, i) => {
                this._indexOf[c] = i;          // 0-based basis index
                basisNames.push(`d${c}`);
            });

            this.algebra = new Algebra(this.dim, 0, 0, basisNames);

            // dα is coordinate-independent; build it once.
            this._dAlpha = this._buildDAlpha();
        }

        // ---- basis helpers -------------------------------------------------

        /** 1-indexed basis vector e(i) for a coordinate's 1-form dx. */
        _e(coord) {
            return this.algebra.e(this._indexOf[coord] + 1);
        }

        /** bitmap of the grade-1 blade for a coordinate's 1-form. */
        _bit(coord) {
            return 1 << this._indexOf[coord];
        }

        /**
         * Build a grade-1 multivector (a 1-form or a tangent vector) from a
         * {coord: value} dictionary.
         */
        oneForm(components) {
            let mv = this.algebra.zero();
            for (const [coord, val] of Object.entries(components)) {
                if (val === undefined || Math.abs(val) < EPSILON) continue;
                if (this._indexOf[coord] === undefined) continue;
                mv = mv.add(this._e(coord).scale(val));
            }
            return mv;
        }

        /** Coefficient of the grade-1 basis blade for `coord` in mv. */
        coeff(mv, coord) {
            return mv.coeffs[this._bit(coord)] || 0;
        }

        // ---- α and dα ------------------------------------------------------

        /**
         * Canonical contact 1-form α = ds − p_a dq^a as a grade-1 Multivector.
         * @param {ContactPoint} pt
         */
        contactForm(pt) {
            const m = this.manifold;
            let alpha = this._e(m.fiberCoord);            // ds
            for (let a = 0; a < this.n; a++) {
                const p = pt.get(m.momentaCoords[a]);
                alpha = alpha.sub(this._e(m.baseCoords[a]).scale(p));  // − p_a dq^a
            }
            return alpha;
        }

        /**
         * Exterior derivative dα = −dp_a ∧ dq^a = Σ_a dq^a ∧ dp_a (grade-2,
         * coordinate-independent).
         */
        _buildDAlpha() {
            const m = this.manifold;
            let dAlpha = this.algebra.zero();
            for (let a = 0; a < this.n; a++) {
                const dq = this._e(m.baseCoords[a]);
                const dp = this._e(m.momentaCoords[a]);
                dAlpha = dAlpha.add(dq.wedge(dp));   // dq^a ∧ dp_a  (= −dp_a∧dq^a)
            }
            return dAlpha;
        }

        contactFormDerivative() {
            return this._dAlpha.clone();
        }

        // ---- genuine non-degeneracy check ---------------------------------

        /**
         * Wedge α ∧ (dα)^n for arbitrary supplied forms. Returns the resulting
         * top-grade multivector. Exposed so tests can feed a *broken* dα for a
         * negative (degeneracy-detection) check.
         * @param {Multivector} alpha  - grade-1
         * @param {Multivector} dAlpha - grade-2
         */
        wedgeAlphaDAlphaPower(alpha, dAlpha) {
            let power = this.algebra.scalar(1);          // (dα)^0
            for (let k = 0; k < this.n; k++) {
                power = power.wedge(dAlpha);              // (dα)^k → (dα)^{k+1}
            }
            return alpha.wedge(power);                   // α ∧ (dα)^n
        }

        /**
         * Genuine non-degeneracy check: computes α ∧ (dα)^n via the GA outer
         * product (replacing the hardcoded n! stub in
         * ContactManifold.verifyContactCondition). For the canonical form the
         * result is ± n! times the pseudoscalar.
         *
         * @param {ContactPoint} pt
         * @returns {{ form: Multivector, pseudoscalarCoeff: number,
         *             magnitude: number, isNonDegenerate: boolean }}
         */
        nonDegeneracy(pt) {
            const alpha = this.contactForm(pt);
            const top = this.wedgeAlphaDAlphaPower(alpha, this._dAlpha);
            const psBitmap = this.algebra.size - 1;      // top-grade blade
            const coeff = top.coeffs[psBitmap] || 0;
            return {
                form: top,
                pseudoscalarCoeff: coeff,
                magnitude: Math.abs(coeff),
                isNonDegenerate: Math.abs(coeff) > EPSILON
            };
        }

        // ---- Reeb field ----------------------------------------------------

        /**
         * Reeb vector field R = ∂/∂s as a grade-1 Multivector (the ds basis
         * vector under the coordinate identity metric).
         */
        reebField() {
            return this._e(this.manifold.fiberCoord);
        }

        /**
         * Verify the Reeb defining equations via GA: α(R) = ip(R, α) and
         * ι_R dα = ip(R, dα).
         * @param {ContactPoint} pt
         * @returns {{ alphaR: number, iRdAlphaIsZero: boolean, iRdAlpha: Multivector }}
         */
        verifyReeb(pt) {
            const R = this.reebField();
            const alpha = this.contactForm(pt);
            const alphaR = R.dot(alpha).scalar();        // ⟨α, R⟩
            const iRdAlpha = R.dot(this._dAlpha);        // ι_R dα
            return {
                alphaR,
                iRdAlpha,
                iRdAlphaIsZero: iRdAlpha.isZero()
            };
        }

        // ---- contact Hamiltonian vector field -----------------------------

        /**
         * Contact Hamiltonian vector field X_H, constructed by *solving* its GA
         * defining equations
         *     ι_{X_H} dα = dH − R(H) α,     α(X_H) = −H
         * rather than by copying the Bravetti EOM.
         *
         * Reading off the grade-1 coefficients of the RHS multivector:
         *   ι_X dα = Σ_a [ (X·dq^a) dp_a − (X·dp_a) dq^a ], hence
         *     X^{q_a} =  coeff_{dp_a}(RHS)
         *     X^{p_a} = −coeff_{dq_a}(RHS)
         *   and the fiber component from α(X)=−H:
         *     X^{s}   = −H + Σ_a p_a X^{q_a}.
         *
         * @param {ContactHamiltonian} hamiltonian
         * @param {ContactPoint} pt
         * @returns {{ vector: Multivector, components: Object }}
         */
        hamiltonianVectorField(hamiltonian, pt) {
            const m = this.manifold;
            const grad = hamiltonian.gradient(pt);
            // ∂H/∂s = 0 for many Hamiltonians; an analytic gradient may simply
            // omit the fiber component, leaving it undefined. Default to 0 so
            // R(H) α does not inject NaN into the vector field.
            const RH = grad[m.fiberCoord] ?? 0;
            const Hval = hamiltonian.evaluate(pt);

            // dH as a grade-1 form, and α at the point.
            const dH = this.oneForm(grad);
            const alpha = this.contactForm(pt);

            // RHS = dH − R(H) α
            const rhs = dH.sub(alpha.scale(RH));

            const components = {};
            let pDotXq = 0;
            for (let a = 0; a < this.n; a++) {
                const qa = m.baseCoords[a];
                const pa = m.momentaCoords[a];
                const Xq = this.coeff(rhs, pa);          //  coeff of dp_a
                const Xp = -this.coeff(rhs, qa);         // −coeff of dq_a
                components[qa] = Xq;
                components[pa] = Xp;
                pDotXq += pt.get(pa) * Xq;
            }
            components[m.fiberCoord] = -Hval + pDotXq;

            return { vector: this.oneForm(components), components };
        }

        /**
         * Verify the X_H defining equations via GA. Returns residuals that
         * should all be ~0 for a correct construction.
         * @param {ContactHamiltonian} hamiltonian
         * @param {ContactPoint} pt
         */
        verifyHamiltonianVectorField(hamiltonian, pt) {
            const m = this.manifold;
            const grad = hamiltonian.gradient(pt);
            // Default missing ∂H/∂s to 0 (see hamiltonianVectorField) so the
            // residual check stays finite rather than NaN.
            const RH = grad[m.fiberCoord] ?? 0;
            const Hval = hamiltonian.evaluate(pt);

            const { vector: X } = this.hamiltonianVectorField(hamiltonian, pt);
            const alpha = this.contactForm(pt);
            const dH = this.oneForm(grad);

            // α(X_H) should equal −H
            const alphaX = X.dot(alpha).scalar();
            const contactResidual = Math.abs(alphaX - (-Hval));

            // ι_{X_H} dα should equal dH − R(H) α
            const lhs = X.dot(this._dAlpha);
            const rhs = dH.sub(alpha.scale(RH));
            const symplecticResidual = lhs.sub(rhs).norm();

            return { alphaX, contactResidual, symplecticResidual };
        }
    }

    const GAContactModule = { GAContactForm };

    if (typeof module !== 'undefined' && module.exports) {
        module.exports = GAContactModule;
    } else if (typeof define === 'function' && define.amd) {
        define([], () => GAContactModule);
    } else {
        global.GAContactForm = GAContactModule;
    }

})(typeof window !== 'undefined' ? window : global);
