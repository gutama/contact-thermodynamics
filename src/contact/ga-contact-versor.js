/**
 * GA-Native Contact Versor Module (Phase 1b)
 *
 * A focused follow-up to Phase 1 (`src/contact/ga-contact-form.js`). Phase 1
 * built α, dα, the Reeb field and X_H in a *Euclidean* carrier Cl(2n+1,0,0) and
 * an external reviewer showed that construction is **signature-invariant**:
 * running it in Cl(5,0,0), Cl(0,5,0), Cl(3,2,0), Cl(2,2,1) gives numerically
 * identical results, so the geometric product / metric does no work there — it
 * is exterior (Grassmann) algebra in Clifford clothing.
 *
 * This module tests the reviewer's follow-up question in the *degenerate*
 * carrier that Phase 1 §3.3 rejected: **Cl(n,n,1)**, where dα is a genuine
 * metric bivector. The question: can contactomorphisms (diffeomorphisms φ with
 * φ*α = f·α, f a nonvanishing conformal factor) be represented as **versors** —
 * exponentials of bivectors applied by a sandwich v ↦ V v V⁻¹, exactly as
 * ordinary rotors generate rotations/boosts?
 *
 * ── Carrier Cl(n,n,1) ─────────────────────────────────────────────────────
 *   eq_1..eq_n   e² = +1     (q-directions)
 *   ep_1..ep_n   e² = −1     (p-directions)
 *   es           e² =  0     (null generator, the ds / Reeb axis)
 *
 * With this signature the contact bivector
 *     B = dα = Σ_a eq_a ∧ ep_a
 * squares to (eq_a ep_a)² = −eq_a² ep_a² = −(+1)(−1) = +1 per pair — a genuine
 * *hyperbolic* metric bivector, unlike the Euclidean carrier where the choice of
 * signature was irrelevant. This is the whole point: here the metric is supposed
 * to do work.
 *
 * ── Darboux ↔ light-cone identification ─────────────────────────────────────
 * The physical Darboux coordinates (q_a, p_a) are identified with a **null
 * (light-cone) basis** of each hyperbolic (eq_a, ep_a) plane:
 *     nq_a = eq_a − ep_a          (⟨nq_a, nq_a⟩ = 0)
 *     np_a = (eq_a + ep_a)/2      (⟨np_a, np_a⟩ = 0,  ⟨nq_a, np_a⟩ = 1)
 * A point (q, p, s) is carried by the grade-1 vector
 *     X = Σ_a ( q_a nq_a + p_a np_a ) + s·es.
 * A hyperbolic versor V = exp(θ B/2) then acts *diagonally* on this basis
 * (nq_a ↦ e^{θ} nq_a, np_a ↦ e^{−θ} np_a), i.e. as the symplectic squeeze
 * q_a ↦ e^{θ} q_a, p_a ↦ e^{−θ} p_a — the natural candidate contactomorphism.
 *
 * Everything reuses the generic `Algebra`/`Multivector` engine in
 * `src/algebra/multivector.js` (which already supports r>0 degenerate
 * generators and `exp`/`sandwich`); this module only adds the contact-specific
 * basis bookkeeping and an **honest, non-tautological contactomorphism test**
 * (`contactResidual`) that pulls α back through the induced coordinate map and
 * checks φ*α = f·α against the mathematical definition.
 *
 * See `docs/PHASE1B_VERSOR_CONTACTOMORPHISMS.md` for the findings (which are
 * mixed: the squeeze IS a genuine versor-generated contactomorphism, but the
 * Reeb translation and the conformal dilation are provably NOT).
 *
 * @module contact/ga-contact-versor
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

    // Guard: the generic Algebra precomputes O(size²) = O(4^dim) product tables.
    // Cl(n,n,1) has dim = 2n+1, so n=5 → dim 11 (size 2048) is fine; n=6 → dim 13
    // (size 8192) needs a sparse backend (Phase 2). Mirrors the guard in
    // `ga-contact-form.js`. Also caps dim < 31 so the `1 << dim` blade sizing
    // stays within 32-bit signed integer range.
    const MAX_ALGEBRA_SIZE = 4096;

    /**
     * Degenerate contact carrier Cl(n,n,1) with the light-cone Darboux
     * identification and a genuine contactomorphism checker.
     */
    class GAContactVersor {
        /**
         * @param {number} n - number of (q,p) pairs (base dimension). n=1 is the
         *   canonical α = ds − p dq case.
         */
        constructor(n = 1) {
            if (!Number.isInteger(n) || n < 1) {
                throw new Error('GAContactVersor: n must be a positive integer');
            }
            this.n = n;
            this.dim = 2 * n + 1;

            const size = 1 << this.dim;
            if (this.dim >= 31 || size > MAX_ALGEBRA_SIZE) {
                throw new Error(
                    `GAContactVersor: base dim n=${n} → Cl(${n},${n},1) has dim ${this.dim} ` +
                    `and ${this.dim >= 31 ? '2^dim' : size} blades, exceeding ` +
                    `MAX_ALGEBRA_SIZE=${MAX_ALGEBRA_SIZE}. The generic Algebra precompute is ` +
                    `O(4^dim); large carriers need a sparse exterior-algebra backend ` +
                    `(deferred to Phase 2).`
                );
            }

            // Signature: n positive (q), n negative (p), 1 null (s).
            const names = [];
            for (let a = 0; a < n; a++) names.push(`eq${a + 1}`);
            for (let a = 0; a < n; a++) names.push(`ep${a + 1}`);
            names.push('es');
            this.algebra = new Algebra(n, n, 1, names);

            // 1-indexed basis-vector positions.
            this._iq = a => 1 + a;             // eq_{a} (a: 0..n-1)
            this._ip = a => 1 + n + a;         // ep_{a}
            this._is = 2 * n + 1;              // es
        }

        // ---- basis vectors -------------------------------------------------

        /** q-generator eq_a (a: 0-based), e² = +1. */
        eq(a) { return this.algebra.e(this._iq(a)); }
        /** p-generator ep_a (a: 0-based), e² = −1. */
        ep(a) { return this.algebra.e(this._ip(a)); }
        /** null generator es (the ds / Reeb axis), e² = 0. */
        es() { return this.algebra.e(this._is); }

        /** Light-cone basis vector nq_a = eq_a − ep_a (⟨nq_a,nq_a⟩ = 0). */
        nq(a) { return this.eq(a).sub(this.ep(a)); }
        /** Light-cone basis vector np_a = (eq_a + ep_a)/2 (⟨nq_a,np_a⟩ = 1). */
        np(a) { return this.eq(a).add(this.ep(a)).scale(0.5); }

        // ---- forms ---------------------------------------------------------

        /**
         * The contact bivector B = dα = Σ_a eq_a ∧ ep_a. In this carrier it is a
         * genuine metric bivector: (eq_a ep_a)² = +1 per pair.
         */
        metricBivector() {
            let B = this.algebra.zero();
            for (let a = 0; a < this.n; a++) {
                B = B.add(this.eq(a).wedge(this.ep(a)));
            }
            return B;
        }

        // ---- point ↔ vector ------------------------------------------------

        /**
         * Grade-1 carrier vector for a Darboux point.
         * @param {number[]|number} q - q-coordinates (array, or scalar for n=1)
         * @param {number[]|number} p - p-coordinates
         * @param {number} s - fiber coordinate
         */
        pointVector(q, p, s) {
            const qa = Array.isArray(q) ? q : [q];
            const pa = Array.isArray(p) ? p : [p];
            let X = this.es().scale(s || 0);
            for (let a = 0; a < this.n; a++) {
                X = X.add(this.nq(a).scale(qa[a] || 0));
                X = X.add(this.np(a).scale(pa[a] || 0));
            }
            return X;
        }

        /**
         * Read Darboux coordinates (q, p, s) back out of a grade-1 vector.
         * Inverts the light-cone identification: for each pair, from the
         * (eq_a, ep_a) coefficients (c_q, c_p) we have q_a = (c_q − c_p)/2 and
         * p_a = c_q + c_p; s is the es coefficient.
         * @returns {{ q: number[], p: number[], s: number }}
         */
        readCoords(mv) {
            const q = [], p = [];
            for (let a = 0; a < this.n; a++) {
                const cq = mv.coeffs[1 << (this._iq(a) - 1)] || 0;
                const cp = mv.coeffs[1 << (this._ip(a) - 1)] || 0;
                q.push((cq - cp) / 2);
                p.push(cq + cp);
            }
            const s = mv.coeffs[1 << (this._is - 1)] || 0;
            return { q, p, s };
        }

        // ---- versors -------------------------------------------------------

        /**
         * Squeeze versor V = exp(Σ_a θ_a eq_a∧ep_a / 2). Because each
         * eq_a∧ep_a squares to +1 (hyperbolic), this is a genuine boost/squeeze
         * rotor: R R̃ = 1. In light-cone coordinates it acts as the symplectic
         * squeeze q_a ↦ e^{θ_a} q_a, p_a ↦ e^{−θ_a} p_a.
         * @param {number[]|number} thetas - rapidity per pair (scalar for n=1)
         */
        squeezeVersor(thetas) {
            const th = Array.isArray(thetas) ? thetas : [thetas];
            // The planes eq_a∧ep_a are orthogonal and their versors commute, but
            // their sum is a *non-simple* bivector whose square is not a scalar,
            // so the generic single-shot `Algebra.exp` (valid only for simple
            // bivectors) is inexact for n>1. Exponentiate each plane and compose.
            let V = this.algebra.scalar(1);
            for (let a = 0; a < this.n; a++) {
                const Ba = this.eq(a).wedge(this.ep(a)).scale((th[a] || 0) * 0.5);
                V = V.mul(this.algebra.exp(Ba));
            }
            return V;
        }

        /**
         * Null-generator "translator" versor exp(t es∧e / 2) where e is a
         * chosen q- or p-generator. These are the candidate versors for the Reeb
         * translation; es∧e is a *null* bivector so exp is parabolic (1 + t/2·B).
         * @param {number} a - pair index (0-based)
         * @param {'q'|'p'} dir - which generator to pair with es
         * @param {number} t - flow parameter
         */
        nullTranslatorVersor(a, dir, t) {
            if (dir !== 'q' && dir !== 'p') {
                throw new Error(
                    `GAContactVersor.nullTranslatorVersor: dir must be 'q' or 'p', got ${JSON.stringify(dir)}`
                );
            }
            if (!Number.isInteger(a) || a < 0 || a >= this.n) {
                throw new Error(
                    `GAContactVersor.nullTranslatorVersor: pair index a must be an integer in [0, ${this.n - 1}], got ${a}`
                );
            }
            const other = dir === 'q' ? this.eq(a) : this.ep(a);
            const B = this.es().wedge(other);
            return this.algebra.exp(B.scale(t / 2));
        }

        /**
         * Versor action v ↦ V v V⁻¹ using the *true* inverse V⁻¹ = Ṽ/(V Ṽ), not
         * merely the reverse Ṽ. This makes the operation scale-invariant — the
         * mathematical definition of the versor action — so that scaling V (e.g.
         * V ↦ 2V) leaves the result unchanged, matching v ↦ V v V⁻¹ exactly. Using
         * the bare reverse (Multivector.sandwich) only coincides when V Ṽ = 1.
         */
        apply(V, v) {
            return V.mul(v).mul(V.inverse());
        }

        /**
         * Induced linear map on Darboux coordinates as a (2n+1)² matrix
         * (row/col order [q_1..q_n, p_1..p_n, s]). Column j is the image of the
         * j-th Darboux basis point under the sandwich, read back via readCoords.
         */
        inducedMatrix(V) {
            const d = this.dim;
            const T = Array.from({ length: d }, () => new Array(d).fill(0));
            const basisPts = [];
            for (let a = 0; a < this.n; a++) basisPts.push(this.nq(a));      // q_a
            for (let a = 0; a < this.n; a++) basisPts.push(this.np(a));      // p_a
            basisPts.push(this.es());                                        // s
            for (let j = 0; j < d; j++) {
                const img = this.readCoords(this.apply(V, basisPts[j]));
                const col = [...img.q, ...img.p, img.s];
                for (let i = 0; i < d; i++) T[i][j] = col[i];
            }
            return T;
        }

        // ---- honest contactomorphism test ---------------------------------

        /**
         * Non-tautological test of the contactomorphism definition φ*α = f·α for
         * an affine coordinate map φ(X) = M·X + b, with α = ds − Σ_a p_a dq^a.
         *
         * The pullback acts as (φ*α)_X(v) = (Mv)_s − Σ_a (M X + b)_{p_a}(Mv)_{q_a}.
         * Requiring this to equal f·(v_s − Σ_a X_{p_a} v_{q_a}) for a *constant* f
         * is a genuine check against the definition (not satisfied by
         * construction). We sample (X, v) from a fixed, seeded pseudo-random
         * sequence, solve f from the ds-coefficient, and return the worst residual
         * over all coefficients and samples. The seeded sequence makes the check
         * fully deterministic and reproducible across runs and in CI.
         *
         * @param {number[][]} M - (2n+1)² linear part in Darboux order
         * @param {number[]} [b] - optional translation part (defaults to 0)
         * @returns {{ f: number|null, isContacto: boolean, residual: number }}
         */
        contactResidual(M, b = null) {
            const d = this.dim;
            const n = this.n;
            const trans = b || new Array(d).fill(0);
            const qIdx = a => a;            // q_a position in coord vector
            const pIdx = a => n + a;        // p_a position
            const sIdx = 2 * n;            // s position

            const matVec = (Mat, x) => {
                const out = new Array(d).fill(0);
                for (let i = 0; i < d; i++) {
                    let acc = 0;
                    for (let j = 0; j < d; j++) acc += Mat[i][j] * x[j];
                    out[i] = acc;
                }
                return out;
            };

            // Deterministic mulberry32 PRNG with a fixed seed, mapped to (−1, 1),
            // so the sampled (X, v) sequence is identical on every run. Replaces
            // Math.random(), which made the check flaky and could leave f=null.
            let seed = 0x9e3779b9 >>> 0;
            const rand = () => {
                seed = (seed + 0x6d2b79f5) >>> 0;
                let x = seed;
                x = Math.imul(x ^ (x >>> 15), x | 1);
                x ^= x + Math.imul(x ^ (x >>> 7), x | 61);
                return (((x ^ (x >>> 14)) >>> 0) / 4294967296) * 2 - 1;
            };

            let f0 = null;
            let maxRes = 0;
            for (let sample = 0; sample < 60; sample++) {
                const X = Array.from({ length: d }, () => rand());
                const v = Array.from({ length: d }, () => rand());
                const Mv = matVec(M, v);
                const phiX = matVec(M, X).map((val, i) => val + trans[i]);

                // LHS = (φ*α)(v)
                let lhs = Mv[sIdx];
                for (let a = 0; a < n; a++) lhs -= phiX[pIdx(a)] * Mv[qIdx(a)];

                // α(v) at X
                let alphaV = v[sIdx];
                for (let a = 0; a < n; a++) alphaV -= X[pIdx(a)] * v[qIdx(a)];

                // Solve f from a sample where α(v) is well away from 0, else use
                // the residual form directly.
                if (Math.abs(alphaV) > 0.1) {
                    const f = lhs / alphaV;
                    if (f0 === null) f0 = f;
                    maxRes = Math.max(maxRes, Math.abs(f - f0));
                }
                if (f0 !== null) {
                    maxRes = Math.max(maxRes, Math.abs(lhs - f0 * alphaV));
                }
            }
            return {
                f: f0,
                isContacto: f0 !== null && Math.abs(f0) > EPSILON && maxRes < 1e-9,
                residual: maxRes
            };
        }

        /**
         * Convenience: is the versor V a contactomorphism, and with what
         * conformal factor f? Runs `contactResidual` on the induced linear map.
         */
        versorContactResidual(V) {
            return this.contactResidual(this.inducedMatrix(V));
        }

        /** Determinant of a (2n+1)² matrix (isometry ⇒ |det| = 1). */
        static det(M) {
            const d = M.length;
            const A = M.map(r => r.slice());
            let det = 1;
            for (let col = 0; col < d; col++) {
                let piv = col;
                for (let r = col + 1; r < d; r++) {
                    if (Math.abs(A[r][col]) > Math.abs(A[piv][col])) piv = r;
                }
                if (Math.abs(A[piv][col]) < 1e-300) return 0;
                if (piv !== col) { [A[piv], A[col]] = [A[col], A[piv]]; det = -det; }
                det *= A[col][col];
                for (let r = col + 1; r < d; r++) {
                    const factor = A[r][col] / A[col][col];
                    for (let c = col; c < d; c++) A[r][c] -= factor * A[col][c];
                }
            }
            return det;
        }
    }

    const GAContactVersorModule = { GAContactVersor };

    if (typeof module !== 'undefined' && module.exports) {
        module.exports = GAContactVersorModule;
    } else if (typeof define === 'function' && define.amd) {
        define([], () => GAContactVersorModule);
    } else {
        global.GAContactVersor = GAContactVersorModule;
    }

})(typeof window !== 'undefined' ? window : global);
