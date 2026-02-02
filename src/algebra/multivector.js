/**
 * Geometric Algebra Module for Contact Thermodynamics
 * 
 * Provides Clifford algebra operations for enhanced geometric computations.
 * Based on patterns from the geometric-algebra skill.
 * 
 * Key classes:
 * - Algebra: Defines Cl(p,q,r) with precomputed product tables
 * - Multivector: Full GA operations (geometric/wedge/inner products, duality)
 * - ContactAlgebra: Factory for contact geometry algebras
 * 
 * @module multivector
 * @license MIT
 */

(function(global) {
    'use strict';

    const EPSILON = 1e-10;
    const abs = Math.abs;
    const sqrt = Math.sqrt;
    const sin = Math.sin;
    const cos = Math.cos;
    const cosh = Math.cosh;
    const sinh = Math.sinh;
    const atan2 = Math.atan2;

    function isZero(x) { return abs(x) < EPSILON; }

    // ============================================================================
    // ALGEBRA CLASS - Defines Cl(p,q,r)
    // ============================================================================

    /**
     * Clifford Algebra Cl(p,q,r)
     * 
     * @param {number} p - Positive signature dimensions (e²=+1)
     * @param {number} q - Negative signature dimensions (e²=-1)
     * @param {number} r - Zero signature dimensions (e²=0, degenerate)
     * @param {string[]} basisNames - Optional custom names for basis vectors
     */
    class Algebra {
        constructor(p, q, r = 0, basisNames = null) {
            this.p = p;
            this.q = q;
            this.r = r;
            this.n = p + q + r;  // total dimensions
            this.size = 1 << this.n;  // 2^n basis blades

            // Metric signature: first p are +1, next q are -1, last r are 0
            this.signature = [];
            for (let i = 0; i < p; i++) this.signature.push(1);
            for (let i = 0; i < q; i++) this.signature.push(-1);
            for (let i = 0; i < r; i++) this.signature.push(0);

            // Basis vector names
            if (basisNames) {
                this.basisNames = basisNames;
            } else {
                this.basisNames = [];
                for (let i = 0; i < this.n; i++) {
                    this.basisNames.push(`e${i + 1}`);
                }
            }

            // Precompute product tables
            this._precomputeProducts();
        }

        /**
         * Popcount - count number of 1 bits (grade of blade)
         */
        _popcount(n) {
            let count = 0;
            while (n) { count += n & 1; n >>= 1; }
            return count;
        }

        /**
         * Convert bitmap to blade name
         */
        _bitmapToName(bitmap) {
            if (bitmap === 0) return '1';
            let name = '';
            for (let i = 0; i < this.n; i++) {
                if (bitmap & (1 << i)) {
                    name += this.basisNames[i];
                }
            }
            return name;
        }

        /**
         * Precompute geometric product table and signs
         */
        _precomputeProducts() {
            // productSign[a][b] = sign of geometric product of blades a and b
            // productResult[a][b] = resulting blade bitmap
            this.productSign = [];
            this.productResult = [];

            for (let a = 0; a < this.size; a++) {
                this.productSign[a] = [];
                this.productResult[a] = [];
                for (let b = 0; b < this.size; b++) {
                    const { sign, result } = this._computeProduct(a, b);
                    this.productSign[a][b] = sign;
                    this.productResult[a][b] = result;
                }
            }
        }

        /**
         * Compute product of two basis blades
         * Returns { sign, result } where result is the bitmap of the product blade
         */
        _computeProduct(a, b) {
            let sign = 1;
            let result = a;

            // For each basis vector in b, move it through result
            for (let i = 0; i < this.n; i++) {
                if (b & (1 << i)) {
                    // Count swaps needed to move e_i to its position
                    let swaps = 0;
                    for (let j = i + 1; j < this.n; j++) {
                        if (result & (1 << j)) swaps++;
                    }
                    if (swaps % 2 === 1) sign *= -1;

                    if (result & (1 << i)) {
                        // e_i * e_i = signature[i]
                        sign *= this.signature[i];
                        result ^= (1 << i); // Remove the basis vector
                    } else {
                        result ^= (1 << i); // Add the basis vector
                    }
                }
            }

            return { sign, result };
        }

        // ----- Factory methods -----

        /**
         * Create scalar multivector
         */
        scalar(s) {
            return new Multivector(this, { 0: s });
        }

        /**
         * Create zero multivector
         */
        zero() {
            return new Multivector(this, {});
        }

        /**
         * Create basis vector e_i (1-indexed)
         */
        e(i) {
            const bitmap = 1 << (i - 1);
            return new Multivector(this, { [bitmap]: 1 });
        }

        /**
         * Create vector from components
         */
        vector(components) {
            const coeffs = {};
            for (let i = 0; i < Math.min(components.length, this.n); i++) {
                if (!isZero(components[i])) {
                    coeffs[1 << i] = components[i];
                }
            }
            return new Multivector(this, coeffs);
        }

        /**
         * Create bivector from component array
         * For 3D, components = [e12, e13, e23] or [e12, e23, e31]
         */
        bivector(components) {
            const coeffs = {};
            let idx = 0;
            for (let i = 0; i < this.n; i++) {
                for (let j = i + 1; j < this.n; j++) {
                    if (idx < components.length && !isZero(components[idx])) {
                        const bitmap = (1 << i) | (1 << j);
                        coeffs[bitmap] = components[idx];
                    }
                    idx++;
                }
            }
            return new Multivector(this, coeffs);
        }

        /**
         * Create pseudoscalar (highest grade element)
         */
        pseudoscalar() {
            const bitmap = this.size - 1;
            return new Multivector(this, { [bitmap]: 1 });
        }

        // ----- Product operations -----

        /**
         * Geometric product of two multivectors
         */
        gp(a, b) {
            const coeffs = {};
            for (const [ka, va] of Object.entries(a.coeffs)) {
                for (const [kb, vb] of Object.entries(b.coeffs)) {
                    const ia = parseInt(ka);
                    const ib = parseInt(kb);
                    const sign = this.productSign[ia][ib];
                    const result = this.productResult[ia][ib];
                    if (sign !== 0) {
                        coeffs[result] = (coeffs[result] || 0) + sign * va * vb;
                    }
                }
            }
            // Clean up near-zero entries
            for (const k of Object.keys(coeffs)) {
                if (isZero(coeffs[k])) delete coeffs[k];
            }
            return new Multivector(this, coeffs);
        }

        /**
         * Outer (wedge) product
         */
        op(a, b) {
            const coeffs = {};
            for (const [ka, va] of Object.entries(a.coeffs)) {
                for (const [kb, vb] of Object.entries(b.coeffs)) {
                    const ia = parseInt(ka);
                    const ib = parseInt(kb);
                    // Only include if blades don't share any basis vectors
                    if ((ia & ib) === 0) {
                        const sign = this.productSign[ia][ib];
                        const result = this.productResult[ia][ib];
                        coeffs[result] = (coeffs[result] || 0) + sign * va * vb;
                    }
                }
            }
            for (const k of Object.keys(coeffs)) {
                if (isZero(coeffs[k])) delete coeffs[k];
            }
            return new Multivector(this, coeffs);
        }

        /**
         * Inner (left contraction) product
         */
        ip(a, b) {
            const coeffs = {};
            for (const [ka, va] of Object.entries(a.coeffs)) {
                for (const [kb, vb] of Object.entries(b.coeffs)) {
                    const ia = parseInt(ka);
                    const ib = parseInt(kb);
                    const gradeA = this._popcount(ia);
                    const gradeB = this._popcount(ib);
                    const sign = this.productSign[ia][ib];
                    const result = this.productResult[ia][ib];
                    const gradeResult = this._popcount(result);
                    // Left contraction: grade drops by grade(a)
                    if (gradeResult === gradeB - gradeA && gradeA <= gradeB) {
                        coeffs[result] = (coeffs[result] || 0) + sign * va * vb;
                    }
                }
            }
            for (const k of Object.keys(coeffs)) {
                if (isZero(coeffs[k])) delete coeffs[k];
            }
            return new Multivector(this, coeffs);
        }

        /**
         * Dual: multiply by pseudoscalar inverse
         */
        dual(a) {
            const I = this.pseudoscalar();
            return this.gp(a, I.inverse());
        }

        /**
         * Undual: multiply by pseudoscalar
         */
        undual(a) {
            const I = this.pseudoscalar();
            return this.gp(a, I);
        }

        /**
         * Regressive (vee) product: a ∨ b = (a* ∧ b*)*
         */
        vee(a, b) {
            return this.undual(this.op(this.dual(a), this.dual(b)));
        }

        /**
         * Create rotor from bivector and angle
         * R = cos(θ/2) + sin(θ/2) * B̂
         */
        rotor(bivector, angle) {
            const halfAngle = angle / 2;
            const Bsq = bivector.mul(bivector).scalar();
            
            if (Bsq < -EPSILON) {
                // Elliptic case (rotation): B² < 0
                const normB = sqrt(-Bsq);
                const Bhat = bivector.scale(1 / normB);
                return this.scalar(cos(halfAngle)).add(Bhat.scale(sin(halfAngle)));
            } else if (Bsq > EPSILON) {
                // Hyperbolic case (boost): B² > 0
                const normB = sqrt(Bsq);
                const Bhat = bivector.scale(1 / normB);
                return this.scalar(cosh(halfAngle)).add(Bhat.scale(sinh(halfAngle)));
            } else {
                // Parabolic case (translation): B² = 0
                return this.scalar(1).add(bivector.scale(halfAngle));
            }
        }

        /**
         * Exponential of bivector (general case with classification)
         * exp(B) adapts formula based on B²
         */
        exp(bivector) {
            const Bsq = bivector.mul(bivector).scalar();
            
            if (Bsq < -EPSILON) {
                // Elliptic (B² < 0): exp(B) = cos(|B|) + sin(|B|)/|B| * B
                const normB = sqrt(-Bsq);
                return this.scalar(cos(normB)).add(bivector.scale(sin(normB) / normB));
            } else if (Bsq > EPSILON) {
                // Hyperbolic (B² > 0): exp(B) = cosh(|B|) + sinh(|B|)/|B| * B
                const normB = sqrt(Bsq);
                return this.scalar(cosh(normB)).add(bivector.scale(sinh(normB) / normB));
            } else {
                // Parabolic (B² = 0): exp(B) = 1 + B
                return this.scalar(1).add(bivector);
            }
        }

        /**
         * Logarithm of rotor (inverse of exp)
         */
        log(rotor) {
            const s = rotor.grade(0).scalar();
            const B = rotor.grade(2);
            const Bsq = B.mul(B).scalar();
            
            if (isZero(Bsq)) {
                // Near identity
                return B;
            }
            
            if (Bsq < 0) {
                // Elliptic
                const normB = sqrt(-Bsq);
                const angle = atan2(normB, s);
                return B.scale(angle / normB);
            } else {
                // Hyperbolic 
                const normB = sqrt(Bsq);
                // acosh implementation
                const angle = Math.log(s + sqrt(s * s - 1));
                return B.scale(angle / normB);
            }
        }

        /**
         * Get info string
         */
        toString() {
            return `Cl(${this.p},${this.q},${this.r}) with ${this.size} basis blades`;
        }
    }

    // ============================================================================
    // MULTIVECTOR CLASS
    // ============================================================================

    /**
     * Multivector in a Clifford algebra
     * Stores coefficients in a sparse object indexed by basis blade bitmap
     */
    class Multivector {
        constructor(algebra, coeffs = {}) {
            this.algebra = algebra;
            this.coeffs = {};
            for (const [k, v] of Object.entries(coeffs)) {
                if (!isZero(v)) {
                    this.coeffs[k] = v;
                }
            }
        }

        clone() {
            return new Multivector(this.algebra, { ...this.coeffs });
        }

        /**
         * Grade extraction - returns multivector with only grade-g components
         */
        grade(g) {
            const coeffs = {};
            for (const [k, v] of Object.entries(this.coeffs)) {
                if (this.algebra._popcount(parseInt(k)) === g) {
                    coeffs[k] = v;
                }
            }
            return new Multivector(this.algebra, coeffs);
        }

        /**
         * Scalar part (grade 0)
         */
        scalar() {
            return this.coeffs[0] || 0;
        }

        /**
         * Reverse (reversion): reverses order of basis vectors in each blade
         * ~(e1∧e2∧...∧ek) = (-1)^(k(k-1)/2) * (e1∧e2∧...∧ek)
         */
        reverse() {
            const coeffs = {};
            for (const [k, v] of Object.entries(this.coeffs)) {
                const grade = this.algebra._popcount(parseInt(k));
                const sign = ((grade * (grade - 1) / 2) % 2 === 0) ? 1 : -1;
                coeffs[k] = sign * v;
            }
            return new Multivector(this.algebra, coeffs);
        }

        /**
         * Grade involution (main involution)
         * Negates odd-grade components
         */
        involute() {
            const coeffs = {};
            for (const [k, v] of Object.entries(this.coeffs)) {
                const grade = this.algebra._popcount(parseInt(k));
                coeffs[k] = (grade % 2 === 0) ? v : -v;
            }
            return new Multivector(this.algebra, coeffs);
        }

        /**
         * Clifford conjugate = involute then reverse
         */
        conjugate() {
            return this.involute().reverse();
        }

        // ----- Arithmetic operations -----

        add(other) {
            const coeffs = { ...this.coeffs };
            for (const [k, v] of Object.entries(other.coeffs)) {
                coeffs[k] = (coeffs[k] || 0) + v;
                if (isZero(coeffs[k])) delete coeffs[k];
            }
            return new Multivector(this.algebra, coeffs);
        }

        sub(other) {
            const coeffs = { ...this.coeffs };
            for (const [k, v] of Object.entries(other.coeffs)) {
                coeffs[k] = (coeffs[k] || 0) - v;
                if (isZero(coeffs[k])) delete coeffs[k];
            }
            return new Multivector(this.algebra, coeffs);
        }

        scale(s) {
            const coeffs = {};
            for (const [k, v] of Object.entries(this.coeffs)) {
                coeffs[k] = s * v;
            }
            return new Multivector(this.algebra, coeffs);
        }

        // ----- Products -----

        mul(other) {
            return this.algebra.gp(this, other);
        }

        wedge(other) {
            return this.algebra.op(this, other);
        }

        dot(other) {
            return this.algebra.ip(this, other);
        }

        vee(other) {
            return this.algebra.vee(this, other);
        }

        /**
         * Commutator product: A × B = ½(AB - BA)
         * This is the "rotation" action - for bivectors, [ω, v] rotates v
         * Essential for curvature computation: Ω = dω + ω ∧ ω uses [ω_a, ω_b]
         */
        commutator(other) {
            const ab = this.mul(other);
            const ba = other.mul(this);
            return ab.sub(ba).scale(0.5);
        }

        /**
         * Anticommutator product: A ○ B = ½(AB + BA)
         */
        anticommutator(other) {
            const ab = this.mul(other);
            const ba = other.mul(this);
            return ab.add(ba).scale(0.5);
        }

        /**
         * Sandwich product: this * other * ~this
         */
        sandwich(other) {
            return this.mul(other).mul(this.reverse());
        }

        // ----- Duality -----

        dual() {
            return this.algebra.dual(this);
        }

        undual() {
            return this.algebra.undual(this);
        }

        // ----- Norms -----

        /**
         * Norm squared: <A ~A>_0
         */
        normSq() {
            const prod = this.mul(this.reverse());
            return abs(prod.scalar());
        }

        norm() {
            return sqrt(this.normSq());
        }

        normalized() {
            const n = this.norm();
            return isZero(n) ? this.clone() : this.scale(1 / n);
        }

        /**
         * Inverse (for versors): ~A / (A * ~A)
         */
        inverse() {
            const rev = this.reverse();
            const normSq = this.mul(rev).scalar();
            if (isZero(normSq)) {
                throw new Error('Cannot invert: norm squared is zero');
            }
            return rev.scale(1 / normSq);
        }

        isZero() {
            return Object.keys(this.coeffs).length === 0;
        }

        /**
         * String representation
         */
        toString() {
            if (this.isZero()) return '0';
            const terms = [];
            for (const [k, v] of Object.entries(this.coeffs)) {
                const name = this.algebra._bitmapToName(parseInt(k));
                if (name === '1') {
                    terms.push(v.toFixed(4));
                } else {
                    terms.push(`${v.toFixed(4)}*${name}`);
                }
            }
            return terms.join(' + ');
        }
    }

    // ============================================================================
    // CONTACT ALGEBRA FACTORY
    // ============================================================================

    /**
     * Factory for creating algebras suitable for contact geometry
     */
    const ContactAlgebra = {
        /**
         * Create algebra for contact manifold of base dimension n
         * Uses Cl(n+1, n, 0) which includes time-like and space-like dimensions
         * 
         * For relativistic contact geometry, use spacetime signature
         */
        create(baseDim, type = 'euclidean') {
            if (type === 'euclidean') {
                // Pure Euclidean: Cl(n, 0, 0)
                return new Algebra(baseDim, 0, 0);
            } else if (type === 'spacetime') {
                // Spacetime: Cl(1, baseDim, 0) with signature (+,-,-,-)
                return new Algebra(1, baseDim, 0, ['e0', 'e1', 'e2', 'e3'].slice(0, baseDim + 1));
            } else if (type === 'pga') {
                // Projective: Cl(baseDim, 0, 1)
                return new Algebra(baseDim, 0, 1);
            } else if (type === 'cga') {
                // Conformal: Cl(baseDim+1, 1, 0)
                return new Algebra(baseDim + 1, 1, 0);
            }
            return new Algebra(baseDim, 0, 0);
        },

        /**
         * Spacetime algebra for relativistic contact mechanics
         * Cl(1,3,0) - the standard spacetime algebra
         */
        spacetime() {
            return new Algebra(1, 3, 0, ['e0', 'e1', 'e2', 'e3']);
        },

        /**
         * 3D Euclidean algebra for spatial computations
         */
        euclidean3D() {
            return new Algebra(3, 0, 0, ['e1', 'e2', 'e3']);
        }
    };

    // ============================================================================
    // BIVECTOR CLASSIFICATION - Key insight from GA skill
    // ============================================================================

    /**
     * Classify bivector type based on B²
     * This determines the correct exponential formula
     * 
     * @param {Multivector} B - Bivector to classify
     * @returns {Object} { type, Bsq } where type is 'elliptic'|'parabolic'|'hyperbolic'
     */
    function classifyBivector(B) {
        const Bsq = B.mul(B).scalar();
        if (Bsq < -EPSILON) {
            return { type: 'elliptic', Bsq, description: 'rotation' };
        } else if (Bsq > EPSILON) {
            return { type: 'hyperbolic', Bsq, description: 'boost' };
        } else {
            return { type: 'parabolic', Bsq, description: 'translation' };
        }
    }

    // ============================================================================
    // EXPORTS
    // ============================================================================

    const GA = {
        Algebra,
        Multivector,
        ContactAlgebra,
        classifyBivector,
        
        // Common algebras
        cl3: () => new Algebra(3, 0, 0),      // Cl(3) - 3D rotors
        pga3: () => new Algebra(3, 0, 1),     // PGA 3D
        cga3: () => new Algebra(4, 1, 0),     // CGA 3D
        sta: () => new Algebra(1, 3, 0),      // Spacetime algebra

        EPSILON
    };

    // Export for different module systems
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = GA;
    }
    if (typeof global !== 'undefined') {
        global.GA = GA;
    }
    if (typeof window !== 'undefined') {
        window.GA = GA;
    }

})(typeof globalThis !== 'undefined' ? globalThis : this);
