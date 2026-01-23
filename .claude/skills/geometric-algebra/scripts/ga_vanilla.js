/**
 * Geometric Algebra - Vanilla JavaScript Implementation
 * 
 * Supports:
 * - Cl(3): 3D Euclidean rotors
 * - PGA 2D (2,0,1): 2D Projective Geometric Algebra
 * - PGA 3D (3,0,1): 3D Projective Geometric Algebra  
 * - CGA 2D (3,1,0): 2D Conformal Geometric Algebra
 * - CGA 3D (4,1,0): 3D Conformal Geometric Algebra
 * 
 * Based on:
 * - "From Invariant Decomposition to Spinors" (Roelfs & De Keninck)
 * - "Graded Symmetry Groups: Plane and Simple" (Roelfs & De Keninck)
 * 
 * @license MIT
 */

(function(global) {
    'use strict';

    // ============================================================================
    // UTILITY FUNCTIONS
    // ============================================================================

    const EPSILON = 1e-10;
    const abs = Math.abs;
    const sqrt = Math.sqrt;
    const sin = Math.sin;
    const cos = Math.cos;
    const atan2 = Math.atan2;
    const acos = Math.acos;
    const exp = Math.exp;
    const log = Math.log;

    function isZero(x) { return abs(x) < EPSILON; }
    
    function normalize(arr) {
        const norm = sqrt(arr.reduce((s, x) => s + x*x, 0));
        return norm < EPSILON ? arr : arr.map(x => x / norm);
    }

    // ============================================================================
    // MULTIVECTOR BASE CLASS
    // ============================================================================

    /**
     * Generic Multivector for any Clifford algebra Cl(p,q,r)
     * Stores coefficients in a flat array indexed by basis blade bitmap
     */
    class Multivector {
        constructor(algebra, coeffs = {}) {
            this.algebra = algebra;
            this.coeffs = new Float64Array(algebra.size);
            for (const [blade, val] of Object.entries(coeffs)) {
                const idx = algebra.bladeIndex[blade];
                if (idx !== undefined) this.coeffs[idx] = val;
            }
        }

        clone() {
            const mv = new Multivector(this.algebra);
            mv.coeffs.set(this.coeffs);
            return mv;
        }

        // Grade extraction
        grade(g) {
            const mv = new Multivector(this.algebra);
            for (let i = 0; i < this.algebra.size; i++) {
                if (this.algebra.grades[i] === g) mv.coeffs[i] = this.coeffs[i];
            }
            return mv;
        }

        // Scalar part
        get scalar() { return this.coeffs[0]; }

        // Reverse (reversion): reverses order of basis vectors in each blade
        reverse() {
            const mv = new Multivector(this.algebra);
            for (let i = 0; i < this.algebra.size; i++) {
                const g = this.algebra.grades[i];
                const sign = (g * (g - 1) / 2) % 2 === 0 ? 1 : -1;
                mv.coeffs[i] = sign * this.coeffs[i];
            }
            return mv;
        }

        // Conjugate (Clifford conjugate): grade involution + reversion
        conjugate() {
            const mv = new Multivector(this.algebra);
            for (let i = 0; i < this.algebra.size; i++) {
                const g = this.algebra.grades[i];
                const sign = (g * (g + 1) / 2) % 2 === 0 ? 1 : -1;
                mv.coeffs[i] = sign * this.coeffs[i];
            }
            return mv;
        }

        // Grade involution (main involution)
        involute() {
            const mv = new Multivector(this.algebra);
            for (let i = 0; i < this.algebra.size; i++) {
                const sign = this.algebra.grades[i] % 2 === 0 ? 1 : -1;
                mv.coeffs[i] = sign * this.coeffs[i];
            }
            return mv;
        }

        // Addition
        add(other) {
            const mv = this.clone();
            for (let i = 0; i < this.algebra.size; i++) {
                mv.coeffs[i] += other.coeffs[i];
            }
            return mv;
        }

        // Subtraction
        sub(other) {
            const mv = this.clone();
            for (let i = 0; i < this.algebra.size; i++) {
                mv.coeffs[i] -= other.coeffs[i];
            }
            return mv;
        }

        // Scalar multiplication
        scale(s) {
            const mv = this.clone();
            for (let i = 0; i < this.algebra.size; i++) {
                mv.coeffs[i] *= s;
            }
            return mv;
        }

        // Geometric product
        mul(other) {
            return this.algebra.geometricProduct(this, other);
        }

        // Outer (wedge) product
        wedge(other) {
            return this.algebra.outerProduct(this, other);
        }

        // Inner (dot) product (left contraction)
        dot(other) {
            return this.algebra.innerProduct(this, other);
        }

        // Regressive (vee) product: a ∨ b = (a* ∧ b*)*
        vee(other) {
            return this.dual().wedge(other.dual()).dual();
        }

        // Dual (Poincaré duality)
        dual() {
            return this.algebra.dual(this);
        }

        // Undual
        undual() {
            return this.algebra.undual(this);
        }

        // Sandwich product: this * other * this.reverse()
        sandwich(other) {
            return this.mul(other).mul(this.reverse());
        }

        // Norm squared
        normSq() {
            const rev = this.reverse();
            return this.mul(rev).scalar;
        }

        // Norm
        norm() {
            const ns = this.normSq();
            return ns >= 0 ? sqrt(ns) : sqrt(-ns);
        }

        // Normalized
        normalized() {
            const n = this.norm();
            return n < EPSILON ? this.clone() : this.scale(1 / n);
        }

        // Inverse (for versors): ~A / (A * ~A)
        inverse() {
            const rev = this.reverse();
            const ns = this.mul(rev).scalar;
            return isZero(ns) ? null : rev.scale(1 / ns);
        }

        // Check if zero
        isZero() {
            return this.coeffs.every(c => isZero(c));
        }

        // String representation
        toString() {
            const terms = [];
            for (let i = 0; i < this.algebra.size; i++) {
                if (!isZero(this.coeffs[i])) {
                    const blade = this.algebra.bladeNames[i];
                    if (blade === '1') {
                        terms.push(this.coeffs[i].toFixed(4));
                    } else {
                        terms.push(`${this.coeffs[i].toFixed(4)}${blade}`);
                    }
                }
            }
            return terms.length ? terms.join(' + ') : '0';
        }
    }

    // ============================================================================
    // ALGEBRA CLASS - Defines Cl(p,q,r)
    // ============================================================================

    class Algebra {
        /**
         * Create Clifford algebra Cl(p,q,r)
         * @param {number} p - positive signature dimensions
         * @param {number} q - negative signature dimensions  
         * @param {number} r - zero signature dimensions (degenerate)
         * @param {string[]} basisNames - optional custom names for basis vectors
         */
        constructor(p, q, r = 0, basisNames = null) {
            this.p = p;
            this.q = q;
            this.r = r;
            this.n = p + q + r;  // total dimensions
            this.size = 1 << this.n;  // 2^n basis blades

            // Metric: e_i^2 = +1 for i < p, -1 for p <= i < p+q, 0 for i >= p+q
            this.metric = [];
            for (let i = 0; i < this.n; i++) {
                if (i < p) this.metric.push(1);
                else if (i < p + q) this.metric.push(-1);
                else this.metric.push(0);
            }

            // Basis vector names
            this.basisNames = basisNames || 
                Array.from({length: this.n}, (_, i) => `e${i + 1}`);

            // Generate blade names and indices
            this.bladeNames = [];
            this.bladeIndex = {};
            this.grades = [];

            for (let i = 0; i < this.size; i++) {
                const name = this._bitmapToName(i);
                this.bladeNames.push(name);
                this.bladeIndex[name] = i;
                this.grades.push(this._popcount(i));
            }

            // Precompute geometric product signs
            this._precomputeProducts();

            // Pseudoscalar
            this.I = this.blade((1 << this.n) - 1);
        }

        _popcount(n) {
            let count = 0;
            while (n) { count += n & 1; n >>= 1; }
            return count;
        }

        _bitmapToName(bitmap) {
            if (bitmap === 0) return '1';
            const indices = [];
            for (let i = 0; i < this.n; i++) {
                if (bitmap & (1 << i)) indices.push(this.basisNames[i]);
            }
            return indices.join('');
        }

        _precomputeProducts() {
            // Precompute sign and result bitmap for geometric product
            this.gpSign = new Int8Array(this.size * this.size);
            this.gpResult = new Uint16Array(this.size * this.size);

            for (let a = 0; a < this.size; a++) {
                for (let b = 0; b < this.size; b++) {
                    const [sign, result] = this._bladeProduct(a, b);
                    const idx = a * this.size + b;
                    this.gpSign[idx] = sign;
                    this.gpResult[idx] = result;
                }
            }
        }

        _bladeProduct(a, b) {
            // Compute e_A * e_B where A and B are bitmaps
            let sign = 1;
            let result = a ^ b;  // XOR gives resulting blade

            // Count swaps needed to reorder
            let swaps = 0;
            let temp = a;
            for (let i = 0; i < this.n; i++) {
                if (b & (1 << i)) {
                    // Count bits in a above position i
                    for (let j = i + 1; j < this.n; j++) {
                        if (temp & (1 << j)) swaps++;
                    }
                }
                temp &= ~(1 << i);  // Clear processed bit
            }
            if (swaps % 2) sign = -sign;

            // Apply metric for contracted basis vectors
            const contracted = a & b;
            for (let i = 0; i < this.n; i++) {
                if (contracted & (1 << i)) {
                    sign *= this.metric[i];
                }
            }

            return [sign, result];
        }

        // Create multivector from single blade
        blade(bitmap, coeff = 1) {
            const mv = new Multivector(this);
            mv.coeffs[bitmap] = coeff;
            return mv;
        }

        // Create scalar
        scalar(s) {
            return this.blade(0, s);
        }

        // Create vector from components
        vector(...components) {
            const mv = new Multivector(this);
            for (let i = 0; i < Math.min(components.length, this.n); i++) {
                mv.coeffs[1 << i] = components[i];
            }
            return mv;
        }

        // Geometric product
        geometricProduct(a, b) {
            const result = new Multivector(this);
            for (let i = 0; i < this.size; i++) {
                if (isZero(a.coeffs[i])) continue;
                for (let j = 0; j < this.size; j++) {
                    if (isZero(b.coeffs[j])) continue;
                    const idx = i * this.size + j;
                    const k = this.gpResult[idx];
                    result.coeffs[k] += this.gpSign[idx] * a.coeffs[i] * b.coeffs[j];
                }
            }
            return result;
        }

        // Outer product
        outerProduct(a, b) {
            const result = new Multivector(this);
            for (let i = 0; i < this.size; i++) {
                if (isZero(a.coeffs[i])) continue;
                for (let j = 0; j < this.size; j++) {
                    if (isZero(b.coeffs[j])) continue;
                    // Only if blades don't share basis vectors
                    if ((i & j) === 0) {
                        const idx = i * this.size + j;
                        const k = this.gpResult[idx];
                        result.coeffs[k] += this.gpSign[idx] * a.coeffs[i] * b.coeffs[j];
                    }
                }
            }
            return result;
        }

        // Left contraction (inner product)
        innerProduct(a, b) {
            const result = new Multivector(this);
            for (let i = 0; i < this.size; i++) {
                if (isZero(a.coeffs[i])) continue;
                const gradeA = this.grades[i];
                for (let j = 0; j < this.size; j++) {
                    if (isZero(b.coeffs[j])) continue;
                    const gradeB = this.grades[j];
                    const idx = i * this.size + j;
                    const k = this.gpResult[idx];
                    // Left contraction: grade(result) = grade(b) - grade(a)
                    if (this.grades[k] === gradeB - gradeA && gradeA <= gradeB) {
                        result.coeffs[k] += this.gpSign[idx] * a.coeffs[i] * b.coeffs[j];
                    }
                }
            }
            return result;
        }

        // Dual: A* = A ⌋ I^(-1) or A * I^(-1) for pseudoscalar
        dual(a) {
            const Iinv = this.I.inverse();
            return Iinv ? a.mul(Iinv) : a.mul(this.I.reverse());
        }

        // Undual: A = A* * I
        undual(a) {
            return a.mul(this.I);
        }

        // Create rotor from bivector and angle
        // R = cos(θ/2) + sin(θ/2) * B̂
        rotor(bivector, angle) {
            const B = bivector.normalized();
            const half = angle / 2;
            return this.scalar(cos(half)).add(B.scale(sin(half)));
        }

        // Exponential of bivector: exp(B) = cos(|B|) + sin(|B|) * B̂
        exp(bivector) {
            const norm = bivector.norm();
            if (isZero(norm)) return this.scalar(1);
            return this.scalar(cos(norm)).add(bivector.scale(sin(norm) / norm));
        }

        // Logarithm of rotor
        log(rotor) {
            const s = rotor.scalar;
            const B = rotor.grade(2);
            const Bnorm = B.norm();
            if (isZero(Bnorm)) return this.scalar(0).grade(2);
            const angle = atan2(Bnorm, s);
            return B.scale(angle / Bnorm);
        }
    }

    // ============================================================================
    // PGA 2D - Projective Geometric Algebra Cl(2,0,1)
    // ============================================================================

    class PGA2D extends Algebra {
        constructor() {
            super(2, 0, 1, ['e1', 'e2', 'e0']);
            
            // Blade shortcuts
            this.e0 = this.blade(0b100);  // e0 (degenerate)
            this.e1 = this.blade(0b001);
            this.e2 = this.blade(0b010);
            this.e01 = this.blade(0b101);
            this.e02 = this.blade(0b110);
            this.e12 = this.blade(0b011);
            this.e012 = this.blade(0b111);  // Pseudoscalar
        }

        // Point at (x, y): P = e12 + x*e02 - y*e01
        point(x, y) {
            return this.e12.add(this.e02.scale(x)).sub(this.e01.scale(y));
        }

        // Ideal point (direction): D = x*e02 - y*e01
        direction(x, y) {
            return this.e02.scale(x).sub(this.e01.scale(y));
        }

        // Line ax + by + c = 0: L = a*e1 + b*e2 + c*e0
        line(a, b, c) {
            return this.e1.scale(a).add(this.e2.scale(b)).add(this.e0.scale(c));
        }

        // Extract (x, y) from point
        pointCoords(p) {
            const w = p.coeffs[this.bladeIndex['e1e2']] || 1;
            const x = p.coeffs[this.bladeIndex['e0e2']] / w;
            const y = -p.coeffs[this.bladeIndex['e0e1']] / w;
            return [x, y];
        }

        // Join: line through two points
        join(p1, p2) {
            return p1.vee(p2);
        }

        // Meet: intersection of two lines
        meet(l1, l2) {
            return l1.wedge(l2);
        }

        // Translator by (dx, dy)
        translator(dx, dy) {
            // T = 1 + (dx*e01 + dy*e02)/2
            return this.scalar(1)
                .add(this.e01.scale(dx / 2))
                .add(this.e02.scale(dy / 2));
        }

        // Rotor around origin by angle
        rotorOrigin(angle) {
            return this.scalar(cos(angle / 2))
                .add(this.e12.scale(sin(angle / 2)));
        }

        // Rotor around point by angle
        rotorPoint(px, py, angle) {
            const T = this.translator(px, py);
            const R = this.rotorOrigin(angle);
            const Tinv = this.translator(-px, -py);
            return T.mul(R).mul(Tinv);
        }

        // Motor (combined rotation + translation)
        motor(angle, dx, dy) {
            const R = this.rotorOrigin(angle);
            const T = this.translator(dx, dy);
            return T.mul(R);
        }

        // Reflect point in line
        reflect(point, line) {
            const L = line.normalized();
            return L.sandwich(point);
        }
    }

    // ============================================================================
    // PGA 3D - Projective Geometric Algebra Cl(3,0,1)
    // ============================================================================

    class PGA3D extends Algebra {
        constructor() {
            super(3, 0, 1, ['e1', 'e2', 'e3', 'e0']);

            // Basis blades
            this.e0 = this.blade(0b1000);
            this.e1 = this.blade(0b0001);
            this.e2 = this.blade(0b0010);
            this.e3 = this.blade(0b0100);
            
            // Bivectors
            this.e01 = this.blade(0b1001);
            this.e02 = this.blade(0b1010);
            this.e03 = this.blade(0b1100);
            this.e12 = this.blade(0b0011);
            this.e31 = this.blade(0b0101);
            this.e23 = this.blade(0b0110);
            
            // Trivectors
            this.e012 = this.blade(0b1011);
            this.e013 = this.blade(0b1101);
            this.e023 = this.blade(0b1110);
            this.e123 = this.blade(0b0111);
            
            // Pseudoscalar
            this.e0123 = this.blade(0b1111);
        }

        // Point at (x, y, z): P = e123 - x*e023 + y*e013 - z*e012
        point(x, y, z) {
            return this.e123
                .sub(this.e023.scale(x))
                .add(this.e013.scale(y))
                .sub(this.e012.scale(z));
        }

        // Ideal point (direction)
        direction(x, y, z) {
            return this.e023.scale(-x)
                .add(this.e013.scale(y))
                .sub(this.e012.scale(z));
        }

        // Plane ax + by + cz + d = 0
        plane(a, b, c, d) {
            return this.e1.scale(a)
                .add(this.e2.scale(b))
                .add(this.e3.scale(c))
                .add(this.e0.scale(d));
        }

        // Extract coordinates from point
        pointCoords(p) {
            const w = p.coeffs[this.bladeIndex['e1e2e3']] || 1;
            const x = -p.coeffs[this.bladeIndex['e0e2e3']] / w;
            const y = p.coeffs[this.bladeIndex['e0e1e3']] / w;
            const z = -p.coeffs[this.bladeIndex['e0e1e2']] / w;
            return [x, y, z];
        }

        // Line through two points
        join(p1, p2) {
            return p1.vee(p2);
        }

        // Line as intersection of two planes
        meet(l1, l2) {
            return l1.wedge(l2);
        }

        // Translator by (dx, dy, dz)
        translator(dx, dy, dz) {
            // T = 1 + (dx*e01 + dy*e02 + dz*e03)/2
            return this.scalar(1)
                .add(this.e01.scale(dx / 2))
                .add(this.e02.scale(dy / 2))
                .add(this.e03.scale(dz / 2));
        }

        // Rotor from axis (ax, ay, az) and angle
        rotorAxis(ax, ay, az, angle) {
            const norm = sqrt(ax*ax + ay*ay + az*az);
            if (norm < EPSILON) return this.scalar(1);
            const half = angle / 2;
            const s = sin(half) / norm;
            // Axis → bivector: ax*e23 + ay*e31 + az*e12
            return this.scalar(cos(half))
                .add(this.e23.scale(ax * s))
                .add(this.e31.scale(ay * s))
                .add(this.e12.scale(az * s));
        }

        // Rotor from bivector plane
        rotorBivector(bivector, angle) {
            const B = bivector.grade(2).normalized();
            const half = angle / 2;
            return this.scalar(cos(half)).add(B.scale(sin(half)));
        }

        // Motor: rotation + translation
        // Following "Graded Symmetry Groups" decomposition
        motor(ax, ay, az, angle, dx, dy, dz) {
            const R = this.rotorAxis(ax, ay, az, angle);
            const T = this.translator(dx, dy, dz);
            return T.mul(R);
        }

        // Apply versor to point: V * P * ~V
        apply(versor, point) {
            return versor.sandwich(point);
        }

        // Euler angles to rotor (ZYX convention)
        rotorEuler(roll, pitch, yaw) {
            const Rx = this.rotorAxis(1, 0, 0, roll);
            const Ry = this.rotorAxis(0, 1, 0, pitch);
            const Rz = this.rotorAxis(0, 0, 1, yaw);
            return Rz.mul(Ry).mul(Rx);
        }

        // Efficient rotation for animation (bypasses full algebra)
        rotatePoint(x, y, z, rotor) {
            // Extract rotor components
            const s = rotor.scalar;
            const xy = rotor.coeffs[this.bladeIndex['e1e2']];
            const xz = rotor.coeffs[this.bladeIndex['e1e3']];  // e31
            const yz = rotor.coeffs[this.bladeIndex['e2e3']];
            
            // Quaternion-equivalent formula
            const i = yz, j = -xz, k = xy;
            const cx = j*z - k*y, cy = k*x - i*z, cz = i*y - j*x;
            const ccx = j*cz - k*cy, ccy = k*cx - i*cz, ccz = i*cy - j*cx;
            
            return [
                x + 2*(s*cx + ccx),
                y + 2*(s*cy + ccy),
                z + 2*(s*cz + ccz)
            ];
        }
    }

    // ============================================================================
    // CGA 2D - Conformal Geometric Algebra Cl(3,1,0)
    // ============================================================================

    class CGA2D extends Algebra {
        constructor() {
            // Signature: e1²=1, e2²=1, e+²=1, e-²=-1
            super(3, 1, 0, ['e1', 'e2', 'ep', 'em']);
            
            // Null basis: n_o (origin), n_∞ (infinity)
            // n_o = (e- - e+)/2, n_∞ = e- + e+
            this.ep = this.blade(0b0100);  // e+
            this.em = this.blade(0b1000);  // e-
            
            this.no = this.em.sub(this.ep).scale(0.5);  // Origin
            this.ni = this.em.add(this.ep);             // Infinity
            
            this.e1 = this.blade(0b0001);
            this.e2 = this.blade(0b0010);
        }

        // Embed Euclidean point into CGA
        // X = x + 0.5*x²*n_∞ + n_o
        point(x, y) {
            const euclidean = this.e1.scale(x).add(this.e2.scale(y));
            const sq = x*x + y*y;
            return euclidean.add(this.ni.scale(0.5 * sq)).add(this.no);
        }

        // Extract (x, y) from CGA point
        pointCoords(X) {
            // X = x*e1 + y*e2 + ..., normalize by -X·n_∞
            const div = -X.dot(this.ni).scalar;
            if (isZero(div)) return [0, 0];
            const x = X.coeffs[this.bladeIndex['e1']] / div;
            const y = X.coeffs[this.bladeIndex['e2']] / div;
            return [x, y];
        }

        // Circle through 3 points
        circle(p1, p2, p3) {
            const P1 = this.point(p1[0], p1[1]);
            const P2 = this.point(p2[0], p2[1]);
            const P3 = this.point(p3[0], p3[1]);
            return P1.wedge(P2).wedge(P3);
        }

        // Circle from center and radius
        circleRadius(cx, cy, r) {
            const C = this.point(cx, cy);
            return C.sub(this.ni.scale(0.5 * r * r));
        }

        // Line through 2 points (includes n_∞)
        line(p1, p2) {
            const P1 = this.point(p1[0], p1[1]);
            const P2 = this.point(p2[0], p2[1]);
            return P1.wedge(P2).wedge(this.ni);
        }

        // Point pair from 2 points
        pointPair(p1, p2) {
            const P1 = this.point(p1[0], p1[1]);
            const P2 = this.point(p2[0], p2[1]);
            return P1.wedge(P2);
        }

        // Translator by (dx, dy)
        translator(dx, dy) {
            const t = this.e1.scale(dx).add(this.e2.scale(dy));
            return this.scalar(1).sub(t.wedge(this.ni).scale(0.5));
        }

        // Rotor (rotation around origin)
        rotorOrigin(angle) {
            const B = this.e1.wedge(this.e2);
            return this.scalar(cos(angle/2)).add(B.scale(sin(angle/2)));
        }

        // Dilator (uniform scaling from origin)
        dilator(factor) {
            const alpha = log(factor) / 2;
            return this.scalar(cos(alpha)).add(this.no.wedge(this.ni).scale(sin(alpha)));
        }

        // Inverter (circle inversion, radius 1)
        inverter() {
            return this.no.sub(this.ni.scale(0.5));
        }
    }

    // ============================================================================
    // CGA 3D - Conformal Geometric Algebra Cl(4,1,0)
    // Based on "From Invariant Decomposition to Spinors"
    // ============================================================================

    class CGA3D extends Algebra {
        constructor() {
            // Signature: e1²=e2²=e3²=e+²=1, e-²=-1
            super(4, 1, 0, ['e1', 'e2', 'e3', 'ep', 'em']);
            
            this.e1 = this.blade(0b00001);
            this.e2 = this.blade(0b00010);
            this.e3 = this.blade(0b00100);
            this.ep = this.blade(0b01000);
            this.em = this.blade(0b10000);
            
            // Null basis
            this.no = this.em.sub(this.ep).scale(0.5);  // Origin
            this.ni = this.em.add(this.ep);             // Infinity
            
            // E3 pseudoscalar
            this.I3 = this.e1.wedge(this.e2).wedge(this.e3);
        }

        // Embed Euclidean point: X = x + 0.5*x²*n_∞ + n_o
        point(x, y, z) {
            const euclidean = this.e1.scale(x).add(this.e2.scale(y)).add(this.e3.scale(z));
            const sq = x*x + y*y + z*z;
            return euclidean.add(this.ni.scale(0.5 * sq)).add(this.no);
        }

        // Extract coordinates
        pointCoords(X) {
            const div = -X.dot(this.ni).scalar;
            if (isZero(div)) return [0, 0, 0];
            return [
                X.coeffs[this.bladeIndex['e1']] / div,
                X.coeffs[this.bladeIndex['e2']] / div,
                X.coeffs[this.bladeIndex['e3']] / div
            ];
        }

        // Sphere through 4 points
        sphere4(p1, p2, p3, p4) {
            const pts = [p1, p2, p3, p4].map(p => this.point(...p));
            return pts[0].wedge(pts[1]).wedge(pts[2]).wedge(pts[3]);
        }

        // Sphere from center and radius
        sphereRadius(cx, cy, cz, r) {
            const C = this.point(cx, cy, cz);
            return C.sub(this.ni.scale(0.5 * r * r));
        }

        // Plane through 3 points (+ n_∞)
        plane3(p1, p2, p3) {
            const pts = [p1, p2, p3].map(p => this.point(...p));
            return pts[0].wedge(pts[1]).wedge(pts[2]).wedge(this.ni);
        }

        // Plane from normal and distance
        planeNd(nx, ny, nz, d) {
            const n = this.e1.scale(nx).add(this.e2.scale(ny)).add(this.e3.scale(nz));
            return n.add(this.ni.scale(d));
        }

        // Circle through 3 points
        circle3(p1, p2, p3) {
            const pts = [p1, p2, p3].map(p => this.point(...p));
            return pts[0].wedge(pts[1]).wedge(pts[2]);
        }

        // Line through 2 points
        line2(p1, p2) {
            const P1 = this.point(...p1);
            const P2 = this.point(...p2);
            return P1.wedge(P2).wedge(this.ni);
        }

        // Point pair
        pointPair(p1, p2) {
            return this.point(...p1).wedge(this.point(...p2));
        }

        // Translator
        translator(dx, dy, dz) {
            const t = this.e1.scale(dx).add(this.e2.scale(dy)).add(this.e3.scale(dz));
            return this.scalar(1).sub(t.wedge(this.ni).scale(0.5));
        }

        // Rotor from axis and angle
        rotorAxis(ax, ay, az, angle) {
            const norm = sqrt(ax*ax + ay*ay + az*az);
            if (norm < EPSILON) return this.scalar(1);
            
            // Bivector for rotation plane
            const half = angle / 2;
            const s = sin(half) / norm;
            const B = this.e2.wedge(this.e3).scale(ax * s)
                .add(this.e3.wedge(this.e1).scale(ay * s))
                .add(this.e1.wedge(this.e2).scale(az * s));
            
            return this.scalar(cos(half)).add(B);
        }

        // Dilator (uniform scaling)
        dilator(factor) {
            const alpha = log(factor) / 2;
            const B = this.no.wedge(this.ni);
            return this.scalar(cos(alpha)).add(B.scale(sin(alpha)));
        }

        // Transversor (special conformal transformation)
        // Composition: translation, inversion, translation, inversion
        transversor(ax, ay, az) {
            const a = this.e1.scale(ax).add(this.e2.scale(ay)).add(this.e3.scale(az));
            return this.scalar(1).add(this.no.wedge(a).scale(0.5));
        }

        // Motor (rotation + translation) per Roelfs & De Keninck
        motor(ax, ay, az, angle, dx, dy, dz) {
            const R = this.rotorAxis(ax, ay, az, angle);
            const T = this.translator(dx, dy, dz);
            return T.mul(R);
        }

        // Spherical linear interpolation for rotors
        slerp(R1, R2, t) {
            const R12 = R1.reverse().mul(R2);
            const logR = this.log(R12);
            const Rt = this.exp(logR.scale(t));
            return R1.mul(Rt);
        }

        // Efficient point rotation (bypasses full algebra)
        rotatePoint(x, y, z, rotor) {
            const s = rotor.scalar;
            const e12 = rotor.coeffs[this.bladeIndex['e1e2']];
            const e13 = rotor.coeffs[this.bladeIndex['e1e3']];
            const e23 = rotor.coeffs[this.bladeIndex['e2e3']];
            
            // Cross product formulation
            const i = e23, j = -e13, k = e12;
            const cx = j*z - k*y, cy = k*x - i*z, cz = i*y - j*x;
            const ccx = j*cz - k*cy, ccy = k*cx - i*cz, ccz = i*cy - j*cx;
            
            return [x + 2*(s*cx + ccx), y + 2*(s*cy + ccy), z + 2*(s*cz + ccz)];
        }
    }

    // ============================================================================
    // SIMPLE ROTOR CLASS (for quick 3D rotations without full algebra overhead)
    // ============================================================================

    class Rotor3D {
        constructor(s = 1, xy = 0, xz = 0, yz = 0) {
            this.s = s;    // scalar
            this.xy = xy;  // e12 bivector
            this.xz = xz;  // e13 bivector
            this.yz = yz;  // e23 bivector
        }

        static fromAxisAngle(ax, ay, az, angle) {
            const norm = sqrt(ax*ax + ay*ay + az*az);
            if (norm < EPSILON) return new Rotor3D(1, 0, 0, 0);
            const half = angle / 2;
            const sinH = sin(half) / norm;
            return new Rotor3D(cos(half), az*sinH, ay*sinH, ax*sinH);
        }

        static fromEuler(roll, pitch, yaw) {
            const Rx = Rotor3D.fromAxisAngle(1, 0, 0, roll);
            const Ry = Rotor3D.fromAxisAngle(0, 1, 0, pitch);
            const Rz = Rotor3D.fromAxisAngle(0, 0, 1, yaw);
            return Rz.mul(Ry.mul(Rx));
        }

        mul(r) {
            return new Rotor3D(
                this.s*r.s - this.xy*r.xy - this.xz*r.xz - this.yz*r.yz,
                this.s*r.xy + this.xy*r.s + this.xz*r.yz - this.yz*r.xz,
                this.s*r.xz + this.xz*r.s - this.xy*r.yz + this.yz*r.xy,
                this.s*r.yz + this.yz*r.s + this.xy*r.xz - this.xz*r.xy
            );
        }

        rev() { return new Rotor3D(this.s, -this.xy, -this.xz, -this.yz); }

        norm() { return sqrt(this.s*this.s + this.xy*this.xy + this.xz*this.xz + this.yz*this.yz); }

        normalized() {
            const n = this.norm();
            return new Rotor3D(this.s/n, this.xy/n, this.xz/n, this.yz/n);
        }

        // Sandwich product: R*v*R̃
        apply(x, y, z) {
            const w = this.s, i = this.yz, j = -this.xz, k = this.xy;
            const cx = j*z - k*y, cy = k*x - i*z, cz = i*y - j*x;
            const ccx = j*cz - k*cy, ccy = k*cx - i*cz, ccz = i*cy - j*cx;
            return [x + 2*(w*cx + ccx), y + 2*(w*cy + ccy), z + 2*(w*cz + ccz)];
        }

        static slerp(r1, r2, t) {
            let dot = r1.s*r2.s + r1.xy*r2.xy + r1.xz*r2.xz + r1.yz*r2.yz;
            if (dot < 0) { r2 = new Rotor3D(-r2.s, -r2.xy, -r2.xz, -r2.yz); dot = -dot; }
            if (dot > 0.9995) {
                return new Rotor3D(
                    r1.s + t*(r2.s - r1.s), r1.xy + t*(r2.xy - r1.xy),
                    r1.xz + t*(r2.xz - r1.xz), r1.yz + t*(r2.yz - r1.yz)
                ).normalized();
            }
            const theta = acos(dot);
            const sinT = sin(theta);
            const a = sin((1-t)*theta) / sinT, b = sin(t*theta) / sinT;
            return new Rotor3D(
                a*r1.s + b*r2.s, a*r1.xy + b*r2.xy,
                a*r1.xz + b*r2.xz, a*r1.yz + b*r2.yz
            );
        }

        toString() {
            return `${this.s.toFixed(3)} + ${this.xy.toFixed(3)}e₁₂ + ${this.xz.toFixed(3)}e₁₃ + ${this.yz.toFixed(3)}e₂₃`;
        }
    }

    // ============================================================================
    // EXPORTS
    // ============================================================================

    const GA = {
        Algebra,
        Multivector,
        PGA2D,
        PGA3D,
        CGA2D,
        CGA3D,
        Rotor3D,
        
        // Factory functions
        pga2d: () => new PGA2D(),
        pga3d: () => new PGA3D(),
        cga2d: () => new CGA2D(),
        cga3d: () => new CGA3D(),
        
        // Quick rotor creation
        rotor: Rotor3D.fromAxisAngle,
        rotorEuler: Rotor3D.fromEuler
    };

    // Export for different module systems
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = GA;
    } else if (typeof define === 'function' && define.amd) {
        define([], () => GA);
    } else {
        global.GA = GA;
    }

})(typeof window !== 'undefined' ? window : global);
