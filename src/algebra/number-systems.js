/**
 * Number Systems & Bivector Classification
 * 
 * Implements the three 2D number systems (Complex, Dual, Hyperbolic)
 * and bivector classification utilities for transformation analysis.
 * 
 * @license MIT
 */

(function(global) {
    'use strict';

    const EPSILON = 1e-10;
    const { sqrt, sin, cos, sinh, cosh, atan2, abs, PI } = Math;

    // ============================================================================
    // THREE NUMBER SYSTEMS
    // ============================================================================

    /**
     * Complex Numbers: z = a + bi, i² = -1
     * Represents Cl(0,1) - elliptic/rotation
     */
    class Complex {
        constructor(a = 0, b = 0) {
            this.a = a;  // Real part
            this.b = b;  // Imaginary part
        }

        // Arithmetic
        add(z) { return new Complex(this.a + z.a, this.b + z.b); }
        sub(z) { return new Complex(this.a - z.a, this.b - z.b); }
        scale(s) { return new Complex(this.a * s, this.b * s); }
        
        mul(z) {
            return new Complex(
                this.a * z.a - this.b * z.b,  // i² = -1
                this.a * z.b + this.b * z.a
            );
        }

        div(z) {
            const d = z.a * z.a + z.b * z.b;
            if (d < EPSILON) return null;
            return new Complex(
                (this.a * z.a + this.b * z.b) / d,
                (this.b * z.a - this.a * z.b) / d
            );
        }

        // Properties
        conjugate() { return new Complex(this.a, -this.b); }
        normSq() { return this.a * this.a + this.b * this.b; }  // |z|² = a² + b²
        norm() { return sqrt(this.normSq()); }
        arg() { return atan2(this.b, this.a); }

        normalized() {
            const n = this.norm();
            return n < EPSILON ? new Complex(1, 0) : this.scale(1 / n);
        }

        // Exponential: exp(iθ) = cos(θ) + i·sin(θ)
        static exp(theta) {
            return new Complex(cos(theta), sin(theta));
        }

        // Logarithm (principal branch)
        log() {
            return new Complex(Math.log(this.norm()), this.arg());
        }

        // Power
        pow(n) {
            const r = this.norm();
            const theta = this.arg();
            const rn = Math.pow(r, n);
            return new Complex(rn * cos(n * theta), rn * sin(n * theta));
        }

        // Unit curve: circle |z| = 1
        static unitCurve(t) {
            return new Complex(cos(t), sin(t));
        }

        toString() {
            const sign = this.b >= 0 ? '+' : '-';
            return `${this.a.toFixed(4)} ${sign} ${abs(this.b).toFixed(4)}i`;
        }

        static get type() { return 'elliptic'; }
        static get iSquared() { return -1; }
    }


    /**
     * Dual Numbers: z = a + bε, ε² = 0
     * Represents Cl(0,0,1) - parabolic/translation
     */
    class Dual {
        constructor(a = 0, b = 0) {
            this.a = a;  // Real part
            this.b = b;  // Dual part (infinitesimal)
        }

        // Arithmetic
        add(z) { return new Dual(this.a + z.a, this.b + z.b); }
        sub(z) { return new Dual(this.a - z.a, this.b - z.b); }
        scale(s) { return new Dual(this.a * s, this.b * s); }
        
        mul(z) {
            return new Dual(
                this.a * z.a,                  // ε² = 0, so b₁b₂ term vanishes
                this.a * z.b + this.b * z.a
            );
        }

        div(z) {
            if (abs(z.a) < EPSILON) return null;
            return new Dual(
                this.a / z.a,
                (this.b * z.a - this.a * z.b) / (z.a * z.a)
            );
        }

        // Properties
        conjugate() { return new Dual(this.a, -this.b); }
        normSq() { return this.a * this.a; }  // |z|² = a² (only real part)
        norm() { return abs(this.a); }

        normalized() {
            if (abs(this.a) < EPSILON) return new Dual(1, 0);
            return new Dual(this.a / abs(this.a), this.b / abs(this.a));
        }

        // Exponential: exp(εθ) = 1 + εθ (series terminates!)
        static exp(theta) {
            return new Dual(1, theta);
        }

        // Logarithm (when a > 0)
        log() {
            if (this.a <= 0) return null;
            return new Dual(Math.log(this.a), this.b / this.a);
        }

        // Power (integer)
        pow(n) {
            // (a + bε)^n = a^n + n·a^(n-1)·b·ε
            const an = Math.pow(this.a, n);
            const an1 = Math.pow(this.a, n - 1);
            return new Dual(an, n * an1 * this.b);
        }

        // Unit curve: |z|² = 1 means a = ±1 (vertical lines)
        static unitCurve(t) {
            return new Dual(1, t);  // Right line a = 1
        }

        /**
         * Automatic differentiation!
         * f(a + bε) = f(a) + b·f'(a)·ε
         * 
         * Example: To compute f(x) and f'(x) simultaneously:
         *   const result = f(new Dual(x, 1));
         *   // result.a = f(x), result.b = f'(x)
         */
        static differentiate(f, x) {
            const result = f(new Dual(x, 1));
            return { value: result.a, derivative: result.b };
        }

        toString() {
            const sign = this.b >= 0 ? '+' : '-';
            return `${this.a.toFixed(4)} ${sign} ${abs(this.b).toFixed(4)}ε`;
        }

        static get type() { return 'parabolic'; }
        static get iSquared() { return 0; }
    }


    /**
     * Hyperbolic (Split-Complex) Numbers: z = a + bj, j² = +1
     * Represents Cl(1,0) - hyperbolic/boost
     */
    class Hyperbolic {
        constructor(a = 0, b = 0) {
            this.a = a;  // Real part (timelike)
            this.b = b;  // Hyperbolic part (spacelike)
        }

        // Arithmetic
        add(z) { return new Hyperbolic(this.a + z.a, this.b + z.b); }
        sub(z) { return new Hyperbolic(this.a - z.a, this.b - z.b); }
        scale(s) { return new Hyperbolic(this.a * s, this.b * s); }
        
        mul(z) {
            return new Hyperbolic(
                this.a * z.a + this.b * z.b,  // j² = +1, so PLUS sign
                this.a * z.b + this.b * z.a
            );
        }

        div(z) {
            const d = z.a * z.a - z.b * z.b;  // Minkowski norm
            if (abs(d) < EPSILON) return null;  // On null cone
            return new Hyperbolic(
                (this.a * z.a - this.b * z.b) / d,
                (this.b * z.a - this.a * z.b) / d
            );
        }

        // Properties
        conjugate() { return new Hyperbolic(this.a, -this.b); }
        
        // MINKOWSKI norm: a² - b² (can be negative!)
        normSq() { return this.a * this.a - this.b * this.b; }
        
        // Magnitude (absolute value of Minkowski norm)
        norm() { return sqrt(abs(this.normSq())); }

        // Causal character
        isTimelike() { return this.normSq() > EPSILON; }
        isSpacelike() { return this.normSq() < -EPSILON; }
        isNull() { return abs(this.normSq()) < EPSILON; }

        // Rapidity (for timelike, a > |b|)
        rapidity() {
            if (!this.isTimelike() || this.a <= 0) return null;
            return Math.atanh(this.b / this.a);
        }

        normalized() {
            const ns = this.normSq();
            if (abs(ns) < EPSILON) return null;  // Null vectors can't be normalized
            const n = sqrt(abs(ns));
            return new Hyperbolic(this.a / n, this.b / n);
        }

        // Exponential: exp(jφ) = cosh(φ) + j·sinh(φ)
        static exp(phi) {
            return new Hyperbolic(cosh(phi), sinh(phi));
        }

        // Logarithm (for timelike, a > |b|)
        log() {
            if (!this.isTimelike() || this.a <= 0) return null;
            const r = sqrt(this.normSq());
            const phi = Math.atanh(this.b / this.a);
            return new Hyperbolic(Math.log(r), phi);
        }

        // Unit curve: a² - b² = 1 (hyperbola)
        static unitCurve(t) {
            return new Hyperbolic(cosh(t), sinh(t));
        }

        // Left branch of hyperbola
        static unitCurveLeft(t) {
            return new Hyperbolic(-cosh(t), sinh(t));
        }

        // Null directions (light cone): a = ±b
        static nullPlus(t) { return new Hyperbolic(t, t); }
        static nullMinus(t) { return new Hyperbolic(t, -t); }

        /**
         * Lorentz boost application
         * Given rapidity φ, transforms (t, x) coordinates
         */
        static boost(t, x, phi) {
            const ch = cosh(phi), sh = sinh(phi);
            return {
                t: ch * t + sh * x,
                x: sh * t + ch * x
            };
        }

        /**
         * Velocity addition (relativistic)
         * Rapidities add: φ_total = φ₁ + φ₂
         * Velocities compose: v_total = (v₁ + v₂)/(1 + v₁v₂/c²)
         */
        static velocityAdd(v1, v2, c = 1) {
            return (v1 + v2) / (1 + v1 * v2 / (c * c));
        }

        toString() {
            const sign = this.b >= 0 ? '+' : '-';
            return `${this.a.toFixed(4)} ${sign} ${abs(this.b).toFixed(4)}j`;
        }

        static get type() { return 'hyperbolic'; }
        static get iSquared() { return 1; }
    }


    // ============================================================================
    // BIVECTOR CLASSIFICATION
    // ============================================================================

    /**
     * Classification of bivectors by their square
     */
    const BivectorType = {
        ELLIPTIC: 'elliptic',     // B² = -1, rotation
        PARABOLIC: 'parabolic',   // B² = 0,  translation
        HYPERBOLIC: 'hyperbolic', // B² = +1, boost/dilation
        MIXED: 'mixed'            // General case
    };

    /**
     * Classify a bivector based on its square
     * @param {number} Bsquared - The value of B²
     * @returns {string} - One of BivectorType values
     */
    function classifyBivectorSquare(Bsquared) {
        if (abs(Bsquared + 1) < EPSILON) return BivectorType.ELLIPTIC;
        if (abs(Bsquared) < EPSILON) return BivectorType.PARABOLIC;
        if (abs(Bsquared - 1) < EPSILON) return BivectorType.HYPERBOLIC;
        return BivectorType.MIXED;
    }

    /**
     * Compute exp(θB) based on bivector type
     * Returns {scalar, bivector_coeff} such that result = scalar + bivector_coeff * B
     * 
     * @param {number} theta - The angle/parameter
     * @param {string} type - BivectorType
     * @returns {object} - {s: scalar_part, b: bivector_coefficient}
     */
    function expBivectorCoeffs(theta, type) {
        const half = theta / 2;
        switch (type) {
            case BivectorType.ELLIPTIC:
                return { s: cos(half), b: sin(half) };
            case BivectorType.PARABOLIC:
                return { s: 1, b: half };
            case BivectorType.HYPERBOLIC:
                return { s: cosh(half), b: sinh(half) };
            default:
                throw new Error(`Cannot exponentiate mixed bivector type`);
        }
    }

    /**
     * Get the number system corresponding to a bivector type
     */
    function getNumberSystem(type) {
        switch (type) {
            case BivectorType.ELLIPTIC: return Complex;
            case BivectorType.PARABOLIC: return Dual;
            case BivectorType.HYPERBOLIC: return Hyperbolic;
            default: return null;
        }
    }


    // ============================================================================
    // CIRCLE PAIR CLASSIFICATION (CGA)
    // ============================================================================

    /**
     * Classify the relationship between two circles based on their inner product
     * This determines what type of transformation maps one to the other
     * 
     * @param {number} dot - C₁·C₂ inner product
     * @param {number} normSq1 - C₁² 
     * @param {number} normSq2 - C₂²
     * @returns {object} - {type, discriminant}
     */
    function classifyCirclePair(dot, normSq1, normSq2) {
        const discriminant = dot * dot - normSq1 * normSq2;
        
        let type;
        if (discriminant < -EPSILON) {
            type = BivectorType.ELLIPTIC;  // Intersecting
        } else if (abs(discriminant) < EPSILON) {
            type = BivectorType.PARABOLIC; // Tangent
        } else {
            type = BivectorType.HYPERBOLIC; // Nested or disjoint
        }
        
        return { type, discriminant };
    }

    /**
     * Classify two circles given their centers and radii (2D)
     * 
     * @param {number} x1, y1, r1 - First circle center and radius
     * @param {number} x2, y2, r2 - Second circle center and radius
     * @returns {object} - {type, discriminant, description}
     */
    function classifyCircles2D(x1, y1, r1, x2, y2, r2) {
        const dx = x2 - x1, dy = y2 - y1;
        const d = sqrt(dx * dx + dy * dy);  // Distance between centers
        
        let type, description;
        
        if (d < abs(r1 - r2) - EPSILON) {
            // One inside the other (not touching)
            type = BivectorType.HYPERBOLIC;
            description = 'nested';
        } else if (abs(d - abs(r1 - r2)) < EPSILON) {
            // Internally tangent
            type = BivectorType.PARABOLIC;
            description = 'internally tangent';
        } else if (d < r1 + r2 - EPSILON) {
            // Intersecting at two points
            type = BivectorType.ELLIPTIC;
            description = 'intersecting';
        } else if (abs(d - (r1 + r2)) < EPSILON) {
            // Externally tangent
            type = BivectorType.PARABOLIC;
            description = 'externally tangent';
        } else {
            // Disjoint (external)
            type = BivectorType.HYPERBOLIC;
            description = 'disjoint';
        }
        
        // Compute inversive distance (for circles)
        // δ = (d² - r₁² - r₂²) / (2r₁r₂)
        const inversiveDistance = (d * d - r1 * r1 - r2 * r2) / (2 * r1 * r2);
        
        return { type, description, inversiveDistance, centerDistance: d };
    }


    // ============================================================================
    // SPACETIME UTILITIES
    // ============================================================================

    /**
     * Classify a spacetime interval
     * @param {number} dt - Time difference
     * @param {number} dx, dy, dz - Spatial differences
     * @param {number} c - Speed of light (default 1)
     */
    function classifyInterval(dt, dx, dy = 0, dz = 0, c = 1) {
        const ds2 = c * c * dt * dt - (dx * dx + dy * dy + dz * dz);
        
        if (ds2 > EPSILON) {
            return { type: BivectorType.HYPERBOLIC, character: 'timelike', ds2 };
        } else if (ds2 < -EPSILON) {
            return { type: BivectorType.ELLIPTIC, character: 'spacelike', ds2 };
        } else {
            return { type: BivectorType.PARABOLIC, character: 'lightlike', ds2 };
        }
    }


    // ============================================================================
    // UNIFIED INTERFACE
    // ============================================================================

    /**
     * Create a number of the appropriate type
     */
    function createNumber(type, a, b) {
        switch (type) {
            case BivectorType.ELLIPTIC: return new Complex(a, b);
            case BivectorType.PARABOLIC: return new Dual(a, b);
            case BivectorType.HYPERBOLIC: return new Hyperbolic(a, b);
            default: throw new Error(`Unknown type: ${type}`);
        }
    }

    /**
     * Compute exponential for any type
     */
    function expByType(theta, type) {
        switch (type) {
            case BivectorType.ELLIPTIC: return Complex.exp(theta);
            case BivectorType.PARABOLIC: return Dual.exp(theta);
            case BivectorType.HYPERBOLIC: return Hyperbolic.exp(theta);
            default: throw new Error(`Unknown type: ${type}`);
        }
    }


    // ============================================================================
    // EXPORTS
    // ============================================================================

    const NumberSystems = {
        // Number classes
        Complex,
        Dual,
        Hyperbolic,
        
        // Bivector classification
        BivectorType,
        classifyBivectorSquare,
        expBivectorCoeffs,
        getNumberSystem,
        
        // Circle classification
        classifyCirclePair,
        classifyCircles2D,
        
        // Spacetime
        classifyInterval,
        
        // Unified interface
        createNumber,
        expByType,
        
        // Constants
        EPSILON
    };

    // Export for different module systems
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = NumberSystems;
    } else if (typeof define === 'function' && define.amd) {
        define([], () => NumberSystems);
    } else {
        global.NumberSystems = NumberSystems;
    }

})(typeof window !== 'undefined' ? window : global);
