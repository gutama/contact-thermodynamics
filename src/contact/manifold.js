/**
 * Contact Manifold Module
 * 
 * Contact manifolds and points for thermodynamic/mechanics applications.
 * Extracted from index.js for modular organization.
 * 
 * @module contact/manifold
 * @license MIT
 */

(function (global) {
    'use strict';

    const EPSILON = 1e-12;

    // ============================================================================
    // CONTACT MANIFOLD BASE CLASS
    // ============================================================================

    /**
     * ContactManifold: Base class for contact manifolds
     * A contact manifold is (M, α) where α ∧ (dα)^n ≠ 0
     * Dimension is always 2n+1 for some n ≥ 1
     */
    class ContactManifold {
        /**
         * @param {string[]} baseCoords - Base configuration coordinates x^a
         * @param {string[]} momentaCoords - Conjugate momenta p_a
         * @param {string} fiberCoord - Fiber coordinate u (action/potential)
         */
        constructor(baseCoords, momentaCoords, fiberCoord = 'A') {
            if (baseCoords.length !== momentaCoords.length) {
                throw new Error('Base and momenta coordinates must have same dimension');
            }
            this.baseCoords = baseCoords;
            this.momentaCoords = momentaCoords;
            this.fiberCoord = fiberCoord;
            this.n = baseCoords.length;  // dim(Q)
            this.dim = 2 * this.n + 1;   // dim(M) = 2n + 1
            this._allCoords = Object.freeze([
                ...this.baseCoords,
                this.fiberCoord,
                ...this.momentaCoords
            ]);
        }

        /**
         * Dimension theorem: dim(M) = 2n + 1
         */
        get dimension() {
            return this.dim;
        }

        /**
         * Get coordinate names in canonical order
         */
        get allCoords() {
            return this._allCoords.slice();
        }

        /**
         * Create a point on the contact manifold
         * @param {Object} coords - Coordinate values { x^a: val, u: val, p_a: val }
         */
        point(coords) {
            return new ContactPoint(this, coords);
        }

        /**
         * Zero point (origin of coordinates)
         */
        get origin() {
            const coords = {};
            this._allCoords.forEach(c => coords[c] = 0);
            return this.point(coords);
        }

        /**
         * Canonical contact 1-form: α = du - p_a dx^a
         * Returns symbolic representation
         */
        contactFormSymbolic() {
            const terms = [`d${this.fiberCoord}`];
            for (let i = 0; i < this.n; i++) {
                terms.push(`-${this.momentaCoords[i]}·d${this.baseCoords[i]}`);
            }
            return terms.join(' ');
        }

        /**
         * Evaluate contact form at a point with tangent vector
         * α(v) = du(v) - p_a dx^a(v)
         * @param {ContactPoint} pt - Point on manifold
         * @param {Object} tangent - Tangent vector components
         */
        evaluateContactForm(pt, tangent) {
            let result = tangent[this.fiberCoord] || 0;
            for (let i = 0; i < this.n; i++) {
                const p = pt.get(this.momentaCoords[i]);
                const dx = tangent[this.baseCoords[i]] || 0;
                result -= p * dx;
            }
            return result;
        }

        /**
         * Verify non-degeneracy: α ∧ (dα)^n ≠ 0
         * Returns the wedge product value at a point
         */
        verifyContactCondition(pt) {
            // For canonical form, this is always satisfied
            // α ∧ (dα)^n = n! du ∧ dx^1 ∧ dp_1 ∧ ... ∧ dx^n ∧ dp_n
            return this.factorial(this.n);
        }

        factorial(n) {
            return n <= 1 ? 1 : n * this.factorial(n - 1);
        }

        /**
         * Reeb vector field R
         * Defined by: α(R) = 1, ι_R dα = 0
         * For canonical α = du - p_a dx^a, R = ∂/∂u
         */
        reebField(pt) {
            const R = {};
            this._allCoords.forEach(c => R[c] = 0);
            R[this.fiberCoord] = 1;  // ∂/∂u
            return R;
        }

        /**
         * Information string
         */
        toString() {
            return `ContactManifold J¹(Q${this.n}): dim=${this.dim}, α=${this.contactFormSymbolic()}`;
        }
    }

    // ============================================================================
    // CONTACT POINT CLASS
    // ============================================================================

    /**
     * ContactPoint: A point on a contact manifold
     */
    class ContactPoint {
        constructor(manifold, coords = {}) {
            this.manifold = manifold;
            this.coords = {};
            manifold._allCoords.forEach(c => {
                this.coords[c] = coords[c] !== undefined ? coords[c] : 0;
            });
        }

        get(coord) {
            return this.coords[coord];
        }

        set(coord, value) {
            if (this.manifold._allCoords.includes(coord)) {
                this.coords[coord] = value;
            }
            return this;
        }

        clone() {
            return new ContactPoint(this.manifold, { ...this.coords });
        }

        /**
         * Add tangent vector (flow)
         */
        add(tangent, dt = 1) {
            const newPt = this.clone();
            for (const c of this.manifold._allCoords) {
                if (tangent[c] !== undefined) {
                    newPt.coords[c] += tangent[c] * dt;
                }
            }
            return newPt;
        }

        toString() {
            const parts = this.manifold._allCoords.map(c =>
                `${c}=${this.coords[c].toFixed(4)}`
            );
            return `(${parts.join(', ')})`;
        }
    }

    // ============================================================================
    // GRAND CONTACT MANIFOLD M₁₃
    // ============================================================================

    /**
     * GrandContactManifold: The "honest" 13-dimensional contact manifold
     * 
     * Base Q₆: x^a = (q¹, q², q³, t, ℓ, S) where ℓ = log(λ)
     * Momenta: p_a = (k₁, k₂, k₃, ω, Δ, T)
     * Fiber: A (action/generating potential)
     * 
     * dim(Q₆) = 6, dim(M₁₃) = 2·6 + 1 = 13
     * 
     * Contact form: α = dA - k_i dq^i - ω dt - Δ dℓ - T dS
     */
    class GrandContactManifold extends ContactManifold {
        constructor() {
            super(
                ['q1', 'q2', 'q3', 't', 'ell', 'S'],  // Base coords
                ['k1', 'k2', 'k3', 'omega', 'Delta', 'T'],  // Momenta
                'A'  // Fiber (action)
            );

            // Physical interpretations
            this.physicalInterpretations = {
                q1: 'spatial coordinate 1',
                q2: 'spatial coordinate 2',
                q3: 'spatial coordinate 3',
                t: 'time',
                ell: 'log(scale factor) = log(λ)',
                S: 'entropy',
                k1: 'wavenumber conjugate to q1',
                k2: 'wavenumber conjugate to q2',
                k3: 'wavenumber conjugate to q3',
                omega: 'frequency conjugate to t',
                Delta: 'dilatation/anomalous dimension conjugate to ℓ',
                T: 'temperature conjugate to S (under entropy closure)',
                A: 'action/generating potential'
            };
        }

        /**
         * Create point with physical naming
         */
        physicalPoint(q1, q2, q3, t, ell, S, k1, k2, k3, omega, Delta, T, A = 0) {
            return this.point({
                q1, q2, q3, t, ell, S,
                k1, k2, k3, omega, Delta, T,
                A
            });
        }

        /**
         * Spatial position vector (q¹, q², q³)
         */
        spatialPosition(pt) {
            return [pt.get('q1'), pt.get('q2'), pt.get('q3')];
        }

        /**
         * Wave vector (k₁, k₂, k₃)
         */
        waveVector(pt) {
            return [pt.get('k1'), pt.get('k2'), pt.get('k3')];
        }

        /**
         * Extended contact form (symbolic with physical meaning)
         */
        contactFormSymbolic() {
            return 'dA - k₁dq¹ - k₂dq² - k₃dq³ - ω dt - Δ dℓ - T dS';
        }

        /**
         * Verify contact condition: α ∧ (dα)⁶ ≠ 0
         */
        verifyContactCondition(pt) {
            // For canonical form, always non-degenerate
            // α ∧ (dα)⁶ = 6! · volume form
            return 720; // 6!
        }

        toString() {
            return 'Grand Contact Manifold M₁₃ = J¹(Q₆): dim=13';
        }
    }

    // ============================================================================
    // HOLOGRAPHIC CONTACT MANIFOLD M₇
    // ============================================================================

    /**
     * HolographicContactManifold: The 7-dimensional "holographic limit"
     * 
     * Base Q₃: x^a = (t, ℓ, S) where ℓ = log(λ)
     * Momenta: p_a = (ω, Δ, T)
     * Fiber: A (action)
     * 
     * dim(Q₃) = 3, dim(M₇) = 2·3 + 1 = 7
     * 
     * Space q^i becomes emergent: q^i = q^i(t, ℓ, S) as scalar fields on Q₃
     * 
     * Contact form: α = dA - ω dt - Δ dℓ - T dS
     */
    class HolographicContactManifold extends ContactManifold {
        constructor() {
            super(
                ['t', 'ell', 'S'],  // Base coords (reduced)
                ['omega', 'Delta', 'T'],  // Momenta
                'A'  // Fiber (action)
            );

            // Emergent spatial fields (not coordinates, but dependent fields)
            this.emergentFields = ['q1', 'q2', 'q3'];

            this.physicalInterpretations = {
                t: 'time',
                ell: 'log(scale factor) = log(λ)',
                S: 'entropy',
                omega: 'frequency',
                Delta: 'dilatation/anomalous dimension',
                T: 'temperature',
                A: 'action',
                q1: 'emergent spatial field q¹(t,ℓ,S)',
                q2: 'emergent spatial field q²(t,ℓ,S)',
                q3: 'emergent spatial field q³(t,ℓ,S)'
            };
        }

        /**
         * Create holographic point
         */
        holographicPoint(t, ell, S, omega, Delta, T, A = 0) {
            return this.point({ t, ell, S, omega, Delta, T, A });
        }

        /**
         * Define emergent spatial configuration
         * q^i = q^i(t, ℓ, S) as scalar fields
         * @param {Function} fieldFunc - (t, ell, S) => [q1, q2, q3]
         */
        createEmergentSpace(pt, fieldFunc) {
            const [q1, q2, q3] = fieldFunc(pt.get('t'), pt.get('ell'), pt.get('S'));
            return { q1, q2, q3 };
        }

        /**
         * Extended point including emergent fields
         */
        extendedPoint(pt, emergent) {
            return {
                ...pt.coords,
                ...emergent
            };
        }

        contactFormSymbolic() {
            return 'dA - ω dt - Δ dℓ - T dS';
        }

        verifyContactCondition(pt) {
            // α ∧ (dα)³ = 3! · volume form
            return 6; // 3!
        }

        toString() {
            return 'Holographic Contact Manifold M₇ = J¹(Q₃): dim=7, space emergent';
        }
    }

    // ============================================================================
    // EXPORTS
    // ============================================================================

    const Manifold = {
        ContactManifold,
        ContactPoint,
        GrandContactManifold,
        HolographicContactManifold,

        // Convenience factory
        grandManifold: () => new GrandContactManifold(),
        holographicManifold: () => new HolographicContactManifold()
    };

    // Export for different module systems
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = Manifold;
    } else if (typeof define === 'function' && define.amd) {
        define([], () => Manifold);
    } else {
        global.ContactManifold = Manifold;
    }

})(typeof window !== 'undefined' ? window : global);
