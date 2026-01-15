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

(function(global) {
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
    // I. CONTACT MANIFOLD BASE CLASS
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
            return [...this.baseCoords, this.fiberCoord, ...this.momentaCoords];
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
            this.allCoords.forEach(c => coords[c] = 0);
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
            this.allCoords.forEach(c => R[c] = 0);
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
            manifold.allCoords.forEach(c => {
                this.coords[c] = coords[c] !== undefined ? coords[c] : 0;
            });
        }

        get(coord) {
            return this.coords[coord];
        }

        set(coord, value) {
            if (this.manifold.allCoords.includes(coord)) {
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
            for (const c of this.manifold.allCoords) {
                if (tangent[c] !== undefined) {
                    newPt.coords[c] += tangent[c] * dt;
                }
            }
            return newPt;
        }

        toString() {
            const parts = this.manifold.allCoords.map(c => 
                `${c}=${this.coords[c].toFixed(4)}`
            );
            return `(${parts.join(', ')})`;
        }
    }

    // ============================================================================
    // II. GRAND CONTACT MANIFOLD M₁₃
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
    // III. HOLOGRAPHIC CONTACT MANIFOLD M₇
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
    // IV. CONTACT HAMILTONIAN DYNAMICS
    // ============================================================================

    /**
     * ContactHamiltonian: Dynamics on contact manifolds
     * 
     * Given H: M → ℝ, the contact Hamiltonian vector field X_H satisfies:
     *   ι_{X_H} α = -H
     *   ι_{X_H} dα = dH - (RH)α
     * 
     * where R is the Reeb field
     */
    class ContactHamiltonian {
        /**
         * @param {ContactManifold} manifold
         * @param {Function} H - Hamiltonian function H(coords) → number
         * @param {Function} [dH] - Gradient of H as {coord: partial_H}
         */
        constructor(manifold, H, dH = null) {
            this.manifold = manifold;
            this.H = H;
            this._dH = dH;
        }

        /**
         * Evaluate Hamiltonian at a point
         */
        evaluate(pt) {
            return this.H(pt.coords);
        }

        /**
         * Numerical gradient of H
         */
        gradient(pt, h = 1e-7) {
            if (this._dH) {
                return this._dH(pt.coords);
            }
            
            const grad = {};
            for (const c of this.manifold.allCoords) {
                const ptPlus = pt.clone();
                const ptMinus = pt.clone();
                ptPlus.coords[c] += h;
                ptMinus.coords[c] -= h;
                grad[c] = (this.H(ptPlus.coords) - this.H(ptMinus.coords)) / (2 * h);
            }
            return grad;
        }

        /**
         * Reeb component: RH = ∂H/∂u
         */
        reebComponent(pt) {
            const grad = this.gradient(pt);
            return grad[this.manifold.fiberCoord];
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
            const Hval = this.H(pt.coords);
            
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
         * Flow the point under contact Hamiltonian dynamics
         * @param {ContactPoint} pt - Initial point
         * @param {number} dt - Time step
         * @param {number} steps - Number of steps
         */
        flow(pt, dt, steps = 1) {
            let current = pt.clone();
            const trajectory = [current.clone()];
            
            for (let i = 0; i < steps; i++) {
                // RK4 integration
                const k1 = this.vectorField(current);
                
                const pt2 = current.add(k1, dt / 2);
                const k2 = this.vectorField(pt2);
                
                const pt3 = current.add(k2, dt / 2);
                const k3 = this.vectorField(pt3);
                
                const pt4 = current.add(k3, dt);
                const k4 = this.vectorField(pt4);
                
                // Combined step
                const combined = {};
                for (const c of this.manifold.allCoords) {
                    combined[c] = (k1[c] + 2*k2[c] + 2*k3[c] + k4[c]) / 6;
                }
                
                current = current.add(combined, dt);
                trajectory.push(current.clone());
            }
            
            return trajectory;
        }

        /**
         * Check if Hamiltonian is preserved (contact case: H changes along Reeb)
         */
        hamiltonianEvolution(trajectory) {
            return trajectory.map(pt => this.evaluate(pt));
        }
    }

    // ============================================================================
    // V. LEGENDRIAN SUBMANIFOLDS & HAMILTON-JACOBI
    // ============================================================================

    /**
     * LegendrianSubmanifold: n-dimensional submanifold L ⊂ M where α|_L = 0
     * 
     * Generated by a function A(x):
     *   u = A(x)
     *   p_a = ∂_a A(x)
     * 
     * Hamilton-Jacobi constraint: H(x, A(x), ∂A(x)) = 0
     */
    class LegendrianSubmanifold {
        /**
         * @param {ContactManifold} manifold
         * @param {Function} generatingFunc - A(x) → number
         * @param {Function} [dA] - Gradient ∂_a A(x)
         */
        constructor(manifold, generatingFunc, dA = null) {
            this.manifold = manifold;
            this.A = generatingFunc;
            this._dA = dA;
        }

        /**
         * Numerical gradient of generating function
         */
        gradient(x, h = 1e-7) {
            if (this._dA) {
                return this._dA(x);
            }
            
            const grad = {};
            for (const coord of this.manifold.baseCoords) {
                const xPlus = { ...x, [coord]: x[coord] + h };
                const xMinus = { ...x, [coord]: x[coord] - h };
                grad[coord] = (this.A(xPlus) - this.A(xMinus)) / (2 * h);
            }
            return grad;
        }

        /**
         * Lift base point to contact manifold via Legendrian embedding
         * @param {Object} x - Base coordinates { x^a: value }
         */
        lift(x) {
            const coords = { ...x };
            
            // u = A(x)
            coords[this.manifold.fiberCoord] = this.A(x);
            
            // p_a = ∂_a A(x)
            const grad = this.gradient(x);
            for (let i = 0; i < this.manifold.n; i++) {
                const xCoord = this.manifold.baseCoords[i];
                const pCoord = this.manifold.momentaCoords[i];
                coords[pCoord] = grad[xCoord];
            }
            
            return this.manifold.point(coords);
        }

        /**
         * Verify Legendrian condition: α|_L = 0
         * For any tangent vector to L, α(v) should vanish
         */
        verifyLegendrianCondition(x) {
            // For generating function construction, this is automatic
            // α = du - p_a dx^a = dA - ∂_a A dx^a = 0 on L
            return true;
        }

        /**
         * Hamilton-Jacobi constraint: H(x, A(x), ∂A(x)) = 0
         * @param {ContactHamiltonian} hamiltonian
         */
        hamiltonJacobiResidual(x, hamiltonian) {
            const pt = this.lift(x);
            return hamiltonian.evaluate(pt);
        }

        /**
         * Sample the Legendrian submanifold
         * @param {Function} sampler - () => base coordinates x
         * @param {number} n - Number of samples
         */
        sample(sampler, n) {
            const points = [];
            for (let i = 0; i < n; i++) {
                const x = sampler();
                points.push(this.lift(x));
            }
            return points;
        }
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
                const f = 1 - 2*M/r;
                const sinTheta2 = sin(theta) ** 2;
                return [
                    [f, 0, 0, 0],
                    [0, -1/f, 0, 0],
                    [0, 0, -r*r, 0],
                    [0, 0, 0, -r*r*sinTheta2]
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
                dHdx.push((this.evaluate(xPlus, p) - this.evaluate(xMinus, p)) / (2*h));
            }
            
            // ∂H/∂p_μ
            for (let mu = 0; mu < 4; mu++) {
                const pPlus = [...p]; pPlus[mu] += h;
                const pMinus = [...p]; pMinus[mu] -= h;
                dHdp.push((this.evaluate(x, pPlus) - this.evaluate(x, pMinus)) / (2*h));
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
                
                x = x.map((xi, mu) => xi + dtau/6 * (k1.xDot[mu] + 2*k2.xDot[mu] + 2*k3.xDot[mu] + k4.xDot[mu]));
                p = p.map((pi, mu) => pi + dtau/6 * (k1.pDot[mu] + 2*k2.pDot[mu] + 2*k3.pDot[mu] + k4.pDot[mu]));
                
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
                    result += 0.5 * mass * (k1*k1 + k2*k2 + k3*k3);
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
                const kSq = k1*k1 + k2*k2 + k3*k3;
                
                if (mass === 0) {
                    return omega - c * sqrt(kSq);
                } else {
                    return omega - sqrt(c*c*kSq + mass*mass*c*c*c*c);
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
        
        // Utilities
        DifferentialForm,
        summaryTable,
        
        // Constants
        EPSILON,
        
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
