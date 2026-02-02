/**
 * Contact Hamiltonian Module
 * 
 * Contact Hamiltonian dynamics for dissipative systems.
 * Extracted from index.js for modular organization.
 * 
 * @module contact/hamiltonian
 * @license MIT
 */

(function (global) {
    'use strict';

    // ============================================================================
    // CONTACT HAMILTONIAN DYNAMICS
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
                    combined[c] = (k1[c] + 2 * k2[c] + 2 * k3[c] + k4[c]) / 6;
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
    // EXPORTS
    // ============================================================================

    const Hamiltonian = {
        ContactHamiltonian
    };

    // Export for different module systems
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = Hamiltonian;
    } else if (typeof define === 'function' && define.amd) {
        define([], () => Hamiltonian);
    } else {
        global.ContactHamiltonian = Hamiltonian;
    }

})(typeof window !== 'undefined' ? window : global);
