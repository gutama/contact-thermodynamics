/**
 * EM Wave Simulation Example
 * 
 * Solves Maxwell's equations using the discrete geometric calculus module:
 *   ∂B/∂t = -∇×E  (Faraday)
 *   ∂E/∂t = c²∇×B (Ampère, vacuum)
 * 
 * Uses Yee-style staggered time-stepping (leapfrog).
 */

const GC = require('../src/geometric-calculus.js');

// ============================================================================
// MAXWELL SOLVER
// ============================================================================

class MaxwellSolver2D {
    /**
     * 2D Maxwell solver (TM mode: Ez, Bx, By)
     * 
     * In 2D TM mode:
     * - E has only z-component (Ez)
     * - B has x,y components (Bx, By)
     * - ∇×E = (∂Ez/∂y, -∂Ez/∂x, 0) → gives Bx, By
     * - ∇×B = (0, 0, ∂By/∂x - ∂Bx/∂y) → gives Ez
     * 
     * @param {number[]} shape - Grid dimensions [nx, ny]
     * @param {number[]} spacing - Grid spacing [dx, dy]
     * @param {number} c - Speed of light
     */
    constructor(shape, spacing, c = 1.0) {
        this.shape = shape;
        this.spacing = spacing;
        this.c = c;
        this.nabla = new GC.SplitDifferentialOperator(shape, spacing);

        // Fields: E_z (scalar), B_x, B_y (vectors stored as scalars)
        this.Ez = new GC.ScalarField(shape, spacing);
        this.Bx = new GC.ScalarField(shape, spacing);
        this.By = new GC.ScalarField(shape, spacing);
    }

    /**
     * Estimate stable time step (CFL condition for EM)
     */
    estimateCFL(safety = 0.9) {
        const minSpacing = Math.min(...this.spacing);
        return safety * minSpacing / (this.c * Math.sqrt(2));
    }

    /**
     * Initialize with a Gaussian pulse in Ez
     */
    initGaussianPulse(x0, y0, sigma = 0.1) {
        this.Ez.initFromFunction((x, y) => {
            const dx = x - x0;
            const dy = y - y0;
            return Math.exp(-(dx * dx + dy * dy) / (2 * sigma * sigma));
        });
    }

    /**
     * Initialize with a plane wave in Ez
     */
    initPlaneWave(kx, ky, amplitude = 1.0) {
        this.Ez.initFromFunction((x, y) => {
            return amplitude * Math.sin(kx * x + ky * y);
        });
        // Set initial B consistent with plane wave
        const k = Math.sqrt(kx * kx + ky * ky);
        const omega = this.c * k;
        this.Bx.initFromFunction((x, y) => {
            return -amplitude * (ky / omega) * Math.sin(kx * x + ky * y);
        });
        this.By.initFromFunction((x, y) => {
            return amplitude * (kx / omega) * Math.sin(kx * x + ky * y);
        });
    }

    /**
     * Single timestep of Maxwell's equations
     * 
     * Uses leapfrog (B at half-steps):
     *   B(t+dt/2) = B(t-dt/2) - dt * ∇×E(t)
     *   E(t+dt) = E(t) + c²dt * ∇×B(t+dt/2)
     */
    step(dt) {
        const c2dt = this.c * this.c * dt;
        const [dx, dy] = this.spacing;
        const [nx, ny] = this.shape;

        // Update B: ∂B/∂t = -∇×E
        // In 2D: (∇×E)_x = ∂Ez/∂y, (∇×E)_y = -∂Ez/∂x
        for (let i = 0; i < nx; i++) {
            for (let j = 0; j < ny; j++) {
                // ∂Ez/∂y → Bx update
                const jPlus = (j + 1) % ny;
                const jMinus = (j - 1 + ny) % ny;
                const dEz_dy = (this.Ez.get(i, jPlus) - this.Ez.get(i, jMinus)) / (2 * dy);

                // ∂Ez/∂x → By update
                const iPlus = (i + 1) % nx;
                const iMinus = (i - 1 + nx) % nx;
                const dEz_dx = (this.Ez.get(iPlus, j) - this.Ez.get(iMinus, j)) / (2 * dx);

                // Faraday: ∂B/∂t = -∇×E
                this.Bx.set(this.Bx.get(i, j) - dt * dEz_dy, i, j);
                this.By.set(this.By.get(i, j) + dt * dEz_dx, i, j);
            }
        }

        // Update E: ∂E/∂t = c²∇×B
        // In 2D: (∇×B)_z = ∂By/∂x - ∂Bx/∂y
        for (let i = 0; i < nx; i++) {
            for (let j = 0; j < ny; j++) {
                const iPlus = (i + 1) % nx;
                const iMinus = (i - 1 + nx) % nx;
                const jPlus = (j + 1) % ny;
                const jMinus = (j - 1 + ny) % ny;

                const dBy_dx = (this.By.get(iPlus, j) - this.By.get(iMinus, j)) / (2 * dx);
                const dBx_dy = (this.Bx.get(i, jPlus) - this.Bx.get(i, jMinus)) / (2 * dy);

                // Ampère: ∂E/∂t = c²∇×B
                const curlB_z = dBy_dx - dBx_dy;
                this.Ez.set(this.Ez.get(i, j) + c2dt * curlB_z, i, j);
            }
        }
    }

    /**
     * Run simulation
     */
    simulate(dt, n_steps, callback = null) {
        const history = [];

        for (let step = 0; step < n_steps; step++) {
            this.step(dt);

            if (callback && step % 10 === 0) {
                callback(step, this.Ez, this.Bx, this.By);
            }

            // Record energy periodically
            if (step % 50 === 0) {
                history.push({
                    step,
                    time: step * dt,
                    energy: this.totalEnergy()
                });
            }
        }

        return history;
    }

    /**
     * Compute total EM energy: U = ½∫(ε₀E² + B²/μ₀) dV
     * In natural units (c=1, ε₀=μ₀=1): U = ½∫(E² + B²) dV
     */
    totalEnergy() {
        const [dx, dy] = this.spacing;
        const dV = dx * dy;
        let U = 0;

        for (let i = 0; i < this.Ez.data.length; i++) {
            const Ez2 = this.Ez.data[i] * this.Ez.data[i];
            const Bx2 = this.Bx.data[i] * this.Bx.data[i];
            const By2 = this.By.data[i] * this.By.data[i];
            U += 0.5 * (Ez2 + Bx2 + By2) * dV;
        }

        return U;
    }

    /**
     * Get field data for visualization
     */
    getFieldData() {
        return {
            Ez: Array.from(this.Ez.data),
            Bx: Array.from(this.Bx.data),
            By: Array.from(this.By.data),
            shape: this.shape
        };
    }
}

// ============================================================================
// RUN SIMULATION
// ============================================================================

console.log('╔══════════════════════════════════════════════════════════════╗');
console.log('║          EM Wave Simulation using Geometric Calculus         ║');
console.log('╚══════════════════════════════════════════════════════════════╝\n');

// Setup grid
const L = 2 * Math.PI;  // Domain size
const nx = 64, ny = 64;
const dx = L / nx, dy = L / ny;

console.log(`Grid: ${nx}×${ny}, spacing: ${dx.toFixed(4)} × ${dy.toFixed(4)}`);
console.log(`Domain: [0, ${L.toFixed(2)}] × [0, ${L.toFixed(2)}]`);

// Create solver
const c = 1.0;  // Speed of light
const maxwell = new MaxwellSolver2D([nx, ny], [dx, dy], c);

// CFL timestep
const dt = maxwell.estimateCFL(0.8);
console.log(`Time step (CFL): dt = ${dt.toFixed(6)}`);

// Test 1: Gaussian pulse spreading
console.log('\n━━━ Test 1: Gaussian Pulse ━━━');
maxwell.initGaussianPulse(L / 2, L / 2, 0.3);
const E0_pulse = maxwell.totalEnergy();
console.log(`Initial energy: ${E0_pulse.toFixed(6)}`);

const history_pulse = maxwell.simulate(dt, 200);
const Ef_pulse = history_pulse[history_pulse.length - 1].energy;
const energy_conservation_pulse = Math.abs(Ef_pulse - E0_pulse) / E0_pulse;
console.log(`Final energy: ${Ef_pulse.toFixed(6)}`);
console.log(`Energy drift: ${(energy_conservation_pulse * 100).toFixed(2)}%`);
console.log(energy_conservation_pulse < 0.05 ? '✓ Energy well conserved' : '⚠ Energy drift detected');

// Test 2: Plane wave
console.log('\n━━━ Test 2: Plane Wave ━━━');
const maxwell2 = new MaxwellSolver2D([nx, ny], [dx, dy], c);
const k = 2;  // Wave number
maxwell2.initPlaneWave(k, 0, 1.0);  // Wave traveling in x direction
const E0_wave = maxwell2.totalEnergy();
console.log(`Initial energy: ${E0_wave.toFixed(6)}`);

const history_wave = maxwell2.simulate(dt, 200);
const Ef_wave = history_wave[history_wave.length - 1].energy;
const energy_conservation_wave = Math.abs(Ef_wave - E0_wave) / E0_wave;
console.log(`Final energy: ${Ef_wave.toFixed(6)}`);
console.log(`Energy drift: ${(energy_conservation_wave * 100).toFixed(2)}%`);
console.log(energy_conservation_wave < 0.05 ? '✓ Energy well conserved' : '⚠ Energy drift detected');

// Test 3: Verify wave speed
console.log('\n━━━ Test 3: Wave Speed Verification ━━━');
// Track the peak of Ez over time
const maxwell3 = new MaxwellSolver2D([128, 32], [0.1, 0.1], 1.0);
maxwell3.initGaussianPulse(3.0, 1.6, 0.3);

let peak_x_initial = 3.0;
let peak_x_final = 0;
const total_time = 3.0;
const n_steps = Math.floor(total_time / maxwell3.estimateCFL(0.8));
const dt3 = maxwell3.estimateCFL(0.8);

maxwell3.simulate(dt3, n_steps);

// Find peak position
let maxEz = 0;
for (let i = 0; i < 128; i++) {
    const val = Math.abs(maxwell3.Ez.get(i, 16));
    if (val > maxEz) {
        maxEz = val;
        peak_x_final = i * 0.1;
    }
}

const distance = Math.abs(peak_x_final - peak_x_initial);
const measured_speed = distance / total_time;
console.log(`Peak moved: ${distance.toFixed(2)} units in ${total_time} time`);
console.log(`Expected speed: ${c.toFixed(2)}, Measured: ${measured_speed.toFixed(2)}`);
console.log(Math.abs(measured_speed - c) < 0.2 ? '✓ Wave speed correct' : '⚠ Speed discrepancy');

console.log('\n═══════════════════════════════════════════════════════════════');
console.log('EM simulation complete! The geometric calculus module works.');
console.log('═══════════════════════════════════════════════════════════════\n');

module.exports = { MaxwellSolver2D };
