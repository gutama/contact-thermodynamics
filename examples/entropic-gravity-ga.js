/**
 * Example: Entropic Gravity with Geometric Algebra (Cl(1,3))
 * 
 * Demonstrates the Coordinate-Free Entropic Gravity framework.
 * 
 * We compare:
 * 1. Standard Geometric Mechanics (Tensor-based)
 * 2. Geometric Algebra Mechanics (Connection Bivector-based)
 * 
 * The physical system is a particle falling into a Schwarzschild black hole
 * perturbed by entropic forces.
 */

const EntropicGravity = require('../src/entropic-gravity.js');
const GMET = require('../src/index.js');
const { abs, sqrt, PI, exp } = Math;

console.log('='.repeat(70));
console.log('ENTROPIC GRAVITY: TENSOR vs GEOMETRIC ALGEBRA');
console.log('='.repeat(70));

// ============================================================================
// 1. Setup Framework
// ============================================================================

// 13D Grand Manifold
const M13 = new GMET.GrandContactManifold();

// Metrics
const M_bh = 1.0;
const g = EntropicGravity.StandardMetrics.schwarzschild(M_bh);

// Matter field (scalar cloud)
const G = new EntropicGravity.MatterInducedMetric({
    scalarField: x => {
        const r = x[1];
        return 0.1 * exp(-r / 5.0);
    }
});

const system = new EntropicGravity.TwoMetricSystem(g, G);

// ============================================================================
// 2. Hamiltonian Construction
// ============================================================================

// Standard Tensor Hamiltonian
console.log('\nInitializing Tensor Hamiltonian...');
const H_tensor = new EntropicGravity.EntropicGravityHamiltonian(M13, system, {
    mass: 1.0,
    entropicCoupling: 0.01,
    useGA: false
});

// GA Hamiltonian
// Internally uses riemannian-spacetime.js to compute Connection Bivector
console.log('Initializing GA Hamiltonian (Cl(1,3))...');
const H_ga = new EntropicGravity.EntropicGravityHamiltonian(M13, system, {
    mass: 1.0,
    entropicCoupling: 0.01,
    useGA: true
});

// ============================================================================
// 3. Compare Force Evaluation
// ============================================================================

// Point: r = 6.0, theta = pi/2 (circular orbitish)
const p0 = M13.physicalPoint(
    6.0, PI / 2, 0, // q
    0, 0, 1.0,    // t, l, S
    0, 0, 4.0,    // p (angular momentum 4)
    1.0, 0, 1.0,  // omega, Delta, T
    0
);

console.log(`\nTest Point: r=${p0.get('q1')}, L=${p0.get('k3')}`);

const val_tensor = H_tensor.evaluate(p0.coords);
const val_ga = H_ga.evaluate(p0.coords);

console.log(`Hamiltonian Value (Tensor): ${val_tensor.toFixed(8)}`);
console.log(`Hamiltonian Value (GA):     ${val_ga.toFixed(8)}`);

// Check Entropic Part specifically
const S_tensor = H_tensor.entropicPart(p0.coords);
const S_ga = H_ga.entropicPartGA(p0.coords);

console.log(`\nEntropic Part S (Tensor ~ Tr(G g^-1)): ${S_tensor.toFixed(8)}`);
console.log(`Entropic Part S (GA ~ <ω ω†>):       ${S_ga.toFixed(8)}`);

// Note: They are defined differently!
// - Tensor S = Information distance S(G||g) (Relative Entropy)
// - GA S = Field Strength Magnitude |ω|^2 (Curvature intensity)
// We expect them to scale differently but both drive the system.

// ============================================================================
// 4. Dynamics Comparison
// ============================================================================

const dt = 0.1;
const steps = 50;

console.log(`\nIntegrating Flow (dt=${dt}, steps=${steps})...`);

const traj_tensor = H_tensor.flow(p0, dt, steps);
const traj_ga = H_ga.flow(p0, dt, steps);

const final_tensor = traj_tensor[steps];
const final_ga = traj_ga[steps];

console.log('\nFinal State Comparison:');
console.log('                 r       |      phi      |      pr');
console.log('-'.repeat(60));
console.log(`TENSOR:   ${final_tensor.get('q1').toFixed(6)} | ${final_tensor.get('q3').toFixed(6)} | ${final_tensor.get('k1').toFixed(6)}`);
console.log(`GA:       ${final_ga.get('q1').toFixed(6)} | ${final_ga.get('q3').toFixed(6)} | ${final_ga.get('k1').toFixed(6)}`);

const diff_r = abs(final_tensor.get('q1') - final_ga.get('q1'));
console.log(`\nDifference in radial drift: ${diff_r.toFixed(6)}`);

if (diff_r > 0.01) {
    console.log("-> Dynamics diverge significantly (expected, different potentials).");
} else {
    console.log("-> Dynamics are similar.");
}

console.log('\nConclusion: GA formalism provides an intrinsic alternative to the tensor formulation.');
console.log('The GA force is driven by the field strength magnitude |ω|^2,');
console.log('while the Tensor force is critical relative entropy S(G||g).');
