/**
 * Example: GMET + Bianconi Entropic Gravity
 *
 * Demonstrates the full integration of:
 * - GMET (Grand/Holographic Contact Manifolds)
 * - Bianconi's Entropic Gravity (two-metric system, relative entropy)
 *
 * The operational recipe:
 * 1. Set up the 13D Grand Manifold M₁₃ = J¹(Q₆)
 * 2. Define the matter-induced metric G via topological forms (φ, A, B)
 * 3. Set the Hamiltonian to minimize relative entropy S(G||g)
 * 4. The resulting flow yields modified GR with emergent Λ
 */

const GMET = require('../src/index.js');
const EntropicGravity = require('../src/entropic-gravity.js');

const { sqrt, exp, sin, cos, PI } = Math;

console.log('='.repeat(70));
console.log('GMET + BIANCONI ENTROPIC GRAVITY DEMONSTRATION');
console.log('='.repeat(70));

// ============================================================================
// Step 1: Set up the Grand Contact Manifold M₁₃
// ============================================================================

console.log('\n1. GRAND CONTACT MANIFOLD M₁₃');
console.log('-'.repeat(50));

const M13 = new GMET.GrandContactManifold();

console.log(`Manifold: ${M13.toString()}`);
console.log(`Dimension: ${M13.dimension}`);
console.log(`Base coordinates: ${M13.baseCoords.join(', ')}`);
console.log(`Momenta: ${M13.momentaCoords.join(', ')}`);
console.log(`Fiber: ${M13.fiberCoord}`);
console.log(`Contact form: ${M13.contactFormSymbolic()}`);

// ============================================================================
// Step 2: Define the Two-Metric System (g, G)
// ============================================================================

console.log('\n2. TWO-METRIC SYSTEM (g, G)');
console.log('-'.repeat(50));

// Spacetime metric g: Schwarzschild black hole
const M_bh = 1.0;  // Black hole mass
const g = EntropicGravity.StandardMetrics.schwarzschild(M_bh);
console.log(`Spacetime metric g: Schwarzschild (M = ${M_bh})`);

// Matter-induced metric G: Scalar field + electromagnetic field
const G = new EntropicGravity.MatterInducedMetric({
    // Scalar field: dilaton-like, falling off with distance
    scalarField: x => {
        const r = sqrt(x[1] * x[1] + x[2] * x[2] + x[3] * x[3]);
        return 0.05 * exp(-r / 5);  // Decays with scale r_s = 5
    },

    // Electromagnetic potential: Coulomb-like
    vectorPotential: x => {
        const r = sqrt(x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + 0.1);
        return [0.1 / r, 0, 0, 0];  // A_0 = φ_em = q/r
    },

    // Ricci tensor: from curvature (simplified)
    ricciTensor: x => {
        // For Schwarzschild vacuum, R_μν = 0
        // But with matter, we get corrections
        const r = sqrt(x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + 0.1);
        const correction = 0.001 / (r * r);
        return [
            [correction, 0, 0, 0],
            [0, -correction, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]
        ];
    }
});

console.log('Matter-induced metric G: scalar field + EM field');

// Create two-metric system
const twoMetricSystem = new EntropicGravity.TwoMetricSystem(g, G);

// Evaluate at a sample point
const testPoint = [0, 5, PI / 2, 0];  // (t, r, θ, φ)
const g_sample = twoMetricSystem.spacetimeMetric(testPoint);
const G_sample = twoMetricSystem.matterMetric(testPoint);

console.log(`\nAt r = 5 (outside horizon):`);
console.log(`  g_00 = ${g_sample[0][0].toFixed(4)} (should be 0.6)`);
console.log(`  g_11 = ${g_sample[1][1].toFixed(4)} (should be -1.667)`);
console.log(`  Tr(G) = ${EntropicGravity.utils.mat4Trace(G_sample).toFixed(6)}`);

// ============================================================================
// Step 3: Compute Relative Entropy and Emergent Λ
// ============================================================================

console.log('\n3. RELATIVE ENTROPY ACTION S(G||g)');
console.log('-'.repeat(50));

const entropyAction = new EntropicGravity.RelativeEntropyAction(twoMetricSystem);
const lambdaCalculator = new EntropicGravity.EmergentCosmologicalConstant(twoMetricSystem);

// Local entropy density at various radii
console.log('\nLocal entropy density vs radius:');
for (const r of [3, 5, 10, 20]) {
    const x = [0, r, PI / 2, 0];
    const density = entropyAction.localDensity(x);
    const Lambda = lambdaCalculator.localValue(x);
    console.log(`  r = ${r.toString().padStart(2)}: S_local = ${density.toFixed(6)}, Λ_G = ${Lambda.toFixed(6)}`);
}

// Classify spacetime type
const classification = lambdaCalculator.classify([0, 10, PI / 2, 0]);
console.log(`\nSpacetime classification at r=10: ${classification}`);

// ============================================================================
// Step 4: Create Entropic Gravity Hamiltonian
// ============================================================================

console.log('\n4. ENTROPIC GRAVITY HAMILTONIAN');
console.log('-'.repeat(50));

const H = new EntropicGravity.EntropicGravityHamiltonian(M13, twoMetricSystem, {
    mass: 0.1,                // Test particle mass
    entropicCoupling: 0.001,  // Coupling strength α (small for numerical stability)
    useGA: false              // Explicitly use Tensor mode for this tutorial example
});

console.log('Hamiltonian: H = H_geo + α·S(G||g)');
console.log(`  mass = 0.1`);
console.log(`  entropic coupling α = 0.001`);

// ============================================================================
// Step 5: Set Initial Conditions
// ============================================================================

console.log('\n5. INITIAL CONDITIONS');
console.log('-'.repeat(50));

// Particle starting at r = 15, slightly inward trajectory
const initialPoint = M13.physicalPoint(
    15,     // q1 = r
    PI / 2,   // q2 = θ (equatorial)
    0,      // q3 = φ
    0,      // t
    0,      // ℓ = log(scale)
    1.0,    // S (initial entropy)
    -0.02,  // k1 (inward radial momentum)
    0,      // k2 (no θ momentum)
    0.005,  // k3 (small angular momentum)
    1.0,    // ω (energy)
    0,      // Δ (dilatation)
    1.0,    // T (temperature)
    0       // A (action)
);

console.log('Initial state:');
console.log(`  r = ${initialPoint.get('q1').toFixed(4)}`);
console.log(`  θ = ${(initialPoint.get('q2') * 180 / PI).toFixed(2)}°`);
console.log(`  k_r = ${initialPoint.get('k1').toFixed(4)}`);
console.log(`  L = ${initialPoint.get('k3').toFixed(4)} (angular momentum)`);
console.log(`  S = ${initialPoint.get('S').toFixed(4)} (entropy)`);

// Evaluate initial Hamiltonian
const H_initial = H.evaluate(initialPoint.coords);
const H_geo_initial = H.geometricPart(initialPoint.coords);
const H_ent_initial = H.entropicPart(initialPoint.coords);

console.log('\nInitial Hamiltonian:');
console.log(`  H_total = ${H_initial.toFixed(6)}`);
console.log(`  H_geo = ${H_geo_initial.toFixed(6)}`);
console.log(`  H_entropy = ${H_ent_initial.toFixed(6)}`);

// ============================================================================
// Step 6: Integrate the Flow
// ============================================================================

console.log('\n6. CONTACT HAMILTONIAN FLOW');
console.log('-'.repeat(50));

const dt = 0.5;
const steps = 100;

console.log(`Integrating with dt = ${dt}, steps = ${steps}...`);

const trajectory = H.flow(initialPoint, dt, steps);

// ============================================================================
// Step 7: Analyze Results
// ============================================================================

console.log('\n7. TRAJECTORY ANALYSIS');
console.log('-'.repeat(50));

// Extract key quantities along trajectory
const data = trajectory.map((pt, i) => ({
    step: i,
    tau: i * dt,
    r: pt.get('q1'),
    phi: pt.get('q3'),
    S: pt.get('S'),
    A: pt.get('A'),
    T: pt.get('T'),
    H: H.evaluate(pt.coords)
}));

// Print summary at key steps
console.log('\nTrajectory summary:');
console.log('  step |    τ    |    r    |    φ    |    S    |    A    |    H');
console.log('  ' + '-'.repeat(65));

for (const i of [0, 25, 50, 75, steps]) {
    const d = data[i];
    console.log(`  ${d.step.toString().padStart(4)} | ${d.tau.toFixed(3).padStart(7)} | ` +
        `${d.r.toFixed(3).padStart(7)} | ${d.phi.toFixed(3).padStart(7)} | ` +
        `${d.S.toFixed(3).padStart(7)} | ${d.A.toFixed(3).padStart(7)} | ` +
        `${d.H.toFixed(4).padStart(7)}`);
}

// Compute changes
const r_change = data[steps].r - data[0].r;
const phi_change = data[steps].phi - data[0].phi;
const S_change = data[steps].S - data[0].S;
const A_change = data[steps].A - data[0].A;
const H_change = data[steps].H - data[0].H;

console.log('\nTotal changes:');
console.log(`  Δr = ${r_change.toFixed(4)}`);
console.log(`  Δφ = ${phi_change.toFixed(4)} rad = ${(phi_change * 180 / PI).toFixed(2)}°`);
console.log(`  ΔS = ${S_change.toFixed(4)} (entropy change)`);
console.log(`  ΔA = ${A_change.toFixed(4)} (action accumulated)`);
console.log(`  ΔH = ${H_change.toFixed(6)} (Hamiltonian drift)`);

// ============================================================================
// Step 8: Physical Interpretation
// ============================================================================

console.log('\n8. PHYSICAL INTERPRETATION');
console.log('-'.repeat(50));

console.log(`
The simulation shows a test particle in Schwarzschild spacetime with
Bianconi's entropic gravity corrections:

1. GEOMETRIC SECTOR (H_geo):
   - Mass-shell constraint: ½g^μν p_μ p_ν - ½m² = 0
   - Particle follows approximate geodesic

2. ENTROPIC SECTOR (H_entropy = α·S(G||g)):
   - Relative entropy creates additional "entropic force"
   - Matter distribution (scalar + EM fields) affects trajectory
   - Emergent Λ contributes to dynamics

3. KEY OBSERVATIONS:
   - Radial motion: ${r_change > 0 ? 'outward' : 'inward'} drift (Δr = ${r_change.toFixed(4)})
   - Angular motion: ${(phi_change * 180 / PI).toFixed(2)}° precession
   - Entropy: ${S_change > 0 ? 'increasing' : 'decreasing'} (thermodynamic arrow)
   - Action: accumulated phase A = ${A_change.toFixed(4)}

4. EMERGENT COSMOLOGICAL CONSTANT:
   - Λ_G varies with position (field-dependent)
   - Classification: ${lambdaCalculator.classify([0, data[steps].r, PI / 2, 0])}
`);

// ============================================================================
// Comparison: Standard GR vs Entropic Gravity
// ============================================================================

console.log('\n9. COMPARISON: STANDARD GR vs ENTROPIC GRAVITY');
console.log('-'.repeat(50));

// Run without entropic coupling
const H_pure_geo = new EntropicGravity.EntropicGravityHamiltonian(M13, twoMetricSystem, {
    mass: 0.1,
    entropicCoupling: 0  // No entropic correction
});

const trajectory_geo = H_pure_geo.flow(initialPoint, dt, steps);
const final_r_geo = trajectory_geo[steps].get('q1');
const final_r_ent = trajectory[steps].get('q1');

console.log(`Final radial position:`);
console.log(`  Standard GR (α=0):         r = ${final_r_geo.toFixed(4)}`);
console.log(`  Entropic gravity (α=0.001): r = ${final_r_ent.toFixed(4)}`);
console.log(`  Difference: Δr = ${(final_r_ent - final_r_geo).toFixed(6)}`);

const correction_percent = Math.abs(final_r_ent - final_r_geo) / Math.abs(final_r_geo - data[0].r + 0.001) * 100;
console.log(`  Entropic correction: ~${correction_percent.toFixed(2)}% of total motion`);

// ============================================================================
// Summary
// ============================================================================

console.log('\n' + '='.repeat(70));
console.log('SUMMARY: GMET + BIANCONI FRAMEWORK');
console.log('='.repeat(70));

console.log(`
┌────────────────────────────────────────────────────────────────────┐
│  COMPONENT              │  ROLE                                   │
├────────────────────────────────────────────────────────────────────┤
│  GrandContactManifold   │  13D phase space (kinematics)          │
│  M₁₃ = J¹(Q₆)          │  Coordinates: (q,t,ℓ,S) + momenta + A   │
├────────────────────────────────────────────────────────────────────┤
│  SpacetimeMetric g      │  Background geometry (Schwarzschild)    │
│                         │  Determines geodesic structure          │
├────────────────────────────────────────────────────────────────────┤
│  MatterInducedMetric G  │  Matter fields (φ, A, B forms)          │
│                         │  Creates "dressed" geometry             │
├────────────────────────────────────────────────────────────────────┤
│  RelativeEntropyAction  │  S(G||g) = information distance         │
│                         │  Generates gravitational dynamics       │
├────────────────────────────────────────────────────────────────────┤
│  EmergentΛ              │  Λ_G from (g,G) interplay               │
│                         │  Dynamically determined, not fixed      │
├────────────────────────────────────────────────────────────────────┤
│  EntropicHamiltonian    │  H = H_geo + α·S(G||g)                  │
│                         │  Contact dynamics on M₁₃                │
└────────────────────────────────────────────────────────────────────┘

The framework cleanly separates:
  - GMET: Universal kinematic container (phase space + thermodynamics)
  - Bianconi: Specific dynamical content (entropic gravity engine)

Together they yield Modified Einstein Equations with emergent Λ.
`);
