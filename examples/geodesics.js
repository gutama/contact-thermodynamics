/**
 * Geodesics and General Relativity Extension
 * 
 * Demonstrates:
 * - Spacetime metrics (Minkowski, Schwarzschild, FLRW)
 * - Relativistic Hamiltonian (mass-shell constraint)
 * - Geodesic integration
 * - Conserved quantities
 */

const CT = require('../src/index.js');

console.log('═══════════════════════════════════════════════════════════════');
console.log('  Geodesics & General Relativity Extension');
console.log('═══════════════════════════════════════════════════════════════\n');

// ============================================================================
// 1. Minkowski Spacetime (Flat)
// ============================================================================

console.log('1. MINKOWSKI SPACETIME\n');

const mink = CT.SpacetimeMetric.minkowski();

// Get metric at origin
const g_mink = mink.covariant([0, 0, 0, 0]);
console.log('Minkowski metric η_μν = diag(+1, -1, -1, -1):');
console.log('   g_tt =', g_mink[0][0]);
console.log('   g_xx =', g_mink[1][1]);
console.log('   g_yy =', g_mink[2][2]);
console.log('   g_zz =', g_mink[3][3]);
console.log();

// Inverse metric
const gInv_mink = mink.contravariant([0, 0, 0, 0]);
console.log('Inverse metric η^μν:');
console.log('   g^tt =', gInv_mink[0][0]);
console.log('   g^xx =', gInv_mink[1][1]);
console.log();

// ============================================================================
// 2. Schwarzschild Metric (Black Hole)
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('2. SCHWARZSCHILD METRIC\n');

const M = 1;  // Black hole mass (geometric units: G = c = 1)
const schw = CT.SpacetimeMetric.schwarzschild(M);

console.log('Schwarzschild metric (M = 1, units: G = c = 1):');
console.log('   ds² = (1-2M/r)dt² - (1-2M/r)⁻¹dr² - r²dΩ²');
console.log();

// Evaluate at different radii
const radii = [3, 5, 10, 100];
console.log('Metric components at various radii (θ = π/2):');
console.log('   r       g_tt      g_rr');
console.log('   ─────────────────────────');

for (const r of radii) {
    const g = schw.covariant([0, r, Math.PI/2, 0]);
    console.log(`   ${r.toString().padEnd(6)}  ${g[0][0].toFixed(4).padStart(8)}  ${g[1][1].toFixed(4).padStart(8)}`);
}
console.log();

console.log('Event horizon at r = 2M =', 2*M);
console.log('Photon sphere at r = 3M =', 3*M);
console.log();

// ============================================================================
// 3. FLRW Metric (Cosmology)
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('3. FLRW METRIC (COSMOLOGY)\n');

// Exponential expansion (de Sitter-like)
const H0 = 0.1;  // Hubble parameter
const a = t => Math.exp(H0 * t);

const flrw = CT.SpacetimeMetric.flrw(a, 0);  // k = 0 (flat universe)

console.log('FLRW metric with a(t) = exp(H₀t), H₀ =', H0);
console.log('   ds² = dt² - a(t)²(dχ² + χ²dΩ²)');
console.log();

const times = [0, 5, 10, 20];
console.log('Scale factor evolution:');
console.log('   t       a(t)');
console.log('   ───────────────');

for (const t of times) {
    console.log(`   ${t.toString().padEnd(6)}  ${a(t).toFixed(4)}`);
}
console.log();

// ============================================================================
// 4. Relativistic Hamiltonian
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('4. RELATIVISTIC HAMILTONIAN (MASS SHELL)\n');

console.log('Mass-shell constraint:');
console.log('   H = ½g^μν p_μ p_ν - ½m² = 0');
console.log();

const relH = new CT.RelativisticHamiltonian(mink, 1);  // mass = 1

// Test on-shell condition in Minkowski
// For m = 1, need p₀² - |p|² = m² = 1
// Try p = (√2, 1, 0, 0): p² = 2 - 1 = 1 ✓

const x = [0, 0, 0, 0];
const p_onshell = [Math.sqrt(2), 1, 0, 0];

console.log('Test in Minkowski (m = 1):');
console.log('   p = [√2, 1, 0, 0]');
console.log('   p₀² - |p|² = 2 - 1 = 1 = m²');
console.log('   H =', relH.evaluate(x, p_onshell).toFixed(6), '(expect 0)');
console.log();

// Off-shell
const p_offshell = [2, 1, 0, 0];
console.log('Off-shell test:');
console.log('   p = [2, 1, 0, 0]');
console.log('   H =', relH.evaluate(x, p_offshell).toFixed(6), '(expect ≠ 0)');
console.log();

// ============================================================================
// 5. Geodesic in Minkowski (Straight Line)
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('5. GEODESIC IN MINKOWSKI SPACETIME\n');

const relH_mink = new CT.RelativisticHamiltonian(mink, 1);

// Initial conditions: particle at origin, moving in x-direction
// For on-shell: p_t = √(m² + |p|²) with our convention
const px = 0.5;
const pt = Math.sqrt(1 + px*px);

const x0_mink = [0, 0, 0, 0];
const p0_mink = [pt, px, 0, 0];

console.log('Initial conditions:');
console.log('   x = [0, 0, 0, 0]');
console.log('   p = [' + p0_mink.map(v => v.toFixed(3)).join(', ') + ']');
console.log('   H =', relH_mink.evaluate(x0_mink, p0_mink).toFixed(6));
console.log();

// Integrate
const geodesic_mink = relH_mink.integrateGeodesic(x0_mink, p0_mink, 0.5, 20);

console.log('Geodesic trajectory:');
console.log('   τ       t         x         p_t       p_x');
console.log('   ─────────────────────────────────────────────────');

for (let i = 0; i <= 20; i += 5) {
    const g = geodesic_mink[i];
    console.log(`   ${g.tau.toFixed(1).padStart(4)}    ${g.x[0].toFixed(3).padStart(8)}  ${g.x[1].toFixed(3).padStart(8)}  ${g.p[0].toFixed(3).padStart(8)}  ${g.p[1].toFixed(3).padStart(8)}`);
}
console.log();

// ============================================================================
// 6. Geodesic in Schwarzschild (Orbit)
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('6. GEODESIC IN SCHWARZSCHILD SPACETIME\n');

const relH_schw = new CT.RelativisticHamiltonian(schw, 1);

// Initial conditions for a bound orbit
// At r = 10M, equatorial plane (θ = π/2)
// Need to choose E (energy) and L (angular momentum) for bound orbit
const r0 = 10;
const E = 0.95;   // Energy per unit mass < 1 for bound orbit
const L = 4.0;    // Angular momentum

// Compute initial p_t and p_φ
// p_t = E (conserved), p_φ = L (conserved)
// p_r determined by H = 0

const f0 = 1 - 2*M/r0;
// H = ½[p_t²/f - p_r²f - p_φ²/r²] - ½m² = 0
// Solve for p_r: p_r² = (E²/f - m² - L²/r²) / f
const pr_sq = (E*E/f0 - 1 - L*L/(r0*r0)) / f0;
const pr_init = pr_sq > 0 ? Math.sqrt(pr_sq) : 0;

const x0_schw = [0, r0, Math.PI/2, 0];
const p0_schw = [E, pr_init, 0, L];

console.log('Initial conditions (equatorial orbit):');
console.log('   r₀ =', r0, '(>> 2M = 2)');
console.log('   E =', E, '(energy per unit mass)');
console.log('   L =', L, '(angular momentum)');
console.log('   p_r(0) =', pr_init.toFixed(4));
console.log();

// Check H
const H_schw_init = relH_schw.evaluate(x0_schw, p0_schw);
console.log('   H =', H_schw_init.toFixed(6), '(should be ≈ 0)');
console.log();

// Integrate
const geodesic_schw = relH_schw.integrateGeodesic(x0_schw, p0_schw, 1, 100);

console.log('Schwarzschild geodesic (first portion):');
console.log('   τ        t         r         φ');
console.log('   ────────────────────────────────────────');

for (let i = 0; i <= 100; i += 20) {
    const g = geodesic_schw[i];
    console.log(`   ${g.tau.toFixed(0).padStart(4)}     ${g.x[0].toFixed(2).padStart(8)}   ${g.x[1].toFixed(3).padStart(8)}   ${g.x[3].toFixed(3).padStart(8)}`);
}
console.log();

// Check conservation
const g_final = geodesic_schw[100];
console.log('Conservation check:');
console.log('   E(0) =', E.toFixed(4), ', E(τ) =', g_final.p[0].toFixed(4));
console.log('   L(0) =', L.toFixed(4), ', L(τ) =', g_final.p[3].toFixed(4));
console.log();

// ============================================================================
// 7. Electromagnetic Coupling
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('7. ELECTROMAGNETIC COUPLING\n');

// Uniform electric field in x-direction: A_μ = (Ex, 0, 0, 0)
const E_field = 0.1;
const gaugePotential = x => [E_field * x[1], 0, 0, 0];  // A_t = E·x

const relH_em = new CT.RelativisticHamiltonian(mink, 1, gaugePotential, 1);

console.log('Charged particle in uniform electric field:');
console.log('   A_μ = (E·x, 0, 0, 0)');
console.log('   E =', E_field);
console.log('   q = 1');
console.log();

console.log('The Hamiltonian becomes:');
console.log('   H = ½g^μν(p_μ - qA_μ)(p_ν - qA_ν) - ½m²');
console.log();

// ============================================================================
// Summary
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('SUMMARY: GR EXTENSION\n');

console.log('The key principle:');
console.log('');
console.log('   ┌─────────────────────────────────────────────────────┐');
console.log('   │  Kinematics from α  ←→  Contact structure (flat)   │');
console.log('   │  Curvature from g_μν ←→ Inside Hamiltonian H       │');
console.log('   └─────────────────────────────────────────────────────┘');
console.log('');
console.log('Available metrics:');
console.log('   • Minkowski (flat spacetime)');
console.log('   • Schwarzschild (black hole)');
console.log('   • FLRW (cosmology)');
console.log('');
console.log('Mass-shell constraint:');
console.log('   H = ½g^μν(p_μ - qA_μ)(p_ν - qA_ν) - ½m² = 0');
console.log();

console.log('═══════════════════════════════════════════════════════════════');
console.log('  Geodesics examples complete!');
console.log('═══════════════════════════════════════════════════════════════');
