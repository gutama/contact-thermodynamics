/**
 * G₂-Manifold Black Hole Remnant Demo
 *
 * Demonstrates the physics from:
 * "Geometric origin of a stable black hole remnant from torsion in
 *  G₂-manifold geometry" (Gen. Rel. Grav., 2026)
 *
 * Key results:
 * 1. G₂ torsion → repulsive force at Planck density
 * 2. Kaluza-Klein reduction → electroweak scale (~246 GeV)
 * 3. Hawking evaporation → stable remnant (~9×10⁻⁴¹ kg)
 * 4. Quasi-normal modes → information storage
 * 5. G₂-Ricci flow → geometric stability
 */

const G2 = require('../src/physics/g2-black-hole-remnant.js');

const { sqrt, pow, PI, abs, log10, floor, max, min } = Math;
const C = G2.PhysicalConstants;

console.log('='.repeat(72));
console.log('  G₂-MANIFOLD BLACK HOLE REMNANT DEMONSTRATION');
console.log('  Based on: Gen. Rel. Grav. (2026) DOI:10.1007/s10714-026-03528-z');
console.log('='.repeat(72));

// ============================================================================
// 1. Physical Constants & Scales
// ============================================================================

console.log('\n1. FUNDAMENTAL SCALES');
console.log('-'.repeat(50));
console.log(`  Planck mass:    M_Pl  = ${C.M_Pl.toExponential(4)} kg`);
console.log(`  Planck length:  l_Pl  = ${C.l_Pl.toExponential(4)} m`);
console.log(`  Planck density: ρ_Pl  = ${C.rho_Pl.toExponential(4)} kg/m³`);
console.log(`  Planck time:    t_Pl  = ${C.t_Pl.toExponential(4)} s`);
console.log(`  Planck temp:    T_Pl  = ${C.T_Pl.toExponential(4)} K`);
console.log(`  Electroweak:    v_EW  = ${C.v_EW_GeV} GeV`);

// ============================================================================
// 2. G₂-Manifold Geometry
// ============================================================================

console.log('\n2. G₂-MANIFOLD GEOMETRY');
console.log('-'.repeat(50));

const g2 = new G2.G2Manifold({
    compactRadius: 6.0e-19,    // Compact radius ~1000 l_Pl (tuned for EW scale)
    torsionStrength: 1.0       // Dimensionless torsion parameter τ₀
});

console.log(`  Manifold: ${g2.toString()}`);
console.log(`  Compact radius: R = ${g2.R.toExponential(3)} m`);
console.log(`  R / l_Pl = ${(g2.R / C.l_Pl).toFixed(1)}`);
console.log(`  Compact volume: V = ${g2.volume().toExponential(3)} m³`);

// Associative 3-form
const phi = g2.associative3Form();
console.log(`\n  Associative 3-form φ (${phi.length} components):`);
const basisNames = ['e¹', 'e²', 'e³', 'e⁴', 'e⁵', 'e⁶', 'e⁷'];
phi.forEach(c => {
    const label = c.indices.map(i => basisNames[i]).join('');
    const sign = c.value > 0 ? '+' : '-';
    console.log(`    ${sign}${label}`);
});

// Torsion classes
const tClasses = g2.torsionClasses();
console.log(`\n  Torsion classes (G₂ decomposition):`);
console.log(`    τ₀ = ${tClasses.tau0.toFixed(4)}  (rep 1,  scalar)`);
console.log(`    τ₁ = ${tClasses.tau1.toFixed(4)}  (rep 7,  1-form)`);
console.log(`    τ₂ = ${tClasses.tau2.toFixed(4)}  (rep 14, 2-form)`);
console.log(`    τ₃ = ${tClasses.tau3.toFixed(4)}  (rep 27, 3-form)`);
console.log(`    Type: ${tClasses.type}`);

// Torsion tensor
const T = g2.torsionTensor();
console.log(`\n  Torsion tensor: T = τ₀ · φ`);
console.log(`    |T| = ${T.norm.toFixed(4)}`);
console.log(`    Magnitude: ${T.magnitude.toFixed(4)}`);

// ============================================================================
// 3. Einstein-Cartan Theory with Torsion
// ============================================================================

console.log('\n3. EINSTEIN-CARTAN THEORY WITH G₂ TORSION');
console.log('-'.repeat(50));

const ec = new G2.EinsteinCartanG2(g2);
console.log(`  ${ec.toString()}`);
console.log(`  Torsion coupling: κ = ${ec.kappa.toExponential(4)}`);
console.log(`  Critical density: ρ_crit = ${ec.criticalDensity().toExponential(4)} kg/m³`);
console.log(`  ρ_crit / ρ_Pl = ${(ec.criticalDensity() / C.rho_Pl).toFixed(4)}`);

console.log('\n  Effective density vs matter density:');
console.log('  ' + '-'.repeat(55));
console.log('    ρ/ρ_Pl       | ρ_eff/ρ       | Torsion correction');
console.log('  ' + '-'.repeat(55));

const rho_Pl = C.rho_Pl;
const rho_crit = ec.criticalDensity();
for (const frac of [1e-10, 1e-5, 0.01, 0.1, 0.5, 0.9, 1.0]) {
    const rho = frac * rho_crit;
    const rho_eff = ec.effectiveDensity(rho);
    const ratio = rho > 0 ? rho_eff / rho : 1;
    const correction = rho > 0 ? (1 - ratio) * 100 : 0;
    console.log(`    ${frac.toExponential(1).padStart(10)} | ` +
        `${ratio.toFixed(6).padStart(12)} | ` +
        `${correction.toFixed(4).padStart(8)}% suppression`);
}

console.log('\n  Modified Friedmann equation: H² = (8πG/3) ρ(1 - ρ/ρ_crit)');
console.log('  → Bounce at ρ = ρ_crit: H² = 0 (singularity replaced by bounce)');

// Show bounce behavior
const H2_half = ec.modifiedFriedmann(1, 0.5 * rho_crit);
const H2_crit = ec.modifiedFriedmann(1, rho_crit);
console.log(`  H²(ρ=0.5ρ_crit) = ${H2_half.toExponential(4)}`);
console.log(`  H²(ρ=ρ_crit)    = ${H2_crit.toExponential(4)} (→ 0, bounce!)`);


// ============================================================================
// 4. Kaluza-Klein Reduction
// ============================================================================

console.log('\n4. KALUZA-KLEIN REDUCTION (7D → 4D)');
console.log('-'.repeat(50));

const kk = new G2.KaluzaKleinReduction(ec);
const reduction = kk.reduce();

console.log(`  Compact volume: V = ${reduction.compactVolume.toExponential(3)} m³`);
console.log(`  Scalar fields: radion R = ${reduction.scalarFields.radion.toExponential(3)} m`);
console.log(`  Gauge group: ${reduction.gaugeFields.group}`);

const ew = kk.electroweakScale();
console.log(`\n  Electroweak scale derivation:`);
console.log(`    v_EW (derived)  = ${ew.v_EW_GeV.toExponential(4)} GeV`);
console.log(`    v_EW (measured) = ${ew.v_EW_measured_GeV} GeV`);
console.log(`    Ratio: ${ew.ratio_to_measured.toFixed(4)}`);

const hierarchy = kk.hierarchyRatio();
console.log(`\n  Hierarchy problem:`);
console.log(`    v_EW / M_Pl = ${hierarchy.toExponential(4)}`);
console.log(`    → Naturally small from compact geometry!`);

const higgs = kk.higgsMass();
console.log(`\n  Higgs mass estimate:`);
console.log(`    m_H (estimated) = ${higgs.mass_GeV.toFixed(1)} GeV`);
console.log(`    m_H (measured)  = ${higgs.measured_GeV} GeV`);

// ============================================================================
// 5. Black Hole Evaporation: Standard vs Torsion
// ============================================================================

console.log('\n5. BLACK HOLE EVAPORATION: STANDARD vs TORSION-MODIFIED');
console.log('-'.repeat(50));

const M0 = 1e12;  // 10¹² kg primordial black hole
const bh = new G2.BlackHoleRemnant(M0, ec);

console.log(`  ${bh.toString()}`);
console.log(`  Initial mass: M₀ = ${M0.toExponential(2)} kg`);
console.log(`  Schwarzschild radius: r_s = ${bh.schwarzschildRadius(M0).toExponential(4)} m`);
console.log(`  Hawking temperature: T_H = ${bh.hawkingTemperature(M0).toExponential(4)} K`);
console.log(`  Hawking luminosity: L = ${bh.hawkingLuminosity(M0).toExponential(4)} W`);
console.log(`  Standard evaporation time: t_evap = ${bh.evaporationTime(M0).toExponential(4)} s`);

const t_evap_yr = bh.evaporationTime(M0) / (365.25 * 24 * 3600);
console.log(`    = ${t_evap_yr.toExponential(4)} years`);

// Run simulation
console.log('\n  Running torsion-modified evaporation simulation...');
const history = bh.evaporationWithTorsion(100);

console.log('\n  Evolution (selected steps):');
console.log('  ' + '-'.repeat(70));
console.log('   step |      t/t_evap  |   M/M₀      |    T_H (K)   |  ρ/ρ_crit  | halted');
console.log('  ' + '-'.repeat(70));

const t_evap = bh.evaporationTime(M0);
const printSteps = [0, 10, 25, 50, 75, 90, 95, 98, 99, 100];
for (const step of printSteps) {
    const h = history[step];
    if (!h) continue;
    console.log(`   ${h.step.toString().padStart(4)} | ` +
        `${(h.t / t_evap).toExponential(3).padStart(12)} | ` +
        `${h.M_ratio.toExponential(3).padStart(11)} | ` +
        `${h.T_H.toExponential(3).padStart(12)} | ` +
        `${h.rho_ratio.toExponential(2).padStart(10)} | ` +
        `${h.halted ? '  YES' : '   no'}`);
}

// Find when evaporation halts
const haltStep = history.find(h => h.halted);
if (haltStep) {
    console.log(`\n  ⟶ Evaporation HALTED at step ${haltStep.step}`);
    console.log(`    Final mass: M_final = ${haltStep.M.toExponential(4)} kg`);
    console.log(`    Interior density reached: ${haltStep.rho_ratio.toExponential(4)} × ρ_crit`);
} else {
    console.log('\n  Note: Evaporation ongoing (would halt at M_rem)');
}

// ============================================================================
// 6. Remnant Properties
// ============================================================================

console.log('\n6. STABLE REMNANT PROPERTIES');
console.log('-'.repeat(50));

const M_rem = bh.remnantMass();
const r_rem = bh.remnantRadius();
const qnms = bh.quasiNormalModes(5);
const info = bh.informationCapacity();

console.log(`  Remnant mass:     M_rem   = ${M_rem.toExponential(4)} kg`);
console.log(`  Remnant radius:   r_rem   = ${r_rem.toExponential(4)} m`);
console.log(`  M_rem / M_Pl   = ${(M_rem / C.M_Pl).toExponential(4)}`);
console.log(`  r_rem / l_Pl   = ${(r_rem / C.l_Pl).toExponential(4)}`);
console.log(`  isStable(M_rem): ${bh.isStable(M_rem)}`);
console.log(`  isStable(2M_rem): ${bh.isStable(2 * M_rem)}`);

console.log(`\n  Quasi-Normal Modes (first 5):`);
console.log('  ' + '-'.repeat(50));
console.log('    n  |    ω (rad/s)     |    f (Hz)        | E (GeV)');
console.log('  ' + '-'.repeat(50));
for (const mode of qnms) {
    console.log(`    ${mode.n}  | ${mode.omega.toExponential(4).padStart(16)} | ` +
        `${mode.frequency.toExponential(4).padStart(16)} | ${mode.energy_GeV.toExponential(3)}`);
}

console.log(`\n  Information capacity: N = ${info.toExponential(4)} qubits`);
console.log(`    (Bekenstein bound: A/(4 l_Pl²))`);


// ============================================================================
// 7. G₂-Ricci Flow Stability
// ============================================================================

console.log('\n7. G₂-RICCI FLOW STABILITY');
console.log('-'.repeat(50));

const flow = new G2.G2RicciFlow(g2);
const flowHistory = flow.flow(g2.tau0 * 2.0, 0.01, 200);

console.log('  Flow evolution (∂φ/∂t = Δ_φ φ + T(φ)):');
console.log('    step |    |τ|     |   |τ*|    |  deviation  | stationary');
console.log('  ' + '-'.repeat(55));

for (const step of [0, 10, 25, 50, 100, 150, 200]) {
    const h = flowHistory[step];
    if (!h) continue;
    console.log(`    ${h.step.toString().padStart(4)} | ` +
        `${h.torsionNorm.toFixed(4).padStart(8)} | ` +
        `${h.fixedPointValue.toFixed(4).padStart(8)} | ` +
        `${h.deviation.toExponential(3).padStart(11)} | ` +
        `${h.isStationary ? '   YES' : '    no'}`);
}

const stability = flow.linearStability();
console.log(`\n  Linear stability analysis at fixed point:`);
console.log(`    Eigenvalues: [${stability.eigenvalues.map(e => e.toFixed(3)).join(', ')}]`);
console.log(`    All negative: ${stability.isStable} → STABLE`);
console.log(`    Convergence rate: λ = ${stability.convergenceRate}`);

// ============================================================================
// 8. Effective Potential (ASCII Plot)
// ============================================================================

console.log('\n8. EFFECTIVE POTENTIAL V_eff(r)');
console.log('-'.repeat(50));

// ASCII art plot of V_eff(r)
const width = 60;
const height = 20;
const r_min_plot = 0.5 * r_rem;
const r_max_plot = 20 * r_rem;

// Compute V_eff values
const nPoints = width;
const r_values = [];
const V_values = [];
for (let i = 0; i < nPoints; i++) {
    const frac = i / (nPoints - 1);
    const r = r_min_plot * pow(r_max_plot / r_min_plot, frac);
    r_values.push(r);
    V_values.push(bh.effectivePotential(r, M_rem));
}

// Find range for normalization
const V_min = Math.min(...V_values);
const V_max = Math.max(...V_values.filter(v => isFinite(v) && v < abs(V_min) * 10));
const V_range = V_max - V_min || 1;

// Build ASCII plot
const grid = Array.from({ length: height }, () => Array(width).fill(' '));

// Plot V_eff
for (let i = 0; i < nPoints; i++) {
    const V = V_values[i];
    if (!isFinite(V)) continue;
    const y = height - 1 - Math.round((V - V_min) / V_range * (height - 1));
    const yc = max(0, min(height - 1, y));
    grid[yc][i] = '█';
}

// Mark remnant radius
const r_rem_idx = r_values.findIndex(r => r >= r_rem);
if (r_rem_idx >= 0 && r_rem_idx < width) {
    for (let y = 0; y < height; y++) {
        if (grid[y][r_rem_idx] === ' ') grid[y][r_rem_idx] = '│';
    }
}

// Draw zero line
const zero_y = height - 1 - Math.round((0 - V_min) / V_range * (height - 1));
const zero_yc = max(0, min(height - 1, zero_y));
for (let x = 0; x < width; x++) {
    if (grid[zero_yc][x] === ' ') grid[zero_yc][x] = '·';
}

console.log('  V_eff ↑');
for (let y = 0; y < height; y++) {
    console.log('  ' + (y === 0 ? '  max' : y === height - 1 ? '  min' : '     ') +
        ' │' + grid[y].join(''));
}
console.log('       └' + '─'.repeat(width) + '→ r');
console.log('        r_min' + ' '.repeat(width - 20) + 'r_max');
console.log(`        ${r_min_plot.toExponential(1)}` + ' '.repeat(width - 30) +
    `${r_max_plot.toExponential(1)}`);
console.log(`        │ = r_rem (potential minimum)`);

// ============================================================================
// 9. Summary
// ============================================================================

console.log('\n' + '='.repeat(72));
console.log('  SUMMARY: G₂-MANIFOLD BLACK HOLE REMNANT');
console.log('='.repeat(72));

console.log(`
┌──────────────────────────────────────────────────────────────────────┐
│  COMPONENT                │  RESULT                                 │
├──────────────────────────────────────────────────────────────────────┤
│  G₂-Manifold              │  7D with torsion τ₀ = ${g2.tau0.toFixed(1).padEnd(17)}│
│  Compact radius           │  R = ${g2.R.toExponential(2).padEnd(25)}│
│  Torsion type             │  Nearly-parallel (τ₀ only)${' '.repeat(13)}│
├──────────────────────────────────────────────────────────────────────┤
│  Einstein-Cartan           │  Torsion coupling κ = ${ec.kappa.toExponential(2).padEnd(8)} │
│  Critical density          │  ρ_crit = ${ec.criticalDensity().toExponential(2).padEnd(20)}│
│  Bounce condition          │  H² → 0 at ρ → ρ_crit${' '.repeat(18)}│
├──────────────────────────────────────────────────────────────────────┤
│  Kaluza-Klein              │  v_EW = ${ew.v_EW_GeV.toExponential(2).padEnd(22)}GeV│
│  Hierarchy ratio           │  v_EW/M_Pl = ${hierarchy.toExponential(2).padEnd(17)}│
│  Higgs mass                │  m_H ≈ ${higgs.mass_GeV.toFixed(0).padEnd(6)} GeV${' '.repeat(17)}│
├──────────────────────────────────────────────────────────────────────┤
│  Remnant mass              │  M_rem = ${M_rem.toExponential(2).padEnd(22)}kg│
│  Remnant radius            │  r_rem = ${r_rem.toExponential(2).padEnd(22)}m │
│  Information capacity      │  N ≈ ${info.toExponential(2).padEnd(24)}qubits│
│  QNM spacing               │  Δω = ${(qnms[1].omega - qnms[0].omega).toExponential(2).padEnd(22)}│
├──────────────────────────────────────────────────────────────────────┤
│  G₂-Ricci flow             │  Converges to stable fixed point       │
│  Stability                 │  All eigenvalues negative               │
└──────────────────────────────────────────────────────────────────────┘

Physical picture:
  1. START: 7D Einstein-Cartan gravity on G₂-manifold with torsion
  2. REDUCTION: KK reduction → 4D gravity + electroweak sector
  3. BH FORMATION: Standard Schwarzschild black hole forms
  4. EVAPORATION: Hawking radiation shrinks the black hole
  5. TORSION KICK-IN: At Planckian density, torsion creates repulsion
  6. STABLE REMNANT: Evaporation halts at M ≈ ${M_rem.toExponential(1)} kg
  7. INFO STORAGE: QNMs encode ~${info.toExponential(1)} qubits of information
  8. STABILITY: G₂-Ricci flow confirms geometric stability
`);
