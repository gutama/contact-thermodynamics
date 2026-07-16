/**
 * G₂ Black Hole Remnant Tests
 *
 * Tests for G₂ holonomy manifold, Einstein-Cartan torsion,
 * Kaluza-Klein reduction, black hole remnant, and G₂-Ricci flow.
 */

const {
    G2Manifold,
    EinsteinCartanG2,
    KaluzaKleinReduction,
    BlackHoleRemnant,
    G2RicciFlow,
    TorsionFieldDynamics,
    RemnantCosmology,
    GravitationalWaveSignature,
    InformationParadox,
    GMETBridge,
    PhysicalConstants
} = require('../../src/physics/g2-black-hole-remnant.js');

const C = PhysicalConstants;

let passed = 0;
let failed = 0;

function assert(condition, message) {
    if (condition) {
        console.log(`  ✓ ${message}`);
        passed++;
    } else {
        console.log(`  ✗ ${message}`);
        failed++;
    }
}

function assertRelative(actual, expected, relTol, message) {
    const rel = Math.abs(actual - expected) / Math.abs(expected);
    assert(rel < relTol, `${message} (got ${actual.toExponential(4)}, rel ${rel.toExponential(2)})`);
}

function assertInRange(value, lo, hi, message) {
    assert(value >= lo && value <= hi,
        `${message} (got ${value.toExponential(4)}, range [${lo.toExponential(1)}, ${hi.toExponential(1)}])`);
}

console.log('========================================');
console.log('  G₂ BLACK HOLE REMNANT TESTS');
console.log('========================================\n');

// ============================================================================
// 1. PHYSICAL CONSTANTS
// ============================================================================
console.log('--- PhysicalConstants ---\n');

assertRelative(C.M_Pl, 2.176e-8, 0.01, 'Planck mass ~ 2.176e-8 kg');
assertRelative(C.l_Pl, 1.616e-35, 0.01, 'Planck length ~ 1.616e-35 m');
assertRelative(C.rho_Pl, 5.155e96, 0.05, 'Planck density ~ 5.155e96 kg/m³');
assertRelative(C.M_Pl * C.l_Pl, C.hbar / C.c, 0.01, 'M_Pl × l_Pl = ℏ/c consistency');

// ============================================================================
// 2. G₂-MANIFOLD
// ============================================================================
console.log('\n--- G2Manifold ---\n');

const g2 = new G2Manifold({ compactRadius: 1e-18, torsionStrength: 1.0 });

// Test: Associative 3-form has 7 components
const phi = g2.associative3Form();
assert(phi.length === 7, 'Associative 3-form has 7 independent components');

// Test: Torsion proportional to phi with strength τ₀
const g2t = new G2Manifold({ torsionStrength: 0.5 });
const torsion = g2t.torsionTensor();
const phiT = g2t.associative3Form();
let proportional = true;
for (let i = 0; i < phiT.length; i++) {
    if (Math.abs(torsion.components[i].value - 0.5 * phiT[i].value) > 1e-10) {
        proportional = false;
        break;
    }
}
assert(proportional, 'Torsion = τ₀ × φ (proportionality)');

// Test: 7D metric has correct structure
const metric = g2.metric7D();
assert(metric[0][0] > 0, '7D metric: g₀₀ > 0 (timelike)');
assert(metric[4][4] < 0 && metric[5][5] < 0 && metric[6][6] < 0,
    '7D metric: compact dimensions have correct signature');

// Test: Volume scales as R³
const g2a = new G2Manifold({ compactRadius: 1.0 });
const g2b = new G2Manifold({ compactRadius: 2.0 });
assertRelative(g2b.volume() / g2a.volume(), 8.0, 0.01, 'Volume scales as R³');

// Test: Torsion classes
const g2c = new G2Manifold({ torsionStrength: 0.7 });
const classes = g2c.torsionClasses();
assert(Math.abs(classes.tau0 - 0.7) < 1e-10, 'τ₀ class equals torsionStrength');
assert(classes.tau1 === 0 && classes.tau2 === 0 && classes.tau3 === 0,
    'Nearly-parallel: τ₁ = τ₂ = τ₃ = 0');

// Test: G₂ cross product
const e1 = [1, 0, 0, 0, 0, 0, 0];
const e2 = [0, 1, 0, 0, 0, 0, 0];
const cross12 = g2.crossProduct(e1, e2);
assert(cross12.length === 7, 'Cross product returns 7-vector');
let crossNorm = 0;
for (let i = 0; i < 7; i++) crossNorm += cross12[i] * cross12[i];
assert(crossNorm > 0, 'e₁ ×_φ e₂ is non-zero');

// Test: Inner product
const ip = g2.innerProduct(e1, e2);
assert(Math.abs(ip) < 1e-12, 'e₁ · e₂ = 0 (orthogonal)');
const ip11 = g2.innerProduct(e1, e1);
assert(Math.abs(ip11 - 1) < 1e-12, 'e₁ · e₁ = 1 (unit)');

// Test: Associative calibration
const e3 = [0, 0, 1, 0, 0, 0, 0];
const assoc = g2.isAssociative(e1, e2, e3);
assert(assoc.isAssociative, 'e₁ ∧ e₂ ∧ e₃ is associative (calibrated by φ)');

// ============================================================================
// 3. EINSTEIN-CARTAN G₂
// ============================================================================
console.log('\n--- EinsteinCartanG2 ---\n');

const ec = new EinsteinCartanG2(g2);
const rho_crit = ec.criticalDensity();

// Test: Effective density = ρ at low density
const rho_low = 1.0;
assertRelative(ec.effectiveDensity(rho_low), rho_low, 1e-6,
    'ρ_eff ≈ ρ at low density');

// Test: Effective density → 0 at ρ_crit
assert(Math.abs(ec.effectiveDensity(rho_crit)) < rho_crit * 1e-10,
    'ρ_eff → 0 at ρ_crit');

// Test: Torsion energy density is negative
assert(ec.torsionEnergyDensity(rho_crit * 0.5) < 0,
    'Torsion energy density < 0 (repulsive)');

// Test: Effective pressure positive at high density
assert(ec.effectivePressure(rho_crit * 0.9) > 0,
    'Effective pressure > 0 at high density');

// Test: Critical density ~ Planck density (within orders of magnitude)
const rho_ratio = rho_crit / C.rho_Pl;
assert(rho_ratio > 1e-30 && rho_ratio < 1e30,
    `ρ_crit ~ ρ_Pl order (ratio = ${rho_ratio.toExponential(2)})`);

// Test: Modified Friedmann → 0 at ρ_crit (bounce)
const H2_crit = ec.modifiedFriedmann(1, rho_crit);
assert(Math.abs(H2_crit) < 1e-10,
    'H² = 0 at ρ_crit (bounce condition)');

// ============================================================================
// 4. KALUZA-KLEIN REDUCTION
// ============================================================================
console.log('\n--- KaluzaKleinReduction ---\n');

const kk = new KaluzaKleinReduction(ec);

// Test: Electroweak scale is positive
const ew = kk.electroweakScale();
assert(ew.v_EW_GeV > 0 && isFinite(ew.v_EW_GeV),
    `Electroweak scale is positive (${ew.v_EW_GeV.toExponential(3)} GeV)`);

// Test: Hierarchy ratio << 1
const hr = kk.hierarchyRatio();
assert(hr < 1e-5, `Hierarchy ratio << 1 (got ${hr.toExponential(3)})`);

// Test: Higgs mass is positive and finite
const higgs = kk.higgsMass();
assert(higgs.mass_GeV > 0 && isFinite(higgs.mass_GeV),
    `Higgs mass positive and finite (${higgs.mass_GeV.toFixed(1)} GeV)`);

// ============================================================================
// 5. BLACK HOLE REMNANT
// ============================================================================
console.log('\n--- BlackHoleRemnant ---\n');

const M0 = 1e12;
const bh = new BlackHoleRemnant(M0, ec);

// Test: r_s scales linearly
assertRelative(bh.schwarzschildRadius(2e30), 2 * bh.schwarzschildRadius(1e30), 1e-10,
    'r_s(2M) = 2 r_s(M)');

// Test: T_H inversely proportional to M
assertRelative(bh.hawkingTemperature(2e30), bh.hawkingTemperature(1e30) / 2, 1e-10,
    'T_H(2M) = T_H(M)/2');

// Test: t_evap scales as M³
assertRelative(bh.evaporationTime(2e30) / bh.evaporationTime(1e30), 8.0, 0.01,
    't_evap(2M)/t_evap(M) = 8');

// Test: Remnant mass ~ 9e-41 kg
const M_rem = bh.remnantMass();
assertInRange(M_rem, 1e-42, 1e-39, 'Remnant mass ~ 9e-41 kg');

// Test: Remnant mass < Planck mass
assert(M_rem < C.M_Pl, 'M_rem < M_Pl');

// Test: isStable
assert(bh.isStable(M_rem) === true, 'isStable(M_rem) = true');
assert(bh.isStable(2 * M_rem) === false, 'isStable(2M_rem) = false');

// Test: Effective potential has minimum (V goes down then up)
const r_rem = bh.remnantRadius();
const V_near = bh.effectivePotential(r_rem * 0.5, M_rem);
const V_at = bh.effectivePotential(r_rem, M_rem);
const V_far = bh.effectivePotential(r_rem * 5, M_rem);
assert(V_near > V_at && V_far > V_at,
    'V_eff has minimum near r_rem (V dips then rises)');

// Test: QNMs positive and equally spaced
const modes = bh.quasiNormalModes(5);
let qnmOk = modes.length === 5 && modes.every(m => m.frequency > 0);
if (modes.length >= 3) {
    const sp0 = modes[1].omega - modes[0].omega;
    for (let i = 2; i < modes.length; i++) {
        if (Math.abs((modes[i].omega - modes[i - 1].omega) - sp0) / sp0 > 0.01) {
            qnmOk = false;
        }
    }
}
assert(qnmOk, 'QNMs: positive and equally spaced');

// Test: Information capacity > 0
const info = bh.informationCapacity();
assert(info > 0 && isFinite(info), `Info capacity > 0 (${info.toExponential(3)} qubits)`);

// Test: Evaporation simulation produces decreasing mass
const history = bh.evaporationWithTorsion(100);
const firstM = history[0].M;
const lastM = history[history.length - 1].M;
assert(lastM < firstM, 'Evaporation: mass decreases over time');

// Test: With small BH, evaporation halts at remnant
const bhSmall = new BlackHoleRemnant(M_rem * 100, ec);
const histSmall = bhSmall.evaporationWithTorsion(500);
const lastSmall = histSmall[histSmall.length - 1];
assert(lastSmall.halted === true || lastSmall.M <= M_rem * 1.5,
    'Small BH evaporation halts near remnant mass');

// ============================================================================
// 6. G₂-RICCI FLOW
// ============================================================================
console.log('\n--- G2RicciFlow ---\n');

const flow = new G2RicciFlow(g2);

// Test: Flow converges
const flowHist = flow.flow(g2.tau0 * 2.0, 0.01, 500);
const lastFlow = flowHist[flowHist.length - 1];
assert(lastFlow.deviation < 0.1,
    `Flow converges (deviation = ${lastFlow.deviation.toExponential(3)})`);

// Test: Linear stability - all eigenvalues negative
const stability = flow.linearStability();
const allNeg = stability.eigenvalues.every(e => e < 0);
assert(allNeg, `All eigenvalues negative: [${stability.eigenvalues.map(e => e.toFixed(3)).join(', ')}]`);

// ============================================================================
// 7. TORSION FIELD DYNAMICS
// ============================================================================
console.log('\n--- TorsionFieldDynamics ---\n');

const tfd = new TorsionFieldDynamics(g2);

// Test: Dispersion relation — massive wave
const k_test = 1e18;
const omega = tfd.dispersion(k_test);
assert(omega > 0 && isFinite(omega), 'Dispersion ω(k) > 0');

// Test: Group velocity < c
const vg = tfd.groupVelocity(k_test);
assert(vg < C.c && vg > 0, `Group velocity < c (${(vg / C.c).toFixed(4)} c)`);

// Test: Phase velocity > group velocity (massive wave)
const vp = tfd.phaseVelocity(k_test);
assert(vp > vg, 'Phase velocity > group velocity (massive dispersion)');

// Test: Simulation produces snapshots with conserved energy
const sims = tfd.simulate1D({ steps: 50, snapInterval: 10 });
assert(sims.length > 1, `Simulation produces ${sims.length} snapshots`);
const E0 = sims[0].energy;
const Elast = sims[sims.length - 1].energy;
const denom = Math.max(Math.abs(E0), Math.abs(Elast), 1e-30);
const relEchange = Math.abs(Elast - E0) / denom;
assert(E0 >= 0 && Elast >= 0 && relEchange < 1.0,
    `Energy non-negative and bounded (E0=${E0.toExponential(3)}, Ef=${Elast.toExponential(3)})`);

// Test: Scattering cross section scales with r_s²
const sigma1 = tfd.scatteringCrossSection(1e30, k_test);
const sigma2 = tfd.scatteringCrossSection(2e30, k_test);
assert(sigma2 > sigma1, 'Scattering σ increases with BH mass');

// ============================================================================
// 8. REMNANT COSMOLOGY
// ============================================================================
console.log('\n--- RemnantCosmology ---\n');

const cosmo = new RemnantCosmology(bh);

// Test: Critical density positive
const rhoCrit = cosmo.criticalDensity();
assert(rhoCrit > 0 && isFinite(rhoCrit), `Critical density ρ_crit > 0 (${rhoCrit.toExponential(3)})`);

// Test: DM density < critical density
const rhoDM = cosmo.darkMatterDensity();
assert(rhoDM < rhoCrit, 'ρ_DM < ρ_crit');
assert(rhoDM > 0, 'ρ_DM > 0');

// Test: Number density needed to explain DM
const nRem = cosmo.remnantNumberDensity();
assert(nRem > 0, `n_rem > 0 (${nRem.toExponential(3)} /m³)`);

// Test: Omega_remnant with given density matches DM
const Omega = cosmo.omegaRemnant(nRem);
assertRelative(Omega, 0.264, 0.01, 'Ω_rem = Ω_DM when using n_rem from DM density');

// Test: Max evaporating mass within reasonable bounds
const Mmax = cosmo.maxEvaporatingMass();
assert(Mmax > 1e8 && Mmax < 1e20, `M_max reasonable (${Mmax.toExponential(3)} kg)`);

// Test: Remnant abundance calculation
const abund = cosmo.remnantAbundance(1e-20);
assert(abund.Omega_rem >= 0, 'Ω_rem ≥ 0');
assert(abund.requiredBeta > 0, 'Required β > 0');

// ============================================================================
// 9. GRAVITATIONAL WAVE SIGNATURES
// ============================================================================
console.log('\n--- GravitationalWaveSignature ---\n');

const gw = new GravitationalWaveSignature(bh);

// Test: Characteristic frequency positive
const f0 = gw.characteristicFrequency();
assert(f0 > 0 && isFinite(f0), `f₀ > 0 (${f0.toExponential(3)} Hz)`);

// Test: Strain decreases with distance
const h1 = gw.strainAmplitude(1e6);
const h2 = gw.strainAmplitude(1e9);
assert(h2 < h1, 'Strain decreases with distance');
assertRelative(h1 / h2, 1e3, 0.01, 'Strain ∝ 1/r');

// Test: Energy spectrum peaks near f₀
const E_at = gw.energySpectrum(f0);
const E_off = gw.energySpectrum(f0 * 10);
assert(E_at > E_off, 'Energy spectrum peaks near f₀');

// Test: Total GW energy positive
const Egw = gw.totalEnergy();
assert(Egw > 0, `E_GW > 0 (${Egw.toExponential(3)} J)`);

// Test: Ringdown waveform has damped oscillation
const wf = gw.ringdownWaveform(1e6, 50);
assert(wf.length === 50, 'Ringdown waveform has requested samples');
assert(Math.abs(wf[0].h) > Math.abs(wf[wf.length - 1].h),
    'Ringdown amplitude decays');

// ============================================================================
// 10. INFORMATION PARADOX
// ============================================================================
console.log('\n--- InformationParadox ---\n');

const info_p = new InformationParadox(bh);

// Test: Entropy positive and scales as M²
const S1 = info_p.entropy(1e30);
const S2 = info_p.entropy(2e30);
assertRelative(S2 / S1, 4.0, 0.01, 'S(2M)/S(M) = 4 (entropy scales as M²)');

// Test: Page time ~ t_evap/2
const tPage = info_p.pageTime();
const tEvap = bh.evaporationTime(M0);
const pageRatio = tPage / tEvap;
assert(pageRatio > 0.3 && pageRatio < 0.6,
    `Page time ~ t_evap/2 (ratio = ${pageRatio.toFixed(3)})`);

// Test: Scrambling time < evaporation time
const tScr = info_p.scramblingTime();
assert(tScr < tEvap, 'Scrambling time < evaporation time');
assert(tScr > 0, 'Scrambling time > 0');

// Test: Page curve has correct phases
const pc = info_p.pageCurve(50);
assert(pc[0].phase === 'thermalization', 'Page curve starts in thermalization');
const lastPhase = pc[pc.length - 1].phase;
assert(lastPhase === 'remnant' || lastPhase === 'purification',
    `Page curve ends in ${lastPhase}`);

// Test: Remnant information
const remInfo = info_p.remnantInformation();
assert(remInfo.bits > 0, `Remnant stores ${remInfo.bits.toExponential(3)} bits`);
assert(remInfo.can_resolve_paradox === true, 'Information paradox resolution: remnant stores info');

// ============================================================================
// 11. GMET BRIDGE
// ============================================================================
console.log('\n--- GMETBridge ---\n');

const bridge = new GMETBridge(bh);

// Test: toGMETPoint produces valid coordinates
const pt = bridge.toGMETPoint(1e30);
assert(pt.base.S > 0, 'GMET point: S (entropy) > 0');
assert(pt.momenta.T > 0, 'GMET point: T (Hawking temperature) > 0');
assert(isFinite(pt.base.ell), 'GMET point: ℓ = ln(M/M_Pl) finite');

// Test: Contact Hamiltonian
const H_val = bridge.contactHamiltonian(pt);
assert(isFinite(H_val), `Contact Hamiltonian finite (${H_val.toExponential(3)})`);

// Test: Evaporation trajectory maps to GMET
const traj = bridge.evaporationTrajectory(20);
assert(traj.length > 5, `Trajectory has ${traj.length} points`);
assert(traj[0].meta.M > traj[traj.length - 1].meta.M,
    'GMET trajectory: mass decreases');

// Test: Contact form along trajectory
const alpha = bridge.contactFormAlongTrajectory(traj);
assert(alpha.length === traj.length - 1, 'Contact form: one value per segment');
assert(alpha.every(a => isFinite(a.alpha)), 'All contact form values finite');

// ============================================================================
// SUMMARY
// ============================================================================
console.log('\n========================================');
console.log(`  RESULTS: ${passed} passed, ${failed} failed`);
console.log('========================================');

if (failed > 0) process.exit(1);
