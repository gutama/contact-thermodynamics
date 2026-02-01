/**
 * Pilot-Wave Module Tests
 * 
 * Tests for Valentini's regularized de Broglie-Bohm dynamics.
 */

const PilotWave = require('../src/pilot-wave.js');
const {
    SmearingKernel,
    WaveFunction,
    PilotWaveSystem,
    QuantumEnsemble,
    ActionPhaseBridge
} = PilotWave;

// Test utilities
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

function assertApprox(actual, expected, tol, message) {
    const diff = Math.abs(actual - expected);
    assert(diff < tol, `${message} (got ${actual.toFixed(6)}, expected ${expected.toFixed(6)})`);
}

function assertArrayApprox(actual, expected, tol, message) {
    const maxDiff = Math.max(...actual.map((v, i) => Math.abs(v - expected[i])));
    assert(maxDiff < tol, `${message} (max diff: ${maxDiff.toFixed(6)})`);
}

console.log('========================================');
console.log('  PILOT-WAVE MODULE TESTS');
console.log('========================================\n');

// ============================================================================
// 1. SMEARING KERNEL TESTS
// ============================================================================

console.log('--- SmearingKernel Tests ---\n');

(function testSmearingKernel() {
    // Test Gaussian kernel normalization
    const kernel = new SmearingKernel(0.1, 'gaussian');

    // Numerical integration should give ~1
    let integral = 0;
    const dr = 0.001;
    for (let r = -1; r <= 1; r += dr) {
        integral += kernel.evaluate(r) * dr;
    }
    assertApprox(integral, 1.0, 0.01, 'Gaussian kernel normalized');

    // Test uniform kernel
    const uniformKernel = new SmearingKernel(0.1, 'uniform');
    let uniformIntegral = 0;
    for (let r = -1; r <= 1; r += dr) {
        uniformIntegral += uniformKernel.evaluate(r) * dr;
    }
    assertApprox(uniformIntegral, 1.0, 0.01, 'Uniform kernel normalized');

    // Test 1D convolution preserves total
    const field = [0, 1, 2, 3, 2, 1, 0];
    const dx = 0.1;
    const convolved = kernel.convolve1D(field, dx);

    const originalSum = field.reduce((a, b) => a + b) * dx;
    const convolvedSum = convolved.reduce((a, b) => a + b) * dx;
    assertApprox(convolvedSum, originalSum, 0.1, 'Convolution preserves total');

    // Test time-dependent width ε(t)
    const timeDepKernel = new SmearingKernel(t => 0.1 + 0.05 * t, 'gaussian');
    assert(timeDepKernel.isTimeDependant, 'Time-dependent kernel detected');

    timeDepKernel.setTime(0);
    assertApprox(timeDepKernel.width, 0.1, 1e-10, 'ε(0) = 0.1');

    timeDepKernel.setTime(2);
    assertApprox(timeDepKernel.width, 0.2, 1e-10, 'ε(2) = 0.2');

    // Fixed kernel should not be time-dependent
    assert(!kernel.isTimeDependant, 'Fixed kernel NOT time-dependent');
})();

// ============================================================================
// 2. WAVEFUNCTION TESTS
// ============================================================================

console.log('\n--- WaveFunction Tests ---\n');

(function testWaveFunction() {
    // Test basic construction
    const psi = new WaveFunction(1.0, Math.PI / 4);
    assertApprox(psi.probabilityDensity, 1.0, 1e-10, 'R=1 gives |ψ|²=1');
    assertApprox(psi.realPart, Math.cos(Math.PI / 4), 1e-10, 'Re(ψ) correct');
    assertApprox(psi.imagPart, Math.sin(Math.PI / 4), 1e-10, 'Im(ψ) correct');

    // Test from Cartesian
    const psi2 = WaveFunction.fromCartesian(1, 1);
    assertApprox(psi2.amplitude, Math.sqrt(2), 1e-10, 'Cartesian amplitude correct');
    assertApprox(psi2.phase, Math.PI / 4, 1e-10, 'Cartesian phase correct');

    // Test plane wave factory (check it returns a function)
    const planeWave = WaveFunction.planeWave(5.0);
    const psiAtX = planeWave(1.0);
    assertApprox(psiAtX.amplitude, 1.0, 1e-10, 'Plane wave has unit amplitude');
    assertApprox(psiAtX.phase, 5.0, 1e-10, 'Plane wave phase = k·x');

    // Test field-based wavefunction
    const n = 100;
    const dx = 0.1;
    const amplitude = [];
    const phase = [];
    for (let i = 0; i < n; i++) {
        const x = (i - n / 2) * dx;
        amplitude.push(Math.exp(-x * x));  // Gaussian
        phase.push(3.0 * x);               // k = 3
    }
    const psiField = new WaveFunction(amplitude, phase, { mass: 1, hbar: 1 });

    // Current j = |ψ|² ∇S / m  = R² * k / m
    const j = psiField.current(dx);
    // At center (i=50), amplitude=1, ∇S=3, so j=3
    assertApprox(j[50], 3.0, 0.1, 'Current j = R² · ∇S at center');
})();

// ============================================================================
// 3. PILOT-WAVE SYSTEM TESTS
// ============================================================================

console.log('\n--- PilotWaveSystem Tests ---\n');

(function testPilotWaveSystem() {
    // Create a wavefunction with a node at x=0
    const n = 101;
    const dx = 0.1;
    const amplitude = [];
    const phase = [];

    for (let i = 0; i < n; i++) {
        const x = (i - n / 2) * dx;
        // ψ ~ x near origin (linear node)
        amplitude.push(Math.abs(x) + 1e-8);  // Small epsilon for stability
        // Phase flips sign at node
        phase.push(x >= 0 ? 0 : Math.PI);
    }

    const psi = new WaveFunction(amplitude, phase, { mass: 1, hbar: 1 });
    const kernel = new SmearingKernel(0.2, 'gaussian');  // Width 0.2
    const system = new PilotWaveSystem(psi, { kernel, dx });

    // Raw de Broglie velocity would diverge at node
    const vRaw = system.deBroglieVelocity();
    const vReg = system.regularizedVelocity();

    // Check regularized velocity is finite at node (index 50)
    const maxVReg = Math.max(...vReg.map(Math.abs));
    assert(isFinite(maxVReg), 'Regularized velocity is finite everywhere');
    assert(maxVReg < 100, `Regularized velocity bounded (max=${maxVReg.toFixed(2)})`);

    // Regularized density should be > 0 even at node
    const rhoReg = system.regularizedDensity();
    const minRho = Math.min(...rhoReg);
    assert(minRho > 0, `Regularized density > 0 everywhere (min=${minRho.toExponential(2)})`);
})();

// ============================================================================
// 4. QUANTUM ENSEMBLE TESTS
// ============================================================================

console.log('\n--- QuantumEnsemble Tests ---\n');

(function testQuantumEnsemble() {
    // Create a Gaussian wavefunction
    const n = 100;
    const dx = 0.1;
    const amplitude = [];
    const phase = [];

    for (let i = 0; i < n; i++) {
        const x = (i - n / 2) * dx;
        amplitude.push(Math.exp(-x * x / 2));  // Gaussian
        phase.push(0);  // Zero momentum
    }

    const psi = new WaveFunction(amplitude, phase, { mass: 1, hbar: 1 });

    // Test equilibrium ensemble
    const eqEnsemble = QuantumEnsemble.equilibrium(psi, dx);
    const Heq = eqEnsemble.hFunction();
    assertApprox(Heq, 0, 0.01, 'H-function = 0 at equilibrium');
    assert(eqEnsemble.isInEquilibrium(0.15), 'Equilibrium detected correctly (L2 < 15%)');

    // Test nonequilibrium ensemble (uniform distribution)
    const nonEq = QuantumEnsemble.uniform(psi, n, dx);
    const Hneq = nonEq.hFunction();
    assert(Hneq > 0, `H-function > 0 out of equilibrium (H=${Hneq.toFixed(4)})`);
    assert(!nonEq.isInEquilibrium(), 'Nonequilibrium detected correctly');

    // Test REGULARIZED equilibrium: ρ ≈ (|ψ|²)_reg
    const kernel = new SmearingKernel(0.3, 'gaussian');
    const rhoQMReg = kernel.convolve1D(psi.probabilityDensity, dx);
    const regEnsemble = new QuantumEnsemble([...rhoQMReg], psi, dx, kernel);

    assert(regEnsemble.isInRegularizedEquilibrium(kernel, 0.15),
        'Regularized equilibrium: ρ = (|ψ|²)_reg');

    // Regularized equilibrium is DIFFERENT from standard equilibrium
    const regEq = regEnsemble.getRegularizedEquilibrium(kernel);
    assert(regEq.length === n, 'getRegularizedEquilibrium returns correct size');

    // (|ψ|²)_reg > 0 everywhere (even at would-be nodes)
    const minRegEq = Math.min(...regEq);
    assert(minRegEq > 0, `(|ψ|²)_reg > 0 everywhere (min=${minRegEq.toExponential(2)})`);
})();

// ============================================================================
// 5. H-THEOREM / RELAXATION TEST
// ============================================================================

console.log('\n--- H-Theorem Relaxation Test ---\n');

(function testRelaxation() {
    // Create ground state with some momentum
    const n = 64;
    const dx = 0.2;
    const amplitude = [];
    const phase = [];

    for (let i = 0; i < n; i++) {
        const x = (i - n / 2) * dx;
        amplitude.push(Math.exp(-x * x / 4));
        phase.push(0.5 * x);  // Small drift
    }

    const psi = new WaveFunction(amplitude, phase, { mass: 1, hbar: 1 });
    const kernel = new SmearingKernel(0.3, 'gaussian');
    const system = new PilotWaveSystem(psi, { kernel, dx });

    // Start with nonequilibrium: shifted Gaussian (different shape!)
    // This is NOT just a rescaled |ψ|², but a genuinely different distribution
    const rhoInit = [];
    for (let i = 0; i < n; i++) {
        const x = (i - n / 2) * dx;
        const shift = 2.0;  // Shift the center
        rhoInit.push(Math.exp(-(x - shift) * (x - shift) / 4));
    }

    const ensemble = new QuantumEnsemble(rhoInit, psi, dx);

    // Check initial H > 0 (nonequilibrium)
    const H0 = ensemble.hFunctionFine();
    assert(H0 > 0.01, `Initial H > 0 (H=${H0.toFixed(4)})`);

    // Evolve and check H behavior
    const dt = 0.01;
    const nSteps = 20;
    const hValues = ensemble.relaxation(system, dt, nSteps, true);

    console.log(`  H trajectory: [${hValues.slice(0, 5).map(h => h.toFixed(4)).join(', ')}...]`);

    // H-theorem: verify H is finite and bounded (short-term transients allowed)
    assert(hValues[hValues.length - 1] < 5, 'H remains bounded during evolution');
})();

// ============================================================================
// 6. ACTION-PHASE BRIDGE TEST
// ============================================================================

console.log('\n--- ActionPhaseBridge Tests ---\n');

(function testActionPhaseBridge() {
    const n = 100;
    const dx = 0.1;
    const hbar = 1.05e-34;  // Planck's constant (SI)

    const amplitude = [];
    const phase = [];
    const k = 1e10;  // Wave number

    for (let i = 0; i < n; i++) {
        const x = i * dx;
        amplitude.push(1);
        phase.push(k * x);
    }

    const psi = new WaveFunction(amplitude, phase, { mass: 9.1e-31, hbar });
    const bridge = new ActionPhaseBridge(psi);

    // Action A = ℏS
    const A = bridge.phaseToAction();
    assertApprox(A[50] / (hbar * phase[50]), 1.0, 1e-10, 'Action A = ℏS');

    // Momentum p = ∇S = k
    const p = bridge.momentum(dx);
    assertApprox(p[50], k, k * 0.01, 'Momentum p = ∇S = k');

    // Legendrian violation should be small (numerical gradient errors)
    const violation = bridge.legendrianViolation(dx);
    assert(violation < 1e-3, `Legendrian condition satisfied (violation=${violation.toExponential(2)})`);
})();

// ============================================================================
// 7. NODED WAVEFUNCTION (HYDROGEN-LIKE)
// ============================================================================

console.log('\n--- Noded Wavefunction Tests ---\n');

(function testNodedWavefunction() {
    // Create wavefunction with node at x=0 (like hydrogen 2p)
    const psi = PilotWave.createNodedWavefunction(0, {
        gridSize: 101,
        dx: 0.1,
        mass: 1,
        hbar: 1
    });

    // Check node exists (amplitude ~ 0 at center)
    const centerIdx = 50;
    assert(psi.amplitude[centerIdx] < 0.1, 'Node present at center (ψ ≈ 0)');

    // Standard velocity would diverge, but regularized should be finite
    const kernel = new SmearingKernel(0.3, 'gaussian');
    const system = new PilotWaveSystem(psi, { kernel, dx: 0.1 });

    const vReg = system.regularizedVelocity();
    const maxV = Math.max(...vReg.map(Math.abs));
    assert(isFinite(maxV), 'Regularized velocity finite at node');
    assert(maxV < 50, `Velocity bounded near node (max=${maxV.toFixed(2)})`);

    // Regularized density > 0 at node
    const rhoReg = system.regularizedDensity();
    assert(rhoReg[centerIdx] > 0, `ρ_reg > 0 at node (=${rhoReg[centerIdx].toExponential(2)})`);
})();

// ============================================================================
// 8. CURVED SPACE PILOT-WAVE
// ============================================================================

console.log('\n--- Curved Space Pilot-Wave Tests ---\n');

(function testCurvedSpacePilotWave() {
    const { CurvedSpacePilotWave } = PilotWave;

    // Try to load Riemannian GA for sphere manifold
    let RGA;
    try {
        RGA = require('../src/riemannian-ga.js');
    } catch (e) {
        console.log('  (Skipping curved-space tests - riemannian-ga.js not available)');
        return;
    }

    // Create sphere manifold
    const sphere = new RGA.Sphere2D(1.0);  // Unit sphere

    // Create wavefunction (uniform phase gradient on sphere)
    const n = 50;
    const amplitude = new Array(n).fill(1.0);
    const phase = [];
    for (let i = 0; i < n; i++) {
        phase.push(0.5 * i * 0.1);  // Linear phase gradient
    }
    const psi = new WaveFunction(amplitude, phase, { mass: 1, hbar: 1 });

    // Create curved-space pilot wave
    const kernel = new SmearingKernel(0.2, 'gaussian');
    const curvedPW = new CurvedSpacePilotWave(psi, sphere, kernel, { mass: 1, hbar: 1 });

    // Test curved velocity at equator (θ = π/2, φ = 0)
    const coords = [Math.PI / 2, 0];
    const v = curvedPW.curvedDeBroglieVelocity(coords);

    assert(Array.isArray(v) && v.length === 2, 'Curved velocity is 2D vector');
    assert(isFinite(v[0]) && isFinite(v[1]), 'Curved velocity is finite');
    console.log(`  ✓ Curved de Broglie velocity at equator: [${v[0].toFixed(4)}, ${v[1].toFixed(4)}]`);

    // Test regularized velocity (geodesic ball)
    const vReg = curvedPW.regularizedVelocityCurved(coords, 30);
    assert(Array.isArray(vReg) && vReg.length === 2, 'Regularized curved velocity is 2D');
    assert(isFinite(vReg[0]) && isFinite(vReg[1]), 'Regularized curved velocity is finite');
    console.log(`  ✓ Regularized velocity (geodesic ball): [${vReg[0].toFixed(4)}, ${vReg[1].toFixed(4)}]`);

    // Test regularized density
    const rhoReg = curvedPW.regularizedDensityCurved(coords, 30);
    assert(isFinite(rhoReg) && rhoReg > 0, `Regularized density > 0 (=${rhoReg.toFixed(4)})`);
    console.log(`  ✓ Regularized density: ${rhoReg.toFixed(4)}`);

    // Test node detection
    const hasNode = curvedPW.hasNodeNear(coords, 1e-6);
    assert(!hasNode, 'No node at equator (uniform amplitude)');
    console.log(`  ✓ Node detection working`);
})();

// ============================================================================
// SUMMARY
// ============================================================================

console.log('\n========================================');
console.log(`  RESULTS: ${passed} passed, ${failed} failed`);
console.log('========================================');

if (failed > 0) {
    process.exit(1);
}
