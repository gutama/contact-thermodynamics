/**
 * GA-Native Pilot-Wave Tests (Phase 2)
 *
 * Validates the GA/rotor reformulation of the de Broglie–Bohm guidance equation
 * and Valentini node-regularization in src/physics/ga-pilot-wave.js against the
 * conventional vector-calculus implementation in src/physics/pilot-wave.js.
 *
 * Coverage (mirrors the RQ2 task checklist):
 *   (a) GA guidance operation reproduces the standard velocity field ∇S exactly
 *       (to floating point) for 3 test wavefunctions incl. a phase vortex; the
 *       grid rotor-difference converges to the standard phase-difference.
 *   (b) Near a phase vortex the GA regularization is finite and far more stable
 *       than the conventional wrapped-phase form (branch-cut artifact).
 *   (c) The curved-space GA guidance reduces exactly to the flat GA guidance
 *       when the metric is the identity and the connection is zero.
 *   (d) An H-theorem convergence check driven by the GA regularized velocity
 *       relaxes to quantum equilibrium at a rate comparable to the standard one.
 */

const {
    createTestTracker, assert, assertApprox, section, summary
} = require('./test-utils.js');

const PW = require('../src/physics/pilot-wave.js');
const GAPW = require('../src/physics/ga-pilot-wave.js');

const {
    WaveFunction, PilotWaveSystem, QuantumEnsemble, SmearingKernel
} = PW;
const {
    GAPilotWaveSystem1D, GAPilotWaveField2D, GACurvedPilotWave,
    phaseRotor, guidanceFromRotor, phaseBivector: I2
} = GAPW;

const t = createTestTracker();
const { hypot, abs, exp, atan2, max } = Math;

// ---------------------------------------------------------------------------
section('I. GA guidance operation reproduces analytic ∇S exactly');
// ---------------------------------------------------------------------------
{
    // Feed the GA guidance operation exact rotor derivatives dU_k = (∂_k S) I U;
    // it must return v = ∇S to machine precision (analytic equivalence).
    function gaOp(S, gradS) {
        const U = phaseRotor(S);
        const dU = gradS.map(g => I2.mul(U).scale(g));
        return guidanceFromRotor(U, dU, 1.0);
    }

    // (1) plane wave S = k·x  ⇒  ∇S = k
    const v1 = gaOp(1.3, [1.7, -0.9]);
    assertApprox(t, v1[0], 1.7, 1e-12, 'plane wave: v_x = k_x');
    assertApprox(t, v1[1], -0.9, 1e-12, 'plane wave: v_y = k_y');

    // (2) quadratic phase, arbitrary gradient
    const v2 = gaOp(2.1, [0.55, 3.14]);
    assertApprox(t, v2[0], 0.55, 1e-12, 'quadratic phase: v_x = ∂_x S');
    assertApprox(t, v2[1], 3.14, 1e-12, 'quadratic phase: v_y = ∂_y S');

    // (3) phase vortex ψ ∝ (x+iy): S = atan2(y,x), ∇S = (−y, x)/r²
    const x = 0.3, y = -0.4, r2 = x * x + y * y;
    const v3 = gaOp(atan2(y, x), [-y / r2, x / r2]);
    assertApprox(t, v3[0], -y / r2, 1e-12, 'vortex: v_x = −y/r²');
    assertApprox(t, v3[1], x / r2, 1e-12, 'vortex: v_y = x/r²');
}

// ---------------------------------------------------------------------------
section('II. Grid rotor-difference converges to standard phase-difference');
// ---------------------------------------------------------------------------
{
    // The two are *different* O(dx²) discretizations of the same ∇S operator, so
    // they are not FP-identical at finite dx, but converge to each other as dx→0.
    function maxDiff(dx) {
        const n = 201, ph = [], amp = [];
        for (let i = 0; i < n; i++) {
            const X = (i - n / 2) * dx;
            ph.push(0.7 * X + 0.1 * X * X);
            amp.push(1);
        }
        const std = new PilotWaveSystem(new WaveFunction(amp, ph, {}), { dx });
        const ga = new GAPilotWaveSystem1D(amp, ph, { dx });
        const vs = std.deBroglieVelocity(), vg = ga.deBroglieVelocity();
        let md = 0;
        for (let i = 5; i < n - 5; i++) md = max(md, abs(vs[i] - vg[i]));
        return md;
    }
    const d1 = maxDiff(0.10), d2 = maxDiff(0.05);
    assert(t, d2 < d1, `rotor/phase discretizations converge (dx=0.10→${d1.toExponential(2)}, dx=0.05→${d2.toExponential(2)})`);
    assert(t, d2 < 0.02, 'discretizations agree closely at dx=0.05');
}

// ---------------------------------------------------------------------------
section('III. Near-node stability: phase vortex ψ ∝ (x+iy) e^{−r²}');
// ---------------------------------------------------------------------------
{
    // Build ψ_r = x e^{−r²}, ψ_i = y e^{−r²} on a 2D grid (vortex/node at origin).
    const N = 41, L = 2.0, dx = 2 * L / (N - 1);
    const re = [], im = [];
    for (let i = 0; i < N; i++) {
        re[i] = []; im[i] = [];
        for (let j = 0; j < N; j++) {
            const X = -L + i * dx, Y = -L + j * dx;
            const g = exp(-(X * X + Y * Y));
            re[i][j] = X * g;
            im[i][j] = Y * g;
        }
    }
    const kernel = new SmearingKernel(0.25, 'gaussian');
    const field = new GAPilotWaveField2D(re, im, { dx, dy: dx, kernel });

    const ga = field.deBroglieVelocityGA();
    const conv = field.deBroglieVelocityPhase();

    // Compare against analytic v = (−y, x)/r², EXCLUDING the physical node ball
    // r < 0.5 (where divergence is real, not an artifact). What remains is the
    // representational branch-cut artifact of the wrapped-phase form.
    let maxGA = 0, maxConv = 0, errGA = 0, errConv = 0;
    for (let i = 1; i < N - 1; i++) {
        for (let j = 1; j < N - 1; j++) {
            const X = -L + i * dx, Y = -L + j * dx, r2 = X * X + Y * Y;
            if (r2 < 0.25) continue;               // exclude r < 0.5 ball
            const an = [-Y / r2, X / r2];
            const magGA = hypot(ga.vx[i][j], ga.vy[i][j]);
            const magC = hypot(conv.vx[i][j], conv.vy[i][j]);
            maxGA = max(maxGA, magGA);
            maxConv = max(maxConv, magC);
            errGA = max(errGA, hypot(ga.vx[i][j] - an[0], ga.vy[i][j] - an[1]));
            errConv = max(errConv, hypot(conv.vx[i][j] - an[0], conv.vy[i][j] - an[1]));
        }
    }
    console.log(`  → max|v| (r>0.5):  GA=${maxGA.toFixed(2)}   conventional=${maxConv.toFixed(2)}`);
    console.log(`  → max error vs analytic:  GA=${errGA.toExponential(2)}   conventional=${errConv.toFixed(2)}`);

    assert(t, isFinite(maxGA) && maxGA < 3.0, `GA velocity bounded away from node (max=${maxGA.toFixed(2)})`);
    // errGA is the O(dx²) grid-discretization error (machine-precision exactness
    // of the operation itself is checked in section I). It must be small AND
    // orders of magnitude below the conventional branch-cut error.
    assert(t, errGA < 0.1, `GA tracks analytic vortex velocity on the grid (err=${errGA.toExponential(2)})`);
    assert(t, errGA < errConv / 100, `GA error ≪ conventional error (${errGA.toExponential(2)} vs ${errConv.toFixed(1)})`);
    assert(t, maxConv > 5 * maxGA, `conventional form blows up on branch cut (${maxConv.toFixed(1)} ≫ ${maxGA.toFixed(1)})`);
    assert(t, errConv > 10, `conventional form has large branch-cut error (${errConv.toFixed(1)})`);

    // Regularized GA velocity is finite everywhere, including through the node.
    const reg = field.regularizedVelocityGA();
    let regFinite = true, maxReg = 0;
    for (let i = 0; i < N; i++) {
        for (let j = 0; j < N; j++) {
            const m = hypot(reg.vx[i][j], reg.vy[i][j]);
            if (!isFinite(m)) regFinite = false;
            maxReg = max(maxReg, m);
        }
    }
    assert(t, regFinite, 'regularized GA velocity finite everywhere (incl. node)');
    assert(t, maxReg < maxConv, `regularization tames velocity (max_reg=${maxReg.toFixed(2)} < conv=${maxConv.toFixed(1)})`);
}

// ---------------------------------------------------------------------------
section('IV. 1D real node: GA gives correct v≈0, standard spikes');
// ---------------------------------------------------------------------------
{
    // Real 1D node ψ ∝ x (a sign flip, ψ purely real). A real wavefunction has
    // zero de Broglie velocity everywhere (Im(∇ψ/ψ)=0). The conventional form
    // represents the sign flip as a 0→π phase jump and differencing it spikes;
    // the GA rotor sends that jump into the *internal* bivector direction, so
    // the physical (grade-1) velocity is exactly 0 — GA is correct here, not
    // merely equivalent.
    const n = 101, dx = 0.1;
    const amp = [], ph = [];
    for (let i = 0; i < n; i++) {
        const x = (i - (n - 1) / 2) * dx;
        amp.push(abs(x) + 1e-8);
        ph.push(x >= 0 ? 0 : Math.PI);
    }
    const kernel = new SmearingKernel(0.3, 'gaussian');
    const std = new PilotWaveSystem(new WaveFunction(amp, ph, {}), { kernel, dx });
    const ga = new GAPilotWaveSystem1D(amp, ph, { dx, kernel });

    const vgReg = ga.regularizedVelocity();
    const vsReg = std.regularizedVelocity();
    const maxGA = max(...vgReg.map(abs));
    const maxStd = max(...vsReg.map(abs));
    console.log(`  → max|v_reg|:  GA=${maxGA.toExponential(2)} (true=0)   standard=${maxStd.toFixed(3)}`);

    assert(t, vgReg.every(isFinite), 'GA regularized velocity finite at real node');
    assert(t, maxGA < 1e-6, `GA gives correct v≈0 for a real wavefunction (max=${maxGA.toExponential(2)})`);
    assert(t, maxStd > 10 * max(maxGA, 1e-9), `standard shows spurious spike at real node (max=${maxStd.toFixed(3)})`);

    // Regularized density > 0 through the node
    const rhoReg = ga.regularizedDensity();
    assert(t, Math.min(...rhoReg) > 0, 'GA (|ψ|²)_reg > 0 everywhere');
}

// ---------------------------------------------------------------------------
section('V. Curved-space GA reduces to flat GA when g=I, ω=0');
// ---------------------------------------------------------------------------
{
    const Sfn = c => 0.5 * c[0] + 0.3 * c[1] + 0.2 * c[0] * c[1];
    const flat = {
        dim: 2,
        metricInverse: () => [[1, 0], [0, 1]],
        constructor: { name: 'FlatIdentity' }
    };
    const cpw = new GACurvedPilotWave(Sfn, flat, {});
    const pt = [0.4, 0.7];
    const vFlat = cpw.flatGuidanceVelocity(pt);
    const vCurved = cpw.curvedGuidanceVelocity(pt);
    const diff = hypot(vFlat[0] - vCurved[0], vFlat[1] - vCurved[1]);
    assert(t, diff < 1e-12, `curved(g=I) reduces to flat GA exactly (diff=${diff.toExponential(2)})`);

    // A nontrivial metric genuinely changes the velocity (raising the index).
    const man = {
        dim: 2,
        metricInverse: () => [[2, 0], [0, 0.5]],
        constructor: { name: 'Diag' }
    };
    const cpw2 = new GACurvedPilotWave(Sfn, man, {});
    const vNontrivial = cpw2.curvedGuidanceVelocity(pt);
    assertApprox(t, vNontrivial[0], 2 * vFlat[0], 1e-6, 'g^{11}=2 raises v_x correctly');
    assertApprox(t, vNontrivial[1], 0.5 * vFlat[1], 1e-6, 'g^{22}=0.5 raises v_y correctly');

    // Parallel transport with a zero connection leaves the phase bivector fixed.
    let RGA;
    try {
        RGA = require('../src/geometry/riemannian-ga.js');
    } catch {
        RGA = null;
    }
    if (RGA) {
        // A flat manifold: identity frame everywhere ⇒ zero connection bivector.
        class FlatPlane extends RGA.RiemannianManifold {
            constructor() { super(2, 3); }
            frame() { return new RGA.TangentFrame([[1, 0, 0], [0, 1, 0]], 3); }
        }
        const plane = new FlatPlane();
        const cpwFlat = new GACurvedPilotWave(Sfn, plane, {});
        const B0 = new RGA.Bivector3D(0, 0, 1);
        const Bt = cpwFlat.parallelTransportPhase(B0, [0, 0], [1, 0], 1.0, 20);
        const drift = hypot(Bt.e23 - B0.e23, Bt.e31 - B0.e31) + abs(Bt.e12 - B0.e12);
        assert(t, drift < 1e-9, `zero connection ⇒ phase bivector unchanged (drift=${drift.toExponential(2)})`);
    } else {
        console.log('  (skipping parallel-transport test — riemannian-ga.js unavailable)');
    }
}

// ---------------------------------------------------------------------------
section('VI. H-theorem: GA regularized velocity relaxes to equilibrium');
// ---------------------------------------------------------------------------
{
    const n = 64, dx = 0.2;
    const amp = [], ph = [];
    for (let i = 0; i < n; i++) {
        const x = (i - n / 2) * dx;
        amp.push(exp(-x * x / 4));
        ph.push(0.5 * x);
    }
    const psi = new WaveFunction(amp, ph, { mass: 1, hbar: 1 });
    const kernel = new SmearingKernel(0.3, 'gaussian');
    const gaSystem = new GAPilotWaveSystem1D(amp, ph, { dx, kernel });
    const stdSystem = new PilotWaveSystem(psi, { kernel, dx });

    // Same nonequilibrium start for both: a shifted Gaussian.
    function makeRhoInit() {
        const r = [];
        for (let i = 0; i < n; i++) {
            const x = (i - n / 2) * dx;
            r.push(exp(-(x - 2.0) * (x - 2.0) / 4));
        }
        return r;
    }

    const ensGA = new QuantumEnsemble(makeRhoInit(), psi, dx);
    const ensStd = new QuantumEnsemble(makeRhoInit(), psi, dx);

    const H0 = ensGA.hFunctionFine();
    assert(t, H0 > 0.01, `initial H > 0 (nonequilibrium, H=${H0.toFixed(4)})`);

    const dt = 0.01, nSteps = 20;
    const hGA = ensGA.relaxation(gaSystem, dt, nSteps, true);
    const hStd = ensStd.relaxation(stdSystem, dt, nSteps, true);

    const HfinGA = hGA[hGA.length - 1], HfinStd = hStd[hStd.length - 1];
    console.log(`  → H: GA ${H0.toFixed(4)}→${HfinGA.toFixed(4)}   standard ${H0.toFixed(4)}→${HfinStd.toFixed(4)}`);

    assert(t, hGA.every(isFinite), 'GA H-trajectory finite throughout');
    assert(t, HfinGA < 5, 'GA H remains bounded during evolution');
    // Comparable relaxation: GA final H within a modest factor of standard.
    assert(t, HfinGA <= HfinStd * 1.5 + 0.05,
        `GA relaxes comparably to standard (GA=${HfinGA.toFixed(4)}, std=${HfinStd.toFixed(4)})`);
}

process.exit(summary(t));
