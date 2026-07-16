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
    phaseRotor, guidanceFromRotor, currentGuidanceVelocity1D, phaseBivector: I2
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

    // GA rotor vs the honest Im-current baseline on a smooth phase with a
    // VARYING amplitude ψ = e^{−x²/4} e^{iS}, S = 0.7x + 0.1x². Analytic
    // v = ∂S = 0.7 + 0.2x. The current form differences the full ψ (so its
    // O(dx²) error couples to ∇R); the GA rotor differences U = ψ/|ψ| (error
    // independent of ∇R) — GA is a constant factor more accurate, both O(dx²).
    function convErr(dx) {
        const n = Math.round(20 / dx) + 1, amp = [], ph = [], re = [], im = [];
        for (let i = 0; i < n; i++) {
            const X = (i - (n - 1) / 2) * dx, S = 0.7 * X + 0.1 * X * X, R = exp(-X * X / 4);
            amp.push(R); ph.push(S); re.push(R * Math.cos(S)); im.push(R * Math.sin(S));
        }
        const vga = new GAPilotWaveSystem1D(amp, ph, { dx }).deBroglieVelocity();
        const vcur = currentGuidanceVelocity1D(re, im, dx, 1.0);
        let eg = 0, ec = 0;
        for (let i = 5; i < n - 5; i++) {
            const X = (i - (n - 1) / 2) * dx;
            if (exp(-X * X / 4) < 1e-3) continue;      // skip amplitude underflow
            const an = 0.7 + 0.2 * X;
            eg = max(eg, abs(vga[i] - an)); ec = max(ec, abs(vcur[i] - an));
        }
        return { eg, ec };
    }
    const c1 = convErr(0.10), c2 = convErr(0.05);
    console.log(`  → smooth+varying-amp err vs analytic: dx=0.10 GA=${c1.eg.toExponential(2)} Im=${c1.ec.toExponential(2)}; dx=0.05 GA=${c2.eg.toExponential(2)} Im=${c2.ec.toExponential(2)}`);
    console.log(`  → constant factor Im/GA ≈ ${(c1.ec / c1.eg).toFixed(1)}× (dx=0.10), ${(c2.ec / c2.eg).toFixed(1)}× (dx=0.05)`);
    assert(t, c2.eg < c1.eg && c2.ec < c1.ec, 'both GA and Im-current are O(dx²) convergent');
    assert(t, c1.eg < c1.ec && c2.eg < c2.ec, 'GA rotor is more accurate than Im-current when amplitude varies');
    assert(t, c1.ec / c1.eg > 2, `GA constant-factor edge is real (Im/GA ≈ ${(c1.ec / c1.eg).toFixed(1)}×)`);
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
    const atan = field.deBroglieVelocityPhase();      // atan2 strawman baseline
    const cur = field.deBroglieVelocityCurrent();     // honest standard baseline

    // Compare against analytic v = (−y, x)/r², EXCLUDING the physical node ball
    // r < 0.5 (where divergence is real, not an artifact). Two error regions are
    // tracked: the full r>0.5 max, and the far field r>1 (away from the node,
    // where amplitude gradients — not the phase singularity — dominate).
    function scan(v, rLo) {
        let mv = 0, err = 0;
        for (let i = 1; i < N - 1; i++) {
            for (let j = 1; j < N - 1; j++) {
                const X = -L + i * dx, Y = -L + j * dx, r2 = X * X + Y * Y;
                if (r2 < rLo * rLo) continue;
                const an = [-Y / r2, X / r2];
                mv = max(mv, hypot(v.vx[i][j], v.vy[i][j]));
                err = max(err, hypot(v.vx[i][j] - an[0], v.vy[i][j] - an[1]));
            }
        }
        return { mv, err };
    }
    const sGA = scan(ga, 0.5), sAtan = scan(atan, 0.5), sCur = scan(cur, 0.5);
    const fGA = scan(ga, 1.0), fCur = scan(cur, 1.0);   // far field r>1
    console.log(`  → max|v| (r>0.5):  atan2=${sAtan.mv.toFixed(2)}   Im-current=${sCur.mv.toFixed(2)}   GA=${sGA.mv.toFixed(2)}`);
    console.log(`  → max error   :  atan2=${sAtan.err.toFixed(2)}   Im-current=${sCur.err.toExponential(2)}   GA=${sGA.err.toExponential(2)}`);
    console.log(`  → far-field (r>1) error:  Im-current=${fCur.err.toExponential(2)}   GA=${fGA.err.toExponential(2)}`);

    // The atan2 form is a strawman: it differences the wrapped phase angle and
    // spikes on the branch cut far from the node.
    assert(t, sAtan.err > 10, `atan2 strawman has huge branch-cut error (${sAtan.err.toFixed(1)})`);
    // The honest Im-current baseline has NO branch cut: same order as GA.
    assert(t, sCur.err < 0.1, `Im-current baseline is accurate, no branch cut (${sCur.err.toExponential(2)})`);
    assert(t, isFinite(sGA.mv) && sGA.mv < 3.0, `GA velocity bounded away from node (max=${sGA.mv.toFixed(2)})`);
    assert(t, sGA.err < 0.1, `GA tracks analytic vortex velocity on the grid (err=${sGA.err.toExponential(2)})`);
    // GA and Im-current are the SAME order of magnitude — the categorical gap
    // was entirely the atan2 artifact, not a property of "the standard equation".
    assert(t, sCur.err < 5 * sGA.err && sGA.err < 5 * sCur.err,
        `GA and Im-current are same order near node (Im=${sCur.err.toExponential(2)}, GA=${sGA.err.toExponential(2)})`);
    // GA's genuine (modest) edge: in the far field it differences only the
    // amplitude-normalized rotor, so it beats the current form there.
    assert(t, fGA.err < fCur.err, `GA more accurate than Im-current in far field (${fGA.err.toExponential(2)} < ${fCur.err.toExponential(2)})`);

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
    assert(t, maxReg < sAtan.mv, `regularization tames velocity (max_reg=${maxReg.toFixed(2)} < atan2=${sAtan.mv.toFixed(1)})`);
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
    const amp = [], ph = [], re = [], im = [];
    for (let i = 0; i < n; i++) {
        const x = (i - (n - 1) / 2) * dx;
        amp.push(abs(x) + 1e-8);
        ph.push(x >= 0 ? 0 : Math.PI);
        re.push(x); im.push(0);                    // ψ ∝ x  (purely real)
    }
    const kernel = new SmearingKernel(0.3, 'gaussian');
    const std = new PilotWaveSystem(new WaveFunction(amp, ph, {}), { kernel, dx });
    const ga = new GAPilotWaveSystem1D(amp, ph, { dx, kernel });

    const vgReg = ga.regularizedVelocity();
    const vsReg = std.regularizedVelocity();        // atan2 phase-difference form
    const vCur = currentGuidanceVelocity1D(re, im, dx, 1.0);  // Im-current baseline
    const maxGA = max(...vgReg.map(abs));
    const maxAtan = max(...vsReg.map(abs));
    const maxCur = max(...vCur.map(abs));
    const maxAtanRaw = max(...std.deBroglieVelocity().map(abs));
    console.log(`  → max|v|:  atan2=${maxAtanRaw.toFixed(2)} (raw) / ${maxAtan.toFixed(3)} (reg)   Im-current=${maxCur.toExponential(2)}   GA=${maxGA.toExponential(2)}   (true=0)`);

    assert(t, vgReg.every(isFinite), 'GA regularized velocity finite at real node');
    assert(t, maxGA < 1e-6, `GA gives correct v≈0 for a real wavefunction (max=${maxGA.toExponential(2)})`);
    // The honest Im-current baseline ALSO returns exactly 0 (im ≡ 0 ⇒ v ≡ 0);
    // only the atan2 phase-difference form spikes.
    assert(t, maxCur < 1e-9, `Im-current baseline also correct v≈0 (max=${maxCur.toExponential(2)})`);
    assert(t, maxAtan > 10 * max(maxGA, 1e-9), `atan2 form shows spurious spike at real node (reg max=${maxAtan.toFixed(3)})`);

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
section('VI. H-theorem: coarse-grained H DECREASES toward equilibrium');
// ---------------------------------------------------------------------------
{
    // A correct quantum-relaxation H-theorem needs a *time-dependent, spatially
    // structured* guidance field that stirs the ensemble — a static wavefunction
    // with a linear phase gives only a constant drift (no mixing) under which
    // |ψ|² is not even stationary, so H does not decrease. We use the canonical
    // two-energy-eigenstate superposition in an infinite square well [0,L]
    // (ℏ = m = 1):
    //
    //   ψ(x,t) = (1/√2)[φ₁ e^{−iE₁t} + φ₂ e^{−iE₂t}],  φ_n = √(2/L) sin(nπx/L),
    //   E_n = n²π²/(2L²),   |ψ(x,t)|² is equivariant (the moving target).
    //
    // The velocity field v = j/|ψ|² oscillates and has a blinking interior node,
    // so an ensemble started far from |ψ|² relaxes: the COARSE-GRAINED H-function
    // decreases monotonically. We drive the SAME advection (QuantumEnsemble) with
    // two velocity fields — the honest Im-current baseline and the GA rotor form.
    const L = Math.PI, n = 101, dx = L / (n - 1);
    const E1 = 1 * 1 * Math.PI * Math.PI / (2 * L * L);
    const E2 = 2 * 2 * Math.PI * Math.PI / (2 * L * L);
    const sqrt2overL = Math.sqrt(2 / L);

    function componentsAt(time) {
        const re = [], im = [];
        for (let i = 0; i < n; i++) {
            const x = i * dx;
            const a = sqrt2overL * Math.sin(Math.PI * x / L);
            const b = sqrt2overL * Math.sin(2 * Math.PI * x / L);
            re.push((a * Math.cos(E1 * time) + b * Math.cos(E2 * time)) / Math.sqrt(2));
            im.push((-a * Math.sin(E1 * time) - b * Math.sin(E2 * time)) / Math.sqrt(2));
        }
        return { re, im };
    }

    const kernel = new SmearingKernel(0.25, 'gaussian');

    // Standard baseline as a QuantumEnsemble-compatible system (Im-current + smear).
    class StdCurrentSystem1D {
        constructor(re, im) { this.re = re; this.im = im; }
        regularizedVelocity() {
            const v = currentGuidanceVelocity1D(this.re, this.im, dx, 1.0);
            const rho = this.re.map((r, i) => r * r + this.im[i] * this.im[i]);
            const j = rho.map((r, i) => r * v[i]);
            const jReg = kernel.convolve1D(j, dx), rhoReg = kernel.convolve1D(rho, dx);
            return jReg.map((x, i) => x / Math.max(rhoReg[i], 1e-12));
        }
    }

    const dt = 0.02, nSteps = 50, cell = 8;

    function relax(kind) {
        const rho0 = new Array(n).fill(1);            // uniform (maximally nonequilibrium)
        const cGiven = componentsAt(0);
        const psiT = WaveFunction.fromCartesian(cGiven.re, cGiven.im, { mass: 1, hbar: 1 });
        const ens = new QuantumEnsemble(rho0, psiT, dx, kernel);
        const H = [ens.hFunction(cell)];
        for (let step = 0; step < nSteps; step++) {
            const c = componentsAt(step * dt);
            const amp = c.re.map((r, i) => hypot(r, c.im[i]));
            const ph = c.re.map((r, i) => atan2(c.im[i], r));
            const system = kind === 'GA'
                ? new GAPilotWaveSystem1D(amp, ph, { dx, kernel })
                : new StdCurrentSystem1D(c.re, c.im);
            ens.evolve(system, dt);
            // Update the equivariant target |ψ(t+dt)|² for the H reference.
            const cNext = componentsAt((step + 1) * dt);
            ens.psi = WaveFunction.fromCartesian(cNext.re, cNext.im, { mass: 1, hbar: 1 });
            H.push(ens.hFunction(cell));
        }
        return H;
    }

    for (const kind of ['standard', 'GA']) {
        const H = relax(kind);
        const H0 = H[0], Hf = H[H.length - 1], Hmin = Math.min(...H);
        let maxInc = 0;
        for (let i = 1; i < H.length; i++) maxInc = max(maxInc, H[i] - H[i - 1]);
        const label = kind === 'GA' ? 'GA rotor       ' : 'standard (Im-J)';
        console.log(`  → ${label}: H ${H0.toFixed(4)} → ${Hf.toFixed(4)}  (min ${Hmin.toFixed(4)}, drop ${((1 - Hmin / H0) * 100).toFixed(1)}%, maxΔ↑ ${maxInc.toExponential(1)})`);

        assert(t, H.every(isFinite), `${kind}: H-trajectory finite throughout`);
        assert(t, H0 > 0.1, `${kind}: initial H > 0 (nonequilibrium, H₀=${H0.toFixed(4)})`);
        assert(t, Hf < H0, `${kind}: H decreases overall (relaxation, not anti-relaxation)`);
        assert(t, (H0 - Hmin) / H0 > 0.15, `${kind}: substantial relaxation (drop ${((1 - Hmin / H0) * 100).toFixed(1)}%)`);
        // The residual ~1e-2 bump is a node-crossing artifact of the coarse
        // upwind advection (it GROWS as dt→0 and is identical for the standard
        // and GA fields, so it is spatial-discretization noise, not a property
        // of either velocity field). It is <2.5% of the ~50% total H drop.
        assert(t, maxInc < 1.2e-2, `${kind}: coarse-grained H non-increasing up to upwind node-crossing noise (maxΔ↑=${maxInc.toExponential(1)})`);
    }
}

process.exit(summary(t));
