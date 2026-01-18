/**
 * Test Suite: Extended Thermodynamics Contact Geometry
 * 
 * Validates:
 * - Contact manifold dimensions
 * - Contact form structure
 * - Hamilton's equations
 * - Legendrian submanifolds
 * - GR extension
 */

const ET = require('../src/index.js');

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

function assertApprox(a, b, tol, message) {
    const diff = Math.abs(a - b);
    assert(diff < tol, `${message} (got ${a.toFixed(6)}, expected ${b.toFixed(6)}, diff ${diff.toExponential(2)})`);
}

function section(name) {
    console.log(`\n━━━ ${name} ━━━`);
}

// ============================================================================
// I. DIMENSION TESTS
// ============================================================================

section('I. Dimensionality Theorem: dim M = 2n + 1');

{
    const grand = ET.grandManifold();
    assert(grand.n === 6, 'Grand manifold base dimension n = 6');
    assert(grand.dim === 13, 'Grand manifold total dimension = 13');

    const holo = ET.holographicManifold();
    assert(holo.n === 3, 'Holographic manifold base dimension n = 3');
    assert(holo.dim === 7, 'Holographic manifold total dimension = 7');

    const gauge = ET.gaugeExtended();
    assert(gauge.n === 7, 'Gauge-extended base dimension n = 7');
    assert(gauge.dim === 15, 'Gauge-extended total dimension = 15');
}

// ============================================================================
// II. COORDINATE STRUCTURE TESTS
// ============================================================================

section('II. Coordinate Structure');

{
    const grand = ET.grandManifold();

    // Base coordinates
    const expectedBase = ['q1', 'q2', 'q3', 't', 'ell', 'S'];
    const expectedMomenta = ['k1', 'k2', 'k3', 'omega', 'Delta', 'T'];

    assert(
        JSON.stringify(grand.baseCoords) === JSON.stringify(expectedBase),
        'Grand base coords: (q¹, q², q³, t, ℓ, S)'
    );

    assert(
        JSON.stringify(grand.momentaCoords) === JSON.stringify(expectedMomenta),
        'Grand momenta: (k₁, k₂, k₃, ω, Δ, T)'
    );

    assert(grand.fiberCoord === 'A', 'Fiber coordinate is A (action)');
}

{
    const holo = ET.holographicManifold();

    assert(
        JSON.stringify(holo.baseCoords) === JSON.stringify(['t', 'ell', 'S']),
        'Holographic base: (t, ℓ, S)'
    );

    assert(
        JSON.stringify(holo.momentaCoords) === JSON.stringify(['omega', 'Delta', 'T']),
        'Holographic momenta: (ω, Δ, T)'
    );
}

// ============================================================================
// III. CONTACT FORM TESTS
// ============================================================================

section('III. Contact Form α = du - p_a dx^a');

{
    const grand = ET.grandManifold();

    // Create a test point
    const pt = grand.physicalPoint(
        1, 0, 0,    // q¹, q², q³
        0,          // t
        0,          // ℓ
        1,          // S
        0.5, 0, 0,  // k₁, k₂, k₃
        1,          // ω
        0,          // Δ
        1,          // T
        0           // A
    );

    // Tangent vector that should give α = 0 on Legendrian
    // If p_a = ∂A/∂x^a, then α(v) = du(v) - p_a dx^a(v) = 0
    // For Legendrian tangent: if v = (δx, δu = p·δx, δp), then α(v) = 0

    // Test with horizontal tangent (δx^a = 1, δu = p_a, δp = 0)
    const legendrianTangent = {
        q1: 1, q2: 0, q3: 0, t: 0, ell: 0, S: 0,
        k1: 0, k2: 0, k3: 0, omega: 0, Delta: 0, T: 0,
        A: 0.5  // = k₁ · δq¹
    };

    const alphaVal = grand.evaluateContactForm(pt, legendrianTangent);
    assertApprox(alphaVal, 0, 1e-10, 'Contact form vanishes on Legendrian tangent');

    // Test with vertical tangent (δu = 1, rest = 0)
    const verticalTangent = {
        q1: 0, q2: 0, q3: 0, t: 0, ell: 0, S: 0,
        k1: 0, k2: 0, k3: 0, omega: 0, Delta: 0, T: 0,
        A: 1
    };

    const alphaVertical = grand.evaluateContactForm(pt, verticalTangent);
    assertApprox(alphaVertical, 1, 1e-10, 'Contact form α(∂/∂u) = 1');
}

// ============================================================================
// IV. CONTACT NON-DEGENERACY: α ∧ (dα)^n ≠ 0
// ============================================================================

section('IV. Contact Non-Degeneracy');

{
    const grand = ET.grandManifold();
    const pt = grand.origin;

    const volume = grand.verifyContactCondition(pt);
    assert(volume === 720, 'Grand manifold: α ∧ (dα)⁶ = 6! = 720');

    const holo = ET.holographicManifold();
    const holoVol = holo.verifyContactCondition(holo.origin);
    assert(holoVol === 6, 'Holographic manifold: α ∧ (dα)³ = 3! = 6');

    const gauge = ET.gaugeExtended();
    const gaugeVol = gauge.verifyContactCondition(gauge.origin);
    assert(gaugeVol === 5040, 'Gauge-extended: α ∧ (dα)⁷ = 7! = 5040');
}

// ============================================================================
// V. REEB VECTOR FIELD
// ============================================================================

section('V. Reeb Vector Field: α(R) = 1, ι_R dα = 0');

{
    const grand = ET.grandManifold();
    const pt = grand.origin;

    const R = grand.reebField(pt);

    // R should be ∂/∂A (fiber direction)
    assert(R['A'] === 1, 'Reeb field R = ∂/∂A');
    assert(R['q1'] === 0, 'Reeb has no q¹ component');
    assert(R['k1'] === 0, 'Reeb has no k₁ component');

    // α(R) = 1
    const alphaR = grand.evaluateContactForm(pt, R);
    assertApprox(alphaR, 1, 1e-10, 'α(R) = 1');
}

// ============================================================================
// VI. CONTACT HAMILTONIAN DYNAMICS
// ============================================================================

section('VI. Contact Hamiltonian Dynamics');

{
    const grand = ET.grandManifold();

    // Simple Hamiltonian: H = ω (frequency as energy)
    const H = coords => coords.omega;
    const ham = new ET.ContactHamiltonian(grand, H);

    const pt = grand.physicalPoint(
        0, 0, 0, 0, 0, 1,
        0.5, 0, 0, 1, 0, 1,
        0
    );

    // Evaluate
    assertApprox(ham.evaluate(pt), 1, 1e-10, 'H(pt) = ω = 1');

    // Vector field
    const X = ham.vectorField(pt);

    // ẋ^a = ∂H/∂p_a = 0 except for t (since ∂H/∂ω = 1)
    assertApprox(X['t'], 1, 1e-6, 'ṫ = ∂H/∂ω = 1');
    assertApprox(X['q1'], 0, 1e-6, 'q̇¹ = ∂H/∂k₁ = 0');

    // ṗ_a = -∂H/∂x^a - p_a · RH
    // Since H = ω, RH = ∂H/∂A = 0
    assertApprox(X['omega'], 0, 1e-6, 'ω̇ = -∂H/∂t - ω·RH = 0');
}

{
    // Test dispersion relation Hamiltonian
    const grand = ET.grandManifold();
    const disp = ET.ThermodynamicHamiltonian.dispersionRelation(grand, 1, 0);

    const pt = grand.physicalPoint(
        0, 0, 0, 0, 0, 0,
        1, 0, 0,  // k = (1, 0, 0), |k| = 1
        1,        // ω = 1
        0, 0,
        0
    );

    // H = ω - c|k| = 1 - 1 = 0 (on shell)
    assertApprox(disp.evaluate(pt), 0, 1e-10, 'Dispersion: H = ω - |k| = 0 on shell');
}

// ============================================================================
// VII. LEGENDRIAN SUBMANIFOLDS
// ============================================================================

section('VII. Legendrian Submanifolds');

{
    const grand = ET.grandManifold();

    // Simple generating function: A(x) = k₀·q - ω₀·t
    // where k₀ = (1, 0, 0), ω₀ = 1 are constants
    const k0 = [1, 0, 0];
    const omega0 = 1;

    const genFunc = x => {
        return k0[0] * (x.q1 || 0) + k0[1] * (x.q2 || 0) + k0[2] * (x.q3 || 0) - omega0 * (x.t || 0);
    };

    const L = new ET.LegendrianSubmanifold(grand, genFunc);

    // Test lift
    const baseX = { q1: 2, q2: 0, q3: 0, t: 1, ell: 0, S: 0 };
    const lifted = L.lift(baseX);

    // A = k₀·q - ω₀·t = 1·2 - 1·1 = 1
    assertApprox(lifted.get('A'), 1, 1e-6, 'Legendrian lift: A = k₀·q - ω₀·t = 1');

    // p_a = ∂_a A
    assertApprox(lifted.get('k1'), k0[0], 1e-6, 'Legendrian lift: k₁ = ∂A/∂q¹ = 1');
    assertApprox(lifted.get('omega'), -omega0, 1e-6, 'Legendrian lift: ω = -∂A/∂t = -1 (sign!)');

    assert(L.verifyLegendrianCondition(baseX), 'Legendrian condition α|_L = 0 verified');
}

// ============================================================================
// VIII. GRAVITATIONAL EXTENSION
// ============================================================================

section('VIII. Gravitational Extension');

{
    // Minkowski metric
    const mink = ET.SpacetimeMetric.minkowski();
    const g = mink.covariant([0, 0, 0, 0]);

    assert(g[0][0] === 1, 'Minkowski: g₀₀ = +1 (timelike)');
    assert(g[1][1] === -1, 'Minkowski: g₁₁ = -1 (spacelike)');
    assert(g[2][2] === -1, 'Minkowski: g₂₂ = -1');
    assert(g[3][3] === -1, 'Minkowski: g₃₃ = -1');

    const gInv = mink.contravariant([0, 0, 0, 0]);
    assertApprox(gInv[0][0], 1, 1e-10, 'Minkowski inverse: g⁰⁰ = +1');
}

{
    // Schwarzschild metric
    const schw = ET.SpacetimeMetric.schwarzschild(1);
    const r = 10;  // Large radius
    const g = schw.covariant([0, r, Math.PI / 2, 0]);

    const f = 1 - 2 / r;  // 1 - 2M/r = 0.8
    assertApprox(g[0][0], f, 1e-10, `Schwarzschild: g₀₀ = 1 - 2M/r = ${f}`);
    assertApprox(g[1][1], -1 / f, 1e-10, `Schwarzschild: g₁₁ = -(1 - 2M/r)⁻¹`);
}

{
    // Relativistic Hamiltonian (mass shell)
    const metric = ET.SpacetimeMetric.minkowski();
    const relH = new ET.RelativisticHamiltonian(metric, 1, null, 0);

    // Test mass shell: H = ½(g^μν p_μ p_ν + m²) = 0
    // For m=1, p=(√2, 1, 0, 0), we have p² = 2 - 1 = 1 = m²
    const x = [0, 0, 0, 0];
    const p = [Math.sqrt(2), 1, 0, 0];

    const Hval = relH.evaluate(x, p);
    assertApprox(Hval, 0, 1e-10, 'Relativistic H = 0 on mass shell');
}

// ============================================================================
// IX. FLOW CONSERVATION
// ============================================================================

section('IX. Flow Integration');

{
    const grand = ET.grandManifold();

    // Free particle: H = ½|k|²
    const H = coords => {
        const k1 = coords.k1 || 0;
        const k2 = coords.k2 || 0;
        const k3 = coords.k3 || 0;
        return 0.5 * (k1 * k1 + k2 * k2 + k3 * k3);
    };

    const ham = new ET.ContactHamiltonian(grand, H);

    const pt = grand.physicalPoint(
        0, 0, 0, 0, 0, 0,
        1, 0, 0,  // k = (1, 0, 0)
        0, 0, 0,
        0
    );

    // Flow for 10 steps
    const trajectory = ham.flow(pt, 0.1, 10);

    assert(trajectory.length === 11, 'Flow produces correct number of points');

    // For free particle, k should be conserved
    const k1_init = trajectory[0].get('k1');
    const k1_final = trajectory[10].get('k1');
    assertApprox(k1_init, k1_final, 1e-6, 'Momentum k₁ conserved under free evolution');

    // q should evolve as q = q₀ + k·t
    const q1_final = trajectory[10].get('q1');
    const expected_q1 = 1.0;  // k₁ · dt · steps = 1 · 0.1 · 10 = 1
    assertApprox(q1_final, expected_q1, 0.1, 'Position q¹ evolves correctly');
}

// ============================================================================
// X. HOLOGRAPHIC MODEL EMERGENT SPACE
// ============================================================================

section('X. Holographic Emergent Space');

{
    const holo = ET.holographicManifold();

    // Create holographic point
    const pt = holo.holographicPoint(
        0,    // t
        0,    // ℓ
        1,    // S
        1,    // ω
        0,    // Δ
        1,    // T
        0     // A
    );

    assert(pt.get('t') === 0, 'Holographic point has time coordinate');
    assert(pt.get('S') === 1, 'Holographic point has entropy coordinate');

    // Emergent space function: q^i = a(t,ℓ,S) · xî
    const emergent = holo.createEmergentSpace(pt, (t, ell, S) => {
        const a = Math.exp(ell) * (1 + 0.1 * S);  // Simple emergence rule
        return [a, 0, 0];
    });

    assert(typeof emergent.q1 === 'number', 'Emergent q¹ is defined');
    assertApprox(emergent.q1, 1.1, 1e-10, 'Emergent q¹ = exp(ℓ)(1 + 0.1S) = 1.1');
}

// ============================================================================
// XI. GEOMETRIC ALGEBRA OPERATIONS
// ============================================================================

section('XI. Geometric Algebra Operations');

{
    const GA = require('../src/multivector.js');

    // Test Cl(3,0,0) - 3D Euclidean
    const cl3 = new GA.Algebra(3, 0, 0);
    assert(cl3.size === 8, 'Cl(3,0,0) has 2³ = 8 basis blades');

    // Test basis vectors
    const e1 = cl3.e(1);
    const e2 = cl3.e(2);
    const e3 = cl3.e(3);

    // e1² = +1 (positive signature)
    const e1_sq = e1.mul(e1).scalar();
    assertApprox(e1_sq, 1, 1e-10, 'e₁² = +1 in Cl(3,0,0)');

    // e1·e2 anticommutes: e1*e2 = -e2*e1
    const e1e2 = e1.mul(e2);
    const e2e1 = e2.mul(e1);
    const anticomm = e1e2.add(e2e1).scalar();
    assertApprox(anticomm, 0, 1e-10, 'e₁e₂ + e₂e₁ = 0 (anticommutation)');

    // Bivector classification
    const B = cl3.bivector([1, 0, 0]); // e1∧e2
    const Bsq = B.mul(B).scalar();
    assert(Bsq < 0, 'e₁∧e₂ has B² < 0 (elliptic)');

    const info = GA.classifyBivector(B);
    assert(info.type === 'elliptic', 'Bivector correctly classified as elliptic');

    // Rotor: R = cos(θ/2) + sin(θ/2)B̂
    const R = cl3.rotor(B, Math.PI / 2); // 90° rotation
    const RRrev = R.mul(R.reverse()).scalar();
    assertApprox(RRrev, 1, 1e-10, 'Rotor satisfies R~R = 1');
}

{
    // Test spacetime algebra Cl(1,3,0)
    const GA = require('../src/multivector.js');
    const sta = new GA.Algebra(1, 3, 0);

    assert(sta.size === 16, 'Cl(1,3,0) has 2⁴ = 16 basis blades');

    // e0² = +1 (timelike), e1² = -1 (spacelike)
    const e0_sq = sta.e(1).mul(sta.e(1)).scalar();  // First basis is e0
    const e1_sq = sta.e(2).mul(sta.e(2)).scalar();  // Second is e1 (spacelike)
    assertApprox(e0_sq, 1, 1e-10, 'e₀² = +1 (timelike)');
    assertApprox(e1_sq, -1, 1e-10, 'e₁² = -1 (spacelike)');
}

// ============================================================================
// XII. CHRISTOFFEL SYMBOLS
// ============================================================================

section('XII. Christoffel Symbols');

{
    // Minkowski: all Christoffel symbols should be zero
    const mink = ET.SpacetimeMetric.minkowski();
    const christoffel = new ET.ChristoffelSymbols(mink);

    const gamma = christoffel.computeAt([0, 0, 0, 0]);

    // Check that all components are zero
    let maxGamma = 0;
    for (let k = 0; k < 4; k++) {
        for (let i = 0; i < 4; i++) {
            for (let j = 0; j < 4; j++) {
                maxGamma = Math.max(maxGamma, Math.abs(gamma[k][i][j]));
            }
        }
    }
    assert(maxGamma < 1e-6, 'Minkowski: All Γᵏᵢⱼ ≈ 0');

    // Ricci scalar should be zero for flat spacetime
    const R = christoffel.ricciScalarAt([0, 0, 0, 0]);
    assertApprox(R, 0, 1e-4, 'Minkowski: Ricci scalar R = 0');
}

{
    // Schwarzschild: verify specific Christoffel symbol
    const schw = ET.SpacetimeMetric.schwarzschild(1);
    const christoffel = new ET.ChristoffelSymbols(schw);

    const r = 10;
    const theta = Math.PI / 2;
    const x = [0, r, theta, 0];

    const gamma = christoffel.computeAt(x);

    // Γ^r_tt = (M/r²)(1 - 2M/r) for M=1
    // At r=10: Γ^r_tt = (1/100)(1 - 0.2) = 0.008
    const expectedGamma_r_tt = (1 / (r * r)) * (1 - 2 / r);
    assertApprox(gamma[1][0][0], expectedGamma_r_tt, 0.001,
        `Schwarzschild: Γʳₜₜ ≈ ${expectedGamma_r_tt.toFixed(4)}`);

    // Γ^r_rr = -M/r²/(1-2M/r)
    const f = 1 - 2 / r;
    const expectedGamma_r_rr = -(1 / (r * r)) / f;
    assertApprox(gamma[1][1][1], expectedGamma_r_rr, 0.001,
        `Schwarzschild: Γʳᵣᵣ ≈ ${expectedGamma_r_rr.toFixed(4)}`);
}

// ============================================================================
// XIII. COVARIANT DERIVATIVE & PARALLEL TRANSPORT
// ============================================================================

section('XIII. Covariant Derivative & Parallel Transport');

{
    // Test covariant derivative on Minkowski
    const mink = ET.SpacetimeMetric.minkowski();
    const cov = new ET.CovariantDerivative(mink);

    // Constant vector field should have zero covariant derivative
    const V_const = x => [1, 0, 0, 0];
    const nablaV = cov.ofVector(V_const, [0, 0, 0, 0]);

    let maxNabla = 0;
    for (let i = 0; i < 4; i++) {
        for (let j = 0; j < 4; j++) {
            maxNabla = Math.max(maxNabla, Math.abs(nablaV[i][j]));
        }
    }
    assert(maxNabla < 1e-5, 'Minkowski: ∇V = 0 for constant vector field');
}

{
    // Parallel transport in flat space should preserve vectors
    const mink = ET.SpacetimeMetric.minkowski();
    const transport = new ET.ParallelTransport(mink);

    // Straight line in flat space
    const line = t => [t, t, 0, 0];  // Lightlike path
    const V0 = [1, 0, 0, 0];  // Initial vector

    const result = transport.transport(line, V0, 0, 1, 50);
    const V_final = result[result.length - 1].V;

    // In flat space, vector should be unchanged
    assertApprox(V_final[0], V0[0], 0.01, 'Flat transport: V⁰ preserved');
    assertApprox(V_final[1], V0[1], 0.01, 'Flat transport: V¹ preserved');
}

// ============================================================================
// XIV. DISCRETE GEOMETRIC CALCULUS
// ============================================================================

section('XIV. Discrete Geometric Calculus');

{
    const GC = require('../src/geometric-calculus.js');

    // Test gradient of f(x,y) = x
    const shape = [32, 32];
    const spacing = [0.1, 0.1];

    const f = new GC.ScalarField(shape, spacing);
    f.initFromFunction((x, y) => x);

    const nabla = new GC.SplitDifferentialOperator(shape, spacing);
    const grad_f = nabla.grad(f);

    // ∂f/∂x = 1, ∂f/∂y = 0
    const grad_at_center = grad_f.get(16, 16);
    assertApprox(grad_at_center[0], 1, 0.01, '∇(x): ∂f/∂x = 1');
    assertApprox(grad_at_center[1], 0, 0.01, '∇(x): ∂f/∂y = 0');
}

{
    const GC = require('../src/geometric-calculus.js');

    // Test Laplacian of f(x,y) = x² + y² should give 4
    const shape = [64, 64];
    const L = 2 * Math.PI;
    const spacing = [L / 64, L / 64];

    const f = new GC.ScalarField(shape, spacing);
    f.initFromFunction((x, y) => x * x + y * y);

    const nabla = new GC.SplitDifferentialOperator(shape, spacing);
    const lap_f = nabla.laplacian(f);

    // ∇²(x² + y²) = 2 + 2 = 4
    const lap_at_center = lap_f.get(32, 32);
    assertApprox(lap_at_center, 4, 0.1, '∇²(x² + y²) = 4');
}

{
    const GC = require('../src/geometric-calculus.js');

    // Test CFL estimation
    const shape = [32, 32];
    const spacing = [0.1, 0.1];

    const leapfrog = new GC.LeapfrogIntegrator(shape, spacing);
    const dt = leapfrog.estimateCFL(1.0);

    // CFL: dt < h / (c * √dim) = 0.1 / (1 * √2) ≈ 0.0707
    const expected_max = 0.1 / Math.sqrt(2);
    assert(dt < expected_max, `CFL: dt = ${dt.toFixed(4)} < ${expected_max.toFixed(4)}`);
    assert(dt > 0, 'CFL: dt > 0');
}

// ============================================================================
// X. FTGC MESH TESTS (Discrete Geometric Calculus on Triangle Meshes)
// ============================================================================

section('X. FTGC Mesh (Triangle Mesh Discrete Geometric Calculus)');

{
    // Test TriangleMesh creation and topology
    const mesh = ET.TriangleMesh.createGrid(8, 8, 1.0, 1.0);
    
    assert(mesh.nVertices === 64, `Grid mesh: ${mesh.nVertices} vertices (expected 64)`);
    assert(mesh.nFaces === 98, `Grid mesh: ${mesh.nFaces} triangles (expected 98)`);
    assert(mesh.nEdges > 0, `Grid mesh: ${mesh.nEdges} edges computed`);
    
    // Boundary vertices: 4 corners + 4 edges × 6 interior = 28
    const nBoundary = mesh.boundaryVertices.reduce((a, b) => a + b, 0);
    assert(nBoundary === 28, `Grid mesh: ${nBoundary} boundary vertices (expected 28)`);
}

{
    // Test icosahedron (closed mesh)
    const ico = ET.TriangleMesh.createIcosahedron(1.0);
    
    assert(ico.nVertices === 12, `Icosahedron: ${ico.nVertices} vertices (expected 12)`);
    assert(ico.nFaces === 20, `Icosahedron: ${ico.nFaces} faces (expected 20)`);
    
    // Euler: V - E + F = 2 → E = 12 - 20 + 2 = 30
    assert(ico.nEdges === 30, `Icosahedron: ${ico.nEdges} edges (expected 30)`);
    
    // Closed mesh: no boundary
    const nBoundary = ico.boundaryVertices.reduce((a, b) => a + b, 0);
    assert(nBoundary === 0, `Icosahedron: ${nBoundary} boundary vertices (expected 0)`);
}

{
    // Test FTGC identity: ∇∧(∇∧f) = 0 (curl of gradient is zero)
    const mesh = ET.TriangleMesh.createGrid(8, 8, 1.0, 1.0);
    const nabla = new ET.MeshGeometricDerivative(mesh);
    
    const f = new Float64Array(mesh.nVertices);
    for (let i = 0; i < mesh.nVertices; i++) {
        const x = mesh.vertices[i * 3];
        const y = mesh.vertices[i * 3 + 1];
        f[i] = Math.sin(Math.PI * x) * Math.cos(Math.PI * y);
    }
    
    const err = nabla.verifyCurlGradZero(f);
    assert(err < 1e-10, `FTGC: ∇∧(∇∧f) = 0, max error = ${err.toExponential(3)}`);
}

{
    // Test FTGC identity: ∇²(constant) = 0
    const mesh = ET.TriangleMesh.createGrid(8, 8, 1.0, 1.0);
    const nabla = new ET.MeshGeometricDerivative(mesh);
    
    const err = nabla.verifyLaplacianConstantZero(1.0);
    assert(err < 1e-10, `FTGC: ∇²(const) = 0, max error = ${err.toExponential(3)}`);
}

{
    // Test cotan weights are computed (most should be non-zero)
    const mesh = ET.TriangleMesh.createGrid(8, 8, 1.0, 1.0);
    const cotanWeights = mesh.cotanWeights();
    
    let nNonZero = 0;
    for (let i = 0; i < cotanWeights.length; i++) {
        if (Math.abs(cotanWeights[i]) > 1e-10) nNonZero++;
    }
    
    // Most edges should have non-zero cotan weights (boundary edges may have less)
    // On an 8x8 grid we expect >60% non-zero
    const ratio = nNonZero / mesh.nEdges;
    assert(ratio > 0.6, `Cotan weights: ${nNonZero}/${mesh.nEdges} non-zero (${(ratio*100).toFixed(0)}%)`);
}

{
    // Test mixed Voronoi dual areas sum to total mesh area
    const mesh = ET.TriangleMesh.createGrid(8, 8, 1.0, 1.0);
    const dualAreas = mesh.vertexDualAreas();
    const faceAreas = mesh.faceAreas();
    
    let totalDual = 0;
    for (let i = 0; i < dualAreas.length; i++) totalDual += dualAreas[i];
    
    let totalFace = 0;
    for (let i = 0; i < faceAreas.length; i++) totalFace += faceAreas[i];
    
    assertApprox(totalDual, totalFace, 1e-10, `Dual areas sum = face areas sum (${totalFace.toFixed(4)})`);
}

{
    // Test LeapfrogGCMesh heat equation conserves mass (with Neumann-like interior)
    const mesh = ET.TriangleMesh.createIcosahedron(1.0);
    const solver = new ET.LeapfrogGCMesh(mesh);
    
    // Initial condition: constant field
    const u0 = new Float64Array(mesh.nVertices).fill(1.0);
    const initialMass = solver.mass(u0);
    
    // Small time step, no boundary (closed mesh)
    const dt = solver.estimateCFLHeat(0.1);
    const u = solver.heatSimulate(u0, dt, 10, 0.1, {});
    
    const finalMass = solver.mass(u);
    assertApprox(finalMass, initialMass, 1e-8, `Heat eq conserves mass on closed mesh`);
}

{
    // Test Gauss-Bonnet on icosahedron: Σ curvature ≈ 4π
    const ico = ET.TriangleMesh.createIcosahedron(1.0);
    const dualAreas = ico.vertexDualAreas();
    
    // Angle defect at each vertex = 2π - Σ angles
    let totalAngleDefect = 0;
    for (let v = 0; v < ico.nVertices; v++) {
        let angleSum = 0;
        const tris = ico.vertexTriangles[v];
        
        for (const t of tris) {
            const i = ico.faces[t * 3];
            const j = ico.faces[t * 3 + 1];
            const k = ico.faces[t * 3 + 2];
            
            // Find angle at vertex v in this triangle
            let v0, v1, v2;
            if (i === v) { v0 = i; v1 = j; v2 = k; }
            else if (j === v) { v0 = j; v1 = k; v2 = i; }
            else { v0 = k; v1 = i; v2 = j; }
            
            const p0 = ico.getVertex(v0);
            const p1 = ico.getVertex(v1);
            const p2 = ico.getVertex(v2);
            
            const a = [p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]];
            const b = [p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]];
            const dot = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
            const normA = Math.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
            const normB = Math.sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
            
            const angle = Math.acos(Math.max(-1, Math.min(1, dot / (normA * normB))));
            angleSum += angle;
        }
        
        totalAngleDefect += (2 * Math.PI - angleSum);
    }
    
    // Gauss-Bonnet: Σ K_i = 2π χ = 2π × 2 = 4π for sphere
    assertApprox(totalAngleDefect, 4 * Math.PI, 0.01, `Gauss-Bonnet: Σ angle defect ≈ 4π (${(totalAngleDefect/Math.PI).toFixed(3)}π)`);
}

// ============================================================================
// SUMMARY
// ============================================================================

console.log('\n' + '═'.repeat(50));
console.log(`Tests completed: ${passed + failed}`);
console.log(`  Passed: ${passed}`);
console.log(`  Failed: ${failed}`);
console.log('═'.repeat(50));

if (failed > 0) {
    console.log('\n⚠ Some tests failed! Review above for details.');
    process.exit(1);
} else {
    console.log('\n✓ All tests passed! Framework validated.');
    process.exit(0);
}