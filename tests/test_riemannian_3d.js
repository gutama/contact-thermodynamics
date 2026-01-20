/**
 * Tests for 3D Riemannian Manifolds
 * 
 * Tests:
 * 1. Bivector4D operations
 * 2. 3-Sphere (S³): Constant positive curvature K = 1/R²
 * 3. 3-Torus (T³): Flat (K = 0)
 * 4. Hyperbolic 3-space (H³): Constant negative curvature K = -1
 * 5. Connection bivector computation for 3D manifolds
 */

const RiemannianGA = require('../src/riemannian-ga.js');

const {
    Bivector4D,
    TangentFrame,
    ConnectionBivector,
    Sphere3D,
    Torus3D,
    HyperbolicSpace3D,
    dot,
    norm,
    EPSILON
} = RiemannianGA;

const PI = Math.PI;

function assertEqual(actual, expected, tolerance, message) {
    const diff = Math.abs(actual - expected);
    if (diff > tolerance) {
        console.error(`FAIL: ${message}`);
        console.error(`  Expected: ${expected}`);
        console.error(`  Actual: ${actual}`);
        console.error(`  Difference: ${diff}`);
        return false;
    }
    console.log(`PASS: ${message}`);
    return true;
}

function assertTrue(condition, message) {
    if (!condition) {
        console.error(`FAIL: ${message}`);
        return false;
    }
    console.log(`PASS: ${message}`);
    return true;
}

// ============================================================================
// TEST 1: BIVECTOR4D OPERATIONS
// ============================================================================

function testBivector4D() {
    console.log('\n=== Test 1: Bivector4D Operations ===\n');

    let passed = true;

    // Test wedge product
    const u = [1, 0, 0, 0];
    const v = [0, 1, 0, 0];
    const B = Bivector4D.fromWedge(u, v);

    console.log(`e1 ∧ e2 = ${B.toString()}`);
    passed &= assertEqual(B.e12, 1, EPSILON, 'Wedge e1∧e2 has e12=1');
    passed &= assertEqual(B.norm(), 1, EPSILON, 'Wedge norm = 1');

    // Test commutator with vector
    const w = [1, 0, 0, 0];  // e1
    const Bxw = B.commutatorWithVector(w);
    console.log(`B × e1 = [${Bxw.map(x => x.toFixed(3)).join(', ')}]`);
    // B = e12 acts on e1 to give e2
    passed &= assertEqual(Bxw[1], 1, EPSILON, 'B × e1 = e2 (index 1)');

    // Test bivector addition/scaling
    const B2 = new Bivector4D(0, 1, 0, 0, 0, 0);  // e13
    const B_sum = B.add(B2);
    passed &= assertEqual(B_sum.e12, 1, EPSILON, 'Sum preserves e12');
    passed &= assertEqual(B_sum.e13, 1, EPSILON, 'Sum has e13');

    return passed;
}

// ============================================================================
// TEST 2: SPHERE3D (3-SPHERE)
// ============================================================================

function testSphere3D() {
    console.log('\n=== Test 2: 3-Sphere (S³) ===\n');

    let passed = true;
    const R = 1.0;
    const sphere = new Sphere3D(R);

    console.log(`Sphere3D: R = ${R}, dim = ${sphere.dim}, ambientDim = ${sphere.ambientDim}`);

    // Test embedding
    const coords = [PI / 4, PI / 3, PI / 6];
    const p = sphere.embedding(coords);
    const dist = norm(p);
    console.log(`Embedding at (π/4, π/3, π/6): [${p.map(x => x.toFixed(4)).join(', ')}]`);
    console.log(`Distance from origin: ${dist.toFixed(6)}`);
    passed &= assertEqual(dist, R, 0.01, 'Point lies on sphere (|p| = R)');

    // Test frame
    const frame = sphere.frame(coords);
    console.log(`Frame dimension: ${frame.dim}, ambient: ${frame.ambientDim}`);
    passed &= assertEqual(frame.dim, 3, EPSILON, 'Frame has 3 tangent vectors');
    passed &= assertEqual(frame.ambientDim, 4, EPSILON, 'Ambient dimension is 4');

    // Verify e^i · e_j = δ^i_j
    const recipTest = dot(frame.reciprocal[0], frame.frame[0]);
    const recipCross = dot(frame.reciprocal[0], frame.frame[1]);
    console.log(`e⁰·e₀ = ${recipTest.toFixed(6)}, e⁰·e₁ = ${recipCross.toFixed(6)}`);
    passed &= assertEqual(recipTest, 1, 0.01, 'Reciprocal frame: e⁰·e₀ = 1');
    passed &= assertEqual(recipCross, 0, 0.1, 'Reciprocal frame: e⁰·e₁ ≈ 0');

    // Test metric
    const g = frame.metric();
    console.log(`Metric diagonal: [${g[0][0].toFixed(4)}, ${g[1][1].toFixed(4)}, ${g[2][2].toFixed(4)}]`);
    // For S³ with R=1, at this point g_00 = R² = 1
    passed &= assertEqual(g[0][0], R * R, 0.1, 'g_χχ = R²');

    // Theoretical curvature
    const K_theory = sphere.theoreticalCurvature();
    console.log(`Theoretical sectional curvature K = ${K_theory.toFixed(6)}`);
    passed &= assertEqual(K_theory, 1 / (R * R), EPSILON, 'K = 1/R² for unit sphere');

    // Test connection bivector
    const connection = new ConnectionBivector(sphere);
    const omegas = connection.computeAt(coords);
    console.log(`Connection bivectors computed: ${omegas.length}`);
    passed &= assertEqual(omegas.length, 3, EPSILON, 'Three connection bivectors (one per coord)');

    // Check bivector type
    passed &= assertTrue(omegas[0] instanceof Bivector4D, 'Connection is Bivector4D');

    return passed;
}

// ============================================================================
// TEST 3: TORUS3D (3-TORUS)
// ============================================================================

function testTorus3D() {
    console.log('\n=== Test 3: 3-Torus (T³) ===\n');

    let passed = true;
    const torus = new Torus3D(1.0, 1.0, 1.0);

    console.log(`Torus3D: r1=r2=r3=1, dim = ${torus.dim}, ambientDim = ${torus.ambientDim}`);

    // Test at a generic point
    const coords = [1.0, 2.0, 3.0];
    const frame = torus.frame(coords);

    console.log(`Frame dimension: ${frame.dim}`);
    passed &= assertEqual(frame.dim, 3, EPSILON, 'Frame has 3 tangent vectors');

    // For flat torus, frame should be orthogonal
    const e1 = frame.frame[0];
    const e2 = frame.frame[1];
    const e12_dot = dot(e1, e2);
    console.log(`e1·e2 = ${e12_dot.toFixed(6)} (should be 0 for orthogonal)`);
    passed &= assertEqual(e12_dot, 0, 0.01, 'Tangent vectors orthogonal (flat)');

    // Theoretical curvature
    const K_theory = torus.theoreticalCurvature();
    console.log(`Theoretical curvature K = ${K_theory}`);
    passed &= assertEqual(K_theory, 0, EPSILON, 'K = 0 (flat)');

    return passed;
}

// ============================================================================
// TEST 4: HYPERBOLIC 3-SPACE
// ============================================================================

function testHyperbolicSpace3D() {
    console.log('\n=== Test 4: Hyperbolic 3-Space (H³) ===\n');

    let passed = true;
    const H3 = new HyperbolicSpace3D();

    console.log(`HyperbolicSpace3D: dim = ${H3.dim}, ambientDim = ${H3.ambientDim}`);

    // Test at point (0, 0, 1) in upper half-space
    const coords = [0, 0, 1];
    const frame = H3.frame(coords);

    console.log(`Frame at z=1: dim = ${frame.dim}`);
    passed &= assertEqual(frame.dim, 3, EPSILON, 'Frame has 3 vectors');

    // Check metric (should be identity at z=1)
    const g = frame.metric();
    console.log(`Metric at z=1: diag = [${g[0][0].toFixed(4)}, ${g[1][1].toFixed(4)}, ${g[2][2].toFixed(4)}]`);
    passed &= assertEqual(g[0][0], 1, 0.01, 'g_xx = 1/z² = 1 at z=1');

    // Test at different z
    const coords2 = [0, 0, 2];
    const frame2 = H3.frame(coords2);
    const g2 = frame2.metric();
    console.log(`Metric at z=2: diag = [${g2[0][0].toFixed(4)}, ${g2[1][1].toFixed(4)}, ${g2[2][2].toFixed(4)}]`);
    passed &= assertEqual(g2[0][0], 0.25, 0.01, 'g_xx = 1/z² = 0.25 at z=2');

    // Theoretical curvature
    const K_theory = H3.theoreticalCurvature();
    console.log(`Theoretical sectional curvature K = ${K_theory}`);
    passed &= assertEqual(K_theory, -1, EPSILON, 'K = -1 (constant negative)');

    const R_theory = H3.theoreticalScalarCurvature();
    console.log(`Theoretical scalar curvature R = ${R_theory}`);
    passed &= assertEqual(R_theory, -6, EPSILON, 'R = -6');

    return passed;
}

// ============================================================================
// TEST 5: CONNECTION BIVECTOR IN 4D
// ============================================================================

function testConnection4D() {
    console.log('\n=== Test 5: Connection Bivector in 4D ===\n');

    let passed = true;
    const sphere = new Sphere3D(1.0);
    const connection = new ConnectionBivector(sphere);

    // Test at a point away from singularities
    const coords = [PI / 3, PI / 4, PI / 5];
    const omegas = connection.computeAt(coords);

    console.log('Connection bivectors at (π/3, π/4, π/5):');
    for (let i = 0; i < omegas.length; i++) {
        const omega = omegas[i];
        console.log(`  ω_${i}: norm = ${omega.norm().toFixed(6)}`);
        passed &= assertTrue(omega instanceof Bivector4D, `ω_${i} is Bivector4D`);
    }

    // Test that connection bivectors are non-zero (sphere has curvature)
    const omega0_norm = omegas[0].norm();
    passed &= assertTrue(omega0_norm > 0.01, 'ω_0 is non-zero (curved manifold)');

    // Test connection along a direction
    const u = [1, 0, 0];  // χ direction
    const omega_u = connection.along(coords, u);
    console.log(`ω(∂χ): norm = ${omega_u.norm().toFixed(6)}`);
    passed &= assertTrue(omega_u instanceof Bivector4D, 'ω(u) is Bivector4D');

    return passed;
}

// ============================================================================
// MAIN
// ============================================================================

function runAllTests() {
    console.log('========================================');
    console.log('  3D Riemannian Manifolds Tests');
    console.log('========================================');

    let passed = 0;
    let failed = 0;

    const tests = [
        ['Bivector4D Operations', testBivector4D],
        ['3-Sphere (S³)', testSphere3D],
        ['3-Torus (T³)', testTorus3D],
        ['Hyperbolic 3-Space (H³)', testHyperbolicSpace3D],
        ['Connection in 4D', testConnection4D]
    ];

    for (const [name, testFn] of tests) {
        try {
            const result = testFn();
            if (result) passed++;
            else failed++;
        } catch (err) {
            console.error(`\nERROR in test "${name}":`, err.message);
            console.error(err.stack);
            failed++;
        }
    }

    console.log('\n========================================');
    console.log(`  Results: ${passed} passed, ${failed} failed`);
    console.log('========================================');

    process.exit(failed > 0 ? 1 : 0);
}

runAllTests();
