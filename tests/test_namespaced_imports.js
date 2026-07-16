/**
 * Namespaced Import Regression Tests
 *
 * Guards against the class of bug where modules under src/<namespace>/ were
 * copied from the flat src/*.js files without adjusting relative require()
 * depth. Those broken requires resolve to non-existent files and are silently
 * swallowed by catch{} fallbacks, disabling functionality at runtime instead
 * of failing loudly.
 *
 * These tests import the namespaced modules DIRECTLY (not the flat re-exports)
 * and exercise the code paths that depend on cross-namespace requires, so a
 * regressed require path surfaces as a test failure in CI.
 */

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

console.log('========================================');
console.log('  NAMESPACED IMPORT REGRESSION TESTS');
console.log('========================================\n');

// ----------------------------------------------------------------------------
// 1. geometry/riemannian-ga.js -> ../utils.js (cross3)
// ----------------------------------------------------------------------------
console.log('--- geometry/riemannian-ga.js: Utils.cross3 wiring ---\n');
{
    const RGA = require('../src/geometry/riemannian-ga.js');

    // Bivector3D.commutatorWithVector calls cross3 from Utils. If the Utils
    // require silently fell back to {}, cross3 is undefined and this throws.
    const B = new RGA.Bivector3D(0, 0, 1); // axial vector (0,0,1)
    let out = null;
    try {
        out = B.commutatorWithVector([1, 0, 0]);
    } catch (e) {
        out = null;
    }
    assert(Array.isArray(out), 'riemannian-ga loaded real Utils (cross3 callable, no {} fallback)');
    assert(out && out[0] === 0 && out[1] === 1 && out[2] === 0,
        'cross3 via Bivector3D gives (0,0,1)×(1,0,0) = (0,1,0)');
}

// ----------------------------------------------------------------------------
// 2. physics/pilot-wave.js -> ../geometry/riemannian-ga.js (curved-space guidance)
// ----------------------------------------------------------------------------
console.log('\n--- physics/pilot-wave.js: curved-space guidance wiring ---\n');
{
    const PW = require('../src/physics/pilot-wave.js');
    const RGA = require('../src/geometry/riemannian-ga.js');

    const manifold = new RGA.RiemannianManifold(2);
    const psi = new PW.WaveFunction(1.0, 0.0);
    const kernel = new PW.SmearingKernel(0.1, 'gaussian');
    const cs = new PW.CurvedSpacePilotWave(psi, manifold, kernel);

    // The `connection` getter requires('../geometry/riemannian-ga.js'). A wrong
    // depth would hit the catch{} and return null (guidance silently disabled).
    const conn = cs.connection;
    assert(conn !== null, 'CurvedSpacePilotWave.connection engaged (ConnectionBivector loaded, not catch fallback)');
    assert(conn && conn.constructor && conn.constructor.name === 'ConnectionBivector',
        'connection is a ConnectionBivector instance');
}

// ----------------------------------------------------------------------------
// 3. Every namespaced module loads and exports something
// ----------------------------------------------------------------------------
console.log('\n--- all namespaced modules load with real exports ---\n');
{
    const modules = [
        'algebra/multivector.js',
        'algebra/number-systems.js',
        'calculus/mesh.js',
        'calculus/grid.js',
        'calculus/mesh-derivative.js',
        'calculus/solvers.js',
        'contact/hamiltonian.js',
        'contact/legendrian.js',
        'contact/manifold.js',
        'geometry/geodesic.js',
        'geometry/riemannian-discrete.js',
        'geometry/riemannian-ga.js',
        'physics/entropic-gravity.js',
        'physics/spacetime.js',
        'physics/pilot-wave.js',
        'physics/g2-black-hole-remnant.js'
    ];
    for (const rel of modules) {
        let ok = false;
        try {
            const mod = require('../src/' + rel);
            ok = mod && Object.keys(mod).length > 0;
        } catch (e) {
            ok = false;
        }
        assert(ok, `src/${rel} loads with >0 exports`);
    }
}

// ----------------------------------------------------------------------------
// 4. physics/spacetime.js -> multivector (GA available, matInv usable)
// ----------------------------------------------------------------------------
console.log('\n--- physics/spacetime.js: multivector require wiring ---\n');
{
    const ST = require('../src/physics/spacetime.js');
    assert(ST && Object.keys(ST).length > 0, 'spacetime.js exports resolved (multivector require OK)');
    assert(typeof ST.SpacetimeManifoldGA === 'function' || typeof ST.SpacetimeManifoldGA === 'object',
        'SpacetimeManifoldGA is exported from spacetime.js');
}

// ----------------------------------------------------------------------------
// 5. SmearingKernel.isTimeDependent (typo fix regression guard)
// ----------------------------------------------------------------------------
console.log('\n--- SmearingKernel.isTimeDependent API ---\n');
{
    const PW = require('../src/physics/pilot-wave.js');
    const timeDep = new PW.SmearingKernel(t => 0.1 + 0.05 * t, 'gaussian');
    const fixed = new PW.SmearingKernel(0.1, 'gaussian');
    assert('isTimeDependent' in timeDep, 'getter is named isTimeDependent (typo fixed)');
    assert(!('isTimeDependant' in timeDep), 'misspelled isTimeDependant no longer present');
    assert(timeDep.isTimeDependent === true, 'time-dependent kernel reports isTimeDependent === true');
    assert(fixed.isTimeDependent === false, 'fixed-width kernel reports isTimeDependent === false');
}

// ----------------------------------------------------------------------------
console.log('\n========================================');
console.log(`  RESULTS: ${passed} passed, ${failed} failed`);
console.log('========================================');

process.exit(failed > 0 ? 1 : 0);
