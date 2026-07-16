/**
 * Root Re-export Shim Regression Tests
 *
 * Guards against the class of bug where each root-level src/*.js file was a
 * byte-for-byte copy of its src/<namespace>/ twin. Duplicated copies drift out
 * of sync (that is how the broken relative-require paths went unnoticed).
 *
 * The root files are now thin re-export shims. These tests assert that:
 *   1. Each root shim re-exports the EXACT same object as its namespaced twin
 *      (identity, not a structural copy) — so there is a single source of truth.
 *   2. index.js exposes the same ContactManifold class object that lives in
 *      src/contact/manifold.js (extracted, not re-implemented inline).
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
console.log('  ROOT RE-EXPORT SHIM REGRESSION TESTS');
console.log('========================================\n');

// ----------------------------------------------------------------------------
// 1. Each root shim === its namespaced implementation (single source of truth)
// ----------------------------------------------------------------------------
console.log('--- root shims re-export namespaced module identity ---\n');
{
    const pairs = [
        ['../src/multivector.js', '../src/algebra/multivector.js'],
        ['../src/geometric-calculus.js', '../src/calculus/grid.js'],
        ['../src/mesh-ftgc.js', '../src/calculus/mesh-derivative.js'],
        ['../src/mesh.js', '../src/calculus/mesh.js'],
        ['../src/mesh-solvers.js', '../src/calculus/solvers.js'],
        ['../src/riemannian-discrete.js', '../src/geometry/riemannian-discrete.js'],
        ['../src/riemannian-ga.js', '../src/geometry/riemannian-ga.js'],
        ['../src/geodesic-ga.js', '../src/geometry/geodesic.js'],
        ['../src/pilot-wave.js', '../src/physics/pilot-wave.js'],
        ['../src/entropic-gravity.js', '../src/physics/entropic-gravity.js'],
        ['../src/riemannian-spacetime.js', '../src/physics/spacetime.js']
    ];

    for (const [root, ns] of pairs) {
        const rootMod = require(root);
        const nsMod = require(ns);
        assert(rootMod === nsMod, `${root} re-exports identical object as ${ns}`);
    }
}

// ----------------------------------------------------------------------------
// 2. index.js consumes the extracted ContactManifold (no inline re-implementation)
// ----------------------------------------------------------------------------
console.log('\n--- index.js ContactManifold is the extracted contact/manifold.js class ---\n');
{
    const CT = require('../src/index.js');
    const Manifold = require('../src/contact/manifold.js');

    assert(CT.ContactManifold === Manifold.ContactManifold,
        'index.ContactManifold === contact/manifold.ContactManifold');
    assert(CT.ContactPoint === Manifold.ContactPoint,
        'index.ContactPoint === contact/manifold.ContactPoint');
    assert(CT.GrandContactManifold === Manifold.GrandContactManifold,
        'index.GrandContactManifold === contact/manifold.GrandContactManifold');
    assert(CT.HolographicContactManifold === Manifold.HolographicContactManifold,
        'index.HolographicContactManifold === contact/manifold.HolographicContactManifold');

    // The extracted class must still be fully functional through the public API.
    const g = CT.grandManifold();
    assert(g.dimension === 13, 'grandManifold() built from extracted class has dim 13');
    const h = CT.holographicManifold();
    assert(h.dimension === 7, 'holographicManifold() built from extracted class has dim 7');
    // GaugeExtendedManifold (defined in index.js) extends the imported class.
    const ge = CT.gaugeExtended();
    assert(ge instanceof CT.GrandContactManifold,
        'GaugeExtendedManifold still extends the imported GrandContactManifold');
}

// ----------------------------------------------------------------------------
console.log('\n========================================');
console.log(`  RESULTS: ${passed} passed, ${failed} failed`);
console.log('========================================');

process.exit(failed > 0 ? 1 : 0);
