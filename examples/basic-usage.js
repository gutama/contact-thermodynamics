/**
 * Basic Usage Examples
 * 
 * Getting started with Contact Thermodynamics
 */

const CT = require('../src/index.js');

console.log('═══════════════════════════════════════════════════════════════');
console.log('  Contact Thermodynamics - Basic Usage Examples');
console.log('═══════════════════════════════════════════════════════════════\n');

// ============================================================================
// 1. Creating Contact Manifolds
// ============================================================================

console.log('1. CREATING CONTACT MANIFOLDS\n');

// Grand Model M₁₃
const M13 = CT.grandManifold();
console.log('Grand Manifold:', M13.toString());
console.log('   Base coordinates:', M13.baseCoords.join(', '));
console.log('   Momenta:', M13.momentaCoords.join(', '));
console.log('   Fiber:', M13.fiberCoord);
console.log('   Dimension:', M13.dim);
console.log();

// Holographic Model M₇
const M7 = CT.holographicManifold();
console.log('Holographic Manifold:', M7.toString());
console.log('   Dimension:', M7.dim);
console.log();

// Gauge-Extended Model M₁₅
const M15 = CT.gaugeExtended();
console.log('Gauge-Extended Manifold:', M15.toString());
console.log('   Dimension:', M15.dim);
console.log();

// ============================================================================
// 2. Creating Points on Manifolds
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('2. CREATING POINTS\n');

// Using physicalPoint for Grand manifold
const pt = M13.physicalPoint(
    1, 0, 0,      // q¹, q², q³ - spatial position
    0,            // t - time
    0,            // ℓ = log(λ) - scale
    1,            // S - entropy
    0.5, 0, 0,    // k₁, k₂, k₃ - wavenumber
    1,            // ω - frequency
    0,            // Δ - dilatation
    1,            // T - temperature
    0             // A - action
);

console.log('Physical point on M₁₃:');
console.log('   Position:', M13.spatialPosition(pt));
console.log('   Wavenumber:', M13.waveVector(pt));
console.log('   Time:', pt.get('t'));
console.log('   Entropy:', pt.get('S'));
console.log('   Temperature:', pt.get('T'));
console.log('   Action:', pt.get('A'));
console.log();

// Using generic point() method
const pt2 = M13.point({
    q1: 2, q2: 1, q3: 0,
    t: 0.5,
    ell: 0.1,
    S: 2,
    k1: 1, k2: 0.5, k3: 0,
    omega: 2,
    Delta: 0,
    T: 0.8,
    A: 0.5
});

console.log('Generic point construction:');
console.log(pt2.toString());
console.log();

// ============================================================================
// 3. Contact Form
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('3. CONTACT FORM\n');

console.log('Canonical contact form:');
console.log('   α =', M13.contactFormSymbolic());
console.log();

// Verify contact condition: α ∧ (dα)⁶ ≠ 0
const volume = M13.verifyContactCondition(pt);
console.log('Contact non-degeneracy:');
console.log('   α ∧ (dα)⁶ =', volume, '(= 6!)');
console.log();

// Evaluate contact form on tangent vectors
// Legendrian tangent: α(v) = 0 when du = p·dx
const legendrianTangent = {
    q1: 1, q2: 0, q3: 0,
    t: 0, ell: 0, S: 0,
    k1: 0, k2: 0, k3: 0,
    omega: 0, Delta: 0, T: 0,
    A: 0.5  // = k₁ · δq¹
};
console.log('Contact form on Legendrian tangent:');
console.log('   α(v) =', M13.evaluateContactForm(pt, legendrianTangent), '(should be 0)');

// Vertical tangent: α(∂/∂A) = 1
const verticalTangent = { A: 1 };
console.log('   α(∂/∂A) =', M13.evaluateContactForm(pt, verticalTangent), '(should be 1)');
console.log();

// ============================================================================
// 4. Reeb Vector Field
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('4. REEB VECTOR FIELD\n');

const R = M13.reebField(pt);
console.log('Reeb field R = ∂/∂A:');
console.log('   R components:', JSON.stringify(R));
console.log('   α(R) =', M13.evaluateContactForm(pt, R), '(should be 1)');
console.log();

// ============================================================================
// 5. Holographic Model with Emergent Space
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('5. HOLOGRAPHIC MODEL - EMERGENT SPACE\n');

const holoPt = M7.holographicPoint(
    0,    // t
    0.5,  // ℓ
    1,    // S
    1,    // ω
    0,    // Δ
    1,    // T
    0     // A
);

console.log('Holographic point (t, ℓ, S, ω, Δ, T, A):');
console.log('   t =', holoPt.get('t'));
console.log('   ℓ =', holoPt.get('ell'));
console.log('   S =', holoPt.get('S'));
console.log();

// Define emergent space: q = a(t,ℓ,S) · x̂
const emergent = M7.createEmergentSpace(holoPt, (t, ell, S) => {
    const scale = Math.exp(ell);
    return [
        scale * Math.cos(t),
        scale * Math.sin(t),
        0.1 * S
    ];
});

console.log('Emergent spatial configuration:');
console.log('   q¹ =', emergent.q1.toFixed(4));
console.log('   q² =', emergent.q2.toFixed(4));
console.log('   q³ =', emergent.q3.toFixed(4));
console.log();

// ============================================================================
// 6. Summary Table
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('6. SUMMARY TABLE\n');

const table = CT.summaryTable();
console.log('Theory Comparison:');
console.log();

// Print header
console.log('  ' + table.columns.map(c => c.padEnd(20)).join(''));
console.log('  ' + '─'.repeat(100));

// Print rows
for (const row of table.rows) {
    console.log('  ' + row.map(c => c.padEnd(20)).join(''));
}
console.log();

console.log('═══════════════════════════════════════════════════════════════');
console.log('  Examples complete!');
console.log('═══════════════════════════════════════════════════════════════');
