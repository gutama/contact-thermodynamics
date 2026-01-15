/**
 * Legendrian Submanifolds and Hamilton-Jacobi Theory
 * 
 * Demonstrates:
 * - Creating Legendrian submanifolds from generating functions
 * - Lifting base points to contact manifold
 * - Hamilton-Jacobi equation
 * - Wave mechanics interpretation
 */

const CT = require('../src/index.js');

console.log('═══════════════════════════════════════════════════════════════');
console.log('  Legendrian Submanifolds & Hamilton-Jacobi Theory');
console.log('═══════════════════════════════════════════════════════════════\n');

const M13 = CT.grandManifold();

// ============================================================================
// 1. Basic Legendrian: Plane Wave A = k₀·q - ω₀·t
// ============================================================================

console.log('1. PLANE WAVE LEGENDRIAN\n');

const k0 = [1, 0, 0];  // Wave vector
const omega0 = 1;       // Frequency

// Generating function: A(x) = k₀·q - ω₀·t (plane wave action)
const planeWaveAction = x => {
    return k0[0] * (x.q1 || 0) + 
           k0[1] * (x.q2 || 0) + 
           k0[2] * (x.q3 || 0) - 
           omega0 * (x.t || 0);
};

const L_plane = new CT.LegendrianSubmanifold(M13, planeWaveAction);

console.log('Generating function: A(x) = k₀·q - ω₀·t');
console.log('   k₀ = [' + k0.join(', ') + ']');
console.log('   ω₀ =', omega0);
console.log();

// Lift a base point
const basePoint = { q1: 2, q2: 0, q3: 0, t: 1, ell: 0, S: 0 };
const liftedPoint = L_plane.lift(basePoint);

console.log('Base point: (q¹, q², q³, t, ℓ, S) = (2, 0, 0, 1, 0, 0)');
console.log('Lifted to contact manifold:');
console.log('   A = k₀·q - ω₀·t =', liftedPoint.get('A').toFixed(4), '(expect 2·1 - 1·1 = 1)');
console.log('   k₁ = ∂A/∂q¹ =', liftedPoint.get('k1').toFixed(4), '(expect 1)');
console.log('   k₂ = ∂A/∂q² =', liftedPoint.get('k2').toFixed(4), '(expect 0)');
console.log('   ω = -∂A/∂t =', (-liftedPoint.get('omega')).toFixed(4), '(expect 1)');
console.log();

// Verify Legendrian condition
console.log('Legendrian condition α|_L = 0:',
    L_plane.verifyLegendrianCondition(basePoint) ? '✓ Satisfied' : '✗ Failed');
console.log();

// ============================================================================
// 2. Hamilton-Jacobi for Dispersion Relation
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('2. HAMILTON-JACOBI EQUATION\n');

// For H = ω - c|k| = 0 (massless wave)
// The HJ equation is: -∂A/∂t - c|∇A| = 0
// Solution: A = k₀·q - c|k₀|·t (plane wave with dispersion)

const c = 1;
const ham_disp = CT.ThermodynamicHamiltonian.dispersionRelation(M13, c, 0);

// Check HJ residual for our plane wave
const hjResidual = L_plane.hamiltonJacobiResidual(basePoint, ham_disp);
console.log('Hamilton-Jacobi equation: H(x, A(x), ∂A(x)) = 0');
console.log('   For plane wave with |k₀| = ω₀/c:');
console.log('   HJ residual =', Math.abs(hjResidual).toFixed(6), '(expect 0)');
console.log();

// ============================================================================
// 3. Spherical Wave
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('3. SPHERICAL WAVE LEGENDRIAN\n');

// A = -√(q₁² + q₂² + q₃²) (inward spherical wave)
const sphericalAction = x => {
    const r2 = (x.q1 || 0)**2 + (x.q2 || 0)**2 + (x.q3 || 0)**2;
    return -Math.sqrt(r2 + 0.01);  // Regularized at origin
};

const L_spherical = new CT.LegendrianSubmanifold(M13, sphericalAction);

// Sample points on the spherical wavefront
console.log('Generating function: A(x) = -|q| (spherical wavefront)');
console.log();
console.log('Wavefront samples:');
console.log('   Point (q¹,q²,q³)     r        A        |k|');
console.log('   ─────────────────────────────────────────────');

const samples = [
    [1, 0, 0],
    [0, 1, 0],
    [0.707, 0.707, 0],
    [2, 0, 0],
    [0, 0, 1]
];

for (const [q1, q2, q3] of samples) {
    const base = { q1, q2, q3, t: 0, ell: 0, S: 0 };
    const lifted = L_spherical.lift(base);
    const r = Math.sqrt(q1*q1 + q2*q2 + q3*q3);
    const kMag = Math.sqrt(
        lifted.get('k1')**2 + 
        lifted.get('k2')**2 + 
        lifted.get('k3')**2
    );
    console.log(`   (${q1.toFixed(2)}, ${q2.toFixed(2)}, ${q3.toFixed(2)})`.padEnd(20) +
                `${r.toFixed(3).padStart(8)}  ${lifted.get('A').toFixed(3).padStart(8)}  ${kMag.toFixed(3).padStart(8)}`);
}
console.log();
console.log('Note: |k| ≈ 1 everywhere (unit wavenumber for eikonal)');
console.log();

// ============================================================================
// 4. Eikonal Equation (Geometric Optics)
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('4. EIKONAL EQUATION (GEOMETRIC OPTICS)\n');

console.log('The eikonal equation |∇S|² = n²(x) relates to contact geometry:');
console.log('');
console.log('   - Legendrian: S = eikonal (optical path length)');
console.log('   - Momenta: k = ∇S (wave vector)');
console.log('   - Hamiltonian: H = ½(|k|² - n²) = 0');
console.log('');

// Eikonal for constant refractive index
const n = 1.5;  // Refractive index
const ham_eikonal = new CT.ContactHamiltonian(M13, coords => {
    const k1 = coords.k1 || 0;
    const k2 = coords.k2 || 0;
    const k3 = coords.k3 || 0;
    return 0.5 * (k1*k1 + k2*k2 + k3*k3 - n*n);
});

// Plane wave in medium: A = n·q₁ (propagating along x)
const eikonalAction = x => n * (x.q1 || 0);
const L_eikonal = new CT.LegendrianSubmanifold(M13, eikonalAction);

const eikonalBase = { q1: 1, q2: 0, q3: 0, t: 0, ell: 0, S: 0 };
const eikonalLifted = L_eikonal.lift(eikonalBase);

console.log('Plane wave in medium with n =', n);
console.log('   A = n·q₁');
console.log('   k₁ = ∂A/∂q₁ =', eikonalLifted.get('k1').toFixed(4), '(expect', n, ')');
console.log('   |k| =', Math.abs(eikonalLifted.get('k1')).toFixed(4));
console.log('   Eikonal residual H =', ham_eikonal.evaluate(eikonalLifted).toFixed(6), '(expect 0)');
console.log();

// ============================================================================
// 5. Thermodynamic Potential as Legendrian
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('5. THERMODYNAMIC POTENTIALS AS LEGENDRIANS\n');

console.log('In thermodynamics, potentials define Legendrian submanifolds:');
console.log('');
console.log('   Internal Energy:  U(S,V) → T = ∂U/∂S, -P = ∂U/∂V');
console.log('   Free Energy:      F(T,V) → S = -∂F/∂T, -P = ∂F/∂V');
console.log('   Gibbs Energy:     G(T,P) → S = -∂G/∂T, V = ∂G/∂P');
console.log();

// Simple ideal gas: U = (3/2)NkT for given S
// Let's use F = U - TS as generating function
// F(T, V) = NkT(3/2 - ln(V/V₀) - (3/2)ln(T/T₀))

const M7 = CT.holographicManifold();

// Simplified Helmholtz free energy
const helmholtzAction = x => {
    const T = x.T || 1;  // Using T as coordinate here is non-standard
    const S = x.S || 1;
    return 1.5 * T - T * S;  // F ≈ (3/2)kT - TS
};

console.log('Simplified Helmholtz-like action: A ≈ (3/2)T - TS');
console.log('(This is a pedagogical example - real thermodynamic geometry');
console.log(' requires careful treatment of extensive/intensive variables)');
console.log();

// ============================================================================
// 6. Sampling a Legendrian
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('6. SAMPLING A LEGENDRIAN SUBMANIFOLD\n');

// Random samples on the spherical wave Legendrian
const randomSampler = () => {
    const theta = Math.random() * Math.PI;
    const phi = Math.random() * 2 * Math.PI;
    const r = 1 + Math.random();  // r ∈ [1, 2]
    return {
        q1: r * Math.sin(theta) * Math.cos(phi),
        q2: r * Math.sin(theta) * Math.sin(phi),
        q3: r * Math.cos(theta),
        t: 0, ell: 0, S: 0
    };
};

const sampledPoints = L_spherical.sample(randomSampler, 5);

console.log('Random samples on spherical wave Legendrian:');
console.log('   #    q¹        q²        q³        A');
console.log('   ───────────────────────────────────────────');

for (let i = 0; i < sampledPoints.length; i++) {
    const p = sampledPoints[i];
    console.log(`   ${i+1}    ${p.get('q1').toFixed(3).padStart(7)}   ${p.get('q2').toFixed(3).padStart(7)}   ${p.get('q3').toFixed(3).padStart(7)}   ${p.get('A').toFixed(3).padStart(7)}`);
}
console.log();

// ============================================================================
// Summary
// ============================================================================

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('SUMMARY: LEGENDRIAN GEOMETRY\n');

console.log('Key relationships:');
console.log('');
console.log('   Generating function A(x)  ←→  Potential / Action / Phase');
console.log('   Momenta p = ∂A/∂x        ←→  Wave vector / Intensive vars');
console.log('   Legendrian α|_L = 0       ←→  dS = p·dq (action principle)');
console.log('   Hamilton-Jacobi H = 0     ←→  Dispersion / Equation of state');
console.log();

console.log('═══════════════════════════════════════════════════════════════');
console.log('  Legendrian examples complete!');
console.log('═══════════════════════════════════════════════════════════════');
