/**
 * Tests for Number Systems Module
 * 
 * Tests Complex, Dual, Hyperbolic number systems and bivector classification.
 * Added as part of project restructure using geometric-algebra skill.
 */

const NumberSystems = require('../src/algebra/number-systems');

const { Complex, Dual, Hyperbolic, BivectorType, classifyBivectorSquare, expBivectorCoeffs } = NumberSystems;

// Test utilities
let passed = 0, failed = 0;
function test(name, condition) {
    if (condition) {
        console.log(`  ✓ ${name}`);
        passed++;
    } else {
        console.log(`  ✗ ${name}`);
        failed++;
    }
}

function approx(a, b, tol = 1e-10) {
    return Math.abs(a - b) < tol;
}

console.log('\n━━━ Number Systems Tests ━━━\n');

// ============================================================================
// COMPLEX NUMBERS
// ============================================================================

console.log('--- Complex Numbers (i² = -1) ---');

const z1 = new Complex(1, 2);
const z2 = new Complex(3, -1);

test('Complex addition', approx(z1.add(z2).a, 4) && approx(z1.add(z2).b, 1));
test('Complex multiplication', approx(z1.mul(z2).a, 5) && approx(z1.mul(z2).b, 5));
test('Complex conjugate', approx(z1.conjugate().b, -2));
test('Complex norm squared', approx(z1.normSq(), 5));

// Euler's formula: e^(iπ) = -1
const eiPi = Complex.exp(Math.PI);
test('Euler formula: e^(iπ) = -1', approx(eiPi.a, -1) && approx(eiPi.b, 0));

// Unit circle
const unit = Complex.unitCurve(Math.PI / 4);
test('Unit circle |e^(iπ/4)| = 1', approx(unit.normSq(), 1));

// ============================================================================
// DUAL NUMBERS
// ============================================================================

console.log('\n--- Dual Numbers (ε² = 0) ---');

const d1 = new Dual(2, 3);
const d2 = new Dual(4, 1);

test('Dual multiplication', approx(d1.mul(d2).a, 8) && approx(d1.mul(d2).b, 14));
test('Dual series terminates: exp(ε) = 1 + ε', approx(Dual.exp(1).a, 1) && approx(Dual.exp(1).b, 1));

// Automatic differentiation: f(x) = x², f'(x) = 2x
const f = x => x.mul(x);
const result = f(new Dual(3, 1));
test('Auto-diff: (x²)\' at x=3 is 6', approx(result.a, 9) && approx(result.b, 6));

// Power rule
const d3 = new Dual(2, 1);
const d3cubed = d3.pow(3);
test('Dual power: (2+ε)³ = 8 + 12ε', approx(d3cubed.a, 8) && approx(d3cubed.b, 12));

// ============================================================================
// HYPERBOLIC NUMBERS
// ============================================================================

console.log('\n--- Hyperbolic Numbers (j² = +1) ---');

const h1 = new Hyperbolic(3, 1);
const h2 = new Hyperbolic(2, 1);

test('Hyperbolic multiplication', approx(h1.mul(h2).a, 7) && approx(h1.mul(h2).b, 5));
test('Minkowski norm: 3² - 1² = 8', approx(h1.normSq(), 8));
test('Timelike classification', h1.isTimelike());

// Hyperbolic exponential: exp(jφ) = cosh(φ) + j·sinh(φ)
const hExp = Hyperbolic.exp(1);
test('Hyperbolic exp', approx(hExp.a, Math.cosh(1)) && approx(hExp.b, Math.sinh(1)));

// Unit hyperbola
const unitH = Hyperbolic.unitCurve(0.5);
test('Unit hyperbola: cosh²-sinh² = 1', approx(unitH.normSq(), 1));

// Null vector
const nullVec = Hyperbolic.nullPlus(2);
test('Null vector on lightcone', nullVec.isNull());

// Velocity addition (relativistic)
const v_sum = Hyperbolic.velocityAdd(0.5, 0.5);
test('Relativistic velocity addition', approx(v_sum, 0.8)); // (0.5+0.5)/(1+0.25) = 0.8

// ============================================================================
// BIVECTOR CLASSIFICATION
// ============================================================================

console.log('\n--- Bivector Classification ---');

test('Elliptic: B² = -1', classifyBivectorSquare(-1) === BivectorType.ELLIPTIC);
test('Parabolic: B² = 0', classifyBivectorSquare(0) === BivectorType.PARABOLIC);
test('Hyperbolic: B² = +1', classifyBivectorSquare(1) === BivectorType.HYPERBOLIC);
test('Mixed: B² = 0.5', classifyBivectorSquare(0.5) === BivectorType.MIXED);

// Exponential coefficients
const ellipticCoeffs = expBivectorCoeffs(Math.PI, BivectorType.ELLIPTIC);
test('Elliptic exp: cos(π/2) = 0', approx(ellipticCoeffs.s, 0));
test('Elliptic exp: sin(π/2) = 1', approx(ellipticCoeffs.b, 1));

const parabolicCoeffs = expBivectorCoeffs(2, BivectorType.PARABOLIC);
test('Parabolic exp: s = 1', approx(parabolicCoeffs.s, 1));
test('Parabolic exp: b = θ/2', approx(parabolicCoeffs.b, 1));

const hyperbolicCoeffs = expBivectorCoeffs(2, BivectorType.HYPERBOLIC);
test('Hyperbolic exp: s = cosh(1)', approx(hyperbolicCoeffs.s, Math.cosh(1)));
test('Hyperbolic exp: b = sinh(1)', approx(hyperbolicCoeffs.b, Math.sinh(1)));

// ============================================================================
// SUMMARY
// ============================================================================

console.log('\n══════════════════════════════════════════════════');
console.log(`Number Systems Tests: ${passed} passed, ${failed} failed`);
console.log('══════════════════════════════════════════════════\n');

process.exit(failed > 0 ? 1 : 0);
