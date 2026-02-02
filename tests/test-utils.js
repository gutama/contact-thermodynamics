/**
 * Unified Test Utilities for Contact Thermodynamics
 *
 * Provides consistent assertion functions and test tracking across all test files.
 *
 * @module test-utils
 * @license MIT
 */

const { abs, sqrt } = Math;

/**
 * Create a test tracker to count passed/failed tests.
 * @returns {{ passed: number, failed: number }}
 */
function createTestTracker() {
    return { passed: 0, failed: 0 };
}

/**
 * Assert a condition is true.
 * @param {object} tracker - Test tracker { passed, failed }
 * @param {boolean} condition - Condition to test
 * @param {string} message - Test description
 */
function assert(tracker, condition, message) {
    if (condition) {
        console.log(`  ✓ ${message}`);
        tracker.passed++;
    } else {
        console.log(`  ✗ ${message}`);
        tracker.failed++;
    }
}

/**
 * Assert two numbers are approximately equal.
 * @param {object} tracker - Test tracker
 * @param {number} actual - Actual value
 * @param {number} expected - Expected value
 * @param {number} tol - Tolerance
 * @param {string} message - Test description
 */
function assertApprox(tracker, actual, expected, tol, message) {
    const diff = abs(actual - expected);
    const passed = diff < tol;
    const details = `(got ${actual.toFixed(6)}, expected ${expected.toFixed(6)}, diff ${diff.toExponential(2)})`;
    assert(tracker, passed, `${message} ${details}`);
}

/**
 * Assert two arrays are approximately equal element-wise.
 * @param {object} tracker - Test tracker
 * @param {number[]} actual - Actual array
 * @param {number[]} expected - Expected array
 * @param {number} tol - Tolerance per element
 * @param {string} message - Test description
 */
function assertArrayApprox(tracker, actual, expected, tol, message) {
    if (actual.length !== expected.length) {
        assert(tracker, false, `${message} (length mismatch: ${actual.length} vs ${expected.length})`);
        return;
    }

    let maxDiff = 0;
    let maxIndex = 0;
    for (let i = 0; i < actual.length; i++) {
        const diff = abs(actual[i] - expected[i]);
        if (diff > maxDiff) {
            maxDiff = diff;
            maxIndex = i;
        }
    }

    const passed = maxDiff < tol;
    const details = passed
        ? `(max diff ${maxDiff.toExponential(2)} at index ${maxIndex})`
        : `(FAILED at index ${maxIndex}: got ${actual[maxIndex]}, expected ${expected[maxIndex]})`;
    assert(tracker, passed, `${message} ${details}`);
}

/**
 * Assert a value is within a range.
 * @param {object} tracker - Test tracker
 * @param {number} value - Value to test
 * @param {number} min - Minimum bound
 * @param {number} max - Maximum bound
 * @param {string} message - Test description
 */
function assertInRange(tracker, value, min, max, message) {
    const passed = value >= min && value <= max;
    const details = `(got ${value}, expected [${min}, ${max}])`;
    assert(tracker, passed, `${message} ${details}`);
}

/**
 * Print a section header.
 * @param {string} name - Section name
 */
function section(name) {
    console.log(`\n━━━ ${name} ━━━`);
}

/**
 * Print test summary and return exit code.
 * @param {object} tracker - Test tracker
 * @returns {number} Exit code (0 = all passed, 1 = some failed)
 */
function summary(tracker) {
    console.log('\n════════════════════════════════════════');
    console.log(`  Total: ${tracker.passed + tracker.failed}`);
    console.log(`  Passed: ${tracker.passed}`);
    console.log(`  Failed: ${tracker.failed}`);
    console.log('════════════════════════════════════════');

    if (tracker.failed > 0) {
        console.log('\n❌ SOME TESTS FAILED\n');
        return 1;
    } else {
        console.log('\n✅ ALL TESTS PASSED\n');
        return 0;
    }
}

/**
 * Helper: Compute L2 norm of an array.
 * @param {number[]} arr
 * @returns {number}
 */
function l2norm(arr) {
    let sum = 0;
    for (const x of arr) sum += x * x;
    return sqrt(sum);
}

/**
 * Helper: Compute relative error.
 * @param {number} actual
 * @param {number} expected
 * @returns {number}
 */
function relativeError(actual, expected) {
    if (abs(expected) < 1e-15) return abs(actual);
    return abs(actual - expected) / abs(expected);
}

module.exports = {
    createTestTracker,
    assert,
    assertApprox,
    assertArrayApprox,
    assertInRange,
    section,
    summary,
    l2norm,
    relativeError
};
