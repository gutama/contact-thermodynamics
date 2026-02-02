/**
 * Shared Utilities for Contact Thermodynamics
 *
 * Common vector, matrix, and numerical operations used across modules.
 * Consolidates duplicate code from riemannian-ga, geodesic-ga, riemannian-discrete, etc.
 *
 * @module utils
 * @license MIT
 */

(function (global) {
    'use strict';

    // ============================================================================
    // CONSTANTS
    // ============================================================================

    const EPSILON = 1e-10;
    const { abs, sqrt, sin, cos, tan, atan2, acos, PI, min, max, exp, log } = Math;

    // ============================================================================
    // SCALAR UTILITIES
    // ============================================================================

    /**
     * Check if a value is effectively zero within tolerance.
     * @param {number} x - Value to check
     * @param {number} [tol=EPSILON] - Tolerance
     * @returns {boolean}
     */
    function isZero(x, tol = EPSILON) {
        return abs(x) < tol;
    }

    /**
     * Clamp a value to a range.
     * @param {number} x - Value to clamp
     * @param {number} lo - Lower bound
     * @param {number} hi - Upper bound
     * @returns {number}
     */
    function clamp(x, lo, hi) {
        return min(max(x, lo), hi);
    }

    // ============================================================================
    // VECTOR UTILITIES (general n-dimensional)
    // ============================================================================

    /**
     * Dot product of two vectors.
     * @param {number[]} u
     * @param {number[]} v
     * @returns {number}
     */
    function dot(u, v) {
        let sum = 0;
        for (let i = 0; i < u.length; i++) sum += u[i] * v[i];
        return sum;
    }

    /**
     * Euclidean norm (length) of a vector.
     * @param {number[]} v
     * @returns {number}
     */
    function norm(v) {
        return sqrt(dot(v, v));
    }

    /**
     * Squared norm of a vector (avoids sqrt).
     * @param {number[]} v
     * @returns {number}
     */
    function normSq(v) {
        return dot(v, v);
    }

    /**
     * Normalize a vector to unit length.
     * @param {number[]} v
     * @returns {number[]} Normalized vector (or copy if near-zero)
     */
    function normalize(v) {
        const n = norm(v);
        return n > EPSILON ? v.map(x => x / n) : v.slice();
    }

    /**
     * Vector addition.
     * @param {number[]} u
     * @param {number[]} v
     * @returns {number[]}
     */
    function vecAdd(u, v) {
        return u.map((x, i) => x + v[i]);
    }

    /**
     * Vector subtraction.
     * @param {number[]} u
     * @param {number[]} v
     * @returns {number[]}
     */
    function vecSub(u, v) {
        return u.map((x, i) => x - v[i]);
    }

    /**
     * Scalar multiplication of a vector.
     * @param {number[]} v
     * @param {number} s - Scalar
     * @returns {number[]}
     */
    function vecScale(v, s) {
        return v.map(x => x * s);
    }

    /**
     * Linear combination: a*u + b*v.
     * @param {number} a
     * @param {number[]} u
     * @param {number} b
     * @param {number[]} v
     * @returns {number[]}
     */
    function vecCombine(a, u, b, v) {
        return u.map((x, i) => a * x + b * v[i]);
    }

    // ============================================================================
    // 3D VECTOR UTILITIES
    // ============================================================================

    /**
     * Cross product (3D only).
     * @param {number[]} u - 3D vector
     * @param {number[]} v - 3D vector
     * @returns {number[]} u × v
     */
    function cross3(u, v) {
        return [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0]
        ];
    }

    /**
     * 3D dot product (optimized, no loop).
     * @param {number[]} u
     * @param {number[]} v
     * @returns {number}
     */
    function dot3(u, v) {
        return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    }

    /**
     * 3D norm (optimized, no loop).
     * @param {number[]} v
     * @returns {number}
     */
    function norm3(v) {
        return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }

    /**
     * 3D normalize (optimized).
     * @param {number[]} v
     * @returns {number[]}
     */
    function normalize3(v) {
        const n = norm3(v);
        return n > EPSILON ? [v[0] / n, v[1] / n, v[2] / n] : [0, 0, 0];
    }

    /**
     * 3D vector addition (optimized).
     * @param {number[]} u
     * @param {number[]} v
     * @returns {number[]}
     */
    function vecAdd3(u, v) {
        return [u[0] + v[0], u[1] + v[1], u[2] + v[2]];
    }

    /**
     * 3D vector subtraction (optimized).
     * @param {number[]} u
     * @param {number[]} v
     * @returns {number[]}
     */
    function vecSub3(u, v) {
        return [u[0] - v[0], u[1] - v[1], u[2] - v[2]];
    }

    /**
     * 3D scalar multiplication (optimized).
     * @param {number[]} v
     * @param {number} s
     * @returns {number[]}
     */
    function vecScale3(v, s) {
        return [v[0] * s, v[1] * s, v[2] * s];
    }

    // ============================================================================
    // MATRIX UTILITIES
    // ============================================================================

    /**
     * Matrix-vector multiplication.
     * @param {number[][]} A - n×m matrix
     * @param {number[]} v - m-vector
     * @returns {number[]} n-vector
     */
    function matVecMul(A, v) {
        const n = A.length;
        const result = new Array(n).fill(0);
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < v.length; j++) {
                result[i] += A[i][j] * v[j];
            }
        }
        return result;
    }

    /**
     * Matrix-matrix multiplication.
     * @param {number[][]} A - n×m matrix
     * @param {number[][]} B - m×p matrix
     * @returns {number[][]} n×p matrix
     */
    function matMul(A, B) {
        const n = A.length;
        const m = B.length;
        const p = B[0].length;
        const C = [];
        for (let i = 0; i < n; i++) {
            C[i] = new Array(p).fill(0);
            for (let j = 0; j < p; j++) {
                for (let k = 0; k < m; k++) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return C;
    }

    /**
     * Matrix transpose.
     * @param {number[][]} A - n×m matrix
     * @returns {number[][]} m×n matrix
     */
    function transpose(A) {
        const n = A.length;
        const m = A[0].length;
        const T = [];
        for (let j = 0; j < m; j++) {
            T[j] = new Array(n);
            for (let i = 0; i < n; i++) {
                T[j][i] = A[i][j];
            }
        }
        return T;
    }

    /**
     * Create identity matrix.
     * @param {number} n - Dimension
     * @returns {number[][]}
     */
    function identity(n) {
        const I = [];
        for (let i = 0; i < n; i++) {
            I[i] = new Array(n).fill(0);
            I[i][i] = 1;
        }
        return I;
    }

    // ============================================================================
    // DETERMINANTS
    // ============================================================================

    /**
     * Determinant of 2×2 matrix.
     * @param {number[][]} A
     * @returns {number}
     */
    function det2x2(A) {
        return A[0][0] * A[1][1] - A[0][1] * A[1][0];
    }

    /**
     * Determinant of 3×3 matrix.
     * @param {number[][]} A
     * @returns {number}
     */
    function det3x3(A) {
        return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
             - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
             + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
    }

    // ============================================================================
    // MATRIX INVERSION
    // ============================================================================

    /**
     * Inverse of 2×2 matrix.
     * @param {number[][]} A
     * @returns {number[][]}
     * @throws {Error} If matrix is singular
     */
    function inv2x2(A) {
        const d = det2x2(A);
        if (abs(d) < EPSILON) throw new Error('Singular matrix');
        return [
            [A[1][1] / d, -A[0][1] / d],
            [-A[1][0] / d, A[0][0] / d]
        ];
    }

    /**
     * Inverse of 3×3 matrix.
     * @param {number[][]} A
     * @returns {number[][]}
     * @throws {Error} If matrix is singular
     */
    function inv3x3(A) {
        const d = det3x3(A);
        if (abs(d) < EPSILON) throw new Error('Singular matrix');

        // Cofactor matrix
        const C = [
            [
                A[1][1] * A[2][2] - A[1][2] * A[2][1],
                -(A[1][0] * A[2][2] - A[1][2] * A[2][0]),
                A[1][0] * A[2][1] - A[1][1] * A[2][0]
            ],
            [
                -(A[0][1] * A[2][2] - A[0][2] * A[2][1]),
                A[0][0] * A[2][2] - A[0][2] * A[2][0],
                -(A[0][0] * A[2][1] - A[0][1] * A[2][0])
            ],
            [
                A[0][1] * A[1][2] - A[0][2] * A[1][1],
                -(A[0][0] * A[1][2] - A[0][2] * A[1][0]),
                A[0][0] * A[1][1] - A[0][1] * A[1][0]
            ]
        ];

        // Transpose of cofactor / determinant
        return [
            [C[0][0] / d, C[1][0] / d, C[2][0] / d],
            [C[0][1] / d, C[1][1] / d, C[2][1] / d],
            [C[0][2] / d, C[1][2] / d, C[2][2] / d]
        ];
    }

    /**
     * General matrix inverse using Gauss-Jordan elimination.
     * Works for any n×n matrix.
     * @param {number[][]} A
     * @returns {number[][]}
     * @throws {Error} If matrix is singular
     */
    function invertMatrix(A) {
        const n = A.length;

        // Create augmented matrix [A | I]
        const aug = [];
        for (let i = 0; i < n; i++) {
            aug[i] = [...A[i]];
            for (let j = 0; j < n; j++) {
                aug[i].push(i === j ? 1 : 0);
            }
        }

        // Forward elimination with partial pivoting
        for (let col = 0; col < n; col++) {
            // Find pivot
            let maxRow = col;
            for (let row = col + 1; row < n; row++) {
                if (abs(aug[row][col]) > abs(aug[maxRow][col])) {
                    maxRow = row;
                }
            }
            [aug[col], aug[maxRow]] = [aug[maxRow], aug[col]];

            if (abs(aug[col][col]) < EPSILON) {
                throw new Error('Singular matrix');
            }

            // Scale pivot row
            const pivot = aug[col][col];
            for (let j = col; j < 2 * n; j++) {
                aug[col][j] /= pivot;
            }

            // Eliminate column
            for (let row = 0; row < n; row++) {
                if (row !== col) {
                    const factor = aug[row][col];
                    for (let j = col; j < 2 * n; j++) {
                        aug[row][j] -= factor * aug[col][j];
                    }
                }
            }
        }

        // Extract inverse from augmented matrix
        const inv = [];
        for (let i = 0; i < n; i++) {
            inv[i] = aug[i].slice(n);
        }
        return inv;
    }

    // ============================================================================
    // EXPORTS
    // ============================================================================

    const Utils = {
        // Constants
        EPSILON,
        PI,

        // Scalar utilities
        isZero,
        clamp,
        abs,
        sqrt,
        sin,
        cos,
        tan,
        atan2,
        acos,
        min,
        max,
        exp,
        log,

        // General vector operations
        dot,
        norm,
        normSq,
        normalize,
        vecAdd,
        vecSub,
        vecScale,
        vecCombine,

        // 3D vector operations (optimized)
        cross3,
        dot3,
        norm3,
        normalize3,
        vecAdd3,
        vecSub3,
        vecScale3,

        // Matrix operations
        matVecMul,
        matMul,
        transpose,
        identity,

        // Determinants
        det2x2,
        det3x3,

        // Matrix inversion
        inv2x2,
        inv3x3,
        invertMatrix
    };

    // Export for different module systems
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = Utils;
    }
    if (typeof define === 'function' && define.amd) {
        define([], () => Utils);
    }
    if (typeof global !== 'undefined') {
        global.ContactThermoUtils = Utils;
    }

})(typeof globalThis !== 'undefined' ? globalThis : (typeof window !== 'undefined' ? window : global));
