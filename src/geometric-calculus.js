/**
 * Backwards-compatible entry point for Geometric Calculus on regular grids.
 *
 * The implementation lives in `calculus/grid.js`; this file preserves the
 * historical `require('./geometric-calculus.js')` path without duplicating the
 * full module body.
 *
 * @module geometric-calculus
 * @license MIT
 */

(function (global) {
    'use strict';

    let GeometricCalculus;
    if (typeof require !== 'undefined') {
        GeometricCalculus = require('./calculus/grid.js');
    } else if (typeof global !== 'undefined') {
        GeometricCalculus = global.GeometricCalculus;
    }

    if (typeof module !== 'undefined' && module.exports) {
        module.exports = GeometricCalculus;
    }
    if (typeof define === 'function' && define.amd) {
        define([], () => GeometricCalculus);
    }
    if (typeof global !== 'undefined') {
        global.GeometricCalculus = GeometricCalculus;
    }

})(typeof globalThis !== 'undefined' ? globalThis : (typeof window !== 'undefined' ? window : this));
