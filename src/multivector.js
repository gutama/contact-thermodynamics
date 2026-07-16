/**
 * Backwards-compatible entry point for Geometric Algebra multivectors.
 *
 * The implementation lives in `algebra/multivector.js`; this file preserves the
 * historical `require('./multivector.js')` path without duplicating the full
 * module body.
 *
 * @module multivector
 * @license MIT
 */

(function (global) {
    'use strict';

    let GA;
    if (typeof require !== 'undefined') {
        GA = require('./algebra/multivector.js');
    } else if (typeof global !== 'undefined') {
        GA = global.GA;
    }

    if (typeof module !== 'undefined' && module.exports) {
        module.exports = GA;
    }
    if (typeof define === 'function' && define.amd) {
        define([], () => GA);
    }
    if (typeof global !== 'undefined') {
        global.GA = GA;
    }

})(typeof globalThis !== 'undefined' ? globalThis : (typeof window !== 'undefined' ? window : this));
