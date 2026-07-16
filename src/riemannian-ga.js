/**
 * Backwards-compatible entry point for Riemannian Geometry via Geometric Algebra.
 *
 * The implementation lives in `geometry/riemannian-ga.js`; this file preserves
 * the historical `require('./riemannian-ga.js')` path without duplicating the
 * full module body.
 *
 * @module riemannian-ga
 * @license MIT
 */

(function (global) {
    'use strict';

    let RiemannianGA;
    if (typeof require !== 'undefined') {
        RiemannianGA = require('./geometry/riemannian-ga.js');
    } else if (typeof global !== 'undefined') {
        RiemannianGA = global.RiemannianGA;
    }

    if (typeof module !== 'undefined' && module.exports) {
        module.exports = RiemannianGA;
    }
    if (typeof define === 'function' && define.amd) {
        define([], () => RiemannianGA);
    }
    if (typeof global !== 'undefined') {
        global.RiemannianGA = RiemannianGA;
    }

})(typeof globalThis !== 'undefined' ? globalThis : (typeof window !== 'undefined' ? window : this));
