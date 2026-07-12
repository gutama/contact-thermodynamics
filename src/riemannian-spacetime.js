/**
 * Backwards-compatible entry point for Spacetime Geometric Algebra.
 *
 * The implementation lives in `physics/spacetime.js`; this file preserves the
 * historical `require('./riemannian-spacetime.js')` path without duplicating the
 * full module body.
 *
 * @module riemannian-spacetime
 * @license MIT
 */

(function (global) {
    'use strict';

    let RiemannianSpacetimeModule;
    if (typeof require !== 'undefined') {
        RiemannianSpacetimeModule = require('./physics/spacetime.js');
    } else if (typeof global !== 'undefined') {
        RiemannianSpacetimeModule = global.RiemannianSpacetime;
    }

    if (typeof module !== 'undefined' && module.exports) {
        module.exports = RiemannianSpacetimeModule;
    }
    if (typeof define === 'function' && define.amd) {
        define([], () => RiemannianSpacetimeModule);
    }
    if (typeof global !== 'undefined') {
        global.RiemannianSpacetime = RiemannianSpacetimeModule;
    }

})(typeof globalThis !== 'undefined' ? globalThis : (typeof window !== 'undefined' ? window : global));
