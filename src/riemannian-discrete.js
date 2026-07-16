/**
 * Backwards-compatible entry point for Discrete Riemannian Geometry on meshes.
 *
 * The implementation lives in `geometry/riemannian-discrete.js`; this file
 * preserves the historical `require('./riemannian-discrete.js')` path without
 * duplicating the full module body.
 *
 * @module riemannian-discrete
 * @license MIT
 */

(function (global) {
    'use strict';

    let RiemannianDiscrete;
    if (typeof require !== 'undefined') {
        RiemannianDiscrete = require('./geometry/riemannian-discrete.js');
    } else if (typeof global !== 'undefined') {
        RiemannianDiscrete = global.RiemannianDiscrete;
    }

    if (typeof module !== 'undefined' && module.exports) {
        module.exports = RiemannianDiscrete;
    }
    if (typeof define === 'function' && define.amd) {
        define([], () => RiemannianDiscrete);
    }
    if (typeof global !== 'undefined') {
        global.RiemannianDiscrete = RiemannianDiscrete;
    }

})(typeof globalThis !== 'undefined' ? globalThis : (typeof window !== 'undefined' ? window : this));
