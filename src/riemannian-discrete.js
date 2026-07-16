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

(function (root, factory) {
    'use strict';

    if (typeof module === 'object' && module.exports) {
        // CommonJS / Node: load the implementation synchronously.
        module.exports = factory(require('./geometry/riemannian-discrete.js'));
    } else if (typeof define === 'function' && define.amd) {
        // AMD (e.g. RequireJS): declare the dependency so the loader resolves it
        // asynchronously — never call a synchronous require() here.
        define(['./geometry/riemannian-discrete.js'], factory);
    } else {
        // Browser global.
        root.RiemannianDiscrete = factory(root.RiemannianDiscrete);
    }

})(typeof globalThis !== 'undefined' ? globalThis : (typeof self !== 'undefined' ? self : this), function (RiemannianDiscrete) {
    'use strict';
    return RiemannianDiscrete;
});
