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

(function (root, factory) {
    'use strict';

    if (typeof module === 'object' && module.exports) {
        // CommonJS / Node: load the implementation synchronously.
        module.exports = factory(require('./calculus/grid.js'));
    } else if (typeof define === 'function' && define.amd) {
        // AMD (e.g. RequireJS): declare the dependency so the loader resolves it
        // asynchronously — never call a synchronous require() here.
        define(['./calculus/grid.js'], factory);
    } else {
        // Browser global.
        root.GeometricCalculus = factory(root.GeometricCalculus);
    }

})(typeof globalThis !== 'undefined' ? globalThis : (typeof self !== 'undefined' ? self : this), function (GeometricCalculus) {
    'use strict';
    return GeometricCalculus;
});
