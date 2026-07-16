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

(function (root, factory) {
    'use strict';

    if (typeof module === 'object' && module.exports) {
        // CommonJS / Node: load the implementation synchronously.
        module.exports = factory(require('./algebra/multivector.js'));
    } else if (typeof define === 'function' && define.amd) {
        // AMD (e.g. RequireJS): declare the dependency so the loader resolves it
        // asynchronously — never call a synchronous require() here.
        define(['./algebra/multivector.js'], factory);
    } else {
        // Browser global.
        root.GA = factory(root.GA);
    }

})(typeof globalThis !== 'undefined' ? globalThis : (typeof self !== 'undefined' ? self : this), function (GA) {
    'use strict';
    return GA;
});
