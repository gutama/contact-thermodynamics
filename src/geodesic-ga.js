/**
 * Backwards-compatible entry point for the coordinate-free geodesic solver.
 *
 * The implementation lives in `geometry/geodesic.js`; this file preserves the
 * historical `require('./geodesic-ga.js')` path without duplicating the full
 * module body.
 *
 * @module geodesic-ga
 * @license MIT
 */

(function (root, factory) {
    'use strict';

    if (typeof module === 'object' && module.exports) {
        // CommonJS / Node: load the implementation synchronously.
        module.exports = factory(require('./geometry/geodesic.js'));
    } else if (typeof define === 'function' && define.amd) {
        // AMD (e.g. RequireJS): declare the dependency so the loader resolves it
        // asynchronously — never call a synchronous require() here.
        define(['./geometry/geodesic.js'], factory);
    } else {
        // Browser global.
        root.GeodesicGA = factory(root.GeodesicGA);
    }

})(typeof globalThis !== 'undefined' ? globalThis : (typeof self !== 'undefined' ? self : this), function (GeodesicGA) {
    'use strict';
    return GeodesicGA;
});
