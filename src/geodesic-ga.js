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

(function (global) {
    'use strict';

    let GeodesicGA;
    if (typeof require !== 'undefined') {
        GeodesicGA = require('./geometry/geodesic.js');
    } else if (typeof global !== 'undefined') {
        GeodesicGA = global.GeodesicGA;
    }

    if (typeof module !== 'undefined' && module.exports) {
        module.exports = GeodesicGA;
    }
    if (typeof define === 'function' && define.amd) {
        define([], () => GeodesicGA);
    }
    if (typeof global !== 'undefined') {
        global.GeodesicGA = GeodesicGA;
    }

})(typeof globalThis !== 'undefined' ? globalThis : (typeof window !== 'undefined' ? window : this));
