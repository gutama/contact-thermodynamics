/**
 * Backwards-compatible entry point for Pilot-Wave theory.
 *
 * The implementation lives in `physics/pilot-wave.js`; this file preserves the
 * historical `require('./pilot-wave.js')` path without duplicating the full
 * module body.
 *
 * @module pilot-wave
 * @license MIT
 */

(function (root, factory) {
    'use strict';

    if (typeof module === 'object' && module.exports) {
        // CommonJS / Node: load the implementation synchronously.
        module.exports = factory(require('./physics/pilot-wave.js'));
    } else if (typeof define === 'function' && define.amd) {
        // AMD (e.g. RequireJS): declare the dependency so the loader resolves it
        // asynchronously — never call a synchronous require() here.
        define(['./physics/pilot-wave.js'], factory);
    } else {
        // Browser global.
        root.PilotWave = factory(root.PilotWave);
    }

})(typeof globalThis !== 'undefined' ? globalThis : (typeof self !== 'undefined' ? self : this), function (PilotWave) {
    'use strict';
    return PilotWave;
});
