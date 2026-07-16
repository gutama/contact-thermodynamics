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

(function (global) {
    'use strict';

    let PilotWave;
    if (typeof require !== 'undefined') {
        PilotWave = require('./physics/pilot-wave.js');
    } else if (typeof global !== 'undefined') {
        PilotWave = global.PilotWave;
    }

    if (typeof module !== 'undefined' && module.exports) {
        module.exports = PilotWave;
    }
    if (typeof define === 'function' && define.amd) {
        define([], () => PilotWave);
    }
    if (typeof global !== 'undefined') {
        global.PilotWave = PilotWave;
    }

})(typeof globalThis !== 'undefined' ? globalThis : (typeof window !== 'undefined' ? window : this));
