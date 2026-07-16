/**
 * Backwards-compatible entry point for the Entropic Gravity module.
 *
 * The implementation lives in `physics/entropic-gravity.js`; this file preserves
 * the historical `require('./entropic-gravity.js')` path without duplicating the
 * full module body.
 *
 * @module entropic-gravity
 * @license MIT
 */

(function (global) {
    'use strict';

    let EntropicGravity;
    if (typeof require !== 'undefined') {
        EntropicGravity = require('./physics/entropic-gravity.js');
    } else if (typeof global !== 'undefined') {
        EntropicGravity = global.EntropicGravity;
    }

    if (typeof module !== 'undefined' && module.exports) {
        module.exports = EntropicGravity;
    }
    if (typeof define === 'function' && define.amd) {
        define([], () => EntropicGravity);
    }
    if (typeof global !== 'undefined') {
        global.EntropicGravity = EntropicGravity;
    }

})(typeof globalThis !== 'undefined' ? globalThis : (typeof window !== 'undefined' ? window : this));
