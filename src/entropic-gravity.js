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

(function (root, factory) {
    'use strict';

    if (typeof module === 'object' && module.exports) {
        // CommonJS / Node: load the implementation synchronously.
        module.exports = factory(require('./physics/entropic-gravity.js'));
    } else if (typeof define === 'function' && define.amd) {
        // AMD (e.g. RequireJS): declare the dependency so the loader resolves it
        // asynchronously — never call a synchronous require() here.
        define(['./physics/entropic-gravity.js'], factory);
    } else {
        // Browser global.
        root.EntropicGravity = factory(root.EntropicGravity);
    }

})(typeof globalThis !== 'undefined' ? globalThis : (typeof self !== 'undefined' ? self : this), function (EntropicGravity) {
    'use strict';
    return EntropicGravity;
});
