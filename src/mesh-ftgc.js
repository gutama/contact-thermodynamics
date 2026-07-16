/**
 * Backwards-compatible entry point for the mesh FTGC geometric derivative.
 *
 * The implementation lives in `calculus/mesh-derivative.js`; this file preserves
 * the historical `require('./mesh-ftgc.js')` path without duplicating the full
 * module body.
 *
 * @module mesh-ftgc
 * @license MIT
 */

(function (root, factory) {
    'use strict';

    if (typeof module === 'object' && module.exports) {
        // CommonJS / Node: load the implementation synchronously.
        module.exports = factory(require('./calculus/mesh-derivative.js'));
    } else if (typeof define === 'function' && define.amd) {
        // AMD (e.g. RequireJS): declare the dependency so the loader resolves it
        // asynchronously — never call a synchronous require() here.
        define(['./calculus/mesh-derivative.js'], factory);
    } else {
        // Browser global: merge the mesh FTGC exports into the shared namespace.
        const MeshFTGCModule = factory(root.ContactThermo);
        if (MeshFTGCModule) {
            root.ContactThermo = Object.assign(root.ContactThermo || {}, MeshFTGCModule);
        }
    }

})(typeof globalThis !== 'undefined' ? globalThis : (typeof self !== 'undefined' ? self : this), function (MeshFTGCModule) {
    'use strict';
    return MeshFTGCModule;
});
