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

(function (global) {
    'use strict';

    let MeshFTGCModule;
    if (typeof require !== 'undefined') {
        MeshFTGCModule = require('./calculus/mesh-derivative.js');
    } else if (typeof global !== 'undefined') {
        MeshFTGCModule = global.ContactThermo;
    }

    if (typeof module !== 'undefined' && module.exports) {
        module.exports = MeshFTGCModule;
    }
    if (typeof define === 'function' && define.amd) {
        define([], () => MeshFTGCModule);
    }
    if (typeof global !== 'undefined' && MeshFTGCModule) {
        global.ContactThermo = Object.assign(global.ContactThermo || {}, MeshFTGCModule);
    }

})(typeof globalThis !== 'undefined' ? globalThis : (typeof window !== 'undefined' ? window : this));
