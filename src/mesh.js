/**
 * Backwards-compatible entry point for triangle-mesh geometry.
 *
 * The implementation lives in `calculus/mesh.js`; this file preserves the
 * historical `require('./mesh.js')` path without duplicating the full module
 * body.
 *
 * @module mesh
 * @license MIT
 */

(function (global) {
    'use strict';

    let MeshModule;
    if (typeof require !== 'undefined') {
        MeshModule = require('./calculus/mesh.js');
    } else if (typeof global !== 'undefined') {
        MeshModule = global.ContactThermo;
    }

    if (typeof module !== 'undefined' && module.exports) {
        module.exports = MeshModule;
    }
    if (typeof define === 'function' && define.amd) {
        define([], () => MeshModule);
    }
    if (typeof global !== 'undefined' && MeshModule) {
        global.ContactThermo = Object.assign(global.ContactThermo || {}, MeshModule);
    }

})(typeof globalThis !== 'undefined' ? globalThis : (typeof window !== 'undefined' ? window : this));
