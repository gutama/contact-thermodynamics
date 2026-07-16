/**
 * Backwards-compatible entry point for mesh PDE solvers.
 *
 * The implementation lives in `calculus/solvers.js`; this file preserves the
 * historical `require('./mesh-solvers.js')` path without duplicating the full
 * module body.
 *
 * @module mesh-solvers
 * @license MIT
 */

(function (global) {
    'use strict';

    let MeshSolversModule;
    if (typeof require !== 'undefined') {
        MeshSolversModule = require('./calculus/solvers.js');
    } else if (typeof global !== 'undefined') {
        MeshSolversModule = global.ContactThermo;
    }

    if (typeof module !== 'undefined' && module.exports) {
        module.exports = MeshSolversModule;
    }
    if (typeof define === 'function' && define.amd) {
        define([], () => MeshSolversModule);
    }
    if (typeof global !== 'undefined' && MeshSolversModule) {
        global.ContactThermo = Object.assign(global.ContactThermo || {}, MeshSolversModule);
    }

})(typeof globalThis !== 'undefined' ? globalThis : (typeof window !== 'undefined' ? window : this));
