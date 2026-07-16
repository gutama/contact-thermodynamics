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

(function (root, factory) {
    'use strict';

    if (typeof module === 'object' && module.exports) {
        // CommonJS / Node: load the implementation synchronously.
        module.exports = factory(require('./calculus/solvers.js'));
    } else if (typeof define === 'function' && define.amd) {
        // AMD (e.g. RequireJS): declare the dependency so the loader resolves it
        // asynchronously — never call a synchronous require() here.
        define(['./calculus/solvers.js'], factory);
    } else {
        // Browser global: merge the solver exports into the shared namespace.
        const MeshSolversModule = factory(root.ContactThermo);
        if (MeshSolversModule) {
            root.ContactThermo = Object.assign(root.ContactThermo || {}, MeshSolversModule);
        }
    }

})(typeof globalThis !== 'undefined' ? globalThis : (typeof self !== 'undefined' ? self : this), function (MeshSolversModule) {
    'use strict';
    return MeshSolversModule;
});
