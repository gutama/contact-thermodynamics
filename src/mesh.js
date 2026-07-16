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

(function (root, factory) {
    'use strict';

    if (typeof module === 'object' && module.exports) {
        // CommonJS / Node: load the implementation synchronously.
        module.exports = factory(require('./calculus/mesh.js'));
    } else if (typeof define === 'function' && define.amd) {
        // AMD (e.g. RequireJS): declare the dependency so the loader resolves it
        // asynchronously — never call a synchronous require() here.
        define(['./calculus/mesh.js'], factory);
    } else {
        // Browser global: merge the mesh exports into the shared namespace.
        const MeshModule = factory(root.ContactThermo);
        if (MeshModule) {
            root.ContactThermo = Object.assign(root.ContactThermo || {}, MeshModule);
        }
    }

})(typeof globalThis !== 'undefined' ? globalThis : (typeof self !== 'undefined' ? self : this), function (MeshModule) {
    'use strict';
    return MeshModule;
});
