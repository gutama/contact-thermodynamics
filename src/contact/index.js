/**
 * Contact Module
 * 
 * Contact geometry: manifolds, Hamiltonians, Legendrian submanifolds.
 * 
 * @module contact
 */

const Manifold = require('./manifold');
const Hamiltonian = require('./hamiltonian');
const Legendrian = require('./legendrian');

module.exports = {
    ...Manifold,
    ...Hamiltonian,
    ...Legendrian
};
