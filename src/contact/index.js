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
const GAContactForm = require('./ga-contact-form');
const GAContactVersor = require('./ga-contact-versor');

module.exports = {
    ...Manifold,
    ...Hamiltonian,
    ...Legendrian,
    ...GAContactForm,
    ...GAContactVersor
};
