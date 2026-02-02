/**
 * Calculus Module
 * 
 * Discrete Geometric Calculus: split differential operator, grid/mesh derivatives,
 * and PDE solvers (wave, heat, Maxwell).
 * 
 * @module calculus
 */

const Grid = require('./grid');
const Mesh = require('./mesh');
const MeshDerivative = require('./mesh-derivative');
const Solvers = require('./solvers');

module.exports = {
    ...Grid,
    ...Mesh,
    ...MeshDerivative,
    ...Solvers
};
