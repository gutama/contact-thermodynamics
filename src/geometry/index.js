/**
 * Geometry Module
 * 
 * Riemannian geometry via Geometric Algebra: connection bivectors, curvature,
 * geodesics, parallel transport.
 * 
 * @module geometry
 */

const RiemannianGA = require('./riemannian-ga');
const RiemannianDiscrete = require('./riemannian-discrete');
const Geodesic = require('./geodesic');

module.exports = {
    ...RiemannianGA,
    ...RiemannianDiscrete,
    ...Geodesic
};
