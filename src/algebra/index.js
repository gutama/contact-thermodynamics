/**
 * Algebra Module
 * 
 * Core Geometric Algebra components: Multivector, Algebra, Number Systems.
 * 
 * @module algebra
 */

const { Algebra, Multivector } = require('./multivector');
const NumberSystems = require('./number-systems');

module.exports = {
    Algebra,
    Multivector,
    ...NumberSystems
};
