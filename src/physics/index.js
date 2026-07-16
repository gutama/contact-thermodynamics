/**
 * Physics Module
 * 
 * Physics applications: pilot-wave theory, entropic gravity, spacetime.
 * 
 * @module physics
 */

const PilotWave = require('./pilot-wave');
const GAPilotWave = require('./ga-pilot-wave');
const EntropicGravity = require('./entropic-gravity');
const Spacetime = require('./spacetime');
const G2BlackHoleRemnant = require('./g2-black-hole-remnant');

module.exports = {
    ...PilotWave,
    ...GAPilotWave,
    ...EntropicGravity,
    ...Spacetime,
    ...G2BlackHoleRemnant
};
