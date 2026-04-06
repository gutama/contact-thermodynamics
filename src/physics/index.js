/**
 * Physics Module
 * 
 * Physics applications: pilot-wave theory, entropic gravity, spacetime.
 * 
 * @module physics
 */

const PilotWave = require('./pilot-wave');
const EntropicGravity = require('./entropic-gravity');
const Spacetime = require('./spacetime');
const G2BlackHoleRemnant = require('./g2-black-hole-remnant');

module.exports = {
    ...PilotWave,
    ...EntropicGravity,
    ...Spacetime,
    ...G2BlackHoleRemnant
};
