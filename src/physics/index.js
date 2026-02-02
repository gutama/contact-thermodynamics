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

module.exports = {
    ...PilotWave,
    ...EntropicGravity,
    ...Spacetime
};
