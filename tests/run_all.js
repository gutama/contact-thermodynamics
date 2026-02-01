/**
 * Master Test Runner
 * 
 * Runs all test suites: geometry, mesh, pde, contact.
 */

const { execSync } = require('child_process');
const path = require('path');

console.log('========================================');
console.log('  CONTACT-THERMODYNAMICS TEST SUITE');
console.log('========================================\n');

const testSuites = [
    { name: 'Geometry', runner: 'tests/run_geometry.js' },
    { name: 'Mesh', runner: 'tests/run_mesh.js' },
    { name: 'PDE', runner: 'tests/run_pde.js' },
    { name: 'Contact', runner: 'tests/run_contact.js' },
    { name: 'Pilot-Wave', tests: ['tests/test_pilot_wave.js', 'tests/test_known_solutions.js'] }
];

let totalPassed = 0;
let totalFailed = 0;

for (const suite of testSuites) {
    console.log(`\n${'='.repeat(40)}`);
    console.log(`  ${suite.name.toUpperCase()} TESTS`);
    console.log('='.repeat(40) + '\n');

    try {
        if (suite.runner) {
            execSync(`node "${suite.runner}"`, {
                stdio: 'inherit',
                cwd: __dirname.replace(/tests$/, '')
            });
        } else if (suite.tests) {
            for (const test of suite.tests) {
                console.log(`--- Running: ${path.basename(test)} ---\n`);
                execSync(`node "${test}"`, {
                    stdio: 'inherit',
                    cwd: __dirname.replace(/tests$/, '')
                });
            }
        }
        totalPassed++;
    } catch (err) {
        totalFailed++;
    }
}

console.log('\n' + '='.repeat(40));
console.log(`  FINAL: ${totalPassed} suites passed, ${totalFailed} failed`);
console.log('='.repeat(40));

process.exit(totalFailed > 0 ? 1 : 0);
