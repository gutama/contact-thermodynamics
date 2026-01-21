/**
 * Test Runner for PDE Tests
 */

const { execSync } = require('child_process');
const path = require('path');
const fs = require('fs');

const testDir = path.join(__dirname, 'pde');

console.log('========================================');
console.log('  PDE TEST SUITE');
console.log('========================================\n');

let passed = 0;
let failed = 0;

const testFiles = fs.readdirSync(testDir).filter(f => f.endsWith('.js'));

for (const file of testFiles) {
    const testPath = path.join(testDir, file);
    console.log(`\n--- Running: ${file} ---\n`);

    try {
        execSync(`node "${testPath}"`, {
            stdio: 'inherit',
            cwd: path.join(__dirname, '..')
        });
        passed++;
    } catch (err) {
        failed++;
    }
}

console.log('\n========================================');
console.log(`  RESULTS: ${passed} passed, ${failed} failed`);
console.log('========================================');

process.exit(failed > 0 ? 1 : 0);
