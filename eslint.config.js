const globals = require('globals');

// ESLint 9+ flat config. Migrated from the legacy .eslintrc.js.
// The codebase is CommonJS/UMD (require + module.exports, with AMD shim
// branches guarded by `typeof define === 'function' && define.amd`), so we
// use sourceType 'commonjs' and declare the node/browser/amd globals that the
// old `env` block provided.
module.exports = [
    {
        files: ['**/*.js'],
        languageOptions: {
            ecmaVersion: 2021,
            sourceType: 'commonjs',
            globals: {
                ...globals.node,
                ...globals.browser,
                ...globals.amd
            }
        },
        rules: {
            'no-unused-vars': ['warn', { args: 'none', varsIgnorePattern: '^_' }],
            'no-undef': 'warn',
            'no-console': 'off',
            'eqeqeq': 'warn',
            'no-var': 'error',
            'prefer-const': 'warn'
        }
    }
];
