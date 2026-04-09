module.exports = {
    env: {
        browser: true,
        node: true,
        es2021: true
    },
    parserOptions: {
        ecmaVersion: 2021
    },
    rules: {
        'no-unused-vars': ['warn', { args: 'none', varsIgnorePattern: '^_' }],
        'no-undef': 'warn',
        'no-console': 'off',
        'eqeqeq': 'warn',
        'no-var': 'error',
        'prefer-const': 'warn'
    }
};
