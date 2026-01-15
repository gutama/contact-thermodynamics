# Contributing to Contact Thermodynamics

Thank you for your interest in contributing! This document provides guidelines for contributing to the project.

## Code of Conduct

Please be respectful and constructive in all interactions. We welcome contributors of all backgrounds and experience levels.

## How to Contribute

### Reporting Issues

1. Check if the issue already exists
2. Use a clear, descriptive title
3. Include:
   - Steps to reproduce
   - Expected behavior
   - Actual behavior
   - Node.js version
   - Any relevant code snippets

### Submitting Pull Requests

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/your-feature`
3. Make your changes
4. Add/update tests as needed
5. Run the test suite: `npm test`
6. Commit with clear messages
7. Push to your fork
8. Open a Pull Request

### Development Setup

```bash
git clone https://github.com/gutama/contact-thermodynamics.git
cd contact-thermodynamics
npm install
npm test
```

### Running Examples

```bash
npm run example:basic
npm run example:dynamics
npm run example:geodesics
```

## Code Style

- Use ES6+ features
- 4-space indentation
- Meaningful variable names
- Comment complex mathematical operations
- Include JSDoc comments for public methods

### Example

```javascript
/**
 * Compute the contact Hamiltonian vector field X_H
 * 
 * For canonical coordinates, this gives:
 *   ẋ^a = ∂H/∂p_a
 *   ṗ_a = -∂H/∂x^a - p_a · RH
 *   u̇ = p_a · ∂H/∂p_a - H
 * 
 * @param {ContactPoint} pt - Point on the manifold
 * @returns {Object} Vector field components
 */
vectorField(pt) {
    // Implementation...
}
```

## Testing

- All new features should have tests
- Tests go in `tests/test.js`
- Use the existing test utilities:

```javascript
function assert(condition, message) { ... }
function assertApprox(a, b, tol, message) { ... }
function section(name) { ... }
```

## Documentation

- Update `docs/API.md` for new public APIs
- Add examples for new features
- Update tutorials if the feature is significant

## Mathematical Notation

When documenting mathematical concepts:

- Use LaTeX in markdown: `$\alpha \wedge (d\alpha)^n \neq 0$`
- Provide both mathematical and code representations
- Include physical interpretations where applicable

## Areas for Contribution

### High Priority

- [ ] TypeScript type definitions (`src/index.d.ts`)
- [ ] Additional spacetime metrics (Kerr, Reissner-Nordström)
- [ ] Symplectic integrators for long-time evolution
- [ ] Browser bundle (webpack/rollup)

### Medium Priority

- [ ] Additional thermodynamic equations of state
- [ ] Contact transformations / Legendrian isotopies
- [ ] Visualization improvements
- [ ] Performance optimizations

### Research Directions

- [ ] Connection to information geometry
- [ ] Non-equilibrium thermodynamics
- [ ] Quantum contact geometry
- [ ] Numerical methods for Legendrian fronts

## Questions?

Open an issue with the "question" label or reach out to the maintainers.

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
