# Tutorial 6: Discrete Geometric Calculus on Meshes

This tutorial covers the **Fundamental Theorem of Geometric Calculus (FTGC)** 
applied to triangle meshes — enabling wave, heat, and Maxwell simulations.

## Prerequisites

Familiarity with:
- Basic mesh concepts (vertices, faces, edges)
- Differential operators (gradient, divergence, Laplacian)
- The wave and heat equations

## 1. Creating a Triangle Mesh

### Grid Mesh

```javascript
const CT = require('contact-thermodynamics');
const { TriangleMesh } = CT;

// Create a 16×16 grid on [0,1]×[0,1]
const mesh = TriangleMesh.createGrid(16, 16, 1.0, 1.0);

console.log(`Vertices: ${mesh.nVertices}`);   // 256
console.log(`Faces: ${mesh.nFaces}`);          // 450
console.log(`Edges: ${mesh.nEdges}`);          // 705
```

### Disk Mesh

```javascript
const disk = TriangleMesh.createDisk(8, 24, 1.0);
// Radial×Angular subdivision, radius 1.0
```

### Custom Mesh

```javascript
// Unit triangle
const verts = new Float64Array([
    0, 0, 0,   // vertex 0
    1, 0, 0,   // vertex 1
    0.5, 0.866, 0  // vertex 2
]);
const faces = new Uint32Array([0, 1, 2]);

const custom = new TriangleMesh(verts, faces);
```

## 2. Mesh Geometry

The mesh provides precomputed geometric quantities:

```javascript
// Face areas (one per triangle)
const areas = mesh.faceAreas();

// Edge lengths
const lengths = mesh.edgeLengths();

// Cotan weights (encode mesh metric)
const weights = mesh.cotanWeights();

// Dual areas per vertex (mixed Voronoi)
const dualAreas = mesh.vertexDualAreas();

// Total area
const totalArea = dualAreas.reduce((a, b) => a + b, 0);
```

### Boundary Detection

```javascript
// Boundary vertex mask (1 = on boundary)
const boundaryV = mesh.boundaryVertices;

// Boundary edge mask
const boundaryE = mesh.boundaryEdges;

// Count boundary vertices
let nBoundary = 0;
for (let i = 0; i < mesh.nVertices; i++) {
    if (mesh.boundaryVertices[i]) nBoundary++;
}
```

## 3. The Geometric Derivative ∇

The `MeshGeometricDerivative` class provides discrete differential operators:

```javascript
const { MeshGeometricDerivative } = CT;

const nabla = new MeshGeometricDerivative(mesh);
```

### Gradient (∇f)

Maps scalars (vertices) to vectors (edges):

```javascript
// Scalar field on vertices
const f = new Float64Array(mesh.nVertices);
for (let i = 0; i < mesh.nVertices; i++) {
    const [x, y, z] = mesh.getVertex(i);
    f[i] = Math.sin(Math.PI * x) * Math.sin(Math.PI * y);
}

// Gradient on edges
const gradF = nabla.grad(f);
console.log(`Gradient has ${gradF.length} values (one per edge)`);
```

### Divergence (∇·V)

Maps vectors (edges) to scalars (vertices):

```javascript
const divV = nabla.div(gradF);
// divV is the Laplacian of f!
```

### Curl (∇×V)

Maps vectors (edges) to bivectors (faces):

```javascript
const curlV = nabla.curl(gradF);
// curl(grad f) = 0 (up to numerical error)
```

### Laplacian (∇²f)

The cotan Laplacian on vertices:

```javascript
const lapF = nabla.laplacian(f);
```

## 4. Verifying FTGC Identities

The library provides verification methods:

```javascript
// curl(grad f) should be zero
const curlGradErr = nabla.verifyCurlGradZero(f);
console.log(`curl(grad f) error: ${curlGradErr.toExponential(2)}`);

// Laplacian of constant should be zero
const lapConstErr = nabla.verifyLaplacianConstantZero(1.0);
console.log(`lap(const) error: ${lapConstErr.toExponential(2)}`);
```

## 5. Heat Equation Simulation

The heat equation ∂u/∂t = α∇²u models diffusion:

```javascript
const { LeapfrogGCMesh, boundaryDirichletMask } = CT;

const solver = new LeapfrogGCMesh(mesh);

// Initial condition: Gaussian bump
const u0 = new Float64Array(mesh.nVertices);
for (let i = 0; i < mesh.nVertices; i++) {
    const [x, y] = mesh.getVertex(i);
    const dx = x - 0.5, dy = y - 0.5;
    u0[i] = Math.exp(-20 * (dx*dx + dy*dy));
}

// Dirichlet BC: fix boundary to zero
const dirichletMask = boundaryDirichletMask(mesh);
const dirichletValues = new Float64Array(mesh.nVertices); // all zeros

// Parameters
const alpha = 0.1;  // diffusion coefficient
const dt = solver.estimateCFLHeat(alpha) * 0.5;  // safe time step
const nSteps = 100;

// Simulate
const uFinal = solver.heatSimulate(u0, dt, nSteps, alpha, {
    implicit: true,  // unconditionally stable
    dirichletMask,
    dirichletValues,
    callback: (step, u) => {
        if (step % 20 === 0) {
            const variance = solver.variance(u);
            console.log(`Step ${step}: variance = ${variance.toExponential(2)}`);
        }
    }
});
```

## 6. Wave Equation Simulation

The wave equation ∂²u/∂t² = c²∇²u models oscillations:

```javascript
// Initial position
const u0 = new Float64Array(mesh.nVertices);
for (let i = 0; i < mesh.nVertices; i++) {
    const [x, y] = mesh.getVertex(i);
    const dx = x - 0.5, dy = y - 0.5;
    u0[i] = Math.exp(-30 * (dx*dx + dy*dy));
}

// Initial velocity (zero)
const v0 = new Float64Array(mesh.nVertices);

// Parameters
const c = 1.0;  // wave speed
const dtWave = solver.estimateCFL(c) * 0.5;
const nStepsWave = 200;

// Simulate
const uWave = solver.waveSimulate(u0, v0, dtWave, nStepsWave, c, {
    dirichletMask,  // fixed boundaries → reflection
    callback: (step, u) => {
        if (step % 50 === 0) {
            const energy = solver.energy(u);
            console.log(`Wave step ${step}: energy = ${energy.toExponential(3)}`);
        }
    }
});
```

## 7. Multivector Fields

For full geometric calculus, use `MeshMultivectorField`:

```javascript
const { MeshMultivectorField } = CT;

// Create fields of each grade
const scalar = MeshMultivectorField.scalarField(mesh, scalarValues);
const vector = MeshMultivectorField.vectorField(mesh, edgeValues);
const bivector = MeshMultivectorField.bivectorField(mesh, faceValues);

// Arithmetic
const sum = scalar.add(scalar);
const scaled = vector.scale(2.0);

// Grade projection
const grade0 = sum.gradeSelect(0);
```

## 8. Sparse Laplacian Matrix

For advanced use, access the Laplacian as a sparse matrix:

```javascript
const L = nabla.laplacianMatrix();

// Matrix-vector product
const Lu = L.matvec(u);

// Matrix properties
console.log(`Laplacian: ${L.nRows}×${L.nCols}, ${L.nnz} non-zeros`);
```

## Summary

| Component | Purpose |
|-----------|---------|
| `TriangleMesh` | Mesh data structure with typed arrays |
| `MeshGeometricDerivative` | Discrete ∇ (grad, div, curl, laplacian) |
| `MeshMultivectorField` | Staggered multivector storage |
| `LeapfrogGCMesh` | Time-stepping solver for PDEs |

## Next Steps

- Explore [mesh-heat-ftgc.js](../../examples/mesh-heat-ftgc.js) for entropy diffusion
- See [API Reference](../API.md) for full method documentation
- Read [THEORY.md](../THEORY.md) for mathematical background
