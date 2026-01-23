# Versor Extensions Reference

This document describes the Versor-style operations added to extend the geometric-algebra skill.

## Quick Start

```javascript
// Load extensions (after ga_vanilla.js)
const cga = new VersorGA.CGA3DVersor();

// Create a booster from point pair
const pp = cga.pointPair([1, 0, 0], [-1, 0, 0]);
const boost = cga.Gen.bst(pp.scale(0.5));

// Apply to a point
const p = cga.point(0, 2, 0);
const boosted = cga.spin(p, boost);
```

---

## 1. Booster Operations

### Gen.bst(pointPair) → Booster
Generate a conformal boost from a point pair.

```javascript
const pp = cga.pointPair([0, -1, 0], [0, 1, 0]);
const boost = cga.Gen.bst(pp.scale(amount));

// Boost type depends on point pair:
// - Real PP (positive sq) → Hyperbolic boost
// - Imaginary PP (negative sq) → Elliptic boost
// - Null PP → Translation-like
```

### Gen.bstPoints(p1, p2, amount) → Booster
Convenience method to create booster from two points.

```javascript
const boost = cga.Gen.bstPoints([1, 0, 0], [-1, 0, 0], 0.5);
```

### Gen.logBoost(booster) → PointPair
Extract the point pair generator from a booster.

```javascript
const pp = cga.Gen.logBoost(boost);
```

### Boosted Surface Generation

```javascript
// Apply booster field to mesh vertices
const vertices = [[0,0,0], [1,0,0], [0,1,0], ...];
const transformed = cga.boostSurface(vertices, boost);

// Create field from multiple generators
const field = cga.boosterField([
    { p1: [-1, 0, 0], p2: [1, 0, 0], weight: 0.3 },
    { p1: [0, -1, 0], p2: [0, 1, 0], weight: -0.2 }
]);
const localBoost = field(x, y, z);
```

---

## 2. Tangent Elements

### cga.tangentVector(px, py, pz, vx, vy, vz) → Tnv
Create a tangent vector at point (px,py,pz) with direction (vx,vy,vz).

```javascript
const tnv = cga.tangentVector(0, 0, 0, 1, 0, 0);
// Represents: Tnv = Point ∧ (Direction ∧ n∞)
```

### cga.tangentVectorAt(point, direction) → Tnv
Create tangent from existing multivectors.

```javascript
const p = cga.point(1, 2, 3);
const dir = cga.e1;  // x direction
const tnv = cga.tangentVectorAt(p, dir);
```

### cga.tangentBivector(px, py, pz, bxy, bxz, byz) → Tnb
Tangent bivector (area element at a point).

```javascript
const tnb = cga.tangentBivector(0, 0, 0, 1, 0, 0);  // xy-plane at origin
```

### cga.tangentTrivector(px, py, pz, scale) → Tnt
Tangent trivector (volume element at a point).

---

## 3. Direction Elements

### cga.directionVector(vx, vy, vz) → Drv
Free vector without position.

```javascript
const drv = cga.directionVector(1, 0, 0);
// Drv = v ∧ n∞
// Direction vectors transform correctly under rotations
// but are unaffected by translations
```

### cga.directionBivector(bxy, bxz, byz) → Drb
Free area element.

### cga.directionTrivector(scale) → Drt
Free volume element.

---

## 4. Extended Round Operations

Access via `cga.Round`:

### Round.null(x, y, z) → Point
Create null point (same as cga.point).

### Round.dls(cx, cy, cz, r) → DualSphere
Dual sphere from center and radius.

```javascript
const sphere = cga.Round.dls(0, 0, 0, 2);  // radius 2 at origin
const imaginary = cga.Round.dls(0, 0, 0, -1);  // imaginary (r² < 0)
```

### Round.dlsFromPoints(center, surface) → DualSphere
Dual sphere from center point and a point on surface.

### Round.carrier(round) → Flat
Get the flat that contains the round.
- carrier(Circle) → Plane
- carrier(PointPair) → Line

```javascript
const circle = cga.circle3([1,0,0], [0,1,0], [-1,0,0]);
const plane = cga.Round.carrier(circle);
```

### Round.surround(round) → Sphere
Get the sphere that surrounds the round.

### Round.location(round) → Point
Center/location of round element (normalized).

### Round.center(round) → Point
Simple center (not normalized).

### Round.direction(round) → Vector
Direction of round element.
- For circle: normal to plane
- For point pair: line direction

### Round.size(round, signed) → number
Radius of round element.

```javascript
const r = cga.Round.size(sphere);
const signedR = cga.Round.size(sphere, true);  // negative for imaginary
```

### Round.isReal(round) → boolean
### Round.isImaginary(round) → boolean

```javascript
if (cga.Round.isImaginary(sphere)) {
    console.log("Imaginary sphere (no real intersection)");
}
```

### Round.split(pointPair) → [Point, Point]
Extract two points from point pair.

```javascript
const pp = cga.pointPair([1, 0, 0], [-1, 0, 0]);
const [p1, p2] = cga.Round.split(pp);
```

### Round.pointOnCircle(circle, t) → Point
Get point on circle at parameter t ∈ [0, 2π].

### Round.weight(round) → number
Weight/magnitude of round element.

---

## 5. Extended Flat Operations

Access via `cga.Flat`:

### Flat.point(x, y, z) → FlatPoint
Flat point (point on a flat).

### Flat.line(p1, p2) → Line
Line through two points.

### Flat.lineDir(px, py, pz, dx, dy, dz) → Line
Line from point and direction.

### Flat.plane(p1, p2, p3) → Plane
Plane through three points.

### Flat.planeNormal(nx, ny, nz, d) → DualPlane
Plane from normal and distance.

### Flat.carrier(flat) → Flat
Normalized flat.

### Flat.direction(line) → Vector
Direction vector of line.

### Flat.location(flat) → Point
Point on flat closest to origin.

### Flat.normal(plane) → [nx, ny, nz]
Normal vector of plane.

### Flat.weight(flat) → number
Weight/magnitude.

---

## 6. Extended Generators

### Gen.rot(bivector) → Rotor
Rotor from bivector (rotation in plane).

### Gen.trs(dx, dy, dz) → Translator
Translator from displacement.

### Gen.trsVec(direction) → Translator
Translator from direction multivector.

### Gen.trsPoint(point) → Translator
Translator to move origin to point.

### Gen.mot(dualLine) → Motor
Motor from dual line (screw motion).

### Gen.dil(point, scale) → Dilator
Dilator centered at point.

### Gen.dilOrigin(scale) → Dilator
Dilator at origin.

### Gen.trv(vx, vy, vz) → Transversor
Transversor from tangent direction.

### Gen.trvVec(tangentVec) → Transversor
Transversor from tangent multivector.

### Gen.logRotor(rotor) → Bivector
Logarithm of rotor.

### Gen.logMotor(motor) → DualLine
Logarithm of motor.

### Gen.ratio(v1, v2) → Versor
Ratio between two versors: v2 * v1⁻¹

### Gen.slerp(v1, v2, t) → Versor
Interpolate versors (0 ≤ t ≤ 1).

---

## 7. Transformation Application

### cga.spin(element, versor) → Element
Sandwich product: V * X * Ṽ

```javascript
const rotated = cga.spin(point, rotor);
const twisted = cga.spin(line, motor);
const bent = cga.spin(circle, booster);
```

### cga.reflect(element, reflector) → Element
Reflection: R * X * R

```javascript
const reflected = cga.reflect(point, plane);
```

### cga.invert(element, sphere) → Element
Circle inversion in sphere.

### cga.meet(a, b) → Element
Intersection (regressive product).

### cga.join(a, b) → Element
Join (outer product).

---

## 8. Arbitrary Algebra Creation

### VersorGA.create(config) → Algebra

```javascript
// Euclidean 4D
const E4 = VersorGA.create({ metric: [1, 1, 1, 1] });

// Conformal 2D (equivalent to CGA2D)
const C2 = VersorGA.create({
    metric: [1, 1, 1, -1],
    conformal: true
});

// Spacetime algebra
const STA = VersorGA.create({ metric: [1, -1, -1, -1] });

// With named types
const custom = VersorGA.create({
    metric: [1, 1, 1],
    types: [
        { name: 'Rot', bases: ['1', 'e12', 'e13', 'e23'] },
        { name: 'Vec', bases: ['e1', 'e2', 'e3'] }
    ]
});

// Use named type
const r = custom.Rot(1, 0, 0, 0);  // Identity rotor
```

---

## Common Patterns

### Surface Deformation with Boosters

```javascript
// Define curvature generators
const generators = [
    { p1: [-1, 0, 0], p2: [1, 0, 0], weight: 0.3 },
    { p1: [0, -1, 0], p2: [0, 1, 0], weight: -0.2 }
];

// Create deformation field
const field = cga.boosterField(generators);

// Deform mesh
const deformed = mesh.vertices.map(v => {
    const boost = field(v[0], v[1], v[2]);
    const p = cga.point(...v);
    const q = cga.spin(p, boost);
    return cga.pointCoords(q);
});
```

### Smooth Interpolation

```javascript
// Interpolate between two rotations
const R1 = cga.rotorAxis(0, 0, 1, 0);
const R2 = cga.rotorAxis(0, 0, 1, Math.PI);

for (let t = 0; t <= 1; t += 0.1) {
    const Rt = cga.Gen.slerp(R1, R2, t);
    const p = cga.spin(point, Rt);
    // Draw p...
}
```

### Handling Imaginary Rounds

```javascript
const sphere = cga.Round.dls(0, 0, 0, radius);

if (cga.Round.isImaginary(sphere)) {
    // This is an imaginary sphere
    // Represents no real intersection
    console.log("No real points on this sphere");
} else {
    const r = cga.Round.size(sphere);
    console.log("Real sphere with radius", r);
}
```

---

## Bivector Classification (Review)

From `ga_number_systems.js`:

| B² | Type | Transformation |
|----|------|----------------|
| -1 | Elliptic | Rotation (circular) |
| 0 | Parabolic | Translation (linear) |
| +1 | Hyperbolic | Boost (hyperbolic) |

Point pairs in CGA follow the same pattern:
- Real PP → Hyperbolic boost
- Imaginary PP → Elliptic boost
- Null PP → Parabolic (translation-like)
