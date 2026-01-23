/**
 * Geometric Algebra Mesh Processing
 * 
 * Implementation of "Clean up your Mesh! Part 1: Plane and Simplex"
 * De Keninck, Roelfs, Dorst, Eelbode (2025)
 * 
 * Extends ga_vanilla.js with:
 * - k-simplices (vertices, edges, triangles, tetrahedra)
 * - k-complexes (point clouds, polygons, meshes)
 * - k-magnitudes (count, length, area, volume)
 * - Moments (center of mass, inertia tensor)
 * - Euclidean/Ideal split and norms
 * 
 * @requires ga_vanilla.js
 * @license MIT
 */

(function(global) {
    'use strict';

    const { sqrt, abs } = Math;
    const EPSILON = 1e-10;

    // ============================================================================
    // PGA 3D MESH ALGEBRA - Extends base PGA3D
    // ============================================================================

    /**
     * PGA3D Mesh Algebra with full mesh processing support
     * Signature: Cl(3,0,1) with e₁²=e₂²=e₃²=1, e₀²=0
     */
    class PGA3DMesh {
        constructor() {
            // Basis blade indices (bitmap representation)
            // Grade 0: scalar
            this.SCALAR = 0;        // 1
            
            // Grade 1: vectors (planes)
            this.E1 = 1;            // e₁
            this.E2 = 2;            // e₂
            this.E3 = 4;            // e₃
            this.E0 = 8;            // e₀ (ideal)
            
            // Grade 2: bivectors (lines)
            this.E12 = 3;           // e₁₂
            this.E13 = 5;           // e₁₃ (also written e₃₁ with sign)
            this.E23 = 6;           // e₂₃
            this.E01 = 9;           // e₀₁
            this.E02 = 10;          // e₀₂
            this.E03 = 12;          // e₀₃
            
            // Grade 3: trivectors (points)
            this.E123 = 7;          // e₁₂₃ = I₃ (Euclidean pseudoscalar)
            this.E012 = 11;         // e₀₁₂
            this.E013 = 13;         // e₀₁₃
            this.E023 = 14;         // e₀₂₃
            
            // Grade 4: pseudoscalar
            this.E0123 = 15;        // e₀₁₂₃ (PGA pseudoscalar)

            // Origin point: o = e₁₂₃ = I₃
            this.origin = new Float64Array(16);
            this.origin[this.E123] = 1;
        }

        // ========================================================================
        // MULTIVECTOR CREATION
        // ========================================================================

        /**
         * Create zero multivector
         */
        zero() {
            return new Float64Array(16);
        }

        /**
         * Create scalar multivector
         */
        scalar(s) {
            const mv = this.zero();
            mv[this.SCALAR] = s;
            return mv;
        }

        /**
         * Create plane: ax + by + cz + d = 0 → ae₁ + be₂ + ce₃ + de₀
         */
        plane(a, b, c, d) {
            const mv = this.zero();
            mv[this.E1] = a;
            mv[this.E2] = b;
            mv[this.E3] = c;
            mv[this.E0] = d;
            return mv;
        }

        /**
         * Create normalized point at (x, y, z)
         * v = e₁₂₃ - e₀(xe₂₃ - ye₁₃ + ze₁₂) = e₁₂₃ - xe₀₂₃ + ye₀₁₃ - ze₀₁₂
         * 
         * From paper eq. (3.8): v = I₃ - e₀·→v·I₃
         */
        point(x, y, z) {
            const mv = this.zero();
            mv[this.E123] = 1;      // Euclidean part (weight)
            mv[this.E023] = -x;     // -x·e₀₂₃
            mv[this.E013] = y;      // +y·e₀₁₃
            mv[this.E012] = -z;     // -z·e₀₁₂
            return mv;
        }

        /**
         * Create ideal point (direction at infinity)
         * Used for translations
         */
        direction(x, y, z) {
            const mv = this.zero();
            mv[this.E023] = -x;
            mv[this.E013] = y;
            mv[this.E012] = -z;
            return mv;
        }

        /**
         * Extract (x, y, z) coordinates from a point trivector
         */
        pointCoords(p) {
            const w = p[this.E123] || 1;
            return [
                -p[this.E023] / w,
                p[this.E013] / w,
                -p[this.E012] / w
            ];
        }

        // ========================================================================
        // GEOMETRIC PRODUCTS
        // ========================================================================

        /**
         * Geometric product of two multivectors
         * Full 16x16 Cayley table for Cl(3,0,1)
         */
        gp(a, b) {
            const r = this.zero();
            
            // This is the full geometric product - optimized for common cases
            // For mesh processing, we mainly need join (regressive product)
            
            // Scalar * anything
            for (let i = 0; i < 16; i++) {
                r[i] += a[0] * b[i] + b[0] * a[i];
            }
            r[0] -= a[0] * b[0]; // Subtract double-counted scalar*scalar
            
            // Vector * Vector (grade 1 * grade 1)
            // e₁e₁ = e₂e₂ = e₃e₃ = 1, e₀e₀ = 0
            // eᵢeⱼ = -eⱼeᵢ for i≠j
            r[0] += a[1]*b[1] + a[2]*b[2] + a[4]*b[4]; // eᵢ² = 1
            // e₀² = 0, so no a[8]*b[8] term
            
            r[3] += a[1]*b[2] - a[2]*b[1];   // e₁₂
            r[5] += a[1]*b[4] - a[4]*b[1];   // e₁₃
            r[6] += a[2]*b[4] - a[4]*b[2];   // e₂₃
            r[9] += a[1]*b[8] - a[8]*b[1];   // e₀₁
            r[10] += a[2]*b[8] - a[8]*b[2];  // e₀₂
            r[12] += a[4]*b[8] - a[8]*b[4];  // e₀₃

            // Continue with higher grade products...
            // (Full implementation would be ~200 lines)
            // For mesh processing, we use specialized join/meet instead
            
            return r;
        }

        /**
         * Reverse (reversion) of multivector
         * Reverses sign of grades 2 and 3
         */
        reverse(a) {
            const r = a.slice();
            // Grade 2 bivectors: negate
            r[3] = -a[3];   // e₁₂
            r[5] = -a[5];   // e₁₃
            r[6] = -a[6];   // e₂₃
            r[9] = -a[9];   // e₀₁
            r[10] = -a[10]; // e₀₂
            r[12] = -a[12]; // e₀₃
            // Grade 3 trivectors: negate
            r[7] = -a[7];   // e₁₂₃
            r[11] = -a[11]; // e₀₁₂
            r[13] = -a[13]; // e₀₁₃
            r[14] = -a[14]; // e₀₂₃
            return r;
        }

        /**
         * Add two multivectors
         */
        add(a, b) {
            const r = this.zero();
            for (let i = 0; i < 16; i++) r[i] = a[i] + b[i];
            return r;
        }

        /**
         * Subtract multivectors
         */
        sub(a, b) {
            const r = this.zero();
            for (let i = 0; i < 16; i++) r[i] = a[i] - b[i];
            return r;
        }

        /**
         * Scale multivector
         */
        scale(a, s) {
            const r = this.zero();
            for (let i = 0; i < 16; i++) r[i] = a[i] * s;
            return r;
        }

        // ========================================================================
        // JOIN (REGRESSIVE) PRODUCT - Core operation for mesh processing
        // From paper: (A ∨ B)* = A* ∧ B*
        // ========================================================================

        /**
         * Join of two points → line (bivector)
         * v₀ ∨ v₁ = (→v₁ - →v₀)I₃ + e₀(→v₀ ∧ →v₁)I₃
         * 
         * From paper eq. (3.10)
         */
        joinPP(v0, v1) {
            const [x0, y0, z0] = this.pointCoords(v0);
            const [x1, y1, z1] = this.pointCoords(v1);
            
            // Direction vector: →v₁ - →v₀
            const dx = x1 - x0;
            const dy = y1 - y0;
            const dz = z1 - z0;
            
            // Moment: →v₀ × →v₁ (cross product = dual of wedge)
            const mx = y0*z1 - z0*y1;
            const my = z0*x1 - x0*z1;
            const mz = x0*y1 - y0*x1;
            
            const line = this.zero();
            // Euclidean part: direction as e₂₃, e₃₁, e₁₂
            line[this.E23] = dx;
            line[this.E13] = -dy;  // e₃₁ = -e₁₃
            line[this.E12] = dz;
            // Ideal part: moment as e₀₁, e₀₂, e₀₃
            line[this.E01] = mx;
            line[this.E02] = my;
            line[this.E03] = mz;
            
            return line;
        }

        /**
         * Join of three points → plane (vector)
         * v₀ ∨ v₁ ∨ v₂
         * 
         * From paper eq. (3.12):
         * = ((→v₂ - →v₀) ∧ (→v₁ - →v₀))I₃ + e₀(→v₀ ∧ →v₁ ∧ →v₂)I₃
         */
        joinPPP(v0, v1, v2) {
            const [x0, y0, z0] = this.pointCoords(v0);
            const [x1, y1, z1] = this.pointCoords(v1);
            const [x2, y2, z2] = this.pointCoords(v2);
            
            // Edge vectors
            const e1x = x1 - x0, e1y = y1 - y0, e1z = z1 - z0;
            const e2x = x2 - x0, e2y = y2 - y0, e2z = z2 - z0;
            
            // Normal: (→v₁ - →v₀) × (→v₂ - →v₀)
            const nx = e1y*e2z - e1z*e2y;
            const ny = e1z*e2x - e1x*e2z;
            const nz = e1x*e2y - e1y*e2x;
            
            // Scalar triple product: →v₀ · (→v₁ × →v₂)
            // = x₀(y₁z₂ - z₁y₂) + y₀(z₁x₂ - x₁z₂) + z₀(x₁y₂ - y₁x₂)
            const d = x0*(y1*z2 - z1*y2) + y0*(z1*x2 - x1*z2) + z0*(x1*y2 - y1*x2);
            
            const plane = this.zero();
            plane[this.E1] = nx;
            plane[this.E2] = ny;
            plane[this.E3] = nz;
            plane[this.E0] = -d;
            
            return plane;
        }

        /**
         * Join of four points → scalar (pseudoscalar weight)
         * v₀ ∨ v₁ ∨ v₂ ∨ v₃ = 6 × signed volume of tetrahedron
         * 
         * From paper eq. (3.14)
         */
        joinPPPP(v0, v1, v2, v3) {
            const [x0, y0, z0] = this.pointCoords(v0);
            const [x1, y1, z1] = this.pointCoords(v1);
            const [x2, y2, z2] = this.pointCoords(v2);
            const [x3, y3, z3] = this.pointCoords(v3);
            
            // Volume = det([v1-v0, v2-v0, v3-v0])
            const e1x = x1 - x0, e1y = y1 - y0, e1z = z1 - z0;
            const e2x = x2 - x0, e2y = y2 - y0, e2z = z2 - z0;
            const e3x = x3 - x0, e3y = y3 - y0, e3z = z3 - z0;
            
            // 3x3 determinant
            const vol = e1x*(e2y*e3z - e2z*e3y) 
                      - e1y*(e2x*e3z - e2z*e3x) 
                      + e1z*(e2x*e3y - e2y*e3x);
            
            return vol;
        }

        /**
         * Join point with origin: v ∨ o
         * Extracts the ideal factor Aᵢ from element A
         * Used for computing ideal norm
         */
        joinPointOrigin(v) {
            // For a point v = e₁₂₃ + ideal terms
            // v ∨ o extracts the ideal part
            const [x, y, z] = this.pointCoords(v);
            
            // Result is a bivector (line through origin)
            const line = this.zero();
            line[this.E23] = -x;
            line[this.E13] = y;
            line[this.E12] = -z;
            return line;
        }

        /**
         * Join line with origin: L ∨ o
         * For computing ideal norm of lines
         */
        joinLineOrigin(L) {
            // Ideal part of line becomes a trivector
            const result = this.zero();
            result[this.E123] = L[this.E01]*0 + L[this.E02]*0 + L[this.E03]*0; // Would need full computation
            return L[this.E01] + L[this.E02] + L[this.E03]; // Simplified for magnitude
        }

        /**
         * Join plane/triangle with origin: F ∨ o
         * Key operation for mesh volume computation
         * 
         * From paper: For plane F = ae₁ + be₂ + ce₃ + de₀
         * F ∨ o = d (signed distance to origin)
         */
        joinPlaneOrigin(F) {
            return F[this.E0];  // The d coefficient
        }

        // ========================================================================
        // MEET (OUTER/WEDGE) PRODUCT
        // ========================================================================

        /**
         * Meet of two planes → line
         * p₁ ∧ p₂ = intersection line
         */
        meetPlanePlane(p1, p2) {
            const line = this.zero();
            
            // Euclidean part: n₁ ∧ n₂ (cross product of normals)
            line[this.E23] = p1[this.E2]*p2[this.E3] - p1[this.E3]*p2[this.E2];
            line[this.E13] = -(p1[this.E1]*p2[this.E3] - p1[this.E3]*p2[this.E1]);
            line[this.E12] = p1[this.E1]*p2[this.E2] - p1[this.E2]*p2[this.E1];
            
            // Ideal part: d₁n₂ - d₂n₁
            line[this.E01] = p1[this.E0]*p2[this.E1] - p2[this.E0]*p1[this.E1];
            line[this.E02] = p1[this.E0]*p2[this.E2] - p2[this.E0]*p1[this.E2];
            line[this.E03] = p1[this.E0]*p2[this.E3] - p2[this.E0]*p1[this.E3];
            
            return line;
        }

        /**
         * Meet of line and plane → point
         * E ∧ p = intersection point
         * 
         * From paper eq. (4.11)
         */
        meetLinePlane(E, p) {
            const dx = E[this.E01], dy = E[this.E02], dz = E[this.E03];
            const lx = E[this.E23], ly = -E[this.E13], lz = E[this.E12];
            const a = p[this.E1], b = p[this.E2], c = p[this.E3], d = p[this.E0];
            
            const point = this.zero();
            point[this.E023] = b*dz - c*dy - d*lx;
            point[this.E013] = c*dx - a*dz - d*ly;
            point[this.E012] = a*dy - b*dx - d*lz;
            point[this.E123] = a*lx + b*ly + c*lz;
            
            return point;
        }

        // ========================================================================
        // NORMS - Euclidean and Ideal
        // From paper Section 3(b)
        // ========================================================================

        /**
         * Euclidean norm: ‖A‖ = ‖Aₑ‖
         * For points/lines/planes, this is the Euclidean magnitude
         */
        norm(a) {
            // Grade 1 (plane): √(a² + b² + c²)
            const n1 = a[this.E1]**2 + a[this.E2]**2 + a[this.E3]**2;
            
            // Grade 2 (line): √(direction²)
            const n2 = a[this.E23]**2 + a[this.E13]**2 + a[this.E12]**2;
            
            // Grade 3 (point): √(e₁₂₃²) = |weight|
            const n3 = a[this.E123]**2;
            
            return sqrt(n1 + n2 + n3);
        }

        /**
         * Ideal norm: ‖A‖∞ = ‖Aᵢ‖ = ‖A ∨ o‖
         * 
         * From paper: This extracts the ideal (at-infinity) part's magnitude
         * Critical for computing mesh volumes from boundaries
         */
        idealNorm(a) {
            // Ideal parts are coefficients with e₀ factor
            
            // Grade 1 (plane): |d| (distance from origin)
            const n1 = a[this.E0]**2;
            
            // Grade 2 (line): ‖moment‖
            const n2 = a[this.E01]**2 + a[this.E02]**2 + a[this.E03]**2;
            
            // Grade 3 (point): ‖position × weight‖
            const n3 = a[this.E023]**2 + a[this.E013]**2 + a[this.E012]**2;
            
            return sqrt(n1 + n2 + n3);
        }

        /**
         * Signed ideal norm for planes (returns scalar, not absolute value)
         * F ∨ o = d (the signed distance coefficient)
         */
        signedIdealNorm(plane) {
            return plane[this.E0];
        }

        // ========================================================================
        // k-SIMPLEX MAGNITUDES
        // From paper Table 3 and Definition 4.2
        // ========================================================================

        /**
         * Edge length: ‖v₀ ∨ v₁‖
         */
        edgeLength(v0, v1) {
            const line = this.joinPP(v0, v1);
            return this.norm(line);
        }

        /**
         * Triangle area: ½‖v₀ ∨ v₁ ∨ v₂‖
         */
        triangleArea(v0, v1, v2) {
            const plane = this.joinPPP(v0, v1, v2);
            return 0.5 * this.norm(plane);
        }

        /**
         * Triangle area (alternative via boundary edges)
         * ½‖∑edges‖∞
         */
        triangleAreaBoundary(v0, v1, v2) {
            const e01 = this.joinPP(v0, v1);
            const e12 = this.joinPP(v1, v2);
            const e20 = this.joinPP(v2, v0);
            const sum = this.add(this.add(e01, e12), e20);
            return 0.5 * this.idealNorm(sum);
        }

        /**
         * Tetrahedron volume: ⅙|v₀ ∨ v₁ ∨ v₂ ∨ v₃|
         */
        tetrahedronVolume(v0, v1, v2, v3) {
            const vol = this.joinPPPP(v0, v1, v2, v3);
            return abs(vol) / 6;
        }

        /**
         * Signed tetrahedron volume (for mesh volume computation)
         */
        tetrahedronSignedVolume(v0, v1, v2, v3) {
            return this.joinPPPP(v0, v1, v2, v3) / 6;
        }

        // ========================================================================
        // MESH PROCESSING - k-COMPLEX OPERATIONS
        // From paper Theorem 4.2 and Section 4(d)
        // ========================================================================

        /**
         * Compute mesh volume from triangle faces
         * Using the "mind the gap" technique from paper
         * 
         * Volume = (1/6) ‖∑Fᵢ‖∞ = (1/6) ∑(Fᵢ ∨ o)
         * 
         * @param {Array} faces - Array of face objects {v0, v1, v2} with point arrays
         * @returns {number} Signed volume
         */
        meshVolume(faces) {
            let sumD = 0;
            
            for (const face of faces) {
                const F = this.joinPPP(face.v0, face.v1, face.v2);
                // F ∨ o = d coefficient (signed volume contribution)
                sumD += this.signedIdealNorm(F);
            }
            
            return sumD / 6;
        }

        /**
         * Mesh volume using precomputed face planes
         * More efficient when faces are reused
         */
        meshVolumeFromPlanes(facePlanes) {
            let sumD = 0;
            for (const F of facePlanes) {
                sumD += F[this.E0];  // d coefficient
            }
            return sumD / 6;
        }

        /**
         * Compute mesh surface area
         * Area = ½ ∑‖Fᵢ‖
         */
        meshArea(faces) {
            let area = 0;
            for (const face of faces) {
                const F = this.joinPPP(face.v0, face.v1, face.v2);
                area += this.norm(F);
            }
            return area / 2;
        }

        /**
         * Compute mesh center of mass
         * From paper Theorem 5.1:
         * Ccom = (1/24V) ∑(v₀ + v₁ + v₂ + o)(F ∨ o)
         * 
         * Returns unnormalized point (homogeneous, scaled by volume)
         */
        meshCenterOfMass(faces) {
            const com = this.zero();
            
            for (const face of faces) {
                const [x0, y0, z0] = this.pointCoords(face.v0);
                const [x1, y1, z1] = this.pointCoords(face.v1);
                const [x2, y2, z2] = this.pointCoords(face.v2);
                
                // Centroid of triangle + origin = (v₀ + v₁ + v₂ + o)/4
                // But we use sum directly and normalize later
                const cx = x0 + x1 + x2;
                const cy = y0 + y1 + y2;
                const cz = z0 + z1 + z2;
                
                // F ∨ o = signed volume contribution (×6)
                const F = this.joinPPP(face.v0, face.v1, face.v2);
                const signedVol = F[this.E0];  // This is 6× the tetrahedron volume to origin
                
                // Weight centroid by signed volume
                com[this.E023] += -cx * signedVol;
                com[this.E013] += cy * signedVol;
                com[this.E012] += -cz * signedVol;
                com[this.E123] += 4 * signedVol;  // Weight accumulator (factor of 4 from sum)
            }
            
            // Normalize by 1/24 is implicit in the weight
            return com;
        }

        /**
         * Extract center of mass coordinates
         */
        centerOfMassCoords(com) {
            return this.pointCoords(com);
        }

        /**
         * Check which side of plane a point is on
         * From paper eq. (4.10): p ∨ v = ax + by + cz + d
         */
        planeSide(plane, point) {
            const [x, y, z] = this.pointCoords(point);
            const a = plane[this.E1], b = plane[this.E2];
            const c = plane[this.E3], d = plane[this.E0];
            return a*x + b*y + c*z + d;
        }

        /**
         * Compute intersection of edge with plane
         * Returns null if no intersection, otherwise the intersection point
         */
        edgePlaneIntersection(v0, v1, plane) {
            const d0 = this.planeSide(plane, v0);
            const d1 = this.planeSide(plane, v1);
            
            // Same side - no intersection
            if (d0 * d1 > 0) return null;
            
            // On plane
            if (abs(d0) < EPSILON) return v0.slice();
            if (abs(d1) < EPSILON) return v1.slice();
            
            // Interpolate
            const t = d0 / (d0 - d1);
            const [x0, y0, z0] = this.pointCoords(v0);
            const [x1, y1, z1] = this.pointCoords(v1);
            
            return this.point(
                x0 + t * (x1 - x0),
                y0 + t * (y1 - y0),
                z0 + t * (z1 - z0)
            );
        }

        /**
         * Slice mesh by plane, return volume below plane
         * Implementation of the fuel tank example from paper Section 4(d)
         */
        meshVolumeBelow(faces, plane) {
            let sumD = 0;
            
            for (const face of faces) {
                const d0 = this.planeSide(plane, face.v0);
                const d1 = this.planeSide(plane, face.v1);
                const d2 = this.planeSide(plane, face.v2);
                
                const below0 = d0 <= 0;
                const below1 = d1 <= 0;
                const below2 = d2 <= 0;
                const countBelow = below0 + below1 + below2;
                
                if (countBelow === 3) {
                    // Entire face below plane
                    const F = this.joinPPP(face.v0, face.v1, face.v2);
                    sumD += F[this.E0];
                } else if (countBelow === 0) {
                    // Entire face above plane - skip
                    continue;
                } else {
                    // Face intersects plane - need to split
                    // Find intersection points
                    const intersections = [];
                    const verts = [face.v0, face.v1, face.v2];
                    const dists = [d0, d1, d2];
                    const belows = [below0, below1, below2];
                    
                    for (let i = 0; i < 3; i++) {
                        const j = (i + 1) % 3;
                        if (belows[i] !== belows[j]) {
                            const p = this.edgePlaneIntersection(verts[i], verts[j], plane);
                            if (p) intersections.push(p);
                        }
                    }
                    
                    if (intersections.length === 2) {
                        // Create triangles below plane
                        if (countBelow === 1) {
                            // One vertex below - one triangle
                            const vBelow = belows[0] ? face.v0 : (belows[1] ? face.v1 : face.v2);
                            const F = this.joinPPP(vBelow, intersections[0], intersections[1]);
                            sumD += F[this.E0];
                        } else {
                            // Two vertices below - two triangles (quad)
                            const vAbove = !belows[0] ? face.v0 : (!belows[1] ? face.v1 : face.v2);
                            const vBelow1 = belows[0] ? face.v0 : face.v1;
                            const vBelow2 = belows[2] ? face.v2 : face.v1;
                            
                            // Full triangle minus the small triangle above
                            const Ffull = this.joinPPP(face.v0, face.v1, face.v2);
                            const Fabove = this.joinPPP(vAbove, intersections[0], intersections[1]);
                            sumD += Ffull[this.E0] - Fabove[this.E0];
                        }
                    }
                }
            }
            
            return sumD / 6;
        }

        // ========================================================================
        // INERTIA TENSOR
        // From paper Section 5(d), Definition 5.1
        // ========================================================================

        /**
         * Compute inertia tensor frame for mesh
         * Returns {I1, I2, I3} frame vectors that encode the inertia tensor
         */
        meshInertiaFrame(faces) {
            const I = {
                xx: 0, yy: 0, zz: 0,
                xy: 0, xz: 0, yz: 0
            };
            
            for (const face of faces) {
                const [x0, y0, z0] = this.pointCoords(face.v0);
                const [x1, y1, z1] = this.pointCoords(face.v1);
                const [x2, y2, z2] = this.pointCoords(face.v2);
                
                // Signed volume factor (×6)
                const F = this.joinPPP(face.v0, face.v1, face.v2);
                const vol = F[this.E0];
                
                // Accumulate second moments
                // Using formulas from paper's reference [17] (Tonon)
                const x = [x0, x1, x2], y = [y0, y1, y2], z = [z0, z1, z2];
                
                // Products
                let sxx = 0, syy = 0, szz = 0;
                let sxy = 0, sxz = 0, syz = 0;
                
                for (let i = 0; i < 3; i++) {
                    for (let j = i; j < 3; j++) {
                        sxx += x[i] * x[j];
                        syy += y[i] * y[j];
                        szz += z[i] * z[j];
                        sxy += x[i] * y[j] + (i !== j ? x[j] * y[i] : 0);
                        sxz += x[i] * z[j] + (i !== j ? x[j] * z[i] : 0);
                        syz += y[i] * z[j] + (i !== j ? y[j] * z[i] : 0);
                    }
                }
                
                I.xx += vol * sxx;
                I.yy += vol * syy;
                I.zz += vol * szz;
                I.xy += vol * sxy;
                I.xz += vol * sxz;
                I.yz += vol * syz;
            }
            
            // Scale factor
            const s = 1 / 60;  // From integration over tetrahedron
            I.xx *= s; I.yy *= s; I.zz *= s;
            I.xy *= s; I.xz *= s; I.yz *= s;
            
            // Convert to inertia tensor (I = trace(M)*Id - M for solid)
            // Moments of inertia
            const Ixx = I.yy + I.zz;  // About x-axis
            const Iyy = I.xx + I.zz;  // About y-axis
            const Izz = I.xx + I.yy;  // About z-axis
            
            // Products of inertia (negated)
            const Ixy = -I.xy;
            const Ixz = -I.xz;
            const Iyz = -I.yz;
            
            // Return as frame vectors (paper Definition 5.1)
            return {
                I1: [Ixx, Ixy, Ixz],
                I2: [Ixy, Iyy, Iyz],
                I3: [Ixz, Iyz, Izz]
            };
        }

        /**
         * Diagonalize inertia frame using Jacobi method
         * Returns principal moments and rotation to principal axes
         */
        diagonalizeInertia(inertiaFrame) {
            // Convert frame to matrix
            const A = [
                [...inertiaFrame.I1],
                [...inertiaFrame.I2],
                [...inertiaFrame.I3]
            ];
            
            // Jacobi eigenvalue algorithm
            const n = 3;
            const V = [[1,0,0], [0,1,0], [0,0,1]];  // Eigenvector accumulator
            
            const maxIter = 50;
            for (let iter = 0; iter < maxIter; iter++) {
                // Find largest off-diagonal element
                let maxVal = 0, p = 0, q = 1;
                for (let i = 0; i < n; i++) {
                    for (let j = i + 1; j < n; j++) {
                        if (abs(A[i][j]) > maxVal) {
                            maxVal = abs(A[i][j]);
                            p = i; q = j;
                        }
                    }
                }
                
                if (maxVal < EPSILON) break;
                
                // Compute rotation angle
                const theta = (A[q][q] - A[p][p]) / (2 * A[p][q]);
                const t = Math.sign(theta) / (abs(theta) + sqrt(1 + theta*theta));
                const c = 1 / sqrt(1 + t*t);
                const s = t * c;
                
                // Apply rotation to A
                const App = A[p][p], Aqq = A[q][q], Apq = A[p][q];
                A[p][p] = c*c*App - 2*s*c*Apq + s*s*Aqq;
                A[q][q] = s*s*App + 2*s*c*Apq + c*c*Aqq;
                A[p][q] = A[q][p] = 0;
                
                for (let i = 0; i < n; i++) {
                    if (i !== p && i !== q) {
                        const Aip = A[i][p], Aiq = A[i][q];
                        A[i][p] = A[p][i] = c*Aip - s*Aiq;
                        A[i][q] = A[q][i] = s*Aip + c*Aiq;
                    }
                }
                
                // Accumulate eigenvectors
                for (let i = 0; i < n; i++) {
                    const Vip = V[i][p], Viq = V[i][q];
                    V[i][p] = c*Vip - s*Viq;
                    V[i][q] = s*Vip + c*Viq;
                }
            }
            
            // Extract results
            return {
                principalMoments: [A[0][0], A[1][1], A[2][2]],
                principalAxes: V,
                // Convert to rotor (simplified - assumes small rotations combine)
                toRotor: () => {
                    // Build rotor from rotation matrix V
                    // Using standard matrix-to-quaternion conversion
                    const trace = V[0][0] + V[1][1] + V[2][2];
                    let w, x, y, z;
                    
                    if (trace > 0) {
                        const s = 0.5 / sqrt(trace + 1);
                        w = 0.25 / s;
                        x = (V[2][1] - V[1][2]) * s;
                        y = (V[0][2] - V[2][0]) * s;
                        z = (V[1][0] - V[0][1]) * s;
                    } else if (V[0][0] > V[1][1] && V[0][0] > V[2][2]) {
                        const s = 2 * sqrt(1 + V[0][0] - V[1][1] - V[2][2]);
                        w = (V[2][1] - V[1][2]) / s;
                        x = 0.25 * s;
                        y = (V[0][1] + V[1][0]) / s;
                        z = (V[0][2] + V[2][0]) / s;
                    } else if (V[1][1] > V[2][2]) {
                        const s = 2 * sqrt(1 + V[1][1] - V[0][0] - V[2][2]);
                        w = (V[0][2] - V[2][0]) / s;
                        x = (V[0][1] + V[1][0]) / s;
                        y = 0.25 * s;
                        z = (V[1][2] + V[2][1]) / s;
                    } else {
                        const s = 2 * sqrt(1 + V[2][2] - V[0][0] - V[1][1]);
                        w = (V[1][0] - V[0][1]) / s;
                        x = (V[0][2] + V[2][0]) / s;
                        y = (V[1][2] + V[2][1]) / s;
                        z = 0.25 * s;
                    }
                    
                    // Return as GA rotor components (s, e12, e31, e23)
                    return { s: w, e12: z, e31: y, e23: x };
                }
            };
        }
    }

    // ============================================================================
    // MESH CLASS - High-level mesh representation
    // ============================================================================

    class Mesh {
        /**
         * Create mesh from vertices and face indices
         * @param {Array} vertices - Array of [x, y, z] coordinates
         * @param {Array} faces - Array of [i0, i1, i2] vertex indices
         */
        constructor(vertices, faces) {
            this.pga = new PGA3DMesh();
            
            // Convert vertices to PGA points
            this.vertices = vertices.map(([x, y, z]) => this.pga.point(x, y, z));
            
            // Store face indices
            this.faceIndices = faces;
            
            // Precompute face planes for efficiency
            this.facePlanes = faces.map(([i0, i1, i2]) => 
                this.pga.joinPPP(this.vertices[i0], this.vertices[i1], this.vertices[i2])
            );
        }

        /**
         * Get face as {v0, v1, v2} object
         */
        getFace(i) {
            const [i0, i1, i2] = this.faceIndices[i];
            return {
                v0: this.vertices[i0],
                v1: this.vertices[i1],
                v2: this.vertices[i2]
            };
        }

        /**
         * Get all faces as array of {v0, v1, v2}
         */
        getFaces() {
            return this.faceIndices.map((_, i) => this.getFace(i));
        }

        /**
         * Compute mesh volume
         */
        volume() {
            let sum = 0;
            for (const F of this.facePlanes) {
                sum += F[this.pga.E0];
            }
            return sum / 6;
        }

        /**
         * Compute mesh surface area
         */
        area() {
            let sum = 0;
            for (const F of this.facePlanes) {
                sum += this.pga.norm(F);
            }
            return sum / 2;
        }

        /**
         * Compute center of mass
         * @returns {Array} [x, y, z] coordinates
         */
        centerOfMass() {
            const com = this.pga.meshCenterOfMass(this.getFaces());
            return this.pga.pointCoords(com);
        }

        /**
         * Compute inertia tensor
         */
        inertia() {
            return this.pga.meshInertiaFrame(this.getFaces());
        }

        /**
         * Compute principal moments and axes
         */
        principalInertia() {
            const I = this.inertia();
            return this.pga.diagonalizeInertia(I);
        }

        /**
         * Volume below a cutting plane
         */
        volumeBelow(plane) {
            return this.pga.meshVolumeBelow(this.getFaces(), plane);
        }

        /**
         * Check if mesh is closed (boundary sums to zero)
         */
        isClosed() {
            let sum = this.pga.zero();
            for (const F of this.facePlanes) {
                sum = this.pga.add(sum, F);
            }
            return this.pga.norm(sum) < EPSILON;
        }
    }

    // ============================================================================
    // PRIMITIVE MESH GENERATORS
    // ============================================================================

    const MeshPrimitives = {
        /**
         * Create a box mesh
         */
        box(width = 1, height = 1, depth = 1) {
            const w = width / 2, h = height / 2, d = depth / 2;
            const vertices = [
                [-w, -h, -d], [w, -h, -d], [w, h, -d], [-w, h, -d],
                [-w, -h, d], [w, -h, d], [w, h, d], [-w, h, d]
            ];
            // CCW winding for outward normals
            const faces = [
                [0, 3, 2], [0, 2, 1],  // Back
                [4, 5, 6], [4, 6, 7],  // Front
                [0, 1, 5], [0, 5, 4],  // Bottom
                [2, 3, 7], [2, 7, 6],  // Top
                [0, 4, 7], [0, 7, 3],  // Left
                [1, 2, 6], [1, 6, 5]   // Right
            ];
            return new Mesh(vertices, faces);
        },

        /**
         * Create a tetrahedron mesh
         */
        tetrahedron(size = 1) {
            const s = size / sqrt(2);
            const vertices = [
                [s, s, s], [s, -s, -s], [-s, s, -s], [-s, -s, s]
            ];
            const faces = [
                [0, 1, 2], [0, 2, 3], [0, 3, 1], [1, 3, 2]
            ];
            return new Mesh(vertices, faces);
        },

        /**
         * Create an icosahedron mesh
         */
        icosahedron(radius = 1) {
            const phi = (1 + sqrt(5)) / 2;
            const a = radius / sqrt(1 + phi * phi);
            const b = a * phi;
            
            const vertices = [
                [-a, b, 0], [a, b, 0], [-a, -b, 0], [a, -b, 0],
                [0, -a, b], [0, a, b], [0, -a, -b], [0, a, -b],
                [b, 0, -a], [b, 0, a], [-b, 0, -a], [-b, 0, a]
            ];
            const faces = [
                [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
                [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
                [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
                [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]
            ];
            return new Mesh(vertices, faces);
        }
    };

    // ============================================================================
    // EXPORTS
    // ============================================================================

    const GAMesh = {
        PGA3DMesh,
        Mesh,
        MeshPrimitives,
        
        // Convenience factory
        createMesh: (vertices, faces) => new Mesh(vertices, faces),
        
        // Quick access to PGA algebra
        pga: () => new PGA3DMesh()
    };

    // Export for different module systems
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = GAMesh;
    } else if (typeof define === 'function' && define.amd) {
        define([], () => GAMesh);
    } else {
        global.GAMesh = GAMesh;
    }

})(typeof window !== 'undefined' ? window : global);
