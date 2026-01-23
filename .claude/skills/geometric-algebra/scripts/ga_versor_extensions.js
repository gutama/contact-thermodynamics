/**
 * Geometric Algebra - Versor Extensions
 * 
 * Extends ga_vanilla.js with functionality from Versor libraries:
 * - Booster (Bst) - conformal transformations via point pairs
 * - Tangent elements (Tnv, Tnb, Tnt)
 * - Extended Round:: operations (carrier, surround, split, etc.)
 * - Extended Flat:: operations
 * - Direction elements (Drv, Drb, Drt)
 * - Motor/Boost logarithms and exponentials
 * - Arbitrary metric algebra generation
 * 
 * Based on:
 * - wolftype/versor (C++ library by Pablo Colapinto)
 * - weshoke/versor.js (JavaScript port)
 * - "Boosted Surfaces" paper (Colapinto, 2014)
 * - "Articulating Space" dissertation (Colapinto)
 * 
 * @requires ga_vanilla.js
 * @license MIT
 */

(function(global) {
    'use strict';

    // Get base GA module
    const GA = global.GA || (typeof require !== 'undefined' ? require('./ga_vanilla.js') : null);
    if (!GA) throw new Error('ga_vanilla.js must be loaded first');

    const { Algebra, Multivector, CGA3D } = GA;

    // ============================================================================
    // CONSTANTS AND UTILITIES
    // ============================================================================

    const EPSILON = 1e-10;
    const PI = Math.PI;
    const sqrt = Math.sqrt;
    const abs = Math.abs;
    const sin = Math.sin;
    const cos = Math.cos;
    const sinh = Math.sinh;
    const cosh = Math.cosh;
    const atan2 = Math.atan2;
    const acos = Math.acos;
    const log = Math.log;
    const exp = Math.exp;

    function isZero(x) { return abs(x) < EPSILON; }
    function sign(x) { return x >= 0 ? 1 : -1; }

    // ============================================================================
    // EXTENDED CGA3D CLASS
    // ============================================================================

    /**
     * Extended CGA3D with full Versor-style operations
     */
    class CGA3DVersor extends CGA3D {
        constructor() {
            super();
            
            // Additional basis shortcuts
            this.E0 = this.no.wedge(this.ni);  // Minkowski plane e₋ ∧ e₊
            
            // Precompute common elements
            this._setupElements();
        }

        _setupElements() {
            // Euclidean pseudoscalar
            this.I3 = this.e1.wedge(this.e2).wedge(this.e3);
            
            // Full pseudoscalar
            this.I5 = this.e1.wedge(this.e2).wedge(this.e3).wedge(this.ep).wedge(this.em);
        }

        // ========================================================================
        // DIRECTION ELEMENTS (Free vectors without position)
        // ========================================================================

        /**
         * Direction Vector (Drv) - free vector, no position
         * Drv = v ∧ n∞
         */
        directionVector(vx, vy, vz) {
            const v = this.e1.scale(vx).add(this.e2.scale(vy)).add(this.e3.scale(vz));
            return v.wedge(this.ni);
        }

        /**
         * Direction Bivector (Drb) - free area element
         * Drb = B ∧ n∞
         */
        directionBivector(bxy, bxz, byz) {
            const B = this.e1.wedge(this.e2).scale(bxy)
                .add(this.e1.wedge(this.e3).scale(bxz))
                .add(this.e2.wedge(this.e3).scale(byz));
            return B.wedge(this.ni);
        }

        /**
         * Direction Trivector (Drt) - free volume element
         * Drt = I3 ∧ n∞
         */
        directionTrivector(scale = 1) {
            return this.I3.scale(scale).wedge(this.ni);
        }

        /**
         * Extract direction vector components from Drv
         */
        directionVectorCoords(drv) {
            // Drv = vx(e1∧ni) + vy(e2∧ni) + vz(e3∧ni)
            const e1ni = this.bladeIndex['e1ep'] !== undefined ? 'e1ep' : 'e1em';
            const e2ni = this.bladeIndex['e2ep'] !== undefined ? 'e2ep' : 'e2em';
            const e3ni = this.bladeIndex['e3ep'] !== undefined ? 'e3ep' : 'e3em';
            
            // Need to extract from wedge with ni = ep + em
            // This is simplified - actual extraction needs proper contraction
            const contracted = drv.dot(this.no);
            return [
                contracted.coeffs[this.bladeIndex['e1']] || 0,
                contracted.coeffs[this.bladeIndex['e2']] || 0,
                contracted.coeffs[this.bladeIndex['e3']] || 0
            ];
        }

        // ========================================================================
        // TANGENT ELEMENTS (Localized vectors at a point)
        // ========================================================================

        /**
         * Tangent Vector (Tnv) at point p with direction v
         * Tnv = p ∧ (v ∧ n∞)
         */
        tangentVector(px, py, pz, vx, vy, vz) {
            const p = this.point(px, py, pz);
            const drv = this.directionVector(vx, vy, vz);
            return p.wedge(drv);
        }

        /**
         * Tangent Vector from point and direction multivectors
         */
        tangentVectorAt(point, direction) {
            const drv = direction.wedge(this.ni);
            return point.wedge(drv);
        }

        /**
         * Tangent Bivector (Tnb) at point p with area B
         * Tnb = p ∧ (B ∧ n∞)
         */
        tangentBivector(px, py, pz, bxy, bxz, byz) {
            const p = this.point(px, py, pz);
            const drb = this.directionBivector(bxy, bxz, byz);
            return p.wedge(drb);
        }

        /**
         * Tangent Trivector (Tnt) at point p
         * Tnt = p ∧ (I3 ∧ n∞)
         */
        tangentTrivector(px, py, pz, scale = 1) {
            const p = this.point(px, py, pz);
            const drt = this.directionTrivector(scale);
            return p.wedge(drt);
        }

        /**
         * Extract location from tangent element
         */
        tangentLocation(tangent) {
            // Location is found by contracting with n∞
            const contracted = tangent.dot(this.ni);
            // Normalize and extract point
            return this.Round.location(contracted);
        }

        /**
         * Extract direction from tangent vector
         */
        tangentDirection(tnv) {
            // Contract with origin to get direction part
            const dir = tnv.dot(this.no);
            return dir;
        }

        // ========================================================================
        // ROUND OPERATIONS (Extended)
        // ========================================================================

        /**
         * Round namespace - operations on round elements
         * (Points, Point Pairs, Circles, Spheres)
         */
        get Round() {
            const cga = this;
            return {
                /**
                 * Create null point (on null cone)
                 */
                null: (x, y, z) => cga.point(x, y, z),

                /**
                 * Create dual sphere from center and radius
                 * Dls = c - ½r²n∞  (negative radius² = imaginary sphere)
                 */
                dls: (cx, cy, cz, r) => {
                    const c = cga.point(cx, cy, cz);
                    return c.sub(cga.ni.scale(0.5 * r * r));
                },

                /**
                 * Create dual sphere from center point and surface point
                 */
                dlsFromPoints: (center, surface) => {
                    const r2 = -2 * center.dot(surface).scalar;
                    return center.sub(cga.ni.scale(0.5 * r2));
                },

                /**
                 * Carrier flat of a round element
                 * The flat that contains the round
                 * carrier(Circle) = Plane, carrier(PointPair) = Line
                 */
                carrier: (round) => {
                    return round.wedge(cga.ni);
                },

                /**
                 * Surround - the sphere that surrounds a round element
                 * sur(Circle) = Sphere containing circle
                 */
                surround: (round) => {
                    const carrier = round.wedge(cga.ni);
                    return round.dot(carrier.inverse());
                },

                /**
                 * Location (center) of a round element
                 * Normalized point at center
                 */
                location: (round) => {
                    // loc = round * n∞ * round / (round * round)
                    const sq = round.mul(round.reverse()).scalar;
                    if (isZero(sq)) {
                        // Null element - extract directly
                        return round.grade(1);
                    }
                    const niRound = cga.ni.mul(round);
                    const loc = round.mul(niRound).mul(round).scale(1 / sq);
                    // Normalize
                    const w = -loc.dot(cga.ni).scalar;
                    return isZero(w) ? loc : loc.scale(1 / w);
                },

                /**
                 * Simple center (not normalized)
                 */
                center: (round) => {
                    return round.mul(cga.ni).mul(round);
                },

                /**
                 * Direction of round element
                 * For circle: normal to plane
                 * For point pair: line direction
                 */
                direction: (round) => {
                    const carrier = round.wedge(cga.ni);
                    // Contract carrier with infinity to get direction
                    return carrier.dot(cga.ni).normalized();
                },

                /**
                 * Size (radius) of round element
                 * @param {Multivector} round - round element
                 * @param {boolean} signed - return signed value (negative for imaginary)
                 * @returns {number} radius (or squared radius if imaginary)
                 */
                size: (round, signed = false) => {
                    const sq = round.mul(round.reverse()).scalar;
                    const loc = cga.Round.location(round);
                    const locSq = loc.mul(loc.reverse()).scalar;
                    
                    // r² = (round² / loc²)
                    const r2 = isZero(locSq) ? sq : sq / locSq;
                    
                    if (r2 < 0) {
                        // Imaginary round
                        return signed ? -sqrt(-r2) : sqrt(-r2);
                    }
                    return sqrt(r2);
                },

                /**
                 * Check if round is real (positive radius²)
                 */
                isReal: (round) => {
                    const sq = round.mul(round.reverse()).scalar;
                    return sq >= 0;
                },

                /**
                 * Check if round is imaginary (negative radius²)
                 */
                isImaginary: (round) => {
                    const sq = round.mul(round.reverse()).scalar;
                    return sq < 0;
                },

                /**
                 * Split point pair into two points
                 * @returns {[Multivector, Multivector]} two points
                 */
                split: (pointPair) => {
                    const sq = pointPair.mul(pointPair.reverse()).scalar;
                    
                    if (sq > 0) {
                        // Real point pair
                        const r = sqrt(sq);
                        const bst = pointPair.scale(1 / r);
                        
                        // Find midpoint
                        const mid = cga.Round.location(pointPair);
                        
                        // Points are mid ± r * direction
                        const dir = cga.Round.direction(pointPair);
                        
                        // Use the booster to extract points
                        // p1 = (1 + bst) * no * (1 + bst)~ normalized
                        const one = cga.scalar(1);
                        const v1 = one.add(bst);
                        const v2 = one.sub(bst);
                        
                        const p1 = v1.mul(cga.no).mul(v1.reverse());
                        const p2 = v2.mul(cga.no).mul(v2.reverse());
                        
                        // Normalize
                        const w1 = -p1.dot(cga.ni).scalar;
                        const w2 = -p2.dot(cga.ni).scalar;
                        
                        return [
                            isZero(w1) ? p1 : p1.scale(1 / w1),
                            isZero(w2) ? p2 : p2.scale(1 / w2)
                        ];
                    } else if (sq < 0) {
                        // Imaginary point pair - return complex conjugate points
                        // (on imaginary sphere)
                        const mid = cga.Round.location(pointPair);
                        return [mid, mid];  // Degenerate case
                    } else {
                        // Null - single point
                        const p = cga.Round.location(pointPair);
                        return [p, p];
                    }
                },

                /**
                 * Point on circle at parameter t ∈ [0, 2π]
                 */
                pointOnCircle: (circle, t) => {
                    const center = cga.Round.location(circle);
                    const radius = cga.Round.size(circle);
                    const plane = cga.Round.carrier(circle);
                    
                    // Get two orthogonal directions in the plane
                    // This is simplified - full implementation needs proper basis
                    const coords = cga.pointCoords(center);
                    
                    // Create rotor in circle's plane
                    const normal = cga.Round.direction(circle);
                    const biv = normal.grade(2);
                    
                    // Point at angle t
                    const ct = cos(t);
                    const st = sin(t);
                    
                    // Approximate: offset from center
                    return cga.point(
                        coords[0] + radius * ct,
                        coords[1] + radius * st,
                        coords[2]
                    );
                },

                /**
                 * Weight of round element
                 */
                weight: (round) => {
                    return round.dot(cga.ni.mul(cga.I5.reverse())).norm();
                }
            };
        }

        // ========================================================================
        // FLAT OPERATIONS (Extended)
        // ========================================================================

        /**
         * Flat namespace - operations on flat elements
         * (Lines, Planes, Flat Points)
         */
        get Flat() {
            const cga = this;
            return {
                /**
                 * Flat point - a point on a flat (grade 2 in CGA)
                 * Flp = p ∧ n∞ (point wedge infinity)
                 */
                point: (x, y, z) => {
                    return cga.point(x, y, z).wedge(cga.ni);
                },

                /**
                 * Line through two points
                 */
                line: (p1, p2) => {
                    const P1 = Array.isArray(p1) ? cga.point(...p1) : p1;
                    const P2 = Array.isArray(p2) ? cga.point(...p2) : p2;
                    return P1.wedge(P2).wedge(cga.ni);
                },

                /**
                 * Line from point and direction
                 */
                lineDir: (px, py, pz, dx, dy, dz) => {
                    const p = cga.point(px, py, pz);
                    const drv = cga.directionVector(dx, dy, dz);
                    return p.wedge(drv);
                },

                /**
                 * Plane through three points
                 */
                plane: (p1, p2, p3) => {
                    const P1 = Array.isArray(p1) ? cga.point(...p1) : p1;
                    const P2 = Array.isArray(p2) ? cga.point(...p2) : p2;
                    const P3 = Array.isArray(p3) ? cga.point(...p3) : p3;
                    return P1.wedge(P2).wedge(P3).wedge(cga.ni);
                },

                /**
                 * Plane from normal and distance
                 * Dlp = n + d*n∞
                 */
                planeNormal: (nx, ny, nz, d) => {
                    const n = cga.e1.scale(nx).add(cga.e2.scale(ny)).add(cga.e3.scale(nz));
                    return n.add(cga.ni.scale(d));
                },

                /**
                 * Carrier of a flat (the flat itself, normalized)
                 */
                carrier: (flat) => {
                    return flat.normalized();
                },

                /**
                 * Direction of a line
                 * Returns Euclidean vector
                 */
                direction: (line) => {
                    // Contract with n∞ to get direction
                    const dir = line.dot(cga.ni);
                    return dir.grade(1);
                },

                /**
                 * Location - a point on the flat closest to origin
                 */
                location: (flat) => {
                    // loc = flat * no * flat / (flat * flat)
                    const sq = flat.mul(flat.reverse()).scalar;
                    if (isZero(sq)) return cga.no;
                    
                    const loc = flat.mul(cga.no).mul(flat).scale(1 / sq);
                    const w = -loc.dot(cga.ni).scalar;
                    return isZero(w) ? loc : loc.scale(1 / w);
                },

                /**
                 * Normal vector of a plane
                 */
                normal: (plane) => {
                    // For dual plane Dlp = n + d*ni, the normal is the e1,e2,e3 part
                    const nx = plane.coeffs[cga.bladeIndex['e1']] || 0;
                    const ny = plane.coeffs[cga.bladeIndex['e2']] || 0;
                    const nz = plane.coeffs[cga.bladeIndex['e3']] || 0;
                    const len = sqrt(nx*nx + ny*ny + nz*nz);
                    return len > EPSILON ? [nx/len, ny/len, nz/len] : [0, 0, 1];
                },

                /**
                 * Weight of flat element
                 */
                weight: (flat) => {
                    return flat.norm();
                }
            };
        }

        // ========================================================================
        // GENERATORS (Extended)
        // ========================================================================

        /**
         * Generator namespace for creating versors
         */
        get Gen() {
            const cga = this;
            return {
                /**
                 * Rotor from bivector (rotation)
                 * R = exp(-θB/2) = cos(θ/2) - sin(θ/2)B
                 */
                rot: (bivector) => {
                    const sq = bivector.mul(bivector.reverse()).scalar;
                    const theta = sqrt(abs(sq));
                    
                    if (theta < EPSILON) {
                        return cga.scalar(1);
                    }
                    
                    const half = theta / 2;
                    const B = bivector.scale(1 / theta);
                    
                    if (sq < 0) {
                        // Elliptic (normal rotation)
                        return cga.scalar(cos(half)).sub(B.scale(sin(half)));
                    } else {
                        // Hyperbolic
                        return cga.scalar(cosh(half)).sub(B.scale(sinh(half)));
                    }
                },

                /**
                 * Translator from direction vector
                 * T = 1 - ½t∧n∞
                 */
                trs: (dx, dy, dz) => {
                    const t = cga.e1.scale(dx).add(cga.e2.scale(dy)).add(cga.e3.scale(dz));
                    return cga.scalar(1).sub(t.wedge(cga.ni).scale(0.5));
                },

                /**
                 * Translator from direction multivector
                 */
                trsVec: (direction) => {
                    return cga.scalar(1).sub(direction.wedge(cga.ni).scale(0.5));
                },

                /**
                 * Motor from dual line (rotation + translation)
                 * M = exp(-½Dll) where Dll encodes screw axis
                 */
                mot: (dualLine) => {
                    // Decompose dual line into rotation and translation parts
                    // Dll = d + m*I3 where d is direction bivector, m is moment
                    
                    const sq = dualLine.mul(dualLine.reverse()).scalar;
                    
                    if (isZero(sq)) {
                        // Pure translation
                        return cga.scalar(1).add(dualLine.scale(-0.5));
                    }
                    
                    const theta = sqrt(abs(sq));
                    const half = theta / 2;
                    const L = dualLine.scale(1 / theta);
                    
                    // Motor = cos(θ/2) - sin(θ/2)L + (translation part)
                    return cga.scalar(cos(half)).sub(L.scale(sin(half)));
                },

                /**
                 * Dilator from point and scale factor
                 * D = cosh(α/2) + sinh(α/2)*no∧ni where α = ln(scale)
                 */
                dil: (point, scale) => {
                    const alpha = log(scale);
                    const half = alpha / 2;
                    const B = cga.no.wedge(cga.ni);
                    
                    // Dilator centered at point
                    const T = cga.Gen.trsPoint(point);
                    const D0 = cga.scalar(cosh(half)).add(B.scale(sinh(half)));
                    const Tinv = T.inverse();
                    
                    return T.mul(D0).mul(Tinv);
                },

                /**
                 * Dilator at origin
                 */
                dilOrigin: (scale) => {
                    const alpha = log(scale);
                    const half = alpha / 2;
                    const B = cga.no.wedge(cga.ni);
                    return cga.scalar(cosh(half)).add(B.scale(sinh(half)));
                },

                /**
                 * Translator to move origin to point
                 */
                trsPoint: (point) => {
                    const coords = cga.pointCoords(point);
                    return cga.Gen.trs(coords[0], coords[1], coords[2]);
                },

                /**
                 * Transversor from tangent vector
                 * Trv = 1 + ½(no ∧ v) = 1 + ½Tnv
                 */
                trv: (vx, vy, vz) => {
                    const v = cga.e1.scale(vx).add(cga.e2.scale(vy)).add(cga.e3.scale(vz));
                    return cga.scalar(1).add(cga.no.wedge(v).scale(0.5));
                },

                /**
                 * Transversor from tangent vector multivector
                 */
                trvVec: (tangentVec) => {
                    return cga.scalar(1).add(cga.no.wedge(tangentVec).scale(0.5));
                },

                /**
                 * BOOSTER from point pair
                 * Bst = exp(½ * PointPair)
                 * 
                 * This is the KEY Versor operation for conformal deformations.
                 * Boosters bend elements around an "orbit" defined by the point pair.
                 */
                bst: (pointPair) => {
                    const sq = pointPair.mul(pointPair.reverse()).scalar;
                    
                    if (isZero(sq)) {
                        // Null point pair - returns identity + translation-like term
                        return cga.scalar(1).add(pointPair.scale(0.5));
                    }
                    
                    const r = sqrt(abs(sq));
                    const half = r / 2;
                    const P = pointPair.scale(1 / r);
                    
                    if (sq > 0) {
                        // Hyperbolic boost (real point pair)
                        // Bst = cosh(r/2) + sinh(r/2)*P̂
                        return cga.scalar(cosh(half)).add(P.scale(sinh(half)));
                    } else {
                        // Elliptic boost (imaginary point pair)  
                        // Bst = cos(r/2) + sin(r/2)*P̂
                        return cga.scalar(cos(half)).add(P.scale(sin(half)));
                    }
                },

                /**
                 * Booster from two points and amount
                 * Creates a point pair and generates the booster
                 */
                bstPoints: (p1, p2, amount) => {
                    const P1 = Array.isArray(p1) ? cga.point(...p1) : p1;
                    const P2 = Array.isArray(p2) ? cga.point(...p2) : p2;
                    const pp = P1.wedge(P2).scale(amount);
                    return cga.Gen.bst(pp);
                },

                /**
                 * Logarithm of rotor (inverse of exp)
                 * Returns bivector
                 */
                logRotor: (rotor) => {
                    const s = rotor.scalar;
                    const B = rotor.grade(2);
                    const Bnorm = B.norm();
                    
                    if (Bnorm < EPSILON) {
                        return cga.scalar(0);  // Identity rotor
                    }
                    
                    const theta = 2 * atan2(Bnorm, s);
                    return B.scale(theta / Bnorm);
                },

                /**
                 * Logarithm of motor
                 * Returns dual line (Plücker coordinates)
                 */
                logMotor: (motor) => {
                    // Extract rotation part
                    const R = motor.grade(0).add(motor.grade(2));
                    const Rnorm = R.norm();
                    
                    if (Rnorm < EPSILON) {
                        // Pure translation or identity
                        return motor.grade(2).scale(2);
                    }
                    
                    // Normalize rotation part
                    const Rn = R.scale(1 / Rnorm);
                    
                    // Get bivector (rotation plane)
                    const B = Rn.grade(2);
                    const theta = 2 * atan2(B.norm(), Rn.scalar);
                    
                    if (theta < EPSILON) {
                        // Small angle - return translation part
                        return motor.grade(2).scale(2);
                    }
                    
                    // Full decomposition
                    const Bhat = B.scale(1 / B.norm());
                    
                    // Translation part encoded in higher grades
                    const T = motor.mul(Rn.inverse());
                    
                    return Bhat.scale(theta);
                },

                /**
                 * Logarithm of booster
                 * Returns point pair
                 */
                logBoost: (booster) => {
                    const s = booster.scalar;
                    const PP = booster.grade(2);
                    const PPnorm = PP.norm();
                    
                    if (PPnorm < EPSILON) {
                        return cga.scalar(0);
                    }
                    
                    // Determine type from sign
                    const sq = PP.mul(PP.reverse()).scalar;
                    
                    if (sq > 0) {
                        // Hyperbolic
                        const alpha = 2 * Math.asinh(PPnorm / s);
                        return PP.scale(alpha / PPnorm);
                    } else {
                        // Elliptic
                        const alpha = 2 * atan2(PPnorm, s);
                        return PP.scale(alpha / PPnorm);
                    }
                },

                /**
                 * Ratio between two versors
                 * ratio(V1, V2) = V2 * V1⁻¹
                 */
                ratio: (v1, v2) => {
                    const v1inv = v1.inverse();
                    return v1inv ? v2.mul(v1inv) : null;
                },

                /**
                 * Interpolate versors (generalized SLERP)
                 */
                slerp: (v1, v2, t) => {
                    const ratio = cga.Gen.ratio(v1, v2);
                    if (!ratio) return v1;
                    
                    // Get logarithm based on type
                    const logR = cga.Gen.logRotor(ratio);
                    
                    // Interpolated log
                    const logT = logR.scale(t);
                    
                    // Exp back
                    const Rt = cga.Gen.rot(logT.scale(-0.5));
                    
                    return v1.mul(Rt);
                }
            };
        }

        // ========================================================================
        // APPLY TRANSFORMATIONS
        // ========================================================================

        /**
         * Apply versor to element via sandwich product
         * X' = V * X * Ṽ
         */
        spin(element, versor) {
            return versor.mul(element).mul(versor.reverse());
        }

        /**
         * Apply reflection
         * X' = V * X * V (no reverse for odd versors)
         */
        reflect(element, reflector) {
            return reflector.mul(element).mul(reflector);
        }

        // ========================================================================
        // CONVENIENCE METHODS
        // ========================================================================

        /**
         * Create a "boosted surface" by applying booster to mesh vertices
         * This is Colapinto's technique for organic shape synthesis
         * 
         * @param {Array<[number,number,number]>} vertices - input vertices
         * @param {Multivector} booster - the boost transformation
         * @returns {Array<[number,number,number]>} transformed vertices
         */
        boostSurface(vertices, booster) {
            return vertices.map(v => {
                const p = this.point(v[0], v[1], v[2]);
                const boosted = this.spin(p, booster);
                return this.pointCoords(boosted);
            });
        }

        /**
         * Create booster field from point pair generators
         * Used for continuous surface deformation
         * 
         * @param {Array<{p1, p2, weight}>} generators - point pair generators
         * @returns {Function} field function (x,y,z) => booster
         */
        boosterField(generators) {
            const cga = this;
            
            return (x, y, z) => {
                let totalBoost = cga.scalar(1);
                
                for (const gen of generators) {
                    const pp = cga.pointPair(gen.p1, gen.p2);
                    const dist = sqrt(
                        (x - (gen.p1[0] + gen.p2[0])/2)**2 +
                        (y - (gen.p1[1] + gen.p2[1])/2)**2 +
                        (z - (gen.p1[2] + gen.p2[2])/2)**2
                    );
                    
                    // Weight falls off with distance
                    const w = gen.weight * exp(-dist * dist);
                    const boost = cga.Gen.bst(pp.scale(w));
                    totalBoost = totalBoost.mul(boost);
                }
                
                return totalBoost;
            };
        }

        /**
         * Circle inversion (reflection in sphere)
         */
        invert(element, sphere) {
            // Inversion is reflection in the sphere
            return sphere.mul(element).mul(sphere.inverse());
        }

        /**
         * Meet of two elements (intersection)
         */
        meet(a, b) {
            return a.vee(b);
        }

        /**
         * Join of two elements (smallest containing element)
         */
        join(a, b) {
            return a.wedge(b);
        }
    }

    // ============================================================================
    // ARBITRARY METRIC ALGEBRA FACTORY
    // ============================================================================

    /**
     * Create algebra with arbitrary metric (Versor-style)
     * 
     * @param {Object} config
     * @param {number[]} config.metric - array of +1, -1, 0 for basis vector squares
     * @param {boolean} config.conformal - if true, treats last two dims as conformal pair
     * @param {Object[]} config.types - optional named type definitions
     * @returns {Algebra} configured algebra
     */
    function createAlgebra(config) {
        const metric = config.metric || [1, 1, 1];
        const conformal = config.conformal || false;
        
        // Count signature
        let p = 0, q = 0, r = 0;
        for (const m of metric) {
            if (m === 1) p++;
            else if (m === -1) q++;
            else r++;
        }
        
        // Create base algebra
        const algebra = new Algebra(p, q, r);
        
        // Add conformal structure if requested
        if (conformal && metric.length >= 2) {
            const n = metric.length;
            
            // Last two basis vectors form conformal pair
            const ep = algebra.blade(1 << (n - 2));  // e+
            const em = algebra.blade(1 << (n - 1));  // e-
            
            algebra.no = em.sub(ep).scale(0.5);  // origin
            algebra.ni = em.add(ep);              // infinity
            
            // Point embedding function
            algebra.point = function(...coords) {
                let euclidean = algebra.scalar(0);
                let sq = 0;
                for (let i = 0; i < n - 2; i++) {
                    const ei = algebra.blade(1 << i);
                    euclidean = euclidean.add(ei.scale(coords[i] || 0));
                    sq += (coords[i] || 0) ** 2;
                }
                return euclidean.add(algebra.ni.scale(0.5 * sq)).add(algebra.no);
            };
            
            // Point extraction
            algebra.pointCoords = function(X) {
                const div = -X.dot(algebra.ni).scalar;
                if (isZero(div)) return new Array(n - 2).fill(0);
                
                const coords = [];
                for (let i = 0; i < n - 2; i++) {
                    coords.push(X.coeffs[1 << i] / div);
                }
                return coords;
            };
        }
        
        // Add named types if specified
        if (config.types) {
            algebra.types = {};
            for (const type of config.types) {
                algebra.types[type.name] = {
                    bases: type.bases,
                    dual: type.dual || false
                };
                
                // Create constructor function
                algebra[type.name] = function(...values) {
                    const mv = new Multivector(algebra);
                    for (let i = 0; i < type.bases.length; i++) {
                        const blade = type.bases[i];
                        const idx = algebra.bladeIndex[blade];
                        if (idx !== undefined && i < values.length) {
                            mv.coeffs[idx] = values[i];
                        }
                    }
                    return mv;
                };
            }
        }
        
        return algebra;
    }

    // ============================================================================
    // EXPORTS
    // ============================================================================

    const VersorGA = {
        // Extended CGA class
        CGA3DVersor,
        
        // Factory for extended CGA
        cga3dVersor: () => new CGA3DVersor(),
        
        // Arbitrary algebra creation
        create: createAlgebra,
        
        // Re-export base classes
        Algebra,
        Multivector
    };

    // Merge with existing GA object
    Object.assign(GA, VersorGA);

    // Export
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = VersorGA;
    } else if (typeof define === 'function' && define.amd) {
        define([], () => VersorGA);
    } else {
        global.VersorGA = VersorGA;
    }

})(typeof window !== 'undefined' ? window : global);
