/**
 * Contact Thermodynamics - TypeScript Definitions
 * 
 * Contact Geometry for Extended Thermodynamics using 1-Jet Bundles
 */

declare module 'contact-thermodynamics' {
    
    /** Numerical tolerance constant */
    export const EPSILON: number;

    /**
     * A point on a contact manifold
     */
    export class ContactPoint {
        /** The parent manifold */
        readonly manifold: ContactManifold;
        /** Coordinate values */
        readonly coords: Record<string, number>;

        constructor(manifold: ContactManifold, coords?: Record<string, number>);

        /** Get coordinate value */
        get(coord: string): number;

        /** Set coordinate value (returns self for chaining) */
        set(coord: string, value: number): this;

        /** Create a copy of this point */
        clone(): ContactPoint;

        /** Add scaled tangent vector */
        add(tangent: Record<string, number>, dt?: number): ContactPoint;

        /** String representation */
        toString(): string;
    }

    /**
     * Base class for contact manifolds J¹(Q)
     */
    export class ContactManifold {
        /** Base configuration coordinate names */
        readonly baseCoords: string[];
        /** Conjugate momenta coordinate names */
        readonly momentaCoords: string[];
        /** Fiber coordinate name */
        readonly fiberCoord: string;
        /** Dimension of base manifold Q */
        readonly n: number;
        /** Total dimension (2n+1) */
        readonly dim: number;

        constructor(baseCoords: string[], momentaCoords: string[], fiberCoord?: string);

        /** All coordinate names in canonical order */
        readonly allCoords: string[];

        /** Origin point (all coordinates zero) */
        readonly origin: ContactPoint;

        /** Create a point on the manifold */
        point(coords: Record<string, number>): ContactPoint;

        /** Symbolic representation of contact form */
        contactFormSymbolic(): string;

        /** Evaluate contact form at point with tangent vector */
        evaluateContactForm(pt: ContactPoint, tangent: Record<string, number>): number;

        /** Verify non-degeneracy α ∧ (dα)^n ≠ 0 */
        verifyContactCondition(pt: ContactPoint): number;

        /** Get Reeb vector field at a point */
        reebField(pt: ContactPoint): Record<string, number>;

        toString(): string;
    }

    /**
     * Grand Contact Manifold M₁₃ = J¹(Q₆)
     * 
     * Base: (q¹, q², q³, t, ℓ, S)
     * Momenta: (k₁, k₂, k₃, ω, Δ, T)
     */
    export class GrandContactManifold extends ContactManifold {
        /** Physical interpretations of coordinates */
        readonly physicalInterpretations: Record<string, string>;

        constructor();

        /** Create point with physical coordinate names */
        physicalPoint(
            q1: number, q2: number, q3: number,
            t: number, ell: number, S: number,
            k1: number, k2: number, k3: number,
            omega: number, Delta: number, T: number,
            A?: number
        ): ContactPoint;

        /** Extract spatial position (q¹, q², q³) */
        spatialPosition(pt: ContactPoint): [number, number, number];

        /** Extract wave vector (k₁, k₂, k₃) */
        waveVector(pt: ContactPoint): [number, number, number];
    }

    /**
     * Holographic Contact Manifold M₇ = J¹(Q₃)
     * 
     * Base: (t, ℓ, S)
     * Momenta: (ω, Δ, T)
     * Emergent: q^i(t, ℓ, S)
     */
    export class HolographicContactManifold extends ContactManifold {
        /** Names of emergent spatial fields */
        readonly emergentFields: string[];
        /** Physical interpretations */
        readonly physicalInterpretations: Record<string, string>;

        constructor();

        /** Create holographic point */
        holographicPoint(
            t: number, ell: number, S: number,
            omega: number, Delta: number, T: number,
            A?: number
        ): ContactPoint;

        /** Define emergent spatial configuration */
        createEmergentSpace(
            pt: ContactPoint,
            fieldFunc: (t: number, ell: number, S: number) => [number, number, number]
        ): { q1: number; q2: number; q3: number };

        /** Extend point with emergent fields */
        extendedPoint(
            pt: ContactPoint,
            emergent: { q1: number; q2: number; q3: number }
        ): Record<string, number>;
    }

    /**
     * Gauge-Extended Manifold M₁₅
     * 
     * Adds (φ, I) gauge pair to Grand model
     */
    export class GaugeExtendedManifold extends GrandContactManifold {
        constructor();
    }

    /**
     * Hamiltonian dynamics on contact manifolds
     */
    export class ContactHamiltonian {
        /** The contact manifold */
        readonly manifold: ContactManifold;
        /** Hamiltonian function */
        readonly H: (coords: Record<string, number>) => number;

        constructor(
            manifold: ContactManifold,
            H: (coords: Record<string, number>) => number,
            dH?: (coords: Record<string, number>) => Record<string, number>
        );

        /** Evaluate Hamiltonian at a point */
        evaluate(pt: ContactPoint): number;

        /** Compute numerical gradient */
        gradient(pt: ContactPoint, h?: number): Record<string, number>;

        /** Get Reeb component RH = ∂H/∂u */
        reebComponent(pt: ContactPoint): number;

        /** Compute contact Hamiltonian vector field X_H */
        vectorField(pt: ContactPoint): Record<string, number>;

        /** Integrate dynamics using RK4 */
        flow(pt: ContactPoint, dt: number, steps?: number): ContactPoint[];

        /** Get Hamiltonian values along trajectory */
        hamiltonianEvolution(trajectory: ContactPoint[]): number[];
    }

    /**
     * Legendrian submanifold defined by generating function
     */
    export class LegendrianSubmanifold {
        /** The ambient contact manifold */
        readonly manifold: ContactManifold;
        /** Generating function A(x) */
        readonly A: (x: Record<string, number>) => number;

        constructor(
            manifold: ContactManifold,
            generatingFunc: (x: Record<string, number>) => number,
            dA?: (x: Record<string, number>) => Record<string, number>
        );

        /** Compute gradient of generating function */
        gradient(x: Record<string, number>, h?: number): Record<string, number>;

        /** Lift base coordinates to contact manifold */
        lift(x: Record<string, number>): ContactPoint;

        /** Verify Legendrian condition α|_L = 0 */
        verifyLegendrianCondition(x: Record<string, number>): boolean;

        /** Compute Hamilton-Jacobi residual H(x, A(x), ∂A(x)) */
        hamiltonJacobiResidual(
            x: Record<string, number>,
            hamiltonian: ContactHamiltonian
        ): number;

        /** Sample points on the Legendrian */
        sample(
            sampler: () => Record<string, number>,
            n: number
        ): ContactPoint[];
    }

    /**
     * Specialized Hamiltonians for thermodynamic systems
     */
    export class ThermodynamicHamiltonian extends ContactHamiltonian {
        /** Physical parameters */
        readonly params: Record<string, number>;

        constructor(
            manifold: ContactManifold,
            params?: {
                mass?: number;
                potential?: (coords: Record<string, number>) => number;
                thermalCoupling?: number;
                scalingDimension?: number;
            }
        );

        /** Create dispersion relation Hamiltonian H = ω - c|k| or H = ω - √(c²|k|² + m²) */
        static dispersionRelation(
            manifold: ContactManifold,
            c?: number,
            mass?: number
        ): ContactHamiltonian;

        /** Create equation of state Hamiltonian */
        static equationOfState(
            manifold: ContactManifold,
            type?: 'ideal' | 'van_der_waals'
        ): ContactHamiltonian;
    }

    /**
     * Spacetime metric g_μν(x)
     */
    export class SpacetimeMetric {
        constructor(
            metricFunc: (x: number[]) => number[][],
            inverseFunc?: (x: number[]) => number[][]
        );

        /** Get covariant metric g_μν at point */
        covariant(x: number[]): number[][];

        /** Get contravariant metric g^μν at point */
        contravariant(x: number[]): number[][];

        /** Minkowski metric η_μν = diag(+1, -1, -1, -1) */
        static minkowski(): SpacetimeMetric;

        /** Schwarzschild black hole metric */
        static schwarzschild(M?: number): SpacetimeMetric;

        /** FLRW cosmological metric */
        static flrw(a: (t: number) => number, k?: number): SpacetimeMetric;
    }

    /**
     * Relativistic Hamiltonian for geodesic motion
     */
    export class RelativisticHamiltonian {
        /** Spacetime metric */
        readonly metric: SpacetimeMetric;
        /** Particle mass */
        readonly m: number;
        /** Gauge potential A_μ */
        readonly A: (x: number[]) => number[];
        /** Electric charge */
        readonly q: number;

        constructor(
            metric: SpacetimeMetric,
            mass?: number,
            gaugePotential?: (x: number[]) => number[],
            charge?: number
        );

        /** Evaluate mass-shell Hamiltonian */
        evaluate(x: number[], p: number[]): number;

        /** Compute Hamilton's equations */
        geodesicEquations(x: number[], p: number[]): {
            xDot: number[];
            pDot: number[];
        };

        /** Integrate geodesic with RK4 */
        integrateGeodesic(
            x0: number[],
            p0: number[],
            dtau: number,
            steps: number
        ): Array<{
            x: number[];
            p: number[];
            tau: number;
        }>;

        /** Hamilton-Jacobi residual */
        hjResidual(
            x: number[],
            actionFunc: (x: number[]) => number,
            gradA: (x: number[]) => number[]
        ): number;
    }

    /**
     * Symbolic differential form
     */
    export class DifferentialForm {
        readonly degree: number;
        readonly terms: Array<{
            coeff: number;
            basis: string[];
        }>;

        constructor(degree: number, terms?: Array<{ coeff: number; basis: string[] }>);

        /** Create 1-form dx^a */
        static oneForm(coord: string): DifferentialForm;

        /** Wedge product */
        wedge(other: DifferentialForm): DifferentialForm;

        toString(): string;
    }

    /**
     * Get theory summary table
     */
    export function summaryTable(): {
        columns: string[];
        rows: string[][];
    };

    // Factory functions
    export function grandManifold(): GrandContactManifold;
    export function holographicManifold(): HolographicContactManifold;
    export function gaugeExtended(): GaugeExtendedManifold;
}
