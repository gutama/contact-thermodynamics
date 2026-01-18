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

    // ============================================================================
    // MESH / FTGC (Discrete Geometric Calculus on Triangle Meshes)
    // ============================================================================

    /**
     * Triangle mesh with typed array storage and precomputed topology.
     */
    export class TriangleMesh {
        /** Vertex positions (nVertices * 3) */
        readonly vertices: Float64Array;
        /** Triangle indices (nFaces * 3) */
        readonly faces: Uint32Array;
        /** Number of vertices */
        readonly nVertices: number;
        /** Number of faces */
        readonly nFaces: number;
        /** Number of edges */
        readonly nEdges: number;
        /** Edge vertex pairs (nEdges * 2) */
        readonly edges: Uint32Array;
        /** Boundary vertex mask (1 = boundary) */
        readonly boundaryVertices: Uint8Array;
        /** Boundary edge mask (1 = boundary) */
        readonly boundaryEdges: Uint8Array;

        constructor(
            vertices: Float64Array,
            faces: Uint32Array,
            opts?: { buildTopology?: boolean }
        );

        /** Build edge list and adjacency structures */
        buildTopology(): void;

        /** Get edge index for vertex pair */
        getEdgeIndex(v0: number, v1: number): number;

        /** Get vertex position as [x, y, z] */
        getVertex(idx: number): [number, number, number];

        /** Compute triangle areas (cached) */
        faceAreas(): Float64Array;

        /** Compute edge lengths (cached) */
        edgeLengths(): Float64Array;

        /** Compute cotan weights per edge (cached) */
        cotanWeights(): Float64Array;

        /** Compute mixed Voronoi dual areas per vertex (cached) */
        vertexDualAreas(): Float64Array;

        /** Invalidate geometry cache after modifying vertices */
        invalidateGeometryCache(): void;

        /** Create a regular grid mesh */
        static createGrid(nx: number, ny: number, Lx?: number, Ly?: number): TriangleMesh;

        /** Create a triangulated disk mesh */
        static createDisk(nRadial: number, nAngular: number, radius?: number): TriangleMesh;

        /** Create an icosahedron (closed mesh) */
        static createIcosahedron(radius?: number): TriangleMesh;
    }

    /**
     * Simple sparse matrix in CSR format.
     */
    export class SparseMatrix {
        readonly rows: number;
        readonly cols: number;

        constructor(rows: number, cols: number);

        /** Set entry (accumulates) */
        set(row: number, col: number, value: number): void;

        /** Get entry */
        get(row: number, col: number): number;

        /** Build CSR format for efficient matvec */
        buildCSR(): void;

        /** Matrix-vector product: y = A * x */
        matvec(x: Float64Array): Float64Array;

        /** Transpose matrix-vector product: y = Aᵀ * x */
        matvecT(x: Float64Array): Float64Array;

        /** Get diagonal entries */
        diagonal(): Float64Array;
    }

    /** Create a diagonal matrix from an array */
    export function diagMatrix(diag: Float64Array): SparseMatrix;

    /**
     * Multivector field on a mesh with staggered storage.
     * Grade 0 → Vertices, Grade 1 → Edges, Grade 2 → Faces
     */
    export class MeshMultivectorField {
        readonly mesh: TriangleMesh;
        readonly grade0: Float64Array;
        readonly grade1: Float64Array;
        readonly grade2: Float64Array;

        constructor(
            mesh: TriangleMesh,
            grade0?: Float64Array | null,
            grade1?: Float64Array | null,
            grade2?: Float64Array | null
        );

        add(other: MeshMultivectorField): MeshMultivectorField;
        sub(other: MeshMultivectorField): MeshMultivectorField;
        scale(s: number): MeshMultivectorField;
        gradeSelect(k: number): MeshMultivectorField;
        clone(): MeshMultivectorField;

        static scalarField(mesh: TriangleMesh, values: Float64Array): MeshMultivectorField;
        static vectorField(mesh: TriangleMesh, values: Float64Array): MeshMultivectorField;
        static bivectorField(mesh: TriangleMesh, values: Float64Array): MeshMultivectorField;
    }

    /**
     * Discrete geometric derivative ∇ on triangle meshes (FTGC).
     * 
     * Implements the split differential operator: ∇ = Σᵢ eⁱ ⊗ ∂/∂xᵢ
     */
    export class MeshGeometricDerivative {
        readonly mesh: TriangleMesh;

        constructor(mesh: TriangleMesh);

        /** Apply outer derivative ∇∧ (grade-raising) */
        wedge(F: MeshMultivectorField): MeshMultivectorField;

        /** Apply inner derivative ∇· (grade-lowering) */
        inner(F: MeshMultivectorField): MeshMultivectorField;

        /** Apply full geometric derivative ∇ = ∇· + ∇∧ */
        apply(F: MeshMultivectorField): MeshMultivectorField;

        /** Gradient of scalar field */
        grad(f: Float64Array): Float64Array;

        /** Divergence of vector field */
        div(V: Float64Array): Float64Array;

        /** Curl of vector field */
        curl(V: Float64Array): Float64Array;

        /** Laplacian of scalar field (cotan Laplacian) */
        laplacian(f: Float64Array, opts?: {
            dirichletMask?: Uint8Array;
            dirichletValues?: Float64Array;
        }): Float64Array;

        /** Build Laplacian as sparse matrix */
        laplacianMatrix(): SparseMatrix;

        /** Build weak Laplacian (without mass matrix inverse) */
        laplacianWeak(): SparseMatrix;

        /** Verify ∇∧(∇∧f) = 0 */
        verifyCurlGradZero(f: Float64Array): number;

        /** Verify ∇²(constant) ≈ 0 */
        verifyLaplacianConstantZero(c?: number): number;
    }

    /** Apply Dirichlet boundary conditions to a field */
    export function applyDirichlet(
        f: Float64Array,
        mask: Uint8Array,
        values?: Float64Array
    ): void;

    /** Create Dirichlet mask from mesh boundary vertices */
    export function boundaryDirichletMask(mesh: TriangleMesh): Uint8Array;

    /**
     * Leapfrog time-stepping solver using FTGC.
     */
    export class LeapfrogGCMesh {
        readonly mesh: TriangleMesh;
        readonly nabla: MeshGeometricDerivative;

        constructor(mesh: TriangleMesh);

        /** Estimate CFL time step for wave equation */
        estimateCFL(c?: number): number;

        /** Estimate CFL time step for heat equation */
        estimateCFLHeat(alpha?: number): number;

        /** One leapfrog step for wave equation */
        waveStep(
            uPrev: Float64Array,
            uCurr: Float64Array,
            dt: number,
            c?: number,
            opts?: {
                dirichletMask?: Uint8Array;
                dirichletValues?: Float64Array;
            }
        ): Float64Array;

        /** Simulate wave equation */
        waveSimulate(
            u0: Float64Array,
            v0: Float64Array,
            dt: number,
            nSteps: number,
            c?: number,
            opts?: {
                dirichletMask?: Uint8Array;
                dirichletValues?: Float64Array;
                callback?: (step: number, u: Float64Array) => void;
            }
        ): Float64Array;

        /** Explicit Euler step for heat equation */
        heatStepExplicit(
            u: Float64Array,
            dt: number,
            alpha?: number,
            opts?: {
                dirichletMask?: Uint8Array;
                dirichletValues?: Float64Array;
            }
        ): Float64Array;

        /** Implicit Euler step for heat equation */
        heatStepImplicit(
            u: Float64Array,
            dt: number,
            alpha?: number,
            opts?: {
                dirichletMask?: Uint8Array;
                dirichletValues?: Float64Array;
                maxIter?: number;
                tol?: number;
            }
        ): Float64Array;

        /** Simulate heat equation */
        heatSimulate(
            u0: Float64Array,
            dt: number,
            nSteps: number,
            alpha?: number,
            opts?: {
                implicit?: boolean;
                dirichletMask?: Uint8Array;
                dirichletValues?: Float64Array;
                callback?: (step: number, u: Float64Array) => void;
            }
        ): Float64Array;

        /** Yee-style Maxwell step */
        maxwellStep(
            E: Float64Array,
            B: Float64Array,
            dt: number,
            c?: number
        ): { E: Float64Array; B: Float64Array };

        /** Compute L2 energy of scalar field */
        energy(u: Float64Array): number;

        /** Compute total mass of scalar field */
        mass(u: Float64Array): number;

        /** Compute variance of scalar field */
        variance(u: Float64Array): number;
    }
}
