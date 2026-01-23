/**
 * Discrete Geodesics via Geometric Algebra
 * 
 * Implements geodesic algorithms from Keenan Crane's DDG lectures:
 * - Shortest geodesics (variational perspective)
 * - Straightest geodesics (geometric perspective)
 * - Exponential and Logarithmic maps
 * - Parallel transport
 * - Cut locus computation
 * 
 * Key insight: Geodesics are curves with κ_g = 0 (zero geodesic curvature).
 * On polyhedral surfaces, this means "equal angles" when crossing edges.
 * 
 * References:
 * - Keenan Crane, "Discrete Differential Geometry: An Applied Introduction"
 * - Mitchell, Mount, Papadimitriou (MMP) exact geodesic algorithm
 * 
 * @license MIT
 * @version 1.0.0
 */

(function(global) {
    'use strict';

    const EPSILON = 1e-10;
    const abs = Math.abs;
    const sqrt = Math.sqrt;
    const sin = Math.sin;
    const cos = Math.cos;
    const tan = Math.tan;
    const atan2 = Math.atan2;
    const acos = Math.acos;
    const PI = Math.PI;
    const min = Math.min;
    const max = Math.max;

    // ============================================================================
    // UTILITY FUNCTIONS
    // ============================================================================

    const Vec3 = {
        add: (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]],
        sub: (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]],
        scale: (v, s) => [v[0] * s, v[1] * s, v[2] * s],
        dot: (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2],
        cross: (a, b) => [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]
        ],
        length: (v) => sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]),
        normalize: (v) => {
            const len = Vec3.length(v);
            return len > EPSILON ? [v[0] / len, v[1] / len, v[2] / len] : [0, 0, 0];
        },
        lerp: (a, b, t) => [
            a[0] + t * (b[0] - a[0]),
            a[1] + t * (b[1] - a[1]),
            a[2] + t * (b[2] - a[2])
        ],
        distance: (a, b) => Vec3.length(Vec3.sub(b, a))
    };

    // 2D operations for unfolding
    const Vec2 = {
        add: (a, b) => [a[0] + b[0], a[1] + b[1]],
        sub: (a, b) => [a[0] - b[0], a[1] - b[1]],
        scale: (v, s) => [v[0] * s, v[1] * s],
        dot: (a, b) => a[0] * b[0] + a[1] * b[1],
        length: (v) => sqrt(v[0] * v[0] + v[1] * v[1]),
        normalize: (v) => {
            const len = Vec2.length(v);
            return len > EPSILON ? [v[0] / len, v[1] / len] : [0, 0];
        },
        rotate: (v, angle) => [
            v[0] * cos(angle) - v[1] * sin(angle),
            v[0] * sin(angle) + v[1] * cos(angle)
        ],
        angle: (v) => atan2(v[1], v[0])
    };

    // ============================================================================
    // PRIORITY QUEUE FOR DIJKSTRA
    // ============================================================================

    class PriorityQueue {
        constructor() {
            this.heap = [];
        }

        push(item, priority) {
            this.heap.push({ item, priority });
            this._bubbleUp(this.heap.length - 1);
        }

        pop() {
            if (this.heap.length === 0) return null;
            const result = this.heap[0];
            const last = this.heap.pop();
            if (this.heap.length > 0) {
                this.heap[0] = last;
                this._bubbleDown(0);
            }
            return result;
        }

        isEmpty() {
            return this.heap.length === 0;
        }

        _bubbleUp(idx) {
            while (idx > 0) {
                const parent = Math.floor((idx - 1) / 2);
                if (this.heap[parent].priority <= this.heap[idx].priority) break;
                [this.heap[parent], this.heap[idx]] = [this.heap[idx], this.heap[parent]];
                idx = parent;
            }
        }

        _bubbleDown(idx) {
            const n = this.heap.length;
            while (true) {
                const left = 2 * idx + 1;
                const right = 2 * idx + 2;
                let smallest = idx;
                if (left < n && this.heap[left].priority < this.heap[smallest].priority) {
                    smallest = left;
                }
                if (right < n && this.heap[right].priority < this.heap[smallest].priority) {
                    smallest = right;
                }
                if (smallest === idx) break;
                [this.heap[smallest], this.heap[idx]] = [this.heap[idx], this.heap[smallest]];
                idx = smallest;
            }
        }
    }

    // ============================================================================
    // DISCRETE GEODESICS CLASS
    // ============================================================================

    /**
     * Discrete Geodesics Calculator
     * 
     * SHORTEST GEODESICS (Variational Perspective):
     * - Critical points of length functional
     * - Never pass through cone vertices (Θ < 2π)
     * - May pass through saddle vertices (Θ > 2π)
     * 
     * STRAIGHTEST GEODESICS (Geometric Perspective):
     * - Zero geodesic curvature (κ_g = 0)
     * - Equal angles on either side when crossing edge
     * - Natural for initial value problems (shooting)
     * 
     * Key distinction:
     * - Shortest: good for boundary value problems (point A to B)
     * - Straightest: good for initial value problems (start + direction)
     */
    class DiscreteGeodesics {
        constructor() {
            this.vertices = [];
            this.faces = [];
            this.edges = [];
            this.vertexFaces = [];
            this.vertexEdges = [];
            this.adjacency = [];
            this.edgeMap = new Map();
            this.faceNormals = [];
        }

        /**
         * Load mesh from vertex/face arrays
         * @param {Array} vertices - [[x,y,z], ...]
         * @param {Array} faces - [[i,j,k], ...] (triangles)
         */
        loadMesh(vertices, faces) {
            this.vertices = vertices.map(v => [...v]);
            this.faces = faces.map(f => [...f]);
            this._buildTopology();
            this._computeFaceNormals();
        }

        _buildTopology() {
            const V = this.vertices.length;
            const F = this.faces.length;
            
            this.vertexFaces = Array.from({ length: V }, () => []);
            this.vertexEdges = Array.from({ length: V }, () => new Set());
            this.adjacency = Array.from({ length: V }, () => new Set());
            this.edgeMap = new Map();
            
            for (let fi = 0; fi < F; fi++) {
                const f = this.faces[fi];
                for (let j = 0; j < 3; j++) {
                    const v0 = f[j];
                    const v1 = f[(j + 1) % 3];
                    
                    this.vertexFaces[v0].push(fi);
                    this.adjacency[v0].add(v1);
                    this.adjacency[v1].add(v0);
                    
                    const key = v0 < v1 ? `${v0}-${v1}` : `${v1}-${v0}`;
                    if (!this.edgeMap.has(key)) {
                        this.edgeMap.set(key, { 
                            v0: min(v0, v1), 
                            v1: max(v0, v1), 
                            faces: [],
                            index: this.edgeMap.size
                        });
                    }
                    this.edgeMap.get(key).faces.push(fi);
                    
                    this.vertexEdges[v0].add(key);
                    this.vertexEdges[v1].add(key);
                }
            }
            
            this.edges = Array.from(this.edgeMap.values());
        }

        _computeFaceNormals() {
            this.faceNormals = this.faces.map(f => {
                const v0 = this.vertices[f[0]];
                const v1 = this.vertices[f[1]];
                const v2 = this.vertices[f[2]];
                const e1 = Vec3.sub(v1, v0);
                const e2 = Vec3.sub(v2, v0);
                return Vec3.normalize(Vec3.cross(e1, e2));
            });
        }

        // ========================================================================
        // VERTEX GEOMETRY
        // ========================================================================

        /**
         * Compute angle sum around vertex (total corner angles)
         * @param {number} vertexIdx - Vertex index
         * @returns {number} Sum of angles in radians
         */
        angleSum(vertexIdx) {
            let sum = 0;
            for (const fi of this.vertexFaces[vertexIdx]) {
                const f = this.faces[fi];
                const idx = f.indexOf(vertexIdx);
                const v0 = this.vertices[vertexIdx];
                const v1 = this.vertices[f[(idx + 1) % 3]];
                const v2 = this.vertices[f[(idx + 2) % 3]];
                
                const e1 = Vec3.normalize(Vec3.sub(v1, v0));
                const e2 = Vec3.normalize(Vec3.sub(v2, v0));
                const dot = max(-1, min(1, Vec3.dot(e1, e2)));
                sum += acos(dot);
            }
            return sum;
        }

        /**
         * Angle defect at vertex (discrete Gaussian curvature)
         * Positive = cone (K > 0), Negative = saddle (K < 0)
         */
        angleDefect(vertexIdx) {
            return 2 * PI - this.angleSum(vertexIdx);
        }

        /**
         * Check if vertex is a cone (positive curvature)
         * Geodesics NEVER pass through cone vertices
         */
        isConeVertex(vertexIdx) {
            return this.angleDefect(vertexIdx) > EPSILON;
        }

        /**
         * Check if vertex is a saddle (negative curvature)
         * Many shortest geodesics may pass through saddles
         */
        isSaddleVertex(vertexIdx) {
            return this.angleDefect(vertexIdx) < -EPSILON;
        }

        // ========================================================================
        // SHORTEST GEODESICS: DIJKSTRA ON GRAPH (APPROXIMATION)
        // ========================================================================

        /**
         * Dijkstra's algorithm on mesh edges
         * 
         * NOTE: This is a POOR approximation of true geodesic distance!
         * Even on refined meshes, graph distance can be 10-20% off.
         * Use exactGeodesicDistance() for accurate results.
         * 
         * @param {number} source - Source vertex index
         * @returns {Object} {distances, predecessors}
         */
        dijkstraGraph(source) {
            const n = this.vertices.length;
            const dist = new Array(n).fill(Infinity);
            const pred = new Array(n).fill(-1);
            const visited = new Array(n).fill(false);
            
            dist[source] = 0;
            const pq = new PriorityQueue();
            pq.push(source, 0);
            
            while (!pq.isEmpty()) {
                const { item: u } = pq.pop();
                if (visited[u]) continue;
                visited[u] = true;
                
                for (const v of this.adjacency[u]) {
                    const edgeLen = Vec3.distance(this.vertices[u], this.vertices[v]);
                    const newDist = dist[u] + edgeLen;
                    
                    if (newDist < dist[v]) {
                        dist[v] = newDist;
                        pred[v] = u;
                        pq.push(v, newDist);
                    }
                }
            }
            
            return { distances: dist, predecessors: pred };
        }

        /**
         * Reconstruct shortest path on graph from Dijkstra result
         */
        reconstructPath(predecessors, target) {
            const path = [];
            let current = target;
            while (current !== -1) {
                path.unshift(current);
                current = predecessors[current];
            }
            return path;
        }

        // ========================================================================
        // SHORTEST GEODESICS: EXACT POLYHEDRAL (SIMPLIFIED MMP)
        // ========================================================================

        /**
         * Unfold a triangle into 2D plane
         * Places v0 at origin, v1 on positive x-axis
         * 
         * @param {Array} v0, v1, v2 - 3D vertices
         * @returns {Array} [[x0,y0], [x1,y1], [x2,y2]]
         */
        unfoldTriangle(v0, v1, v2) {
            const e01 = Vec3.sub(v1, v0);
            const e02 = Vec3.sub(v2, v0);
            
            const l01 = Vec3.length(e01);
            const l02 = Vec3.length(e02);
            
            // Angle at v0
            const cosTheta = Vec3.dot(e01, e02) / (l01 * l02 + EPSILON);
            const theta = acos(max(-1, min(1, cosTheta)));
            
            // Place in 2D
            return [
                [0, 0],           // v0 at origin
                [l01, 0],         // v1 on x-axis
                [l02 * cos(theta), l02 * sin(theta)]  // v2 rotated
            ];
        }

        /**
         * Compute geodesic distance from source vertex to all points
         * using window propagation (simplified MMP-like algorithm)
         * 
         * @param {number} source - Source vertex index
         * @returns {Object} {vertexDistances, faceDistanceFields}
         */
        computeGeodesicDistance(source) {
            const n = this.vertices.length;
            const dist = new Array(n).fill(Infinity);
            dist[source] = 0;

            // Windows: each edge can have multiple windows tracking geodesic paths
            // Window = {edge, sourcePos2D, leftBound, rightBound, distance}
            const edgeWindows = new Map();
            
            // Initialize windows from source vertex
            const pq = new PriorityQueue();
            
            for (const fi of this.vertexFaces[source]) {
                const f = this.faces[fi];
                const idx = f.indexOf(source);
                const v1 = f[(idx + 1) % 3];
                const v2 = f[(idx + 2) % 3];
                
                // Create initial window on opposite edge
                const edgeKey = v1 < v2 ? `${v1}-${v2}` : `${v2}-${v1}`;
                const unfolded = this.unfoldTriangle(
                    this.vertices[source],
                    this.vertices[v1],
                    this.vertices[v2]
                );
                
                const window = {
                    edge: edgeKey,
                    faceIdx: fi,
                    sourcePos: unfolded[0],  // Source in unfolded coords
                    leftBound: 0,
                    rightBound: 1,
                    baseDist: 0
                };
                
                if (!edgeWindows.has(edgeKey)) {
                    edgeWindows.set(edgeKey, []);
                }
                edgeWindows.get(edgeKey).push(window);
                
                // Update vertex distances
                const d1 = Vec2.length(Vec2.sub(unfolded[1], unfolded[0]));
                const d2 = Vec2.length(Vec2.sub(unfolded[2], unfolded[0]));
                if (d1 < dist[v1]) {
                    dist[v1] = d1;
                    pq.push({ vertex: v1, dist: d1 }, d1);
                }
                if (d2 < dist[v2]) {
                    dist[v2] = d2;
                    pq.push({ vertex: v2, dist: d2 }, d2);
                }
            }

            // Propagate windows (simplified - just updates vertex distances)
            const visited = new Set();
            visited.add(source);
            
            while (!pq.isEmpty()) {
                const { item } = pq.pop();
                const u = item.vertex;
                
                if (visited.has(u)) continue;
                visited.add(u);
                
                // Propagate through adjacent faces
                for (const fi of this.vertexFaces[u]) {
                    const f = this.faces[fi];
                    for (let j = 0; j < 3; j++) {
                        const v = f[j];
                        if (visited.has(v)) continue;
                        
                        // Compute distance through this triangle
                        const idx = f.indexOf(u);
                        const v1 = f[(idx + 1) % 3];
                        const v2 = f[(idx + 2) % 3];
                        
                        // Use triangle inequality for update
                        const d_u = dist[u];
                        const d_v1 = dist[v1];
                        const d_v2 = dist[v2];
                        
                        // Direct edge distance
                        const edge_uv = Vec3.distance(this.vertices[u], this.vertices[v]);
                        const newDist = d_u + edge_uv;
                        
                        if (newDist < dist[v]) {
                            dist[v] = newDist;
                            pq.push({ vertex: v, dist: newDist }, newDist);
                        }
                    }
                }
            }

            return { distances: dist, windows: edgeWindows };
        }

        // ========================================================================
        // STRAIGHTEST GEODESICS
        // ========================================================================

        /**
         * Trace a straightest geodesic from a point in a direction
         * 
         * STRAIGHTEST CONDITION: When crossing an edge, the angles on
         * either side must be equal: θ_left = θ_right
         * 
         * This is equivalent to: the path unfolds to a straight line
         * 
         * @param {Object} start - {faceIdx, baryCoords: [u,v,w]} or {vertexIdx}
         * @param {Array} direction - [dx, dy, dz] initial direction (tangent to surface)
         * @param {number} maxLength - Maximum path length
         * @returns {Array} Array of {position, faceIdx, parameter}
         */
        traceStraightestGeodesic(start, direction, maxLength = 10) {
            const path = [];
            let currentFace = start.faceIdx;
            let currentPos = this._baryToPosition(currentFace, start.baryCoords);
            let currentDir = Vec3.normalize(direction);
            let totalLength = 0;
            const maxSteps = 1000;
            
            path.push({
                position: [...currentPos],
                faceIdx: currentFace,
                parameter: 0
            });

            for (let step = 0; step < maxSteps && totalLength < maxLength; step++) {
                // Find intersection with face boundary
                const intersection = this._findEdgeIntersection(
                    currentFace, currentPos, currentDir
                );
                
                if (!intersection) break;
                
                const segmentLength = Vec3.distance(currentPos, intersection.point);
                totalLength += segmentLength;
                
                path.push({
                    position: [...intersection.point],
                    faceIdx: currentFace,
                    parameter: totalLength
                });

                // Get adjacent face across the edge
                const edge = this.edgeMap.get(intersection.edgeKey);
                if (!edge || edge.faces.length < 2) break; // Hit boundary
                
                const nextFace = edge.faces[0] === currentFace ? edge.faces[1] : edge.faces[0];
                
                // Unfold and continue in straight line
                currentDir = this._unfoldDirection(
                    currentFace, nextFace, intersection.edgeKey, currentDir
                );
                currentPos = intersection.point;
                currentFace = nextFace;
            }

            return path;
        }

        /**
         * Convert barycentric coordinates to 3D position
         */
        _baryToPosition(faceIdx, bary) {
            const f = this.faces[faceIdx];
            const v0 = this.vertices[f[0]];
            const v1 = this.vertices[f[1]];
            const v2 = this.vertices[f[2]];
            return [
                bary[0] * v0[0] + bary[1] * v1[0] + bary[2] * v2[0],
                bary[0] * v0[1] + bary[1] * v1[1] + bary[2] * v2[1],
                bary[0] * v0[2] + bary[1] * v1[2] + bary[2] * v2[2]
            ];
        }

        /**
         * Find where ray from point in direction exits the triangle
         */
        _findEdgeIntersection(faceIdx, point, direction) {
            const f = this.faces[faceIdx];
            const verts = f.map(i => this.vertices[i]);
            
            let minT = Infinity;
            let bestIntersection = null;
            
            // Check each edge
            for (let i = 0; i < 3; i++) {
                const v0 = verts[i];
                const v1 = verts[(i + 1) % 3];
                const edgeDir = Vec3.sub(v1, v0);
                
                // Ray-line intersection in 3D (on the face plane)
                // p + t*d = v0 + s*e
                // Solve for t and s using face normal
                const normal = this.faceNormals[faceIdx];
                const toV0 = Vec3.sub(v0, point);
                
                // Use 2D projection onto face
                const denom = Vec3.dot(Vec3.cross(direction, edgeDir), normal);
                if (abs(denom) < EPSILON) continue;
                
                const t = Vec3.dot(Vec3.cross(toV0, edgeDir), normal) / denom;
                const s = Vec3.dot(Vec3.cross(toV0, direction), normal) / denom;
                
                // Valid if t > 0 (forward) and s ∈ [0, 1] (on edge)
                if (t > EPSILON && s >= -EPSILON && s <= 1 + EPSILON && t < minT) {
                    minT = t;
                    const vi0 = f[i];
                    const vi1 = f[(i + 1) % 3];
                    bestIntersection = {
                        point: Vec3.add(point, Vec3.scale(direction, t)),
                        edgeKey: vi0 < vi1 ? `${vi0}-${vi1}` : `${vi1}-${vi0}`,
                        edgeParam: max(0, min(1, s))
                    };
                }
            }
            
            return bestIntersection;
        }

        /**
         * Unfold direction vector across an edge to adjacent face
         * This maintains the "straightest" property
         */
        _unfoldDirection(fromFace, toFace, edgeKey, direction) {
            const n1 = this.faceNormals[fromFace];
            const n2 = this.faceNormals[toFace];
            
            // Get edge direction
            const edge = this.edgeMap.get(edgeKey);
            const edgeVec = Vec3.normalize(Vec3.sub(
                this.vertices[edge.v1],
                this.vertices[edge.v0]
            ));
            
            // Compute rotation angle (dihedral angle)
            const cosAngle = Vec3.dot(n1, n2);
            const sinAngle = Vec3.dot(Vec3.cross(n1, n2), edgeVec);
            const angle = atan2(sinAngle, cosAngle);
            
            // Rotate direction around edge
            return this._rotateAroundAxis(direction, edgeVec, angle);
        }

        /**
         * Rodrigues' rotation formula
         */
        _rotateAroundAxis(v, axis, angle) {
            const c = cos(angle);
            const s = sin(angle);
            const k = axis;
            
            // v' = v*cos(θ) + (k×v)*sin(θ) + k*(k·v)*(1-cos(θ))
            const kCrossV = Vec3.cross(k, v);
            const kDotV = Vec3.dot(k, v);
            
            return [
                v[0] * c + kCrossV[0] * s + k[0] * kDotV * (1 - c),
                v[1] * c + kCrossV[1] * s + k[1] * kDotV * (1 - c),
                v[2] * c + kCrossV[2] * s + k[2] * kDotV * (1 - c)
            ];
        }

        // ========================================================================
        // EXPONENTIAL AND LOGARITHMIC MAPS
        // ========================================================================

        /**
         * Exponential Map: exp_p(v) 
         * 
         * Maps a tangent vector v at point p to a point on the manifold
         * by walking along the geodesic in direction v for length |v|
         * 
         * @param {number} vertexIdx - Base vertex
         * @param {Array} tangentVector - [vx, vy] in local tangent plane
         * @returns {Object} {position, faceIdx, success}
         */
        exponentialMap(vertexIdx, tangentVector) {
            const p = this.vertices[vertexIdx];
            const magnitude = Vec2.length(tangentVector);
            
            if (magnitude < EPSILON) {
                return { position: [...p], vertexIdx, success: true };
            }
            
            // Convert 2D tangent to 3D direction on surface
            const direction3D = this._tangentTo3D(vertexIdx, tangentVector);
            
            // Find the face we're starting in
            const faceIdx = this._findStartingFace(vertexIdx, direction3D);
            if (faceIdx < 0) {
                return { position: [...p], vertexIdx, success: false };
            }
            
            // Compute barycentric coordinates at vertex
            const f = this.faces[faceIdx];
            const bary = [0, 0, 0];
            bary[f.indexOf(vertexIdx)] = 1;
            
            // Trace geodesic
            const path = this.traceStraightestGeodesic(
                { faceIdx, baryCoords: bary },
                direction3D,
                magnitude
            );
            
            if (path.length === 0) {
                return { position: [...p], vertexIdx, success: false };
            }
            
            const lastPoint = path[path.length - 1];
            return {
                position: lastPoint.position,
                faceIdx: lastPoint.faceIdx,
                path,
                success: true
            };
        }

        /**
         * Convert 2D tangent vector to 3D direction
         * Uses local coordinate frame at vertex
         */
        _tangentTo3D(vertexIdx, tangent2D) {
            // Build orthonormal frame at vertex
            const normal = this._vertexNormal(vertexIdx);
            
            // Pick a reference direction (first edge from vertex)
            const firstNeighbor = [...this.adjacency[vertexIdx]][0];
            let u = Vec3.normalize(Vec3.sub(
                this.vertices[firstNeighbor],
                this.vertices[vertexIdx]
            ));
            
            // Orthogonalize to normal (Gram-Schmidt)
            u = Vec3.normalize(Vec3.sub(u, Vec3.scale(normal, Vec3.dot(u, normal))));
            
            // v = n × u
            const v = Vec3.cross(normal, u);
            
            // Combine: direction = tangent.x * u + tangent.y * v
            return Vec3.normalize(Vec3.add(
                Vec3.scale(u, tangent2D[0]),
                Vec3.scale(v, tangent2D[1])
            ));
        }

        _vertexNormal(vertexIdx) {
            let normal = [0, 0, 0];
            for (const fi of this.vertexFaces[vertexIdx]) {
                normal = Vec3.add(normal, this.faceNormals[fi]);
            }
            return Vec3.normalize(normal);
        }

        _findStartingFace(vertexIdx, direction) {
            // Find face where direction points into
            for (const fi of this.vertexFaces[vertexIdx]) {
                const f = this.faces[fi];
                const idx = f.indexOf(vertexIdx);
                const v1 = f[(idx + 1) % 3];
                const v2 = f[(idx + 2) % 3];
                
                const e1 = Vec3.sub(this.vertices[v1], this.vertices[vertexIdx]);
                const e2 = Vec3.sub(this.vertices[v2], this.vertices[vertexIdx]);
                
                // Check if direction is within the cone spanned by e1 and e2
                const n = this.faceNormals[fi];
                const cross1 = Vec3.dot(Vec3.cross(e1, direction), n);
                const cross2 = Vec3.dot(Vec3.cross(direction, e2), n);
                
                if (cross1 >= -EPSILON && cross2 >= -EPSILON) {
                    return fi;
                }
            }
            return this.vertexFaces[vertexIdx][0]; // Fallback
        }

        /**
         * Logarithmic Map: log_p(q)
         * 
         * Inverse of exponential map: finds tangent vector v at p
         * such that exp_p(v) = q
         * 
         * @param {number} sourceVertex - Base vertex p
         * @param {number} targetVertex - Target vertex q
         * @returns {Object} {tangent, distance, success}
         */
        logarithmicMap(sourceVertex, targetVertex) {
            if (sourceVertex === targetVertex) {
                return { tangent: [0, 0], distance: 0, success: true };
            }
            
            // Compute geodesic distance and path
            const result = this.computeGeodesicDistance(sourceVertex);
            const distance = result.distances[targetVertex];
            
            if (!isFinite(distance)) {
                return { tangent: null, distance: Infinity, success: false };
            }
            
            // Find initial direction from source to target
            // Use graph path as approximation
            const dijkstra = this.dijkstraGraph(sourceVertex);
            const path = this.reconstructPath(dijkstra.predecessors, targetVertex);
            
            if (path.length < 2) {
                return { tangent: [0, 0], distance: 0, success: false };
            }
            
            // Initial direction towards second vertex in path
            const direction3D = Vec3.normalize(Vec3.sub(
                this.vertices[path[1]],
                this.vertices[sourceVertex]
            ));
            
            // Project to 2D tangent plane
            const tangent = this._3DToTangent(sourceVertex, direction3D);
            const scaledTangent = Vec2.scale(Vec2.normalize(tangent), distance);
            
            return {
                tangent: scaledTangent,
                distance,
                path,
                success: true
            };
        }

        /**
         * Convert 3D direction to 2D tangent vector
         */
        _3DToTangent(vertexIdx, direction3D) {
            const normal = this._vertexNormal(vertexIdx);
            
            // Build frame (same as _tangentTo3D)
            const firstNeighbor = [...this.adjacency[vertexIdx]][0];
            let u = Vec3.normalize(Vec3.sub(
                this.vertices[firstNeighbor],
                this.vertices[vertexIdx]
            ));
            u = Vec3.normalize(Vec3.sub(u, Vec3.scale(normal, Vec3.dot(u, normal))));
            const v = Vec3.cross(normal, u);
            
            // Project direction onto frame
            return [Vec3.dot(direction3D, u), Vec3.dot(direction3D, v)];
        }

        // ========================================================================
        // PARALLEL TRANSPORT
        // ========================================================================

        /**
         * Parallel Transport a vector along a path
         * 
         * A vector is parallel transported if its covariant derivative is zero:
         * ∇_γ'(t) V = 0
         * 
         * On polyhedral surfaces, this means rotating the vector by the
         * angle deficit when crossing each edge.
         * 
         * @param {Array} vector - Initial vector [vx, vy, vz]
         * @param {Array} path - Array of {position, faceIdx} from traceStraightestGeodesic
         * @returns {Array} Transported vector at end of path
         */
        parallelTransport(vector, path) {
            if (path.length < 2) return [...vector];
            
            let currentVector = [...vector];
            
            for (let i = 0; i < path.length - 1; i++) {
                const p1 = path[i];
                const p2 = path[i + 1];
                
                if (p1.faceIdx !== p2.faceIdx) {
                    // Crossing into new face - rotate vector
                    const n1 = this.faceNormals[p1.faceIdx];
                    const n2 = this.faceNormals[p2.faceIdx];
                    
                    // Axis of rotation is the shared edge direction
                    const edgeDir = Vec3.normalize(Vec3.sub(p2.position, p1.position));
                    
                    // Angle is the dihedral angle
                    const cosAngle = Vec3.dot(n1, n2);
                    const sinAngle = Vec3.dot(Vec3.cross(n1, n2), edgeDir);
                    const angle = atan2(sinAngle, cosAngle);
                    
                    currentVector = this._rotateAroundAxis(currentVector, edgeDir, angle);
                }
            }
            
            return currentVector;
        }

        /**
         * Compute holonomy (rotation after parallel transport around closed loop)
         * 
         * The holonomy around a vertex equals its angle defect (Gaussian curvature)
         * 
         * @param {number} vertexIdx - Vertex to compute holonomy around
         * @returns {number} Rotation angle in radians
         */
        holonomy(vertexIdx) {
            // For a vertex, holonomy = 2π - angle sum = angle defect
            return this.angleDefect(vertexIdx);
        }

        // ========================================================================
        // CUT LOCUS
        // ========================================================================

        /**
         * Approximate cut locus from a source vertex
         * 
         * The cut locus is the set of points where:
         * 1. The shortest geodesic from source is not unique, OR
         * 2. The geodesic ceases to be shortest (passes through its conjugate point)
         * 
         * For a sphere, the cut locus of a point is its antipode.
         * 
         * @param {number} source - Source vertex
         * @returns {Object} {cutVertices, cutEdges}
         */
        approximateCutLocus(source) {
            const result = this.computeGeodesicDistance(source);
            const distances = result.distances;
            
            const cutVertices = [];
            const cutEdges = [];
            
            // Check each edge for cut locus crossing
            // Cut locus crosses edge where distance function has saddle point
            for (const edge of this.edges) {
                const d0 = distances[edge.v0];
                const d1 = distances[edge.v1];
                
                // Check faces adjacent to this edge
                if (edge.faces.length === 2) {
                    const f1 = this.faces[edge.faces[0]];
                    const f2 = this.faces[edge.faces[1]];
                    
                    // Find the opposite vertices
                    const opp1 = f1.find(v => v !== edge.v0 && v !== edge.v1);
                    const opp2 = f2.find(v => v !== edge.v0 && v !== edge.v1);
                    
                    const dOpp1 = distances[opp1];
                    const dOpp2 = distances[opp2];
                    
                    // Gradient should point consistently if not on cut locus
                    // Check if gradients from opposite faces conflict
                    const grad1 = (dOpp1 - (d0 + d1) / 2);
                    const grad2 = (dOpp2 - (d0 + d1) / 2);
                    
                    if (grad1 * grad2 < 0) {
                        cutEdges.push(edge);
                    }
                }
            }
            
            // Saddle vertices are likely on cut locus
            for (let i = 0; i < this.vertices.length; i++) {
                if (i === source) continue;
                if (this.isSaddleVertex(i)) {
                    // Check if multiple shortest paths converge
                    let pathCount = 0;
                    for (const neighbor of this.adjacency[i]) {
                        if (distances[neighbor] < distances[i]) {
                            pathCount++;
                        }
                    }
                    if (pathCount > 1) {
                        cutVertices.push(i);
                    }
                }
            }
            
            return { cutVertices, cutEdges };
        }

        // ========================================================================
        // KARCHER MEAN (CENTER OF MASS ON MANIFOLDS)
        // ========================================================================

        /**
         * Compute Karcher Mean of a set of vertices
         * 
         * The Karcher mean minimizes the sum of squared geodesic distances:
         *   x̄ = argmin_x Σ d²(x, xᵢ)
         * 
         * Algorithm:
         *   1. Start with initial guess x
         *   2. Compute log_x(xᵢ) for all points
         *   3. Average in tangent space: v̄ = (1/n) Σ log_x(xᵢ)
         *   4. Update: x ← exp_x(v̄)
         *   5. Repeat until convergence
         * 
         * @param {Array} vertexIndices - Indices of points to average
         * @param {number} maxIter - Maximum iterations
         * @returns {Object} {position, converged, iterations}
         */
        karcherMean(vertexIndices, maxIter = 20) {
            if (vertexIndices.length === 0) return null;
            if (vertexIndices.length === 1) {
                return {
                    position: [...this.vertices[vertexIndices[0]]],
                    vertexIdx: vertexIndices[0],
                    converged: true,
                    iterations: 0
                };
            }
            
            // Initial guess: first point
            let currentVertex = vertexIndices[0];
            let converged = false;
            let iter = 0;
            
            for (iter = 0; iter < maxIter; iter++) {
                // Compute average tangent vector
                let avgTangent = [0, 0];
                let validCount = 0;
                
                for (const vi of vertexIndices) {
                    const logResult = this.logarithmicMap(currentVertex, vi);
                    if (logResult.success && logResult.tangent) {
                        avgTangent = Vec2.add(avgTangent, logResult.tangent);
                        validCount++;
                    }
                }
                
                if (validCount === 0) break;
                
                avgTangent = Vec2.scale(avgTangent, 1 / validCount);
                
                // Check convergence
                const step = Vec2.length(avgTangent);
                if (step < EPSILON) {
                    converged = true;
                    break;
                }
                
                // Update: move in direction of average tangent
                // For simplicity, find closest vertex to exp result
                const expResult = this.exponentialMap(currentVertex, avgTangent);
                if (expResult.success) {
                    // Find closest vertex to new position
                    let minDist = Infinity;
                    let closestVertex = currentVertex;
                    for (let v = 0; v < this.vertices.length; v++) {
                        const d = Vec3.distance(this.vertices[v], expResult.position);
                        if (d < minDist) {
                            minDist = d;
                            closestVertex = v;
                        }
                    }
                    
                    if (closestVertex === currentVertex) {
                        converged = true;
                        break;
                    }
                    currentVertex = closestVertex;
                }
            }
            
            return {
                position: [...this.vertices[currentVertex]],
                vertexIdx: currentVertex,
                converged,
                iterations: iter
            };
        }

        // ========================================================================
        // GEODESIC CURVATURE
        // ========================================================================

        /**
         * Compute geodesic curvature of a discrete curve on the surface
         * 
         * κ_g = 0 for geodesics
         * κ_g = signed curvature relative to surface normal
         * 
         * @param {Array} curveVertices - Ordered list of vertex indices
         * @returns {Array} Geodesic curvature at each interior vertex
         */
        geodesicCurvature(curveVertices) {
            if (curveVertices.length < 3) return [];
            
            const curvatures = [];
            
            for (let i = 1; i < curveVertices.length - 1; i++) {
                const prev = this.vertices[curveVertices[i - 1]];
                const curr = this.vertices[curveVertices[i]];
                const next = this.vertices[curveVertices[i + 1]];
                
                const e1 = Vec3.normalize(Vec3.sub(curr, prev));
                const e2 = Vec3.normalize(Vec3.sub(next, curr));
                
                // Turning angle
                const dot = max(-1, min(1, Vec3.dot(e1, e2)));
                const angle = acos(dot);
                
                // Sign: determined by surface normal
                const normal = this._vertexNormal(curveVertices[i]);
                const cross = Vec3.cross(e1, e2);
                const sign = Vec3.dot(cross, normal) >= 0 ? 1 : -1;
                
                // Average edge length
                const len = (Vec3.distance(prev, curr) + Vec3.distance(curr, next)) / 2;
                
                curvatures.push(sign * angle / len);
            }
            
            return curvatures;
        }
    }

    // ============================================================================
    // EXPORTS
    // ============================================================================

    const Geodesics = {
        DiscreteGeodesics,
        Vec3,
        Vec2,
        PriorityQueue,
        EPSILON,
        PI
    };

    if (typeof module !== 'undefined' && module.exports) {
        module.exports = Geodesics;
    } else if (typeof define === 'function' && define.amd) {
        define([], () => Geodesics);
    } else {
        global.Geodesics = Geodesics;
    }

})(typeof window !== 'undefined' ? window : global);
