# Polygons, tetrahedra and face gluings

The API keeps the geometric vertices, exact coefficient fields and face maps
separate from Makie/GeometryBasics rendering. This implementation supports
finite and ideal vertices. Hyperideal vertices and polar truncation are not
part of this (3,3,3) implementation.

```julia
using Kleinian, Poincare
v = (IdealPoint(1,0,0,1), IdealPoint(1,0,0,-1),
     IdealPoint(1,1,0,0), IdealPoint(1,0,1,0))
t = HyperbolicTetrahedron(v)
faceindices(t)                 # face i is opposite vertex i
face(t,1)                     # HyperbolicPolygon
faceplane(t,1)                # outward Lorentz normal
edges(t)                      # vertex-index pairs
dihedralangle(t,(1,2))        # radians
dihedralcos(t,(1,2))          # exact when the coefficient backend permits
grammatrix(t)                 # outward unit face-normal Gram matrix
idealshape(t)                 # crossratio(infinity,0,1,z) = z
edgeparameter(t,(1,2))        # positive-imaginary edge shape
polygonarea(face(t,1))        # π for an ideal triangle
interiorpoint(t) in t
```

`HyperbolicPolygon(vertices)` accepts a convex cyclic list of coplanar finite
or ideal vertices. `HyperbolicPolyhedron(vertices, faces)` accepts explicit
cyclic face-index lists, including quadrilateral and larger faces. It orients
faces outward and checks supporting half-spaces and a closed spherical
boundary. It does not compute a convex hull from an unordered point cloud.
`HyperbolicTetrahedron` is the four-vertex specialization with a fixed opposite
face convention. Nonconvex unions are represented by multiple cells.

`vertices`, `faceindices`, `face`, `faceplane`, `edges`, `interiorpoint`,
`contains`, `edgelength`, `interiorangles`, `polygonarea`, `grammatrix`,
`dihedralcos` and `dihedralangle` expose the geometry. An edge with an ideal
endpoint has infinite length. Areas use the Gauss–Bonnet angle defect.
Angles and areas accept `T=BigFloat`; they are numerical evaluations.
`dihedralcos`, cross-ratios and supported matrix operations can remain exact.
Isometries act on polygons, polyhedra, tetrahedra, and `GeodesicSegment(p,q)`;
`mapcoefficients` reconstructs them in another coefficient backend.

Constructors use exact predicates on exact coefficients. Floating input has
a default allowance of 128 machine epsilons on normalized chart coordinates;
an explicit nonnegative `atol` overrides it. This is a numerical validation
policy, not interval arithmetic or certification of nearly degenerate input.

## Face-pairing complexes

```julia
pair = FacePairing((source_cell,source_face), (target_cell,target_face), g;
                   label="a")
domain = TetrahedralComplex(cells,pairings)
cycles = edgecycles(domain)
report = verify_edge_gluings(domain)
report.valid
```

Each `FacePairing` supplies a source-to-target isometry. Its inverse is
automatic, so list each paired pair once. Every face must occur exactly once.
The constructor obtains the three vertex correspondences by applying the
map and checks that a source interior point goes outside the target face.
Use `internal=true` for an identity gluing across a shared geometric face.
`boundaryfaces(domain)` excludes those internal faces. `contains` treats the
displayed region as the union of its cells, so a reflex boundary edge is kept.

`edgecycles` follows directed face crossings around each tetrahedral edge,
switching to the other incident face after crossing. Its sectors partition
all local tetrahedral edges. A physical edge on an internal face contributes
one sector from each adjacent tetrahedron. The returned `EdgeCycle` stores
the sectors, signed crossings, full holonomy, angle sum, complex shape product
when all cells are ideal, and whether the circuit reverses the edge.
`cycleword(domain,cycle)` writes matrix-composition order: rightmost acts first.

`verify_edge_gluings` checks, for each cycle of cone order m:

- the angle sum is 2π/m, with the correct winding number;
- the holonomy fixes both endpoints individually;
- its m-th power is projectively the identity, with no smaller divisor order;
- for ideal tetrahedra, the product Z of edge shapes satisfies Z^m=1.

The matrix check also detects translation along the edge. The complex shape
check includes scaling, which an angle-only check misses. Edge shapes use
inversion when the vertex ordering must be reversed; replacing that inversion
with conjugation would produce incorrect moduli.

Supply `orders=[...]` to prescribe the orders in the returned cycle order.
Otherwise they are inferred from the angles, and a nonintegral result fails.
`atol=0` checks matrix, endpoint and shape equations exactly. Angle checks
remain numerical, using `angle_atol` (default 1e-9 radians). Each entry of
`report.edge_cycles` records its word, order, residuals and individual checks.
These are local edge equations. They do not establish global tiling, cusp
completeness, or faithfulness of an abstract presentation.

## Makie rendering

With a Makie backend loaded, `plotgeometry` and `hyperbolicplot!` also accept
segments, polygons, polyhedra, tetrahedra and tetrahedral complexes. Complexes
render their exterior faces only. `plotgeometry` returns a standard
GeometryBasics mesh, preserving its mesh views. Add segments separately for
edge colors, or use native Makie wireframes/materials on the returned mesh.

GeometryBasics triangulates each face and supplies the tessellation grid;
the adapter maps its coordinates into the desired hyperbolic model. Faces
incident to the ideal boundary are cut off at Klein radius `1-cutoff`.
Increase `resolution` to resolve curved faces, or reduce `cutoff` to approach
ideal vertices. Finite segment endpoints are retained; ideal segment ends
are shortened by the requested fraction of their Klein chord.

See [the trefoil construction and exact verification](trefoil333.md) and
[the general model/plotting API](visualization.md).
