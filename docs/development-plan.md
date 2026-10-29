# Kleinian.jl development plan

Prepared 2026-10-28. This is the working implementation plan; unchecked items
are proposed work, not implemented functionality.

## Goal and source context

Develop a Julia library for studying Kleinian groups, especially fundamental
domains in hyperbolic 3-space. Use Lorentzian geometric algebra for geometric
computations, the Klein model for convex constructions, and Klein, Poincare
ball, and upper half-space models for visualization. Later work includes
moduli spaces and pleating coordinates.

This plan draws on the accessible requests in:

- [Hyperbolic Geometry in Cl(3,1)](chatgpt-conversation://6ae1176e-7ba0-83ec-b1f1-35106dd0cee9):
  geometric algebra descriptions of reflections, geodesics, planes, and other objects.
- [Julia Kleinian Library Plan](chatgpt-conversation://6ae017b9-4234-83eb-8316-e10256110b42):
  existing free-word code, fundamental domains, multiple models, intersections,
  moduli and pleating coordinates, with Bella and Opti as feature references.
- [Map codebase and plan cleanup](thread://01a131ab-3d3f-7f71-af2a-c884c6028097?hostId=local):
  the existing Farey package, its preferred Monoids implementation, and cleanup direction.

The two ChatGPT conversations were read through the chat tool, but their
assistant answers were returned as content-reference placeholders. Their full
proposed mathematics and architecture were unavailable. This document is a
new synthesis of the visible requirements, not a transcription of those answers.

The local Farey planning files were also inspected at
`/home/tcsiz/Code/Math-Code/farey-words/docs/`: `codebase-map.md` and
`cleanup-tracker.md`. Historical validation recorded there does not validate
this package or establish that the planned Farey migration has happened.

## Architecture and conventions

Keep three layers distinct: Clifford arithmetic, validated hyperbolic objects,
and group/domain algorithms. Plotting and symbolic-word adapters consume the
core through explicit interfaces.

Use Cl(3,1) with ordered basis `(e0,e1,e2,e3)` and metric
`diag(-1,1,1,1)`. Thus `ei*ej + ej*ei = 2gij`. Coordinates are `(t,x,y,z)`;
future means `t > 0`. Record this order in every matrix/model adapter.

Proposed mathematical contracts, to be established by tests:

- A multivector has 16 coefficients indexed by bit masks `0:15`, with bit `i`
  selecting `ei` for `i=0:3`; products use ascending basis order. Grade is the mask popcount.
- Start with an immutable `Multivector{T}` backed by an `SVector{16,T}` (StaticArrays).
  Keep scalar promotion explicit and support exact coefficients where operations permit.
  Benchmark before adding specialized storage or a second backend.
- Geometric product is primary. For homogeneous grades `r,s`, the outer product
  is the grade `r+s` part of the product. Define left contraction as the grade
  `s-r` part when `r <= s`, zero otherwise, extended bilinearly.
  The separate symmetric dot product is explicitly `a ⋅ b = (a*b+b*a)/2`,
  including for mixed grades, as requested during implementation.
- Reversion on grade `r` has sign `(-1)^(r(r-1)/2)`; grade involution has sign
  `(-1)^r`; Clifford conjugation composes them.
- With `I=e0*e1*e2*e3`, `I^2=-1` and `I^-1=-I`. Choose right dual
  `dual(A)=A*I^-1`. Specify the inverse dual separately; dual applied twice
  is not the identity under this convention. Derive meet signs from this choice.
- Do not implement general multivector inverse as `reverse(A)/(A*reverse(A))`
  without checking that the denominator is a nonzero scalar. Start with vectors
  and validated versors; arbitrary inverses are a later decision.

Store finite points as projective future timelike vectors (`X^2<0`, `t>0`),
without normalizing on construction. Normalize only for calculations that need
it, such as distance. Unit representatives have `X^2=-1`. Ideal points are future null
rays, not unit timelike points. A plane is `X dot n=0` with spacelike normal
`n`; its oriented half-space uses `X dot n <= 0`. Scaling a normal by a positive
factor preserves the half-space, while negating it changes the chosen side.

A plane reflection is `X -> -n*X*inv(n)`, equivalently
`X - 2*(X dot n)/(n dot n)*n`. An even versor acts by `R*X*inv(R)`;
reversion replaces inversion only when the normalization justifies it.
Geodesics are intersections with timelike 2-dimensional linear subspaces,
represented by simple bivectors. Plane intersections can be represented by
wedges of normals and dualized into a spanning subspace; classification must
precede normalization or construction of a hyperbolic object.

Finite-point model maps for `X=(t,v)` are Klein `k=v/t` and Poincare ball
`p=v/(t+1)`. Their inverses are `(1,k)/sqrt(1-|k|^2)` and
`(1+|p|^2,2p)/(1-|p|^2)`. Choose upper half-space coordinates explicitly
before implementing that conversion. Boundary conversions need separate contracts.
Distance obeys `cosh(d(X,Y))=-X dot Y`; handle floating-point roundoff near
coincidence without hiding genuinely invalid input.

For background on the Lorentz model and its relation to projective models,
see [Caroline Series' hyperbolic geometry notes](https://warwick.ac.uk/fac/sci/maths/people/staff/caroline_series/hyperbolic_geometry_ma448_lecture_notes.pdf).
The specific Clifford storage, contraction, and dual conventions above are
local design choices to verify, rather than conventions attributed to that source.

## Milestone 1 — geometric algebra kernel (implemented; API review pending)

- [x] Add `src/GeometricAlgebra.jl` with scalar and basis-blade constructors,
  coefficient access, grade projection, arithmetic, and scalar promotion.
- [x] Implement blade multiplication from swaps and metric signs; use it for
  the full geometric product. Keep the implementation small enough to audit.
- [x] Add outer product, explicit left contraction, reverse, grade involution,
  Clifford conjugation, dual/undual, and checked vector/versor inversion.
- [x] Document API names and operator choices; export only the stable subset.
- [x] Add exact tests for all 16-by-16 basis products and associativity on all
  basis triples, with an independent ordered-generator reference calculation.
- [x] Test metric relations, wedge antisymmetry, grade selection, involution
  identities, dual signs, zero divisors, null-vector inversion rejection,
  scalar promotion, and representative mixed multivectors.
- [x] Write `examples/cl31_basics.jl` demonstrating algebra identities and a
  reflection, checked against the independent coordinate reflection formula.

Exit: the exact algebra identities pass, conventions are documented, and the
example runs. Do this before implementing polyhedra or visualizations.

## Milestone 2 — hyperbolic objects and intersections

- [x] Add projective `HyperbolicPoint` in `src/HyperbolicGeometry.jl`, with
  validation, projective equality, positive rescaling, scale-independent
  distance, and explicitly requested unit representatives. See
  [point conventions and API](hyperbolic-points.md).
- [x] Add GA ideal points, oriented planes, oriented geodesics, and horospheres
  retaining their null-vector scale, with homogeneous incidence predicates.
- [x] Add half-turns, translations, finite-order elliptics, parabolics, products,
  inverses, and powers using generic coefficient operations. See
  [geometry and isometries](isometries.md).
- [ ] Add half-spaces, segments, and orientation-reversing plane-reflection objects.
  Geodesic segments are now implemented; half-space and reflection objects remain.
- [ ] Add join/meet operations with explicit primal/dual representation contracts.
- [ ] Classify plane-plane and line-plane intersections, including intersecting,
  asymptotic, disjoint, coincident, and degenerate cases as applicable.
- [ ] For three planes, classify the common homogeneous solution as finite,
  ideal, hyperideal, or degenerate. Return structured results rather than NaNs.
- [ ] Define tolerance policy for floating-point predicates; preserve exact
  algebra paths and allow BigFloat. Near-degenerate results may be indeterminate.
- [ ] Test reflection involutivity, fixed-plane incidence, future-sheet and
  Lorentz-product preservation, with coordinate fixtures and degenerate examples.

Exit: geometric constructions agree with independent coordinate linear algebra;
intersection results retain their classification and orientation information.

## Milestone 3 — models and group representations

- [x] Implement hyperboloid/Klein/ball/half-space conversions, with round-trip
  tests and separate treatment of ideal points and excluded chart points.
- [ ] Define generator-to-isometry maps, evaluation, inverses, and word labels.
- [x] Fix action and multiplication order with a noncommuting example: g*h
  acts by g(h(x)). Preserve source factors and powers for reevaluation.
- [x] Bridge projective complex 2x2 matrices and Clifford versors, including
  explicit SL2 lifts, Hermitian point action, and the scale correction for horospheres.
- [ ] Add an explicit Lorentz 4x4 matrix API.
- [ ] Test boundary Mobius action against interior/model actions.
- [ ] Add fixed points, axes, traces, isometry classification, and complex
  lengths with explicit branch conventions and numerical uncertainty handling.

## Reuse the existing Farey library

Keep Farey as a separate package. Do not copy its old FreeGroup implementation
or make core geometry inherit its group hierarchy. Preferred integration:
`Farey.Monoids` generators/words plus an explicit image map into isometries.
The adapter must distinguish generator blocks (`gens`) from signed letters
(`eachgen`) and preserve multiplication order and identity evaluation.

The existing cleanup plan prefers Monoids, relocates generic defaults into
`GroupInterface.jl`, and aims to remove KnuthBendix and GroupsCore. That cleanup
is a separate task, not a prerequisite for the Clifford kernel. Prior chat
validation found cyclic equality, cross-representation equality, and custom-word
DOT-export issues; check current behavior before depending on those paths.

Use existing `examples/research_*.jl` and `farey-testing/working_setup.jl` as
sources of fixtures for word evaluation, substitutions, Farey/Christoffel words,
continued fractions, and commutators. Freeze the intended commutator and
conjugation conventions explicitly. Do not execute the setup script wholesale:
it mixes representations, top-level research state, and missing includes.

## Milestone 4 — convex domains and side pairings

- [x] Add explicit convex finite/ideal polygons and polyhedra, tetrahedra,
  face-pairing complexes, dihedral/shape calculations, and directed edge-cycle
  verification. See [polyhedra](polyhedra.md) and the
  [exact (3,3,3) trefoil example](trefoil333.md). Global tiling and cusp
  completeness are not inferred from these local checks.

- [ ] Build a Klein half-space intersection layer retaining plane labels,
  vertices, edges, faces, incidence, and finite/ideal/hyperideal status.
- [ ] Start with explicit finite sets of planes and known polyhedra.
- [ ] Construct Dirichlet bisectors from a basepoint and a bounded word search.
- [ ] Track which generators/words produced each candidate face.
- [ ] Validate side-pairing incidence, cycles, and relevant Poincare polyhedron
  theorem hypotheses before claiming a fundamental domain or presentation.
- [ ] Distinguish truncated candidate domains from certified results; expose
  search limits and unresolved faces. Finite enumeration alone proves neither
  discreteness nor completeness of a domain.

## Milestone 5 — visualization and research workflows

- [x] Render the same geometric objects in Klein, ball, and half-space models;
  keep geometry independent of meshes and plotting state.
  See [Makie visualizations](visualization.md); GeometryBasics supplies meshes.
- [ ] Show face labels, side pairings, axes, ideal boundaries, and uncertainty.
- [ ] Add parameterized representations, reproducible parameter sweeps, and
  cached geometric invariants for moduli-space experiments.
- [ ] Specify pleating-coordinate conventions and supported representation
  families before implementing pleating loci and continuation algorithms.
- [ ] Compare concrete workflows with Bella and Opti after identifying their
  intended versions and authoritative sources; feature parity is not assumed.

## Proposed layout and dependency policy

Only the package entry point and scaffold test exist initially. Grow toward:

```text
src/
  Kleinian.jl
  GeometricAlgebra.jl
  HyperbolicGeometry.jl
  Models.jl
  Isometries.jl
  Representations.jl
  Domains.jl
  Deformations.jl
ext/                  # optional Farey, exact-arithmetic, and plotting adapters
examples/             # reproducible mathematical examples
test/                # kernel, geometry, model, and group regression tests
docs/                # conventions, API notes, development plan
```

Use StaticArrays for immutable fixed-size coefficient storage. Evaluate optional numerical,
plotting, polyhedral, and exact-arithmetic packages when their milestones need
them. Hecke and Oscar were mentioned in the source request; they are candidates
for specialized exact computations, not selected dependencies. Assess required
number types, startup costs, licensing, and interoperability first. Decide on
an external Clifford backend only after comparing its conventions and scalar
types with the tested kernel. Avoid installing a broad research stack now.
Nemo 0.56 is now an optional weak dependency; its QQBarFieldElem adapter is
tested separately without adding Nemo to the required runtime dependencies.

## Current status and immediate next session

- [x] Inspect accessible reference chats and Farey planning documents.
- [x] Establish a minimal loadable Julia package and scaffold test.
- [x] Record architecture, conventions, milestones, and reuse boundaries.
- [x] Implement Milestone 1 with `SVector{16,T}` and the time-first basis.
- [x] Implement the GA geometric objects and orientation-preserving isometries,
  with native rational/BigFloat and optional exact Nemo validation.
- [ ] Review storage and API choices with the user before adding compact grade types.

See [the geometric algebra API](geometric-algebra.md) and the runnable
`examples/cl31_basics.jl`, `examples/hyperbolic_isometries.jl`, and optional
`examples/nemo_isometries.jl`. [Model conversions and Makie plotting](visualization.md)
are implemented, with a comparison in `examples/visualizations.jl`.
Geometric intersections and group-word integration remain subsequent milestones.
A full inverse for arbitrary multivectors is outside this initial kernel.
