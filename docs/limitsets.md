# Complex boundary and limit sets

Ideal points remain future null rays in the existing geometric algebra.
`complexidealpoint(z)` lifts a complex number; `complexcoordinate(p)` projects
it back. `nothing` represents infinity in both directions. The display chart
is `(x+iy)/(t-z₃)`, matching half-space coordinates and `crossratio`.
The existing matrix convention acts on its conjugate, so the usual complex
fractional-linear formula uses `conj.(matrix(g))` in this chart.

`limitset(generators; depth=8, seed=complexidealpoint(0.0), maxpoints=1_000_000)`
returns ideal-point samples for words of exactly the requested depth. Inverses
are included automatically and adjacent inverse pairs are excluded. Supply
each generator once. With n generators there are `2n*(2n-1)^(depth-1)` samples
at positive depth; the budget is checked before allocation. This samples a
boundary orbit, with no discreteness or convergence certification. For suitable
Schottky groups, deeper levels approach the limit set. Relations can produce
duplicate points; the plot's density describes the sampling measure.

Load a Makie backend and call `limitsetplot(points)` for a two-dimensional
Datashader view. Keywords such as `colormap`, `operation`, `local_operation`,
`agg`, and `binsize` pass through to Makie. The default is `colormap=:inferno`
and `async=false` so saved images are fully aggregated. Use `async=true` for
interactive exploration or `renderer=:scatter, markersize=1` for dots.
Infinity is omitted, and unrepresentable Float64 values raise an error.
Axis options can set a viewport, for example `axis=(limits=(-3,3,-3,3),)`.
Zooming reaggregates existing samples; it does not generate finer orbit points.

See `examples/limitsets.jl` for a classical two-generator Schottky example.
This follows the finite word-tree perspective in
[Indra's Pearls, chapter 4](https://www.cambridge.org/core/books/abs/indras-pearls/schottky-dance-pages-107-to-120/C797CFAF11D6F5BDBACE00A31320D1BC).
Datashader options are documented in [Makie's API](https://docs.makie.org/stable/api.html#Makie.datashader).

## Implementation stages

1. Implemented: coefficient-aware complex boundary conversion, bounded
   reduced-word sampling, and Datashader/scatter adapters using existing points
   and isometries.
2. Implemented: depth-first streaming and reproducible random reduced-word
   walks, direct stream-to-Datashader projection, and Riley cusp examples.
3. Next: viewport-aware adaptive subdivision with a spherical error criterion,
   retaining branches that cross infinity; compare convergence across depths.
4. Later: parameter families and circle-pairing overlays inspired by the book,
   with interactive controls and high-resolution exports.

## Streaming and random sampling

`limitsetstream(generators; depth, seed, maxpoints)` traverses the reduced-word
tree depth-first, retaining only O(depth) ideal points. `limitset` now collects
this iterator, preserving its previous order and point-count guard.

`randomlimitset(generators; samples=100_000, burnin=100, thin=1, seed, rng)`
returns a lazy random walk retaining O(1) ideal points. It chooses each allowed
letter uniformly, forbidding only the inverse of the previous letter. Burn-in
is counted in actions; each output follows another `thin` actions. The result
is a correlated random-walk measure, not uniform geometric density. At cusps,
long parabolic excursions and rare branches can converge slowly; more samples
or thinning alone do not certify coverage. With one generator the walk cannot
reverse direction, so use deterministic sampling to see both ends.

```julia
using Kleinian, Riley, Poincare, Random, GLMakie
stream = randomlimitset(rileygroup(2.0im);
    samples=2_000_000, burnin=1000, rng=Xoshiro(2026))
view = limitsetplot(stream;
    axis=(limits=(-1.5,1.5,-1.6,0.6),), colormap=:inferno)
save("apollonian.png", view.figure; px_per_unit=2)
```

An explicitly supplied RNG is copied without advancing it. Each traversal
replays from that snapshot; the default draws a fresh seed from Julia's default
RNG. Reproducibility is within the same Julia/RNG implementation, not promised
across Julia versions. `seed` is the initial geometric point, not the RNG seed.
Iterators can be consumed with `for`, `Iterators.take`, or any online accumulator.
Their stated memory bounds count retained points: exact coefficient bit sizes
can grow during long walks, so use floating coefficients for large pictures.

Makie accepts either stream directly, converting each point to `Point2d` once.
It still retains O(samples) compact coordinates so Datashader can reaggregate
on zoom; this is not a constant-memory raster renderer. Calling `collect(stream)`
also intentionally retains all ideal points.

## Riley family and cusp examples

`rileygroup(ρ)` constructs the pair

```math
X=\begin{pmatrix}1&1\\0&1\end{pmatrix},\qquad
Y_\rho=\begin{pmatrix}1&0\\\rho&1\end{pmatrix}
```

in the **displayed complex chart**. These are the usual Riley-family
normalization in [Elzenaar, Martin and Schillewaert, Approximations of the Riley
slice](https://aelzenaar.github.io/farey/SchottkyGroups.pdf). Internally the
matrices are conjugated to match this library's existing convention. Both
parabolic generators remain ordinary `MobiusIsometry` objects. The constructor
accepts any finite nonzero parameter and does not decide slice membership.

The example `examples/riley_limitsets.jl` compares three parameters:

- `ρ=4`: the real boundary cusp; `tr(X*inv(Y)) = 2-ρ = -2`, with a
  nonidentity parabolic word. Its Fuchsian limit set lies on the extended real
  line, which appears as a line in this chart.
- `ρ=2.2im`: an approach toward the imaginary cusp along the imaginary axis.
- `ρ=2im`: the Apollonian gasket group, as recorded in
  [Elzenaar's group zoo](https://aelzenaar.github.io/zoo/index.html).
  `tr(X*Y*inv(X)*inv(Y)) = 2+ρ² = -2`, and the word is nonidentity.

The tests check the word traces exactly at these cusp values and the displayed
Möbius actions. A parabolic trace equation by itself does **not** identify a
valid boundary cusp for an arbitrary word or parameter. General rational cusp
finding needs Farey words and branch selection; it is not implemented here.

For headless export, load CairoMakie in the plotting environment, include
`examples/riley_data.jl`, and call `save("riley.png", rileyfigure())`. The example
uses two million samples per panel, a fixed viewport, and seeded RNGs.
It conjugates the displayed groups by `w=1/(z+i)` to show bounded circle
packings, including a circle for the real cusp, rather than an unbounded strip.
Use a seed constructed at the same precision as the generators for BigFloat
walks; promoting an already rounded Float64 null vector cannot restore its
null norm at higher precision.

## Sampled convex hulls in the Klein ball

`limitsethull(points; maxpoints=4096, atol=1e-10)` computes a numerical hull of
an ideal-point iterable. Since geodesics are straight in the Klein model,
the Euclidean convex hull there represents the sampled hyperbolic hull.
The implementation uses an incremental hull algorithm with quadratic worst-case
work, intended for hundreds or a few thousand points. No new dependency is
required. It retains original ideal-point representatives alongside Float64
Klein coordinates and triangulated faces.

```julia
stream = randomlimitset(rileygroup(2.0im); samples=30_000, rng=Xoshiro(2026))
hull = limitsethull(Iterators.take(stream,256))
view = hyperbolicplot(hull; color=(:steelblue,0.25), transparency=true)
hyperbolicplot!(view.axis,collect(stream);model=:klein,
                color=(:darkorange,0.5),markersize=2)
```

`vertices(hull)`, `faceindices(hull)`, and `edges(hull)` expose the mesh.
`hull.dimension` distinguishes a point, segment, planar polygon, and 3D hull.
Planar hulls are triangulated filled polygons; their `edges` omit internal fan
diagonals. A real Fuchsian limit set produces a planar hull, including after
Möbius conjugation. Three-dimensional edges include triangulation edges.
`hull.samplecount` records the input count before merging nearby vertices.

The absolute Klein-distance tolerance controls merging and numerical affine
dimension. Nearly coplanar samples are treated as planar. Small caps and close
vertices may be discarded at this tolerance; predicates are not exact. The
3D output is checked for closed edge incidence and supporting face planes.
Empty input and samples exceeding `maxpoints` raise errors. Use
`Iterators.take(stream,n)` to set a deliberate budget; the function does not
silently subsample. Distinct sphere points are extreme, so a dense limit-set
sample is generally unsuitable as the hull budget.

Rendering supports `model=:klein`, `:ball` (Poincaré), and `:halfspace`.
The default remains Klein, where full triangles and straight edges are exact.
For the curved models, `plotgeometry(hull; model, resolution=16, cutoff=1e-3)`
cuts each ideal corner at fraction `cutoff` along its two incident Klein edges,
then subdivides the remaining face polygon in its original supporting plane
before mapping to the chosen model. `resolution` is the number of subdivisions
per finite triangle edge. The resulting mesh approximates the curved surface;
its vertices lie on the face. Adjacent faces use matching trimmed boundaries.
Require `0 < cutoff < 1/2`; increasing cutoff removes more of the ideal ends.
The open ends are not capped, and their heights in half-space are not bounded
by a user-specified height. Very small faces can exceed Float64 precision;
reduce the hull sample budget if chart conversion fails.

`hulledgesgeometry(hull; model, resolution=64, cutoff=1e-3)` returns edge curves
with NaN separators for `lines!`, using the same corner trimming. For example:

```julia
view = hyperbolicplot(hull;model=:ball,boundary=false,
                      color=(:steelblue,0.25),transparency=true)
lines!(view.axis,hulledgesgeometry(hull;model=:ball);color=:black)
```

`hyperbolicplot` also renders dimensions 0 and 1 as a point or geodesic segment.
A singleton at half-space infinity raises `DomainError`, as for other ideal
point plots. Planar hulls are supported in all three models.
This remains the hull of a finite sample, not the entire limit set or the
convex core quotient. Adaptive convergence and exact predicates remain future work.

## Dragging the Riley parameter

Run:

```sh
julia --project=experiments examples/riley_interactive.jl
```

The GLMakie example links a parameter plane, a complex boundary orbit, and a
transparent Klein hull. Click anywhere in the parameter plane and drag to move
ρ. The parameter panel is restricted to `Re(ρ) ≥ 0`, `Im(ρ) ≥ 0`,
with a fixed view of `[0,4.5] × [0,3.5]`. Parameter-axis panning and zooming
are locked so the panel stays in this quadrant; the geometric views retain
their navigation controls. Scripted viewer updates reject negative components.
Preview updates use at most 3,000 orbit points and 64 hull samples, throttled
to at most eight updates per second; releasing refines to 30,000 points and the
selected hull budget. The controls offer 32–1024 hull samples and independent
face/pleat visibility. The hull view omits the sphere at infinity and shows
only bending edges in black by default. Coplanar triangulation diagonals are
hidden. The Minimum bend slider (0–180 degrees) filters small sampled bends
without resampling or rebuilding the hull. Planar hulls have no pleat edges.
Select Klein, Poincaré ball, or upper half-space with the model menu.
The selector reuses the current sampled hull without resampling the orbit.
All views omit the sphere/plane at infinity and keep black pleats enabled.
The upper half-space viewport initially shows x∈[-1.2,1.2], y∈[-2.2,0.2],
and height∈[0,2.4]; geometry beyond it is clipped only by the view.
Drag the 3D view to rotate it.

The fixed RNG seed keeps the sampled words consistent across parameter changes.
The display uses `w=1/(z+i)` to make the cusp examples compact, with a fixed
complex viewport; points outside that viewport remain in the spherical overlay.
This plane explores the **two-parabolic parameter family**, not a certified
picture of the Riley slice. Arbitrary parameters may give nondiscrete groups.
The mathematical `rileygroup` API remains unrestricted by this viewer symmetry
choice. ρ=0 is excluded; numerical failures hide the previous geometry and show a
message rather than leaving an apparently current result on screen.

For headless output, load CairoMakie and include `examples/riley_hull_explorer.jl`.
`view = riley_hull_explorer()` returns `view.figure`, controls, and
`view.set_parameter!(ρ; preview=false)` for scripted updates. GLMakie provides
the intended interactive transparency; CairoMakie is suitable for static previews.

## Fractal-only viewer for larger samples

```sh
julia --project=experiments examples/riley_fractal_interactive.jl
```

This dedicated viewer devotes most of its window to the fractal, with a smaller
parameter plane and no hull calculation or 3D point storage. Like the hull
explorer, it restricts the parameter to the first quadrant; negative real or
imaginary components are rejected before changing the current view. The sample menu
ranges from 100,000 to 10,000,000 points (one million by default). Click/drag
the parameter to request a 10,000-point preview; release to refine to the menu's
budget. The same seeded words are used across parameters and budgets. Changing
the budget also requests a new refined sample.

Sampling runs in a cooperative task, yielding every 10,000 steps. New requests
cancel older work at the next batch boundary; obsolete results do not replace
the current parameter's picture. Closing the window also cancels generation.
Julia compilation and the final Datashader aggregation can still cause pauses.
The status line reports progress and whether the visible result is a preview
or a refined sample. Invalid parameters clear the picture and show a message.

Scroll to zoom and drag within the fractal panel to pan. The **Reset view**
button restores the initial framing without regenerating points. Colour choices
are Inferno, Magma, Viridis, and Gray. Zooming reaggregates retained samples;
it does not automatically generate more points in the zoomed region.

Coordinates use Float64 `Point2d` storage: 10 million points take approximately
160 MB for one coordinate buffer, plus plotting overhead. During refinement,
the previous and replacement buffers can coexist, so peak memory is higher.
No array of ideal points or convex-hull vertices is retained.

For scripted or headless use, load a backend and include
`examples/riley_fractal_viewer.jl`. `riley_fractal_viewer(samples=2_000_000,
background=false)` finishes sampling synchronously before returning. The
interactive default uses `background=true`; wait for `view.task[]` before
saving a completed render. `view.set_parameter!(ρ; preview=false)` returns the
new task in background mode or a success Boolean in synchronous mode.


## Pleat edges and bending angles

`hullbendingangles(h; planetol=h.atol)` returns `(edge, angle)` records for the
three-dimensional sampled hull. Angles are exterior hyperbolic bending angles
in radians, computed from outward spacelike Lorentzian face normals:

```math
\theta=\arccos\left(\frac{\langle n_1,n_2\rangle}
{\sqrt{\langle n_1,n_1\rangle\langle n_2,n_2\rangle}}\right).
```

`pleatedges(h; minangle=0, planetol=h.atol)` returns edges with bending angle
strictly above `minangle`. Coplanar face pairs have angle zero and are always
excluded. Plane comparison uses Euclidean-unit plane coefficients in the Klein
chart and the numerical tolerance `planetol`. Angles use 256-bit working
arithmetic to reduce cancellation near ideal vertices; the underlying hull
coordinates are still Float64. This does not recover precision already lost
in the sample or certify the limit hull's pleating structure.

```julia
lines!(view.axis, hulledgesgeometry(hull; model=:ball,
        pleats=true, minangle=deg2rad(2)); color=:black)
```

The viewer uses this filter by default; the general `hulledgesgeometry` API
retains all edges unless `pleats=true` is requested. For dimensions 0–2 the
bending list is empty: a planar polygon's perimeter is available through
`edges(h)` or unfiltered `hulledgesgeometry`. Small nonzero bends can be sampling
artifacts; compare increasing hull budgets when interpreting the visible folds.
