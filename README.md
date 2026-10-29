# Kleinian.jl

This repository also contains **Poincare**, **Riley**, and **Farey** as separate
Julia packages under `packages/`, with shared research work in `experiments/`.
See the [repository and development workflow](docs/repository.md) for package
boundaries, setup, Git conventions, and automated checks.

A Julia library in development for studying Kleinian groups, hyperbolic
geometry, fundamental domains, and deformation spaces.

Start with the [development plan](docs/development-plan.md) and the
[geometric algebra API](docs/geometric-algebra.md). The kernel uses immutable
`SVector{16,T}` storage with generic coefficients, basis `(e0,e1,e2,e3)`,
and metric `(-,+,+,+)`. [Hyperbolic points](docs/hyperbolic-points.md) use
future timelike projective representatives without automatic normalization.
The [geometry and isometry API](docs/isometries.md) adds GA geodesics, planes,
horospheres, half-turns, translations, elliptics, parabolics, and projective
PSL₂(C) matrix representations. Nemo algebraic numbers are supported through
an optional extension.

[Makie visualizations](docs/visualization.md) render points, ideal points,
geodesics, planes, and horospheres in Klein, Poincaré ball, and upper half-space
models. Loading a Makie backend enables the optional GeometryBasics adapters.
See `examples/visualizations.jl` for a three-model comparison.

[Complex boundary and limit sets](docs/limitsets.md) adds complex ideal-point
coordinates, reduced-word orbit sampling, and 2D Makie Datashader rendering.
Run `examples/limitsets.jl` for a classical Schottky fractal, or
`examples/riley_limitsets.jl` for million-sample Riley cusp comparisons.
`limitsetstream` and seeded `randomlimitset` walks feed Datashader directly.
`limitsethull` builds a sampled Klein convex hull, including planar Fuchsian
cases. Run `examples/riley_interactive.jl` to drag the Riley parameter and
explore linked limit-set and hull views with a hull-sample budget control
and a selector for Klein, Poincaré ball, and upper half-space rendering.
Black pleat edges exclude coplanar diagonals; a minimum bend-angle slider
filters small bends of the sampled hull. Both interactive viewers restrict
the parameter panel to nonnegative real and imaginary parts.
For a larger fractal-only view with up to 10 million samples, run
`examples/riley_fractal_interactive.jl`.

[Polygons and tetrahedral complexes](docs/polyhedra.md) add convex faces,
dihedral angles, ideal shapes, and computed edge-cycle holonomy. The
[(3,3,3) trefoil example](docs/trefoil333.md) reconstructs the two-tetrahedron
domain, verifies its gluing equations, and provides Makie figures and an
interactive model selector.

From this directory:

```sh
julia --project=. -e 'using Kleinian'
julia --project=. test/runtests.jl
julia --project=. examples/cl31_basics.jl
julia --project=. examples/hyperbolic_points.jl
julia --project=. examples/hyperbolic_isometries.jl
```

The existing Farey library will supply symbolic words through a separate adapter
once its retained interface is established. StaticArrays is the core storage dependency.
