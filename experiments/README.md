# Experiments

This is the shared workspace for composing the local libraries. From the repository root:

```sh
julia experiments/setup.jl
julia --project=experiments experiments/farey_riley.jl
julia --project=experiments experiments/fractals.jl
julia --project=experiments experiments/hull_explorer.jl
julia --project=experiments experiments/hulls.jl
```

`setup.jl` develops Kleinian, Poincare, Riley and Farey together and installs the
plotting dependencies. Run it again after relocating the checkout. Interactive
GLMakie examples need a graphical display; `farey_riley.jl` is a non-graphical
word-evaluation and sampled-hull example.

- `fractals.jl`: draggable first-quadrant Riley parameter and high-density Datashader viewer.
- `hull_explorer.jl`: linked fractal and sampled convex hull in three models, with black pleats and bending-angle filtering.
- `hulls.jl`: small static hull example without the sphere at infinity.
- `farey_riley.jl`: evaluate a Farey word in Riley generators and build a hull using Poincare.
- [sakuma_weeks](sakuma_weeks/README.md): extracted Chesebro polynomials, developed tetrahedra and a model/transformation viewer.
- [farey-testing](farey-testing/README.md): imported research notebooks, topic summary, limitations and separate optional algebra environment.

Validation:

```sh
julia --project=experiments experiments/check.jl
julia --project=experiments test/makie.jl
julia --project=experiments test/trefoil_exact.jl
```

The package checks run in separate Julia processes. The Makie suite uses CairoMakie
and does not require opening a desktop window. The archived research scripts are
not included in these checks. Local manifests are generated, not tracked.
