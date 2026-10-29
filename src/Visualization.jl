"""
    plotgeometry(object; model=:ball, resolution=64, cutoff=1e-3)

Adapt a point, geodesic, segment, plane, horosphere, polygon, polyhedron,
or tetrahedral complex to GeometryBasics points or a
triangle mesh in the chosen model. Requires Makie and GeometryBasics to be
loaded (loading a Makie backend also loads both).

`resolution` controls curve samples and primitive tessellation. `cutoff` in
`(0,1)` trims infinite ends: geodesics/planes stop short of their ideal
boundary; horosphere triangles near the unit sphere are omitted. Surface
vertices lie on the object; triangle interiors approximate curved surfaces.
Polygon meshes omit triangles reaching Klein radius `1-cutoff`.
Rendering uses Float64 coordinates, independently of the algebra's coefficients.
"""
function plotgeometry end

"""
    hyperbolicplot(object; model=:ball, resolution=64, cutoff=1e-3,
                   boundary=true, figure=(;), axis=(;), kwargs...)

Create a Makie Figure and Axis3; return `(figure, axis, plot)` as a named
tuple. `figure` and `axis` are keyword named tuples for Makie constructors.
Other keywords go to scatter!, lines!, or mesh! as appropriate.
Load a Makie backend first, e.g. `using GLMakie` or `using CairoMakie`.
"""
function hyperbolicplot end

"""
    hyperbolicplot!(axis, object; model=:ball, resolution=64, cutoff=1e-3, kwargs...)

Add geometry to a Makie axis. Points use scatter!, geodesics lines!, planes
and horospheres mesh!. A vector of points uses a single scatter plot; a mixed
vector of objects returns a vector of plot handles.
Use the same model for all objects in an axis. Makie styling keywords pass
through directly; transparent surfaces can use `color=(:steelblue,0.25)`.
"""
function hyperbolicplot! end

"""
    modelboundary!(axis; model=:ball, resolution=64, extent=3, kwargs...)

Draw the unit sphere with Makie wireframe! for ball/Klein models, or the
height-zero square `[-extent,extent]²` with mesh! for upper half-space.
Only the finite square of the unbounded ideal plane is shown.
"""
function modelboundary! end
