"""
    LimitSetHull

Numerical convex hull of a finite ideal-point sample in the Klein ball.
`vertices(h)` retains original ideal points, `h.klein` contains Float64 working
coordinates, `faceindices(h)` contains oriented triangles, and `h.dimension`
is the numerical affine dimension (0–3). Planar hulls are filled polygons
triangulated as a fan; dimensions 0 and 1 have no faces.
"""
struct LimitSetHull
    vertices::Vector{IdealPoint}
    klein::Vector{SVector{3,Float64}}
    faces::Vector{NTuple{3,Int}}
    dimension::Int
    atol::Float64
    samplecount::Int
end
vertices(h::LimitSetHull) = h.vertices
faceindices(h::LimitSetHull) = h.faces
function edges(h::LimitSetHull)
    h.dimension == 0 && return Tuple{Int,Int}[]
    h.dimension == 1 && return [(1,2)]
    h.dimension == 2 && return [(i,mod1(i+1,length(h.vertices))) for i in eachindex(h.vertices)]
    sort!(collect(Set(minmax(f[i],f[mod1(i+1,3)]) for f in h.faces for i in 1:3)))
end
_hullnorm(v) = sqrt(sum(abs2,v))
function _hullnormal(ks,f)
    n = _cross3(ks[f[2]]-ks[f[1]],ks[f[3]]-ks[f[1]])
    n/_hullnorm(n)
end
function _hullface(ks,a,b,c,center)
    n = _cross3(ks[b]-ks[a],ks[c]-ks[a])
    _hullnorm(n) > 0 || throw(DomainError((a,b,c),"degenerate hull face; change sample or tolerance"))
    _dot3(n,center-ks[a]) > 0 ? (a,c,b) : (a,b,c)
end

# Monotone-chain hull in a basis of the numerical affine plane.
function _planarhull(ks,origin,u,n,tol)
    v = _cross3(n,u)
    xy = [(_dot3(p-origin,u),_dot3(p-origin,v)) for p in ks]
    ordered = sortperm(xy)
    turn(a,b,c) = (xy[b][1]-xy[a][1])*(xy[c][2]-xy[a][2]) -
                  (xy[b][2]-xy[a][2])*(xy[c][1]-xy[a][1])
    function chain(order)
        out = Int[]
        for i in order
            while length(out) >= 2
                a,b = out[end-1],out[end]
                # Distance tolerance, not a scale-dependent area tolerance.
                turn(a,b,i) > tol*hypot(xy[b][1]-xy[a][1],xy[b][2]-xy[a][2]) && break
                pop!(out)
            end
            push!(out,i)
        end
        out
    end
    lo,hi = chain(ordered),chain(reverse(ordered))
    vcat(lo[1:end-1],hi[1:end-1])
end

"""
    limitsethull(points; maxpoints=4096, atol=1e-10)

Compute a sampled hyperbolic convex hull using straight convexity in the Klein
ball. Accepts any finite iterable of `IdealPoint`s, including limit-set streams.
Input count is bounded by `maxpoints` (throws rather than silently truncating).
Use `Iterators.take(stream,n)` to choose a smaller hull sample explicitly.

Predicates use Float64 and absolute Klein-distance tolerance `atol`. Nearby
points are merged and nearly planar data is treated as planar; this is a
visual approximation, not an exact or certified hull of the entire limit set.
Original representatives are retained for surviving vertices. The incremental
algorithm uses quadratic worst-case work and is intended for modest samples,
separate from the dense limit-set overlay. Empty input is rejected.
"""
function limitsethull(points;maxpoints::Integer=4096,atol::Real=1e-10)
    maxpoints >= 1 || throw(ArgumentError("maxpoints must be positive"))
    isfinite(atol) && 0 < atol < 1 || throw(ArgumentError("atol must lie strictly between 0 and 1"))
    tol = Float64(atol)
    tol > 0 || throw(ArgumentError("atol is below Float64 range"))
    ps = IdealPoint[]
    ks = SVector{3,Float64}[]
    count = 0
    for p in points
        count += 1
        count <= maxpoints || throw(ArgumentError("hull sample exceeds maxpoints; use Iterators.take or increase the budget"))
        p isa IdealPoint || throw(ArgumentError("limitsethull expects ideal points"))
        k = SVector{3,Float64}(modelcoordinates(p;model=:klein))
        all(isfinite,k) || throw(DomainError(p,"hull needs finite Float64 Klein coordinates"))
        any(q -> _hullnorm(k-q) <= tol,ks) && continue
        push!(ps,p); push!(ks,k)
    end
    isempty(ps) && throw(ArgumentError("hull sample is empty"))
    finish(ids,fs,dim) = LimitSetHull(ps[ids],ks[ids],fs,dim,tol,count)
    length(ps) == 1 && return finish([1],NTuple{3,Int}[],0)
    a = 1
    b = argmax([sum(abs2,k-ks[a]) for k in ks])
    u = (ks[b]-ks[a])/_hullnorm(ks[b]-ks[a])
    distances = [_hullnorm(_cross3(k-ks[a],u)) for k in ks]
    c = argmax(distances)
    if distances[c] <= tol
        t = [_dot3(k-ks[a],u) for k in ks]
        return finish([argmin(t),argmax(t)],NTuple{3,Int}[],1)
    end
    n = _hullnormal(ks,(a,b,c))
    heights = [abs(_dot3(n,k-ks[a])) for k in ks]
    d = argmax(heights)
    if heights[d] <= tol
        ids = _planarhull(ks,ks[a],u,n,tol)
        fs = [(1,i,i+1) for i in 2:length(ids)-1]
        return finish(ids,fs,length(ids) >= 3 ? 2 : length(ids)-1)
    end
    center = (ks[a]+ks[b]+ks[c]+ks[d])/4
    fs = [_hullface(ks,f...,center) for f in ((a,b,c),(a,b,d),(a,c,d),(b,c,d))]
    for i in eachindex(ks)
        i in (a,b,c,d) && continue
        visible = [j for j in eachindex(fs) if _dot3(_hullnormal(ks,fs[j]),ks[i]-ks[fs[j][1]]) > tol]
        isempty(visible) && continue
        horizon = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
        for j in visible
            f = fs[j]
            for e in ((f[1],f[2]),(f[2],f[3]),(f[3],f[1]))
                key = minmax(e...)
                haskey(horizon,key) ? delete!(horizon,key) : (horizon[key] = e)
            end
        end
        deleteat!(fs,visible)
        for key in sort!(collect(keys(horizon)))
            e = horizon[key]
            push!(fs,_hullface(ks,e...,i,center))
        end
    end
    # Reject a broken numerical mesh rather than display plausible wrong geometry.
    incidences = Dict{Tuple{Int,Int},Int}()
    for f in fs
        normal = _hullnormal(ks,f)
        all(k -> _dot3(normal,k-ks[f[1]]) <= 8tol,ks) ||
            throw(DomainError(tol,"unstable hull support planes; change sample or tolerance"))
        for j in 1:3
            key = minmax(f[j],f[mod1(j+1,3)])
            incidences[key] = get(incidences,key,0)+1
        end
    end
    all(==(2),values(incidences)) || throw(DomainError(tol,"hull boundary is not closed"))
    ids = sort!(unique!([i for f in fs for i in f]))
    lookup = Dict(i=>j for (j,i) in enumerate(ids))
    finish(ids,[map(i -> lookup[i],f) for f in fs],3)
end

"""
    hulledgesgeometry(h; model=:klein, resolution=64, cutoff=1e-3,
                      pleats=false, minangle=0, planetol=h.atol)

Makie/GeometryBasics adapter returning sampled hull edges with NaN separators,
for `lines!`. Set `pleats=true` to exclude coplanar diagonals and bends at or
below `minangle` radians, using `pleatedges`. The default retains all edges.
Ball and half-space edges trim a fraction `cutoff` from each end
in Klein edge parameters. Requires a loaded Makie backend. Use `0 < cutoff < 1/2`.
`plotgeometry(h; model, resolution=16, cutoff)` tessellates each face in its
supporting plane, clipping its ideal corners by the same edge fraction before
conversion. Klein meshes retain the full original triangles. Curved triangles
are approximations; increase resolution for smoother surfaces. Trimming leaves
open ends, not artificial caps or a closed finite-volume core.
"""
function hulledgesgeometry end

"""
    hullbendingangles(h; planetol=h.atol)

Return `(edge=(i,j), angle=θ)` records for the edges of a three-dimensional
sampled hull. Angles are exterior hyperbolic bending angles in radians, from
outward Lorentzian supporting-plane normals. Coplanar triangulation edges
have angle zero. `planetol` compares Euclidean-unit plane coefficients in the
Klein chart before computing angles. BigFloat working arithmetic reduces
cancellation near ideal vertices but does not recover precision lost in the
sample. Hulls of dimension less than 3 have no two-sided bending edges.
"""
function hullbendingangles(h::LimitSetHull;planetol::Real=h.atol)
    isfinite(planetol) && planetol >= 0 || throw(ArgumentError("planetol must be finite and nonnegative"))
    result = NamedTuple{(:edge,:angle),Tuple{Tuple{Int,Int},Float64}}[]
    h.dimension < 3 && return result
    adjacency = Dict{Tuple{Int,Int},Vector{Int}}()
    for (j,f) in enumerate(h.faces), i in 1:3
        push!(get!(adjacency,minmax(f[i],f[mod1(i+1,3)]),Int[]),j)
    end
    setprecision(BigFloat,256) do
        planes = map(h.faces) do f
            a,b,c = (BigFloat.(h.klein[i]) for i in f)
            n = _cross3(b-a,c-a)
            n /= _hullnorm(n)
            SVector(_dot3(n,a),n...)
        end
        for edge in sort!(collect(keys(adjacency)))
            incident = adjacency[edge]
            length(incident) == 2 || throw(ArgumentError("bending requires two incident faces per edge"))
            a,b = planes[incident[1]],planes[incident[2]]
            if maximum(abs,a-b) <= planetol
                push!(result,(;edge,angle=0.0))
                continue
            end
            aa,bb = _lorentz(a,a),_lorentz(b,b)
            aa > 0 && bb > 0 || throw(DomainError(edge,"supporting plane lost at sample precision"))
            cosine = _lorentz(a,b)/sqrt(aa*bb)
            abs(cosine) <= 1+big"1e-8" || throw(DomainError(edge,"bending angle lost at sample precision"))
            push!(result,(;edge,angle=Float64(acos(clamp(cosine,-one(cosine),one(cosine))))))
        end
    end
    result
end

"""
    pleatedges(h; minangle=0, planetol=h.atol)

Select edges whose exterior bending angle is strictly greater than `minangle`
(in radians, 0–π). Coplanar triangulation diagonals are always excluded. A
planar hull has no pleat edges: its perimeter remains available via `edges`.
These are bends of the finite sampled hull, not certified limit-set pleating.
"""
function pleatedges(h::LimitSetHull;minangle::Real=0,planetol::Real=h.atol)
    isfinite(minangle) && 0 <= minangle <= pi || throw(ArgumentError("minangle must lie between 0 and π"))
    [r.edge for r in hullbendingangles(h;planetol) if r.angle > minangle]
end
