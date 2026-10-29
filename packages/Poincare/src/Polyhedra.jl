const PolyhedralPoint = Union{HyperbolicPoint,IdealPoint}
import Base: contains

"""A bounded geodesic segment, allowing ideal endpoints."""
struct GeodesicSegment{P<:PolyhedralPoint,Q<:PolyhedralPoint}
    start::P
    stop::Q
    function GeodesicSegment(p::P,q::Q) where {P<:PolyhedralPoint,Q<:PolyhedralPoint}
        Geodesic(p,q) # rejects coincident endpoints
        new{P,Q}(p,q)
    end
end
vertices(s::GeodesicSegment) = (s.start,s.stop)
edgelength(s::GeodesicSegment) = s.start isa IdealPoint || s.stop isa IdealPoint ? Inf : distance(s.start,s.stop)
(g::AbstractIsometry)(s::GeodesicSegment) = GeodesicSegment(g(s.start),g(s.stop))
mapcoefficients(f,s::GeodesicSegment) = GeodesicSegment(mapcoefficients(f,s.start),mapcoefficients(f,s.stop))

_spatial(c) = SVector(c[2],c[3],c[4])
_cross3(a,b) = SVector(a[2]*b[3]-a[3]*b[2],a[3]*b[1]-a[1]*b[3],a[1]*b[2]-a[2]*b[1])
_dot3(a,b) = sum(a[i]*b[i] for i in 1:3)
function _poly_atol(vs,atol)
    if isnothing(atol)
        return maximum(vs) do p
            T = typeof(real(scalarpart(representative(p))))
            T <: AbstractFloat ? 128eps(T) : 0
        end
    end
    isreal(atol) && isfinite(atol) && atol >= 0 || throw(ArgumentError("atol must be finite and nonnegative"))
    atol
end
_normalchart(p::HyperbolicPlane) = begin
    c = _pointcoordinates(representative(p))
    map(x -> _fielddiv(x,maximum(abs,c)),c)
end
_planeside(p::HyperbolicPlane,q::PolyhedralPoint) = _lorentz(_normalchart(p),_chartvector(q))
_samevertex(p,q,atol) = iszero(atol) ? p == q :
    all(abs.(_chartvector(p)-_chartvector(q)) .<= atol)

"""
    HyperbolicPolygon(vertices; atol=nothing)

A convex geodesic polygon with finite or ideal vertices in cyclic order.
Rejects repeated vertices, noncoplanar input, collinear corners, and nonconvex
or self-crossing orderings. Exact inputs use exact predicates by default;
floating input uses 128 machine epsilons on normalized chart coordinates.
"""
struct HyperbolicPolygon{V<:Tuple,P<:HyperbolicPlane}
    vertices::V
    plane::P
end
function HyperbolicPolygon(vertices;atol=nothing)
    vs = Tuple(vertices)
    length(vs) >= 3 && all(p -> p isa PolyhedralPoint,vs) ||
        throw(ArgumentError("a polygon needs at least three finite/ideal vertices"))
    tol = _poly_atol(vs,atol)
    cs = map(_chartvector,vs)
    n = dual(clvector(cs[1]) ∧ clvector(cs[2]) ∧ clvector(cs[3]))
    iszero(n) && throw(ArgumentError("polygon's first three vertices are collinear"))
    plane = HyperbolicPlane(n)
    all(p -> abs(_planeside(plane,p)) <= tol,vs) ||
        throw(ArgumentError("polygon vertices must be coplanar"))
    normal = _spatial(_normalchart(plane))
    ks = map(_spatial,cs)
    turn = _dot3(_cross3(ks[2]-ks[1],ks[3]-ks[1]),normal)
    sign = turn > 0 ? 1 : -1
    for i in eachindex(vs)
        j = mod1(i+1,length(vs))
        for k in eachindex(vs)
            k in (i,j) && continue
            side = sign*_dot3(_cross3(ks[j]-ks[i],ks[k]-ks[i]),normal)
            side > tol || throw(ArgumentError("vertices must form a strictly convex polygon in cyclic order"))
        end
    end
    HyperbolicPolygon(vs,plane)
end

"""Return the stored finite/ideal vertices (a tuple)."""
vertices(p::HyperbolicPolygon) = p.vertices
"""Return the supporting plane of a polygon, or an outward face plane of a polyhedron."""
faceplane(p::HyperbolicPolygon) = p.plane
"""Return vertex-index pairs for the edges."""
edges(p::HyperbolicPolygon) = [(i,mod1(i+1,length(vertices(p)))) for i in eachindex(vertices(p))]
"""A finite interior point obtained by averaging vertices in the Klein chart."""
interiorpoint(p::HyperbolicPolygon) = HyperbolicPoint(clvector(sum(_chartvector,vertices(p))))

function Base.in(q::PolyhedralPoint,p::HyperbolicPolygon)
    contains(p,q)
end
"""Test membership in the closed polygon/polyhedron, using normalized half-space residuals."""
function contains(p::HyperbolicPolygon,q::PolyhedralPoint;atol=nothing)
    tol = _poly_atol(vertices(p),atol)
    abs(_planeside(faceplane(p),q)) <= tol || return false
    ks = map(x -> _spatial(_chartvector(x)),vertices(p))
    x = _spatial(_chartvector(q))
    n = _spatial(_normalchart(faceplane(p)))
    s = _dot3(_cross3(ks[2]-ks[1],ks[3]-ks[1]),n) > 0 ? 1 : -1
    all(s*_dot3(_cross3(ks[j]-ks[i],x-ks[i]),n) >= -tol for (i,j) in edges(p))
end

"""Interior polygon angles, in radians; an ideal corner has angle zero."""
function interiorangles(p::HyperbolicPolygon;T=Float64)
    vs = vertices(p)
    map(eachindex(vs)) do i
        vs[i] isa IdealPoint && return zero(T)
        x = T.(_chartvector(vs[i])); x /= sqrt(-_lorentz(x,x))
        a,b = (T.(_chartvector(vs[mod1(j,length(vs))])) for j in (i-1,i+1))
        a += _lorentz(a,x)*x; b += _lorentz(b,x)*x
        c = _lorentz(a,b)/sqrt(_lorentz(a,a)*_lorentz(b,b))
        acos(clamp(c,-one(T),one(T)))
    end
end
"""Hyperbolic polygon area from the Gauss–Bonnet angle defect (curvature -1)."""
polygonarea(p::HyperbolicPolygon;T=Float64) = (length(vertices(p))-2)*T(pi)-sum(interiorangles(p;T))

abstract type AbstractHyperbolicPolyhedron end

"""
    HyperbolicPolyhedron(vertices, faces; atol=nothing)

A convex, closed polyhedron specified by finite/ideal vertices and cyclic
vertex-index lists for its faces. Faces are automatically oriented outwards.
The constructor checks planarity, convexity, supporting half-spaces, and a
closed spherical boundary. Nonconvex unions belong in `TetrahedralComplex`.
"""
struct HyperbolicPolyhedron{V<:Tuple,F<:Tuple,P<:Tuple} <: AbstractHyperbolicPolyhedron
    vertices::V
    faceindices::F
    planes::P
end
function HyperbolicPolyhedron(vertices,faces;atol=nothing)
    vs = Tuple(vertices)
    length(vs) >= 4 && all(p -> p isa PolyhedralPoint,vs) ||
        throw(ArgumentError("a polyhedron needs at least four finite/ideal vertices"))
    tol = _poly_atol(vs,atol)
    for i in eachindex(vs),j in 1:i-1
        _samevertex(vs[i],vs[j],tol) && throw(ArgumentError("polyhedron vertices must be distinct"))
    end
    center = HyperbolicPoint(clvector(sum(_chartvector,vs)))
    fs = Tuple[]; ps = HyperbolicPlane[]
    for indices in faces
        f = Tuple(indices)
        all(i -> i isa Integer && 1 <= i <= length(vs),f) || throw(ArgumentError("invalid face vertex index"))
        length(unique(f)) == length(f) || throw(ArgumentError("repeated face vertex"))
        polygon = HyperbolicPolygon(map(i -> vs[i],f);atol=tol)
        plane = faceplane(polygon)
        side = _planeside(plane,center)
        abs(side) > tol || throw(ArgumentError("polyhedron has no strict three-dimensional interior"))
        if side > 0
            f = reverse(f)
            plane = HyperbolicPlane(-representative(plane))
        end
        all(q -> _planeside(plane,q) <= tol,vs) || throw(ArgumentError("face is not a supporting face of a convex polyhedron"))
        push!(fs,f); push!(ps,plane)
    end
    length(fs) >= 4 || throw(ArgumentError("a polyhedron needs at least four faces"))
    incidence = Dict{Tuple{Int,Int},Vector{Tuple{Int,Int}}}()
    for f in fs,i in eachindex(f)
        a,b = f[i],f[mod1(i+1,length(f))]
        push!(get!(incidence,minmax(a,b),Tuple{Int,Int}[]),(a,b))
    end
    all(v -> length(v)==2 && v[1]==reverse(v[2]),values(incidence)) ||
        throw(ArgumentError("every edge must occur in two oppositely oriented faces"))
    length(vs)-length(incidence)+length(fs)==2 || throw(ArgumentError("boundary must be a sphere with every vertex used"))
    HyperbolicPolyhedron(vs,Tuple(fs),Tuple(ps))
end

"""
    HyperbolicTetrahedron(p, q, r, s; atol=nothing)
    HyperbolicTetrahedron(vertices; atol=nothing)

A finite/ideal tetrahedron. Face i is opposite vertex i, independent of the
input orientation. The four original vertices and coefficient field are kept.
"""
struct HyperbolicTetrahedron{P<:HyperbolicPolyhedron} <: AbstractHyperbolicPolyhedron
    polyhedron::P
    function HyperbolicTetrahedron(p::P) where {P<:HyperbolicPolyhedron}
        length(p.vertices)==4 && length(p.faceindices)==4 &&
            all(i -> sort(collect(p.faceindices[i]))==[j for j in 1:4 if j!=i],1:4) ||
            throw(ArgumentError("tetrahedral face i must be opposite vertex i"))
        new{P}(p)
    end
end
function HyperbolicTetrahedron(vs;atol=nothing)
    length(vs)==4 || throw(ArgumentError("a tetrahedron needs exactly four vertices"))
    HyperbolicTetrahedron(HyperbolicPolyhedron(vs,[Tuple(j for j in 1:4 if j!=i) for i in 1:4];atol))
end
HyperbolicTetrahedron(p::PolyhedralPoint,q::PolyhedralPoint,r::PolyhedralPoint,s::PolyhedralPoint;kwargs...) =
    HyperbolicTetrahedron((p,q,r,s);kwargs...)
_polyhedron(p::HyperbolicPolyhedron) = p
_polyhedron(p::HyperbolicTetrahedron) = p.polyhedron
vertices(p::AbstractHyperbolicPolyhedron) = _polyhedron(p).vertices
"""Cyclic vertex-index tuples, oriented outwards. Tetrahedral face i excludes vertex i."""
faceindices(p::AbstractHyperbolicPolyhedron) = _polyhedron(p).faceindices
faceplane(p::AbstractHyperbolicPolyhedron,i::Integer) = _polyhedron(p).planes[i]
"""Return face i as a geodesic polygon."""
face(p::AbstractHyperbolicPolyhedron,i::Integer) = HyperbolicPolygon(map(j -> vertices(p)[j],faceindices(p)[i]))
function edges(p::AbstractHyperbolicPolyhedron)
    sort!(unique!([minmax(f[i],f[mod1(i+1,length(f))]) for f in faceindices(p) for i in eachindex(f)]))
end
interiorpoint(p::AbstractHyperbolicPolyhedron) = HyperbolicPoint(clvector(sum(_chartvector,vertices(p))))
function contains(p::AbstractHyperbolicPolyhedron,q::PolyhedralPoint;atol=nothing)
    tol = _poly_atol(vertices(p),atol)
    all(i -> _planeside(faceplane(p,i),q) <= tol,eachindex(faceindices(p)))
end
Base.in(q::PolyhedralPoint,p::AbstractHyperbolicPolyhedron) = contains(p,q)
"""Finite edge length; any edge with an ideal endpoint has infinite length."""
function edgelength(p::Union{HyperbolicPolygon,AbstractHyperbolicPolyhedron},edge::Tuple{Int,Int})
    any(e -> minmax(e...)==minmax(edge...),edges(p)) || throw(ArgumentError("not an edge"))
    a,b = vertices(p)[edge[1]],vertices(p)[edge[2]]
    a isa IdealPoint || b isa IdealPoint ? Inf : distance(a,b)
end

"""Gram matrix of outward unit spacelike face normals; exact when the backend supports the square roots."""
function grammatrix(p::AbstractHyperbolicPolyhedron)
    ns = map(i -> _normalchart(faceplane(p,i)),eachindex(faceindices(p)))
    [_fielddiv(_lorentz(a,b),_sqrt_preserving(_lorentz(a,a)*_lorentz(b,b))) for a in ns,b in ns]
end
"""Cosine of the interior dihedral angle at a vertex-index edge."""
function dihedralcos(p::AbstractHyperbolicPolyhedron,edge::Tuple{Int,Int})
    fs = findall(f -> all(i -> i in f,edge),faceindices(p))
    length(fs)==2 && edge[1]!=edge[2] || throw(ArgumentError("not an edge"))
    a,b = (_normalchart(faceplane(p,i)) for i in fs)
    -_fielddiv(_lorentz(a,b),_sqrt_preserving(_lorentz(a,a)*_lorentz(b,b)))
end
"""Interior dihedral angle in radians, evaluated in floating type T (default Float64)."""
dihedralangle(p::AbstractHyperbolicPolyhedron,edge::Tuple{Int,Int};T=Float64) =
    acos(clamp(T(dihedralcos(p,edge)),-one(T),one(T)))

function (g::AbstractIsometry)(p::HyperbolicPolygon)
    HyperbolicPolygon(map(g,vertices(p)))
end
function (g::AbstractIsometry)(p::HyperbolicPolyhedron)
    HyperbolicPolyhedron(map(g,vertices(p)),faceindices(p))
end
(g::AbstractIsometry)(p::HyperbolicTetrahedron) = HyperbolicTetrahedron(map(g,vertices(p)))
mapcoefficients(f,p::HyperbolicPolygon) = HyperbolicPolygon(map(x -> mapcoefficients(f,x),vertices(p)))
mapcoefficients(f,p::HyperbolicPolyhedron) = HyperbolicPolyhedron(map(x -> mapcoefficients(f,x),vertices(p)),faceindices(p))
mapcoefficients(f,p::HyperbolicTetrahedron) = HyperbolicTetrahedron(map(x -> mapcoefficients(f,x),vertices(p)))

"""Ordered ideal tetrahedron shape; its sign records the vertex orientation."""
function idealshape(t::HyperbolicTetrahedron)
    all(p -> p isa IdealPoint,vertices(t)) || throw(ArgumentError("shape requires four ideal vertices"))
    crossratio(vertices(t)...)
end
"""Ideal edge shape with positive imaginary part; its argument is the interior dihedral angle."""
function edgeparameter(t::HyperbolicTetrahedron,edge::Tuple{Int,Int})
    edge in edges(t) || reverse(edge) in edges(t) || throw(ArgumentError("not an edge"))
    all(p -> p isa IdealPoint,vertices(t)) || throw(ArgumentError("edge parameter requires an ideal tetrahedron"))
    rest = Tuple(i for i in 1:4 if !(i in edge))
    z = crossratio(map(i -> vertices(t)[i],(edge...,rest...))...)
    # Swapping the other two vertices inverts the cross-ratio. Conjugating
    # instead would keep the wrong modulus and break complex edge equations.
    imag(z) > 0 ? z : _fielddiv(one(z),z)
end
