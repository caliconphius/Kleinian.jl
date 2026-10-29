module PoincareMakieExt
import Kleinian as K
import Poincare as P
import Makie
import GeometryBasics as GB
using LinearAlgebra: norm, cross, dot
using StaticArrays: SVector
const SurfaceObject = Union{P.HyperbolicPolygon,P.AbstractHyperbolicPolyhedron,P.TetrahedralComplex}
function _settings(model,resolution,cutoff)
    K._checkmodel(model)
    resolution isa Integer && resolution >= 8 ||
        throw(ArgumentError("resolution must be an integer ≥ 8"))
    isfinite(cutoff) && 0 < cutoff < 1 ||
        throw(ArgumentError("cutoff must be finite and between zero and one"))
end

function _floatcoords(v)
    c = map(x -> Float64(real(x)),K._pointcoordinates(v))
    all(isfinite,c) || throw(DomainError(c,"rendering needs finite Float64 coefficients"))
    c
end

function _frame(n)
    # Choose the coordinate axis least parallel to n for a stable tangent frame.
    a = SVector(ntuple(i -> i == argmin(abs.(n)) ? 1.0 : 0.0,3))
    u = cross(n,a)
    u /= norm(u)
    u,cross(n,u)
end

_mapklein(v,model) = GB.Point3d(K.modelcoordinates(K.modelpoint(v;model=:klein);model))
_mapball(v,model) = GB.Point3d(K.modelcoordinates(K.modelpoint(v;model=:ball);model))

function K.plotgeometry(s::P.GeodesicSegment;model=:ball,resolution=64,cutoff=1e-3)
    _settings(model,resolution,cutoff)
    a,b = map(p -> SVector{3,Float64}(K.modelcoordinates(p;model=:klein)),P.vertices(s))
    lo = s.start isa K.IdealPoint ? cutoff : 0.0
    hi = s.stop isa K.IdealPoint ? 1-cutoff : 1.0
    lo < hi || throw(ArgumentError("cutoff removes the whole segment"))
    [_mapklein((1-t)*a+t*b,model) for t in range(lo,hi;length=resolution)]
end

# Keep GeometryBasics' triangulation/connectivity. Only coordinates and the
# small cutoff cap at ideal vertices are specific to hyperbolic geometry.
function _map_patch(vertices,faces,model,cutoff)
    valid = [sum(abs2,p)<(1-cutoff)^2 for p in vertices]
    kept = filter(f -> all(i -> valid[i],f) &&
                  norm(cross(vertices[f[2]]-vertices[f[1]],vertices[f[3]]-vertices[f[1]]))>eps(Float64)^2,faces)
    isempty(kept) && throw(ArgumentError("polygon patch is empty; reduce cutoff or increase resolution"))
    used = sort!(unique!(Int[i for f in kept for i in f]))
    indices = zeros(Int,length(vertices)); indices[used] = eachindex(used)
    fs = [GB.GLTriangleFace(ntuple(j -> indices[f[j]],3)) for f in kept]
    GB.Mesh([_mapklein(vertices[i],model) for i in used],fs)
end

function K.plotgeometry(p::P.HyperbolicPolygon;model=:ball,resolution=32,cutoff=1e-4)
    _settings(model,resolution,cutoff)
    ks = [SVector{3,Float64}(K.modelcoordinates(v;model=:klein)) for v in P.vertices(p)]
    n = _floatcoords(K.representative(P.faceplane(p)))[2:4]; n /= norm(n)
    u,v = _frame(n); origin = first(ks)
    outline = [GB.Point2d(dot(k-origin,u),dot(k-origin,v)) for k in ks]
    triangulation = GB.triangle_mesh(GB.Polygon(outline);pointtype=GB.Point2d)
    lifted = [origin+q[1]*u+q[2]*v for q in GB.coordinates(triangulation)]
    grid = GB.triangle_mesh(GB.Tessellation(GB.Rect2d(0,0,1,1),(resolution,resolution));pointtype=GB.Point2d)
    patches = map(GB.faces(triangulation)) do f
        a,b,c = (lifted[i] for i in f)
        points = [(1-q[1])*a+q[1]*((1-q[2])*b+q[2]*c) for q in GB.coordinates(grid)]
        _map_patch(points,GB.faces(grid),model,cutoff)
    end
    merge(patches)
end
function K.plotgeometry(p::P.AbstractHyperbolicPolyhedron;model=:ball,resolution=32,cutoff=1e-4)
    merge([K.plotgeometry(P.face(p,i);model,resolution,cutoff) for i in eachindex(P.faceindices(p))])
end
function K.plotgeometry(d::P.TetrahedralComplex;model=:ball,resolution=32,cutoff=1e-4)
    merge([K.plotgeometry(P.face(d.cells[c],f);model,resolution,cutoff) for (c,f) in P.boundaryfaces(d)])
end

# Triangular barycentric grid, shared by all finite face patches.
function _hullgrid(n)
    weights = Tuple{Float64,Float64}[]
    ids = Dict{Tuple{Int,Int},Int}()
    for i in 0:n, j in 0:n-i
        push!(weights,(i/n,j/n)); ids[(i,j)] = length(weights)
    end
    faces = GB.GLTriangleFace[]
    for i in 0:n-1, j in 0:n-1-i
        push!(faces,GB.GLTriangleFace(ids[(i,j)],ids[(i+1,j)],ids[(i,j+1)]))
        i+j < n-1 && push!(faces,GB.GLTriangleFace(ids[(i+1,j)],ids[(i+1,j+1)],ids[(i,j+1)]))
    end
    weights,faces
end
function _hullsettings(model,resolution,cutoff)
    _settings(model,resolution,cutoff)
    cutoff < .5 || throw(ArgumentError("hull cutoff must be less than 1/2"))
end
function K.plotgeometry(h::P.LimitSetHull;model=:klein,resolution=16,cutoff=1e-3)
    _hullsettings(model,resolution,cutoff)
    h.dimension >= 2 || throw(ArgumentError("a hull of dimension < 2 has no surface mesh"))
    model === :klein && return GB.Mesh(GB.Point3d.(h.klein),[GB.GLTriangleFace(f) for f in P.faceindices(h)])
    weights,triangles = _hullgrid(resolution)
    points = GB.Point3d[]
    faces = GB.GLTriangleFace[]
    for f in P.faceindices(h)
        a,b,c = (h.klein[i] for i in f)
        center = (a+b+c)/3
        # Cut only ideal corners, within the original supporting plane.
        # Adjacent faces share the same trimmed edge and edge subdivision.
        rim = ((1-cutoff)*a+cutoff*b,(1-cutoff)*b+cutoff*a,
               (1-cutoff)*b+cutoff*c,(1-cutoff)*c+cutoff*b,
               (1-cutoff)*c+cutoff*a,(1-cutoff)*a+cutoff*c)
        for k in 1:6
            u,v = rim[k],rim[mod1(k+1,6)]
            offset = length(points)
            for (α,β) in weights
                q = _mapklein((1-α-β)*center+α*u+β*v,model)
                all(isfinite,q) || throw(DomainError(q,"hull chart overflow; increase cutoff"))
                push!(points,q)
            end
            append!(faces,[GB.GLTriangleFace(ntuple(i -> convert(Int,t[i])+offset,3)) for t in triangles])
        end
    end
    GB.Mesh(points,faces)
end
function P.hulledgesgeometry(h::P.LimitSetHull;model=:klein,resolution=64,cutoff=1e-3,
                              pleats=false,minangle=0,planetol=h.atol)
    _hullsettings(model,resolution,cutoff)
    result = GB.Point3d[]
    selected = pleats ? P.pleatedges(h;minangle,planetol) : P.edges(h)
    for (i,j) in selected
        a,b = h.klein[i],h.klein[j]
        if model === :klein
            append!(result,(GB.Point3d(a),GB.Point3d(b)))
        else
            for t in range(cutoff,1-cutoff;length=resolution+1)
                push!(result,_mapklein((1-t)*a+t*b,model))
            end
        end
        push!(result,GB.Point3d(NaN))
    end
    result
end
function K.hyperbolicplot!(axis,h::P.LimitSetHull;model=:klein,resolution=16,cutoff=1e-3,kwargs...)
    _hullsettings(model,resolution,cutoff)
    h.dimension == 0 && return Makie.scatter!(axis,[K.plotgeometry(only(P.vertices(h));model)];kwargs...)
    h.dimension == 1 && return Makie.lines!(axis,P.hulledgesgeometry(h;model,resolution,cutoff);kwargs...)
    Makie.mesh!(axis,K.plotgeometry(h;model,resolution,cutoff);kwargs...)
end
function K.hyperbolicplot(h::P.LimitSetHull;model=:klein,resolution=16,cutoff=1e-3,
                          boundary=true,figure=(;),axis=(;),kwargs...)
    _hullsettings(model,resolution,cutoff)
    fig = Makie.Figure(;figure...)
    opts = merge((;aspect=:data,xlabel="x",ylabel="y",zlabel="z"),axis)
    ax = Makie.Axis3(fig[1,1];opts...)
    boundary && K.modelboundary!(ax;model,resolution)
    plot = K.hyperbolicplot!(ax,h;model,resolution,cutoff,kwargs...)
    (figure=fig,axis=ax,plot=plot)
end

function K.hyperbolicplot!(axis,s::P.GeodesicSegment;model=:ball,resolution=64,cutoff=1e-3,kwargs...)
    Makie.lines!(axis,K.plotgeometry(s;model,resolution,cutoff);kwargs...)
end
function K.hyperbolicplot!(axis,p::SurfaceObject;model=:ball,resolution=64,cutoff=1e-3,kwargs...)
    Makie.mesh!(axis,K.plotgeometry(p;model,resolution,cutoff);kwargs...)
end
function K.hyperbolicplot(object::Union{P.GeodesicSegment,SurfaceObject};
                          model=:ball,resolution=64,cutoff=1e-3,boundary=true,
                          figure=(;),axis=(;),kwargs...)
    _settings(model,resolution,cutoff)
    fig = Makie.Figure(;figure...)
    ax = Makie.Axis3(fig[1,1];aspect=:data,xlabel="x",ylabel="y",zlabel="z",axis...)
    boundary && K.modelboundary!(ax;model,resolution)
    p = K.hyperbolicplot!(ax,object;model,resolution,cutoff,kwargs...)
    (figure=fig,axis=ax,plot=p)
end


end
