module KleinianMakieExt

import Kleinian as K
import Makie
import GeometryBasics as GB
using LinearAlgebra: norm, cross, dot
using StaticArrays: SVector

const PointObject = Union{K.HyperbolicPoint,K.IdealPoint}
const SurfaceObject = Union{K.HyperbolicPlane,K.Horosphere}
const GeometryObject = Union{PointObject,K.Geodesic,SurfaceObject}

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

function K.plotgeometry(p::PointObject;model=:ball,resolution=64,cutoff=1e-3)
    _settings(model,resolution,cutoff)
    v = GB.Point3d(K.modelcoordinates(p;model))
    all(isfinite,v) || throw(DomainError(v,"point is outside rendering precision"))
    v
end

function K.plotgeometry(g::K.Geodesic;model=:ball,resolution=64,cutoff=1e-3)
    _settings(model,resolution,cutoff)
    L = K.representative(g)
    coeffs = Float64.(real.(K.coefficients(L)))
    all(isfinite,coeffs) || throw(DomainError(coeffs,"rendering needs finite Float64 coefficients"))
    coeffs /= maximum(abs,coeffs)
    # L=(1,c)∧(0,d), c⋅d=0: recover the closest Klein point and direction.
    d = SVector(ntuple(i -> coeffs[(1 | (1 << i))+1],3))
    A = [i == j ? 0.0 : (i < j ? 1.0 : -1.0)*coeffs[((1 << i) | (1 << j))+1]
         for i in 1:3,j in 1:3]
    c = A*d/dot(d,d)
    d /= norm(d)
    radius = sqrt(max(0.0,1-dot(c,c)))
    radius > 0 || throw(DomainError(g,"geodesic lost at rendering precision"))
    [_mapklein(c+(s*radius)*d,model) for s in range(-1+cutoff,1-cutoff;length=resolution)]
end

function K.plotgeometry(p::K.HyperbolicPlane;model=:ball,resolution=64,cutoff=1e-3)
    _settings(model,resolution,cutoff)
    c = _floatcoords(K.representative(p))
    c /= maximum(abs,c)
    n = SVector(c[2],c[3],c[4])
    h = c[1]/norm(n)
    n /= norm(n)
    center = h*n
    radius = sqrt(max(0.0,1-h*h))*(1-cutoff)
    radius > 0 || throw(DomainError(p,"plane lost at rendering precision"))
    u,v = _frame(n)
    # GeometryBasics supplies tessellation and triangle connectivity. Map its
    # square to a disk before lifting into the Klein plane and changing chart.
    grid = GB.triangle_mesh(GB.Tessellation(GB.Rect2d(-1,-1,2,2),(resolution,resolution));pointtype=GB.Point2d)
    vertices = map(GB.coordinates(grid)) do q
        a,b = q
        x,y = a*sqrt(1-b*b/2), b*sqrt(1-a*a/2)
        _mapklein(center+radius*(x*u+y*v),model)
    end
    GB.Mesh(vertices,GB.faces(grid))
end

function K.plotgeometry(h::K.Horosphere;model=:ball,resolution=64,cutoff=1e-3)
    _settings(model,resolution,cutoff)
    c = _floatcoords(K.representative(h))
    # In the ball: center = h_spatial/(h_time+1), radius = 1/(h_time+1).
    center = GB.Point3d(c[2:4]/(c[1]+1))
    radius = 1/(c[1]+1)
    sphere = GB.triangle_mesh(GB.Tessellation(GB.Sphere(center,radius),resolution);pointtype=GB.Point3d)
    vertices = GB.coordinates(sphere)
    inside = [sum(abs2,v) < (1-cutoff)^2 for v in vertices]
    faces = filter(f -> all(i -> inside[i],f),GB.faces(sphere))
    isempty(faces) && throw(ArgumentError("horosphere mesh is empty; reduce cutoff or increase resolution"))
    # Compact the retained GeometryBasics mesh before chart conversion; the
    # omitted cap contains its ideal tangency, where the half-space chart may diverge.
    used = sort!(unique!(Int[i for f in faces for i in f]))
    indices = zeros(Int,length(vertices))
    indices[used] = eachindex(used)
    mapped = [_mapball(vertices[i],model) for i in used]
    mappedfaces = [GB.GLTriangleFace(ntuple(j -> indices[f[j]],3)) for f in faces]
    GB.Mesh(mapped,mappedfaces)
end

function K.hyperbolicplot!(axis,p::PointObject;model=:ball,resolution=64,cutoff=1e-3,kwargs...)
    Makie.scatter!(axis,[K.plotgeometry(p;model,resolution,cutoff)];kwargs...)
end
function K.hyperbolicplot!(axis,points::AbstractVector{<:PointObject};
                          model=:ball,resolution=64,cutoff=1e-3,kwargs...)
    _settings(model,resolution,cutoff)
    vertices = GB.Point3d[K.plotgeometry(p;model,resolution,cutoff) for p in points]
    Makie.scatter!(axis,vertices;kwargs...)
end
function K.hyperbolicplot!(axis,g::K.Geodesic;model=:ball,resolution=64,cutoff=1e-3,kwargs...)
    Makie.lines!(axis,K.plotgeometry(g;model,resolution,cutoff);kwargs...)
end
function K.hyperbolicplot!(axis,p::SurfaceObject;model=:ball,resolution=64,cutoff=1e-3,kwargs...)
    Makie.mesh!(axis,K.plotgeometry(p;model,resolution,cutoff);kwargs...)
end
function K.hyperbolicplot!(axis,objects::AbstractVector;kwargs...)
    [K.hyperbolicplot!(axis,object;kwargs...) for object in objects]
end

function K.modelboundary!(axis;model=:ball,resolution=64,extent=3,kwargs...)
    _settings(model,resolution,1e-3)
    isfinite(extent) && extent > 0 || throw(ArgumentError("extent must be finite and positive"))
    if model === :halfspace
        grid = GB.triangle_mesh(GB.Rect2d(-extent,-extent,2extent,2extent))
        vertices = [GB.Point3d(p[1],p[2],0) for p in GB.coordinates(grid)]
        Makie.mesh!(axis,GB.Mesh(vertices,GB.faces(grid));color=(:gray,0.12),kwargs...)
    else
        sphere = GB.Tessellation(GB.Sphere(GB.Point3d(0),1.0),resolution)
        Makie.wireframe!(axis,sphere;color=(:gray,0.15),kwargs...)
    end
end

function K.hyperbolicplot(object::Union{GeometryObject,AbstractVector};
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
