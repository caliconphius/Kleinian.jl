# Optional integration suite: run with the test/plotting environment.
using Test
using Random
using Kleinian, Poincare, Riley
using CairoMakie
import CairoMakie: Makie
import GeometryBasics as GB

@testset "Makie geometry adapters" begin
    @test Base.get_extension(Kleinian,:KleinianMakieExt) !== nothing
    origin = HyperbolicPoint(1,0,0,0)
    q = HyperbolicPoint(5,1,2,3)
    geodesics = [Geodesic(origin,q), Geodesic(IdealPoint(1,0,0,-1),IdealPoint(1,0,0,1)),
                 Geodesic(q,HyperbolicPoint(3,-1,1,0))]
    planes = [HyperbolicPlane(clvector(0,0,1,0)),HyperbolicPlane(clvector(1,2,1,0))]
    horospheres = [Horosphere(clvector(1,0,0,1)),Horosphere(clvector(2,2,0,0))]
    for m in (:klein,:ball,:halfspace)
        @test plotgeometry(q;model=m) ≈ modelcoordinates(q;model=m)
        for g in geodesics
            vertices = plotgeometry(g;model=m,resolution=32)
            @test length(vertices) == 32
            @test all(v -> isincident(modelpoint(v;model=m),g;atol=1e-9),vertices)
            @test plotgeometry(reverse(g);model=m,resolution=32) ≈ reverse(vertices)
        end
        for object in [planes;horospheres]
            geometry = plotgeometry(object;model=m,resolution=24)
            @test geometry isa GB.Mesh
            @test !isempty(GB.faces(geometry))
            @test all(v -> all(isfinite,v),GB.coordinates(geometry))
            @test all(v -> isincident(modelpoint(v;model=m),object;atol=1e-7),GB.coordinates(geometry))
        end
        f = Figure(size=(600,450))
        ax = Axis3(f[1,1];aspect=:data)
        @test hyperbolicplot!(ax,origin;model=m,color=:black) isa Makie.Scatter
        @test hyperbolicplot!(ax,geodesics[1];model=m,color=:red) isa Makie.Lines
        @test hyperbolicplot!(ax,planes[1];model=m,resolution=16,color=(:blue,0.2)) isa Makie.Mesh
        @test hyperbolicplot!(ax,horospheres[1];model=m,resolution=16,color=(:orange,0.3)) isa Makie.Mesh
        @test modelboundary!(ax;model=m,resolution=16) isa Makie.Plot
        @test hyperbolicplot!(ax,[origin,q];model=m) isa Makie.Scatter
        @test length(hyperbolicplot!(ax,[origin,geodesics[1]];model=m)) == 2
        mktempdir() do dir
            path = joinpath(dir,"$m.png")
            save(path,f)
            @test filesize(path) > 1000
        end
    end
    # -X⋅(a,0,0,a)=a/height=1, so this is the horizontal plane at height a.
    h = Horosphere(clvector(2,0,0,2))
    @test all(v -> v[3] ≈ 2,GB.coordinates(plotgeometry(h;model=:halfspace)))
    # Klein geodesics are straight; a central ball geodesic is a diameter.
    g = Geodesic(origin,HyperbolicPoint(2,1,0,0))
    @test all(v -> iszero(v[2]) && iszero(v[3]),plotgeometry(g;model=:ball))
    # Independent Euclidean shape equations in conformal models.
    arc = Geodesic(IdealPoint(1,1,0,0),IdealPoint(1,0,1,0))
    @test all(v -> sum(abs2,v-GB.Point3d(1,1,0)) ≈ 1,plotgeometry(arc;model=:ball))
    cap = plotgeometry(HyperbolicPlane(clvector(1,2,0,0));model=:ball)
    @test all(v -> sum(abs2,v-GB.Point3d(2,0,0)) ≈ 3,GB.coordinates(cap))
    hemisphere = plotgeometry(HyperbolicPlane(clvector(0,0,0,1));model=:halfspace)
    @test all(v -> sum(abs2,v) ≈ 1,GB.coordinates(hemisphere))
    sphere = plotgeometry(horospheres[2];model=:ball)
    @test all(v -> sum(abs2,v-GB.Point3d(2/3,0,0)) ≈ 1/9,GB.coordinates(sphere))
    @test_throws DomainError plotgeometry(IdealPoint(1,0,0,1);model=:halfspace)
    @test_throws ArgumentError plotgeometry(origin;resolution=3)
    @test_throws ArgumentError plotgeometry(origin;cutoff=0)
    @test_throws ArgumentError plotgeometry(origin;cutoff=NaN)
    @test_throws ArgumentError plotgeometry(origin;model=:typo)
    result = hyperbolicplot(origin;boundary=false,figure=(size=(500,400),),axis=(title="Origin",))
    @test result.figure isa Figure
    @test result.axis isa Axis3
    @test result.plot isa Makie.Scatter
end

include("makie_polyhedra.jl")

@testset "Complex limit-set rendering" begin
    pair(a,b) = MobiusIsometry(conj.([b -a*b-1; 1 -a]))
    points = limitset([pair(-2.0,2.0),pair(-2.0im,2.0im)];depth=6)
    push!(points,complexidealpoint(nothing))
    for renderer in (:scatter,:datashader)
        view = limitsetplot(points;renderer,figure=(size=(600,600),))
        @test view.axis isa Makie.Axis
        @test view.plot isa (renderer === :scatter ? Makie.Scatter : Makie.DataShader)
        mktempdir() do dir
            path = joinpath(dir,"limitset.png")
            save(path,view.figure)
            @test filesize(path) > 1000
        end
    end
    @test_throws ArgumentError limitsetplot(points;renderer=:invalid)
    @test_throws ArgumentError limitsetplot([complexidealpoint(nothing)])
end

@testset "Streaming limit-set rendering" begin
    gs = rileygroup(2.0im)
    for stream in (limitsetstream(gs;depth=5),
                   randomlimitset(gs;samples=2000,rng=Xoshiro(25)))
        view = limitsetplot(stream;axis=(limits=(-1.5,1.5,-1.6,0.6),))
        @test view.plot isa Makie.DataShader
        mktempdir() do dir
            path = joinpath(dir,"stream.png")
            save(path,view.figure)
            @test filesize(path) > 1000
        end
    end
end

include("makie_hulls.jl")
include("makie_fractal_viewer.jl")
include("makie_hull_models.jl")
