using Test, Kleinian, CairoMakie
import CairoMakie: Makie
import GeometryBasics as GB
include(joinpath(@__DIR__,"..","examples","trefoil333_data.jl"))

@testset "Makie polygon and tetrahedral domain adapters" begin
    d = trefoil333().domain
    tetrahedron = d.cells[1]
    polygon = face(tetrahedron,1)
    segment = GeodesicSegment(HyperbolicPoint(1,0,0,0),HyperbolicPoint(2,1,0,0))
    for model in (:klein,:ball,:halfspace)
        curve = plotgeometry(segment;model,resolution=12)
        @test first(curve)≈modelcoordinates(segment.start;model)
        @test last(curve)≈modelcoordinates(segment.stop;model)
        mesh = plotgeometry(polygon;model,resolution=12)
        @test mesh isa GB.Mesh
        @test !isempty(GB.faces(mesh))
        @test all(p -> all(isfinite,p),GB.coordinates(mesh))
        @test all(p -> contains(polygon,modelpoint(p;model);atol=1e-8),GB.coordinates(mesh))
        for object in (tetrahedron,d)
            geometry = plotgeometry(object;model,resolution=12)
            @test !isempty(GB.faces(geometry))
            @test all(p -> all(isfinite,p),GB.coordinates(geometry))
        end
        fig = Figure(size=(500,450))
        ax = Axis3(fig[1,1];aspect=:data)
        @test hyperbolicplot!(ax,d;model,resolution=12,color=(:steelblue,.3)) isa Makie.Mesh
        @test hyperbolicplot!(ax,segment;model) isa Makie.Lines
        mktempdir() do dir
            save(joinpath(dir,"domain.png"),fig)
            @test filesize(joinpath(dir,"domain.png"))>1000
        end
    end
    square = HyperbolicPolygon([modelpoint((x,y,0);model=:klein) for (x,y) in ((-.2,-.2),(.2,-.2),(.2,.2),(-.2,.2))])
    @test !isempty(GB.faces(plotgeometry(square;resolution=8)))
    @test hyperbolicplot(tetrahedron;resolution=8,boundary=false).plot isa Makie.Mesh
    @test_throws ArgumentError plotgeometry(segment;cutoff=NaN)
end

include(joinpath(@__DIR__,"..","examples","trefoil333.jl"))
@testset "Makie trefoil model selector" begin
    view = trefoil333_explorer(resolution=8)
    @test view.menu.selection[]==:ball
    view.menu.selection[]=:klein
    @test !isempty(view.axis.scene.plots)
    view.menu.selection[]=:halfspace
    @test !isempty(view.axis.scene.plots)
    view.menu.selection[]=:ball
    @test !isempty(view.axis.scene.plots)
end
