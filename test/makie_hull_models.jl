@testset "Curved hull surfaces" begin
    # An octahedron includes the half-space point at infinity and faces on
    # independently known planes ±x ±y ±z = 1 in Klein coordinates.
    points = [IdealPoint(1,(i==1 ? s : 0),(i==2 ? s : 0),(i==3 ? s : 0)) for i in 1:3 for s in (-1,1)]
    h = limitsethull(points)
    for model in (:ball,:halfspace)
        mesh = plotgeometry(h;model,resolution=8,cutoff=.02)
        @test all(p -> all(isfinite,p),GB.coordinates(mesh))
        for p in GB.coordinates(mesh)
            k = modelcoordinates(modelpoint(p;model);model=:klein)
            @test sum(abs,k) ≈ 1 atol=1e-10
        end
        model === :ball && @test all(p -> sum(abs2,p)<1,GB.coordinates(mesh))
        model === :halfspace && @test all(p -> p[3]>0,GB.coordinates(mesh))
        edgepoints = hulledgesgeometry(h;model,resolution=8,cutoff=.02)
        @test count(p -> all(isnan,p),edgepoints) == length(edges(h))
        @test all(p -> all(isnan,p) || all(isfinite,p),edgepoints)
        # All trimmed edge samples also lie on the original hull boundary.
        for p in filter(p -> all(isfinite,p),edgepoints)
            k = modelcoordinates(modelpoint(p;model);model=:klein)
            @test sum(abs,k) ≈ 1 atol=1e-10
        end
        planar = limitsethull(points[1:4])
        @test planar.dimension == 2
        planar_mesh = plotgeometry(planar;model,resolution=8)
        @test all(p -> abs(modelcoordinates(modelpoint(p;model);model=:klein)[3])<1e-10,GB.coordinates(planar_mesh))
        @test hyperbolicplot(h;model,resolution=8,boundary=false).plot isa Makie.Mesh
    end
    for model in (:klein,:ball,:halfspace)
        @test isempty(hulledgesgeometry(h;model,pleats=true,minangle=pi))
        @test count(p -> all(isnan,p),hulledgesgeometry(h;model,pleats=true)) == length(pleatedges(h))
    end
    @test_throws ArgumentError plotgeometry(h;cutoff=.5)
    @test_throws ArgumentError hulledgesgeometry(h;model=:unknown)
    @test_throws DomainError hyperbolicplot(limitsethull([complexidealpoint(nothing)]);model=:halfspace,boundary=false)
end
