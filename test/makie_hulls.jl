@testset "Sampled hull rendering and Riley parameter explorer" begin
    tetra = [IdealPoint(3,x,y,z) for (x,y,z) in ((1,2,2),(1,-2,-2),(-1,2,-2),(-1,-2,2))]
    h = limitsethull(tetra)
    @test length(GB.faces(plotgeometry(h))) == 4
    @test length(GB.faces(plotgeometry(h;model=:ball,resolution=8))) > 4
    @test hyperbolicplot(h;boundary=false).plot isa Makie.Mesh
    @test hyperbolicplot(limitsethull(tetra[1:1]);boundary=false).plot isa Makie.Scatter
    @test hyperbolicplot(limitsethull(tetra[1:2]);boundary=false).plot isa Makie.Lines
    include(joinpath(@__DIR__,"..","examples","riley_hull_explorer.jl"))
    view = riley_hull_explorer(samples=1000,hullsamples=32)
    @test_throws ArgumentError view.set_parameter!(-1+2im)
    @test_throws ArgumentError view.set_parameter!(1-2im)
    @test view.parameter[] == 2im
    @test view.parameteraxis.xpanlock[] && view.parameteraxis.ypanlock[]
    @test view.set_parameter!(4.0)
    @test view.state[].hull.dimension == 2
    @test view.set_parameter!(2.2im;preview=true)
    @test view.state[].hull.dimension == 3
    @test !view.set_parameter!(0)
    @test !view.valid[]
    @test view.set_parameter!(2im)
    for model in (:ball,:halfspace,:klein)
        view.modelmenu.selection[] = model
        @test view.valid[]
        @test all(p -> all(isfinite,p),GB.coordinates(view.plots[3][1][]))
    end
    Makie.set_close_to!(view.slider,64)
    @test view.state[].hull.samplecount == 64
    view.faces.active[] = false
    view.wire.active[] = true
    @test !view.plots[3].visible[]
    @test view.plots[4].visible[]
    before = view.state[]
    Makie.set_close_to!(view.bendslider,180)
    @test isempty(view.plots[4][1][])
    @test view.state[] === before
    Makie.set_close_to!(view.bendslider,0)
    @test !isempty(view.plots[4][1][])
    mktempdir() do dir
        path = joinpath(dir,"riley-hull.png")
        save(path,view.figure)
        @test filesize(path) > 1000
    end
    # Exercise the actual mouse callbacks in data-to-screen coordinates.
    ev = events(view.figure)
    function moveto(z)
        vp = Makie.viewport(view.parameteraxis.scene)[]
        lim = view.parameteraxis.finallimits[]
        lo,wh = minimum(lim),GB.widths(lim)
        fraction = ((real(z)-lo[1])/wh[1],(imag(z)-lo[2])/wh[2])
        ev.mouseposition[] = Tuple(minimum(vp) + GB.Point2d(GB.widths(vp) .* fraction))
    end
    moveto(2.1im)
    ev.mousebutton[] = Makie.MouseButtonEvent(Mouse.left,Mouse.press)
    @test view.parameter[] ≈ 2.1im atol=1e-5
    @test occursin("preview",view.status[])
    sleep(.15) # exceed the intentional drag throttle
    moveto(0.5+2.0im)
    @test view.parameter[] ≈ 0.5+2.0im atol=1e-5
    ev.mousebutton[] = Makie.MouseButtonEvent(Mouse.left,Mouse.release)
    @test view.valid[]
    @test occursin("refined",view.status[])
    @test view.state[].hull.samplecount == 64

end
