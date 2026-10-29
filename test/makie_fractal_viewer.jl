@testset "Fractal-only Riley viewer" begin
    include(joinpath(@__DIR__,"..","examples","riley_fractal_viewer.jl"))
    v = riley_fractal_viewer(samples=2000,previewpoints=100,background=false)
    @test length(v.points[]) == 2000
    @test v.valid[]
    @test_throws ArgumentError v.set_parameter!(-1+2im)
    @test_throws ArgumentError v.set_parameter!(1-2im)
    @test v.parameter[] == 2im
    @test v.parameteraxis.xpanlock[] && v.parameteraxis.ypanlock[]
    @test v.set_parameter!(2.2im;preview=true)
    @test v.rendered[].samples == 100
    @test v.rendered[].preview
    @test v.set_parameter!(2im)
    reference = copy(v.points[])
    @test v.set_parameter!(2im)
    @test v.points[] == reference
    @test !v.set_parameter!(0)
    @test !v.valid[]
    @test isempty(v.points[])
    @test v.set_parameter!(2im)
    v.countmenu.selection[] = 100_000
    @test v.rendered[].samples == 100_000
    v.colors.selection[] = :magma
    @test v.colors.selection[] == :magma
    limits!(v.axis,-.2,.2,-1.2,-.8)
    v.reset.clicks[] += 1
    @test v.axis.limits[] == ((-1.1,1.1),(-2.1,.1))
    mktempdir() do dir
        path = joinpath(dir,"fractal.png")
        save(path,v.figure)
        @test filesize(path) > 1000
    end
    ev = events(v.figure)
    function moveto(z)
        vp = Makie.viewport(v.parameteraxis.scene)[]
        lim = v.parameteraxis.finallimits[]
        lo,wh = minimum(lim),GB.widths(lim)
        fraction = ((real(z)-lo[1])/wh[1],(imag(z)-lo[2])/wh[2])
        ev.mouseposition[] = Tuple(minimum(vp)+GB.Point2d(GB.widths(vp).*fraction))
    end
    moveto(2.1im)
    ev.mousebutton[] = Makie.MouseButtonEvent(Mouse.left,Mouse.press)
    @test v.parameter[] ≈ 2.1im atol=1e-5
    @test v.rendered[].preview
    moveto(.5+2im)
    ev.mousebutton[] = Makie.MouseButtonEvent(Mouse.left,Mouse.release)
    @test v.parameter[] ≈ .5+2im atol=1e-5
    @test !v.rendered[].preview
    @test v.rendered[].samples == 100_000
    # A superseded job must not overwrite the current parameter or report success.
    asyncview = riley_fractal_viewer(samples=100_000,previewpoints=100,background=true)
    old = asyncview.task[]
    yield() # let the old request start, then cancel at a batch boundary
    latest = asyncview.set_parameter!(2.2im;preview=true)
    wait(latest); wait(old)
    @test fetch(old) == false
    @test asyncview.rendered[].parameter == 2.2im
    @test asyncview.rendered[].samples == 100
    @test length(asyncview.points[]) == 100
end
