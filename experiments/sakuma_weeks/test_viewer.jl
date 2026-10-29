using Test, CairoMakie
import CairoMakie: Makie
include("viewer.jl")
@testset "Sakuma–Weeks model and transformation controls" begin
    v = sakuma_weeks_viewer(3//8;resolution=12)
    @test length(v.data.tetrahedra)==4
    for model in (:ball,:klein,:halfspace)
        v.modelmenu.selection[]=model
        for transform in (:identity,:U,:W,:S)
            v.transformmenu.selection[]=transform
            @test !isempty(v.axis.scene.plots)
        end
    end
    Makie.set_close_to!(v.layers,1)
    @test v.layers.value[]==1
    @test !isempty(v.axis.scene.plots)
    @test_throws ArgumentError developed_domain(PolynomialSystem(),2//5;parameter=1.0+0im)
    d = developed_domain(PolynomialSystem(),10//33;parameter=.45+.28im)
    @test !d.complete
    @test length(d.tetrahedra)==12
end
