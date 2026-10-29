using Test
using Kleinian
using StaticArrays: SVector

@testset "Package scaffold" begin
    @test nameof(Kleinian) == :Kleinian
end

include("geometric_algebra.jl")
include("dot_product.jl")
include("hyperbolic_points.jl")

include("hyperbolic_objects.jl")
include("isometries.jl")
include("coefficient_protocol.jl")
include("printing.jl")
include("models.jl")
