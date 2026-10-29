using CairoMakie
import CairoMakie: Makie
include("viewer.jl")
view = sakuma_weeks_viewer(10//33)
save(joinpath(@__DIR__,"sakuma_weeks_10_33.png"),view.figure)
