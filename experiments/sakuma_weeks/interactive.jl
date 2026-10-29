using GLMakie
import GLMakie: Makie
include("viewer.jl")
# Replace 2//5 with 3//8, 5//12 or 10//33 to explore other domains.
view = sakuma_weeks_viewer(2//5)
display(view.figure)
