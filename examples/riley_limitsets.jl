using GLMakie
include("riley_data.jl")

fig = rileyfigure()
display(fig)
# save("riley-cusps.png",fig;px_per_unit=2)
