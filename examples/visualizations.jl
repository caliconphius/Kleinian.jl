# First load a backend, e.g. `using GLMakie` or `using CairoMakie`, then include
# this file. Only the backend package needs installing in your environment.
using Kleinian
using GLMakie
"""The same H³ objects, viewed in the three standard models."""
function visualization_comparison(;resolution=48)
    origin = HyperbolicPoint(1,0,0,0)
    point = modelpoint((0.28,-0.2,0.1);model=:ball)
    axis = Geodesic(IdealPoint(5,-4,0,-3),IdealPoint(5,4,3,0))
    plane = HyperbolicPlane(clvector(1,0,2,0))
    horosphere = Horosphere(clvector(2,2,0,0))
    fig = Figure(size=(1500,550),figure_padding=30)
    for (i,(model,title)) in enumerate(zip((:klein,:ball,:halfspace),
                                         ("Klein","Poincaré ball","Upper half-space")))
        ax = Axis3(fig[1,i];title,aspect=:data,azimuth=0.9,elevation=0.45,
                        xlabel="x",ylabel="y",zlabel="z")
        modelboundary!(ax;model,resolution=24,extent=2)
        hyperbolicplot!(ax,plane;model,resolution,color=(:steelblue,0.25))
        hyperbolicplot!(ax,horosphere;model,resolution,color=(:orange,0.5))
        hyperbolicplot!(ax,axis;model,resolution=200,color=:crimson,linewidth=3)
        hyperbolicplot!(ax,[origin,point];model,color=:black,markersize=12)
        if model === :halfspace
            limits!(ax,-2,2,-2,2,0,3)
        else
            limits!(ax,-1,1,-1,1,-1,1)
        end
    end
    Label(fig[2,1:3],"Black: points   ·   Red: geodesic   ·   Blue: plane   ·   Orange: horosphere")
    fig
end
visualization_comparison()
