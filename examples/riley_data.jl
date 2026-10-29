using Kleinian, Riley
using Random

"""
    rileyfigure(; samples=2_000_000, rngseed=2026)

Compare the real cusp, an approach to the imaginary cusp, and the Apollonian
cusp. Sampling is lazy; Makie retains only a compact array of 2D coordinates.
Increase samples for finer density detail; viewport limits exclude long cusp
excursions from the framing without changing the walk.
"""
function rileyfigure(;samples=2_000_000,rngseed=2026)
    fig = Figure(size=(1800,700),backgroundcolor=:black)
    # Move the strip into a bounded circle packing: w = 1/(z+i).
    chart = MobiusIsometry(conj.([0.0 1.0; 1.0 im]))
    cases = ((4.0,"ρ = 4 · real cusp"),
             (2.2im,"ρ = 2.2i · approaching the cusp"),
             (2.0im,"ρ = 2i · Apollonian cusp"))
    for (j,(ρ,title)) in enumerate(cases)
        ax = Axis(fig[1,j];title,titlecolor=:white,backgroundcolor=:black,
                  aspect=DataAspect(),limits=(-1.1,1.1,-2.1,0.1))
        hidedecorations!(ax)
        hidespines!(ax)
        generators = map(g -> MobiusIsometry(matrix(chart*g*inv(chart))),rileygroup(ρ))
        seed = chart(complexidealpoint(0.0))
        stream = randomlimitset(generators;seed,samples,burnin=1000,rng=Xoshiro(rngseed+j))
        limitsetplot!(ax,stream;colormap=:inferno,local_operation=log1p,async=false)
    end
    Label(fig[2,1:3],"Chart w = 1/(z+i) · reduced random walks · $(samples) samples per panel";
          color=:gray80)
    fig
end
