using Kleinian, Poincare, Riley, Random, GLMakie
stream = randomlimitset(rileygroup(2.0im);samples=256,rng=Xoshiro(2026))
hull = limitsethull(stream)
view = hyperbolicplot(hull;model=:ball,boundary=false,color=(:steelblue,.25),transparency=true)
lines!(view.axis,hulledgesgeometry(hull;model=:ball,pleats=true);color=:black)
display(view.figure)
