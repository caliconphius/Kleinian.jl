using Kleinian, Riley
using GLMakie

# Pair disjoint unit circles centered at -2 and 2, then -2im and 2im.
# In the displayed complex chart: w = b - 1/(z-a).
# matrix(g) acts on the conjugate coordinate, hence conj below.
paircircles(a,b) = MobiusIsometry(conj.([b -a*b-1; 1 -a]))
generators = [paircircles(-2.0,2.0),paircircles(-2.0im,2.0im)]
points = limitset(generators;depth=10) # 78,732 ideal points
view = limitsetplot(points;figure=(size=(1000,1000),),colormap=:inferno)
display(view.figure)
# save("limitset.png",view.figure)
# For interactive aggregation: limitsetplot(points;async=true)
# For individual dots: limitsetplot(points;renderer=:scatter,markersize=1)
