using Kleinian

# Time comes first; representatives need not have square -1.
origin = HyperbolicPoint(1,0,0,0)
p = HyperbolicPoint(clvector(5,3,0,0))
@assert coordinates(p) == [5,3,0,0]
@assert scalarpart(representative(p)^2) == -16
@assert p == 7p
@assert distance(origin,p) ≈ log(2)
@assert distance(3origin,7p) ≈ distance(origin,p)

# Only request a unit vector when a calculation needs one.
u = unitvector(p)
@assert scalarpart(u*u) ≈ -1
@assert coordinates(p) == [5,3,0,0]

println("Stored homogeneous coordinates: ",coordinates(p))
println("Equivalent scaled point: ",7p)
println("Distance from origin: ",distance(origin,p))
println("Unit representative: ",u)
