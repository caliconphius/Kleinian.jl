using Kleinian

# Coordinates and basis indices put time first.
e0, e1, e2, e3 = (basisvector(i, 1//1) for i in 0:3)
@assert e0^2 == scalar(-1)
@assert e1^2 == e2^2 == e3^2 == scalar(1)
@assert e0*e1 == -(e1*e0)
@assert e0 ⋅ e0 == scalar(-1)
@assert e0 ⋅ e1 == scalar(0)
println("Lorentz dot product e0 ⋅ e0: ", e0 ⋅ e0)

B = e1 ∧ e2
println("Bivector: ", B)
@assert grade(B, 2) == B
@assert reverse(B) == -B
@assert undual(dual(B)) == B

# A future unit timelike point and a spacelike plane normal.
X = clvector(5//3, 4//3, 0//1, 0//1)
n = e1 + 2*e2
@assert X*X == scalar(-1)

Y = -n*X*inv(n)
println("Reflected multivector: ", Y)
coordinate_reflection = X - (2*scalarpart(leftcontract(X,n))/scalarpart(n*n))*n
@assert Y == coordinate_reflection
@assert Y*Y == scalar(-1)
@assert -n*Y*inv(n) == X

println("Cl(3,1), metric (-,+,+,+), basis (e0,e1,e2,e3)")
println("Reflected point (t,x,y,z): ", Tuple(coefficient(Y, 1 << i) for i in 0:3))
