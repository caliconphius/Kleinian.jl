using Kleinian

# Exact rational input, with time first. Write 2*e1 explicitly: 2e1 is a
# floating-point scientific-notation literal in Julia!
e0,e1,e2,e3 = (basisvector(i,big(1)//1) for i in 0:3)
origin = HyperbolicPoint(e0)
p = HyperbolicPoint(clvector(big(5)//1,1,2,3))

# Geodesics are simple timelike bivectors. Half-turns act through GA.
axis = Geodesic(e0∧e3)
halfturn = HalfTurn(axis)
@assert halfturn(halfturn(p)) == p
@assert halfturn(origin) == origin

# Two skew half-turn axes give a translation with a twist.
h1 = HalfTurn(Geodesic(e0∧e1))
h2 = HalfTurn(Geodesic((5*e0+3*e3)∧(e1+e2)))
g = h2*h1                      # h1 acts first
@assert g(p) == h2(h1(p))
@assert distance(origin,g(origin)) ≈ log(4)

# Positive multiplier 2 means signed length log(2), without evaluating it.
t = Translation(axis;multiplier=2//1)
@assert t(origin) == HyperbolicPoint(5,0,0,3)

# A future null vector specifies a horosphere, including its scale.
horosphere = Horosphere(e0+e3)
parabolic = Parabolic(horosphere,e1)
@assert origin in horosphere
@assert parabolic(origin) in horosphere
@assert parabolic(horosphere) == horosphere
@assert t(horosphere) == Horosphere(2*(e0+e3))

# Preserve the order as an integer; choose numerical precision at evaluation.
rotation = Elliptic(axis,7)
setprecision(BigFloat,256) do
    highprecision_rotation = mapcoefficients(BigFloat,rotation)
    @assert isidentity(highprecision_rotation^7;atol=BigFloat(2)^(-230))
end

# Projective matrices are an optional representation. Determinant-one lifts
# are available separately through sl2matrix(g).
A = matrix(g)
@assert MobiusIsometry(A)(p) == g(p)
@assert inv(g)(g(p)) == p

println("Two-half-turn matrix: ",A)
println("Image of the point: ",g(p))
println("Parabolic image of origin: ",parabolic(origin))
println("Exact GA and matrix actions agree.")
