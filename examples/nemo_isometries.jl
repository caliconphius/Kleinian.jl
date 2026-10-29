# Optional: run in an environment where both Kleinian and Nemo are available.
using Kleinian
import Nemo

K = Nemo.algebraic_closure(Nemo.QQ)
e0,e1,e2,e3 = (basisvector(i,K(1)) for i in 0:3)
axis = Geodesic(e0∧(e1+e2))
rotation = Elliptic(axis,7)

# sinpi/cospi and the axis normalization are algebraic and exact here.
@assert isidentity(rotation^7)
@assert !isidentity(rotation)
p = HyperbolicPoint(clvector(K(5),K(1),K(2),K(3)))
@assert MobiusIsometry(Kleinian.matrix(rotation))(p) == rotation(p)
@assert inv(rotation)(rotation(p)) == p
println("Order-seven rotation verified exactly over Nemo's algebraic closure.")
