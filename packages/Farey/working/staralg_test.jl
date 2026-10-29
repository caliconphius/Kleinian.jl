using Farey
import StarAlgebras as SA
using StarAlgebras
include("autF2.jl")
F2

    B = SA.DiracBasis{UInt32}(BFS(F2))
    RG = SA.StarAlgebra(F2, B)
    SA.star(x::Farey.FreeGroupElement) = inv(x)
SA.star(RG(a))

SA.object(B)

B = SA.DiracBasis(F2)

KB.nwords

