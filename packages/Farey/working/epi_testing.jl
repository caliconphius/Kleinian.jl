using Farey
using LinearAlgebra
# computing in FreeGroups/Monoids!

F2ab = FreeGroup(:a, :b)
F2 = FreeGroup(:f, :g)
F2n = FreeGroup(:u, :v)
id = one(F2)
idab = one(F2ab)
f, g = GPC.gens(F2)
F, G = GPC.gens(F2) .|> inv
a, b = GPC.gens(F2ab)
A, B = GPC.gens(F2ab) .|> inv
farey(q) = farey_word(q, f, g)
Γ(p, q) = christoffel(p, q, a, b)
Γ(Q) = christoffel(Q, a, b)
f\g
N1 = Farey.Hom(F2ab, F2ab, (
    a=> a*b,
    b=> a,
))

N2 = Farey.Hom(F2ab, F2ab, (
    a=> b,
    b=> B*a,
))



↑(ϕ, n) = foldl(|>, ITR.take(ITR.cycle([ϕ]), n))

det2(A::Matrix) = A[1,1] * A[2,2] - A[1,2] * A[2,1]

function Γ_hom(A::Matrix{T}, l::S=one(a)) where {T<:Real, S<: Farey.AbstractGroupElement}
    (size(A)==(2,2) && abs(det2(A))==1) || 
        error("$A ∉ PSL₂ℤ")

    x, y = det2(A)==1 ? (a,b) : (b,a)
    A = Int.(A)
    Farey.Hom(F2ab, F2ab, (
        x=>Γ(A[1,1], A[2,1])^l, 
        y=>Γ(A[1,2], A[2,2])^l
    ))
end

function inv2(A::Matrix{T}) where T<:Integer
    1//det2(A) * T[
        A[2,2] -A[1,2];
        -A[2,1] A[1,1]
    ]
end

ab(g) = mapreduce(
    x->
    Int(x)==1 ? [1,  0] :
    Int(x)==2 ? [0,  1] :
    Int(x)==3 ? [-1, 0] :
                [0, -1]
    , + , g.word
)

ab(Γ(213, 314))

(N1↑4)(f)  |> x-> (ab(x), x)
(N1↑4)(g)  |> x-> (ab(x), x)



N2↑10

M = (N1↑4)

M.([f,g]) .|> 
        x-> x ^ (palindrome(a, b, 3))   



Q = 5//8|>ContinuedFraction
dd = denominator(Q)
nn = numerator(Q)
# Q = 1//3
(q1, q2, q3) = farey_neighbours(Q)
Mat = Int[
    numerator(q1|>Rational) numerator(q2|>Rational);
    denominator(q1|>Rational) denominator(q2|>Rational)

]



N1↑4

iMat = [1 0; 0 1]*inv2(Mat)

NQ = Γ_hom(Mat) |> N2
NQi = Γ_hom(iMat)

N1 |> NQ

(NQi∘NQ).([a,b]) .|> x->x^(b^2*a*b)

farey(1010//3436778);