using Farey
F2 = FreeGroup(:f, :g)
id = one(F2)
f, g = GPC.gens(F2)
F, G = GPC.gens(F2) .|> inv

F2ab = FreeGroup(:a, :b)
idab = one(F2ab)
a, b = GPC.gens(F2ab)
A, B = GPC.gens(F2ab) .|> inv

farey(q) = farey_word(q, f, g)
Γ(q) = christoffel(q, a, b)

# Q = @cf0 [1,4]
# qs = [farey_neighbours(Q)..., Q]


# symbs = [
#     "q1 = ",
#     "q2 = ",
#     "q1⊖q2 = ",
#     "q1⊕q2 = ",
# ]


# ["Christoffel: "*sy*"$(Γ(q))" for (sy, q) in zip(symbs, qs)] .|> println;
# ["Farey: "*sy*"$(farey(q))" for (sy, q) in zip(symbs, qs)].|> println;

# rr = numerator(Rational(Q))
# ss = denominator(Rational(Q))
# Q
# W = farey(Q)
# s_seq(Q)[1:rr]'
# W.word[1:ss] |> F2
# q2

# ϕ = Farey.Hom(F2, F2, (
#     f=>f*g,
#     g=>F
# )) 

ϕ = Farey.Hom(F2, F2, (
    f=>f,
    g=>f*g)) 

ψ = Farey.Hom(F2, F2, (
    f=>f,
    g=>F*g
)) 




# u, v = ψ.([f,g])

# farey_word(2//3, u, v)<<2

# u1 = christoffel(q1, u, v)<<2
# u1 = christoffel(q1, f, g)
# # u2 = christoffel(q2, f, g)
# u2 = christoffel(q2, f*g, G*F*G)
# (farey_word(Q, u, v))

# s_seq(Q)'

# N = 10
# Qn = [ContinuedFraction(-n//(3n - n^2 -1)) for n in 1:N]


A0 = ContinuedFraction[0 1 ; -1 3]
A1 = ContinuedFraction[1 1 ; 1 2]



function matrix_eval(w::Farey.AbstractGroupElement, a::Dict{K, AbstractMatrix{T}}) where {T<:BigFloat, K}
    
    add_inverses = Dict([(a|>pairs)..., [Farey.inv(x, w.parent.alphabet)=>inv(y) for (x,y) in (a|>pairs)]...])
    matrices = [add_inverses[Int(i)] for i in w.word]
    prod(matrices)
    
end


function matrix_eval(w::Farey.AbstractGroupElement, a::Dict) 
    
    add_inverses = Dict([(a|>pairs)..., [Farey.inv(x, w.parent.alphabet)=>inv(y) for (x,y) in (a|>pairs)]...])
    matrices = [add_inverses[Int(i)] for i in w.word]
    prod(matrices)
    
end

using LinearAlgebra:det

function inv(A::Matrix{ContinuedFraction})

    
    d = 1/det(A)
    # (size(A)!=(2,2) || abs(d)!=1) && error("$A ∉ SL₂ℤ")

    return d * BigInt[A[2,2] -A[1,2]; -A[2,1] A[1,1]]
end


A0

inv(A0)

integer_inverse(A1)

mat_val = Dict([
    1=>A0,
    2=>A1
])



mat_val2 = Dict([
    1=>A0,
    2=>A0*A1
])

(farey(1//2) |> ϕ) << 1

matrix_eval(farey(545//29355), mat_val2)
AAA = matrix_eval(Γ(3//203), mat_val)
AAA1 = matrix_eval(Γ(3//203), mat_val)
integer_inverse(AAA) - AAA1
AAA * AAA1
A0 * integer_inverse(A0)

farey(234//1007)

@cf 34//71
@cf 3//203

@cf 1/3

