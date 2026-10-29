# α = Kx(ζ^Int(l / a))
# β = Kx(ζ^Int(l / b))
# γ = Kx(ζ^Int(l / c))
# trW = γ + 1 / γ
# trW2 = 2 + trW^2

# printstyled(styled"{yellow:tr(W)} = $?= $(pretty_printing(trW^2-2))\n")
include("working_setup.jl")
using Nemo
Qb = algebraic_closure(QQ)
Rx, x = polynomial_ring(Qb, "z")
Kx = fraction_field(Rx)
bb = one(x)#root_of_unity(Qb, 3, 1)
aa = one(x)#oot_of_unity(Qb, 4, 1)
B = Kx[
  bb BigInt(0);
  x^2 1/bb
]


A = Kx[
  aa BigInt(1);
  0 1/aa
]

function inv2(A::T) where T
  A.base_ring[
    A[2, 2] -A[1, 2];
    -A[2, 1] A[1, 1]
  ]
end
inv2(A)

using Farey
using Farey.Monoids

f, g = FreeGen("f g")
a, b = FreeGen("a b")
gens = [f, g]
ms = [A, B]
function ρ(Q; gens=gens, ms=ms)
  mats = Dict([l => r for (l, r) in zip(gens, ms)])
  W = farey_word(Q, -f, g)
  rW = (mateval(W, mats)) |> trace
  # pretty_printing(rW) |> printstyled
  # println("\n")
end

for i in 1:(k+10)
  Po = ρ(1 // i)
  Po = change_coefficient_ring(R, Po)
  @show i
  # pretty_printing(Po) |> printstyled
  # println("")
  println((Po(z)))
end

using Base.Iterators
Iterators.accumulate(Iterators.cycle([next_calkin_wilf], 10), init=ZZ(0) // 1) do y, x
  x(y)
end |> x -> Iterators.filter(≤(1), x) |> collect

using DataStructures

function f_seq(n::Int)
  rst = SortedSet{Rational}(Rational[0, 1])
  for den in 1:n, num in 1:den-1
    push!(rst, Rational(num, den))
  end
  return rst
end

for n in 1:11
  print("F_$n: ")
  for frac in f_seq(n)
    print(numerator(frac), "/", denominator(frac), " ")
  end
  println()
end


k = 13
R = GF(k)
z = R(7)
for q in f_seq(10)
  Po = ρ(q)
  Po = map_coefficients(BigInt, numerator(Po))
  Po = change_coefficient_ring(R, (Po))
  printstyled(styled"{red:$q}\n")
  # pretty_printing(Po) |> printstyled
  # println("")
  println((Po(z)))
end
Qb = algebraic_closure(QQ)
Rx, x = polynomial_ring(Qb, "z")
Kx = fraction_field(Rx)
bb = one(x)#root_of_unity(Qb, 3, 1)
aa = one(x)#oot_of_unity(Qb, 4, 1)
q = @cf0 [1,1,2]
Po = ρ(q)


B = Kx[
  bb BigInt(0);
  x^2 1/bb
]


A = Kx[
  aa BigInt(1);
  0 1/aa
]

J = Kx[
  Qb(im)*x 0;
  0 -Qb(im)*1//x;
  ]

Θ(Q) = -conj_prefix(f\farey_word(Q, f, g))[1]
Fw(Q) = farey_word(Q, f, g)
Rw(Q)  = farey_word(Q, f, g)
A1 = J * A * inv2(J)
B1 = J * B * inv2(J)
ms = [A, B]
Θ(1//3)
ZZx, zx = polynomial_ring(ZZ,:zx)

function ρw(W; gens=gens, ms=ms, Fk = QQ)
  mats = Dict([l => r for (l, r) in zip(gens, ms)])
  rW = (mateval(W, mats)) |>
    # z->map_coefficients(BigInt, numerator(z))
    z->z
  Zx, x = polynomial_ring(Fk,:z)
  FF = z->map_coefficients(BigInt, numerator(z))
  # FF = z->z
  rW = [el|>FF|>l->Zx(l.coeffs) for el in rW]
  # pretty_printing(rW) |> printstyled
  # println("\n")
  return rW
end
# J *
import Nemo:tr
fg_trace = ρw(Fw(3//5))|>tr

# pol = tr(ρw(Fw(3//5); Fk=Qb))-2
# faccc = pol|>factor
# faccc.fac
# facccc
# C = CalciumField()
# C(zz0)

γ = roots(Qb, cyclotomic(3, zx))
# rls = roots(Qb, pol-γ[1])

# zz0 = rls[1]

# KΩ, lΩ = cyclotomic_field(3, "l")
# Kx  = Nemo.polynomial_ring_only(KΩ, :t)

# facs = factor(trace(ρw(Fw(3//5); Fk = Qb))-γ[1])
# facs.fac
# facs2 = factor(trace(ρw(Fw(3//8); Fk = Qb))-2)
# pws = [first(l) for l in facs.fac|>pairs]
# ρw(Fw(3//5))
# pws[4]

# using JLD2
# saved = Dict([
#   "facc"=>("factor(trace(ρw(Fw(3//5); Fk = Qb))-γ[1])",facs),
#   "facc2"=>("factor(trace(ρw(Fw(3//8)))-2)",facs2),
#   ])
# # JLD2.save("saved.jld2", saved)
# import Hecke
# map_coefficients(Qb, pol)-γ[1]
