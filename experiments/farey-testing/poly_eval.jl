# Pkg.add("Hecke")
# include("working_setup.jl")
# using Farey.Monoids
using Farey
import Hecke
using Hecke: trace
using Nemo
using StyledStrings
include(joinpath(@__DIR__, "../../packages/Farey/src/MonoidHoms.jl"))
import .modhoms
using .modhoms: MonoidHom
const CF = ContinuedFraction

function qpair(qq)
  qq = CF(qq)
  qr = CF(0, qq[2:end-1])
  qt = qq[end] % 2 == 0 ? CF(0, qq[2:end-2]) : CF(0, [qq[2:end-1]..., 1])
  k = qq[end] ÷ 2

  qr, qt, k
end

function qnbs(Q::T) where T<:Number
  Q = CF(Q)
  qs = CF.([1 // (Q[2])])
  for (k, l) in enumerate(Q.L[2:end])
    append!(qs, [positive_form(CF(0, [Q[2:k+1]..., x])) for x in 1:l])
  end
  qs
end

function farey_expr(qq::Number)
  qq = CF(qq)

  idxr = 1 + sum(qq[3:end-1])
  idxt = qq[end] % 2 == 1 ? idxr + 1 : 1 + sum(qq[3:end-2])
  k = qq[end] ÷ 2

  w1, w2, f = MonoidWord("W$idxr W$idxt f")

  # length(qq) % 2 == 0 ? w2 * inv(w1^f)^k*w1^k : w2 * w1^k*inv(w1^-f)^k
  length(qq) % 2 == 0 ? w1^k * w2 * inv(w1^f)^k : inv(w1^-f)^k * w2 * w1^k
end


# begin
#   A = i->CyclicGen("x"; order = i)
#   B = i->CyclicGen("y"; order = i)
#   a, b = MonoidWord("a b")

#   hom(i,j) = modhoms.MonoidHom([
#     a => A(i),
#     b => B(j)
#   ])
#   hom(2,5)
#   fin_words = [
#     1=>(element("a b a b a -b a -b") => (2,3,2)),
#     2=>(element("a b -a b a -b") => (2,3,3)),
#     3=>(element("a b a -b") => (3,3,2)),
#     4=>(element("a b -a -b") => (3,3,2)),
#     5=>(element("a b a b 2") => (2,5,2)),
#     6=>(element("a b a b a b 4") => (2,5,2)),
#     7=>(element("a b a b 2 a b 4") => (2,4,2)),
#     8=>(element("a b a b a -b") => (2,4,2)),
#     9=>(element("a b a -b") => (2,3,2)),
#     10=>(element("a b a b a -b") => (2,3,2)),
#     11=>(element("a b a b a b a -b") => (2,3,2)),
#     12=>(element("a b a b a -b a b a -b") => (2,3,2)),
#     13=>(element("a b a b a b a -b a -b a b a -b a -b") => (2,3,2))
#   ]


# ws =  last.(fin_words) .|> x->(first(x)|hom(last(x)[1], last(x)[2]), last(x)[3]
# )
# wps = last.(fin_words) .|> x->(first(x)^last(x)[3]
# )|hom(last(x)[1], last(x)[2])
# end
function expandhom(q1)
  qns = qnbs(q1)

  f, g, w, ws... = MonoidWord("f g W " * join(["W$i " for i in eachindex(qns)]))
  t0, t1 = MonoidWord("t$(q1[2]) t$(q1[2]+1)")
  u, v = q1[2] % 2 == 0 ? (-f, -g) : (-g, -f)
  MonoidHom([[
      ws[1] => f * t0 * u / t0,
      ws[2] => f * t1 * v / t1,
    ]..., w => farey_word(0, f, g)]),
  MonoidHom([
    ws[i] => farey_expr(q) for (i, q) in Iterators.filter(x -> length(x[2]) > 2, enumerate(qns))
  ])
end


function pretty_printing(x::T; k=5) where T
  if T <: MatElem
    s = "[\n" * join(["$(Tuple(l)):\t" * repr(x[l]) * "\n" for l in eachindex(x)]) * "]"
  else
    s = repr(x)
  end


  s = replace(s,

    # r"[\,|]"=>"\n]",
    r"x\^?-?\d+" => l -> "{green:$("𝒵"*to_superscript(l[3:end]))}",
    r"x" => l -> "{green:$("𝒵")}"
    # r"z(?<subscr>\_\d+)(?<power>\^-?\d*)"=>s"z\g<power> z\g<subscr>"
  )
  s = replace(s,
    r"\^-?\d+" => (l -> "{bright_magenta:$(to_superscript(l[2:end]))}"),
    r"z\_\d+" => (l -> "{bright_blue:ζ$(to_subscript(l[3:end]))}"),
    r"\*" => (l -> ""),
    r"\d+" => (l -> "{code:$l}")
  )

  s = replace(s,
    Regex("\\d{$(2*k),}") => l -> "($(l[1:k]).($(length(l))).$(l[end-k+1:end]))"
  )
  StyledStrings.styled(s)

end


function mateval(x::Farey.AbstractMonoidElement, mats::Dict{K,V}; IFUNC=inv2)::V where {K,V}
  di = merge(mats, Dict([inv(l) => IFUNC(r) for (l, r) in pairs(mats)]))

  mapreduce(*, eachgen(x)) do x
    di[x]
  end
end

# pretty_printing(big(2)^100)|>printstyled
# # printstyled(styled"{red:𝒳}")
# # replace("x^-34", r"z"=>l->to_superscript(l[2:end]))

rW, rW1, ζ, x, α, β, rW2, A, B =
  let (a, b, c) = (2, 3, 5)
    l = lcm(a, b, c)

    # K, ζ = cyclotomic_field_as_cm_extension(2*l);
    # L, b = cyclotomic_field(2*l);
    K, ζ = Nemo.cyclotomic_field(2 * l, "z_$(2*l)")
    # Kl, ξ =  Nemo.cyclotomic_real_subfield(2*l, "z_0")

    Kx, x = polynomial_ring(K, "x")
    α = Kx(ζ^Int(l / a))
    β = Kx(ζ^Int(l / b))
    γ = Kx(ζ^Int(l / c))
    trW = γ + 1 / γ
    trW2 = 2 + trW^2

    printstyled(styled"{yellow:tr(W)} = $(pretty_printing(trW))\n")
    printstyled(styled"{yellow:tr(W)²-2} = $(pretty_printing(trW^2-2))\n")

    B = Kx[
      β BigInt(0);
      x inv(β)
    ]


    A = Kx[
      α BigInt(1);
      0 inv(α)
    ]

    Ai = A^(a - 1)
    Bi = B^(b - 1)

    f = CyclicGen("f")
    g = CyclicGen("g")
    L1 = [2]
    L = [L1..., 2 * c - 1, -reverse(L1)...]
    Q = CF(0, L) |> positive_form
    Q1 = CF(0, L1) |> positive_form
    println(A^a)

    printstyled(styled"Q = {bright_blue:$L} = {code:$Q}\nQ1 = {green:$Q1}\n")

    printstyled(styled"({red:a},{blue:b},{green:c}) = ({red:$a},{blue:$b},{green:$c})\n")

    printstyled(styled"{bright_magenta:lcm(a,b,c)} = {magenta:$l}\n")

    # printstyled(styled"({red:α},{blue:β},{green:γ}) = ({red:$α},{blue:$β},{green:$γ})\n")
    printstyled(styled"({red:α},{blue:β},{green:γ}) = ({red:$(pretty_printing(α))},{blue:$(pretty_printing(β))},{green:$(pretty_printing(γ))})\n")
    W = farey_word(1 // 3, f, g)
    W1 = farey_word(1 // 4, f, g)
    W2 = farey_word(1 // 5, f, g)



    gens = [f, g]
    ms = [A, B]
    mats = Dict([l => r for (l, r) in zip(gens, ms)])
    rW = (mateval(W, mats) |> trace)
    rW1 = (mateval(W1, mats) |> trace)
    rW2 = (mateval(W2, mats) |> trace)



    # printstyled((rW)*"\n")

    pretty_printing((trW)) |> printstyled
    println("\n")

    rW, rW1, ζ, x, α, β, rW2, A, B
  end;

# ms = [A, B]
# f = CyclicGen("f")
# g = CyclicGen("g")
# gens = [f, g]
# mats = Dict([l=>r for (l,r) in zip(gens, ms)])
# rWl(Q) = (mateval(farey_word(Q,f,g), mats)|>trace)

# pretty_printing(rWl(1//2)*rW)|>printstyled;println("\n")
# pretty_printing(rW2)|>printstyled;println("\n")
# polys = k->pretty_printing(gcdx(rWl(CF(0,[k])),rWl(CF(0,[k+1]))))|>printstyled;println("\n")
# # Φ


f = CyclicGen("f")
g = CyclicGen("g")
L1 = [2]
# L = [L1..., 2, -reverse(L1)...]
# Q = CF(0, L) |> positive_form
# Q1 = CF(0, L1) |> positive_form
# println(A^a)

printstyled(styled"Q = {bright_blue:$L} = {code:$Q}\nQ1 = {green:$Q1}\n")

printstyled(styled"({red:a},{blue:b},{green:c}) = ({red:$a},{blue:$b},{green:$c})\n")

printstyled(styled"{bright_magenta:lcm(a,b,c)} = {magenta:$l}\n")

# printstyled(styled"({red:α},{blue:β},{green:γ}) = ({red:$α},{blue:$β},{green:$γ})\n")
printstyled(styled"({red:α},{blue:β},{green:γ}) = ({red:$(pretty_printing(α))},{blue:$(pretty_printing(β))},{green:$(pretty_printing(γ))})\n")
W = farey_word(1 // 3, f, g)
W1 = farey_word(1 // 4, f, g)
W2 = farey_word(1 // 5, f, g)



gens = [f, g]
ms = [A, B]
mats = Dict([l => r for (l, r) in zip(gens, ms)])
rW = (mateval(W, mats) |> trace)
rW1 = (mateval(W1, mats) |> trace)
rW2 = (mateval(W2, mats) |> trace)
