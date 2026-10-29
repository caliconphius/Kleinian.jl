# using Pkg
# Pkg.resolve()
# Pkg.instantiate()
include("working_setup.jl")

# function Base.show(io::IO, t::BigInt)
#     if get(io, :compact, false)
#         print(io, string(t))
#     else
#         s = string(t)
#         rep = length(s) < 15 ? s : s[1:6] * "..." * s[end-5:end]
#         print(io, rep)
#     end
# end


# Pkg.add("Hecke")
using Hecke: cyclotomic_field, QQ, cyclotomic_field_as_cm_extension, trace
using Farey.Monoids
using Nemo
using StyledStrings
include(joinpath(@__DIR__, "../../packages/Farey/src/MonoidHoms.jl"))
# import .modhoms
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
let Q = @cf0 [2, 2, 2, 2]
  println.(qnbs(Q))
  FareyDecomp!(Q)
  farey_expr(Q)
end


Q = @cf0 [1, 1, 1, 1, 1, 1, 2, -1, -1, -1, -1, -1, -1]
farey_expr(Q)
qnbs(Q)

# Rational(Q)|>Float64, 3/5

using .modhoms: MonoidHom
begin
  A = i -> CyclicGen("x"; order=i)
  B = i -> CyclicGen("y"; order=i)
  a, b = MonoidWord("a b")

  hom(i, j) = MonoidHom([
    a => A(i),
    b => B(j)
  ])
  hom(2, 5)
  fin_words = [
    1 => (element("a b a b a -b a -b") => (2, 3, 2)),
    2 => (element("a b -a b a -b") => (2, 3, 3)),
    3 => (element("a b a -b") => (3, 3, 2)),
    4 => (element("a b -a -b") => (3, 3, 2)),
    5 => (element("a b a b 2") => (2, 5, 2)),
    6 => (element("a b a b a b 4") => (2, 5, 2)),
    7 => (element("a b a b 2 a b 4") => (2, 4, 2)),
    8 => (element("a b a b a -b") => (2, 4, 2)),
    9 => (element("a b a -b") => (2, 3, 2)),
    10 => (element("a b a b a -b") => (2, 3, 2)),
    11 => (element("a b a b a b a -b") => (2, 3, 2)),
    12 => (element("a b a b a -b a b a -b") => (2, 3, 2)),
    13 => (element("a b a b a b a -b a -b a b a -b a -b") => (2, 3, 2))
  ]


  ws = last.(fin_words) .|> x -> (first(x) | hom(last(x)[1], last(x)[2]), last(x)[3]
  )
  wps = last.(fin_words) .|> x -> (first(x)^last(x)[3]
  ) | hom(last(x)[1], last(x)[2])
end
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
Q = CF(0, [2, 2, 2, 2, 2])
fn, fw = expandhom(Q)
qs = qnbs(Q)
println.(enumerate(qs))
println(fn)
f, g, w, ws... = MonoidWord("f g W " * join(["W$i " for i in eachindex(qns)]))
fw
[(fn * fw^3)(w) for w in ws]

function pretty_printing(x::T; k=5) where T
  if T <: MatElem
    s = "[\n" * join(["$(Tuple(l)):\t" * repr(x[l]) * "\n" for l in eachindex(x)]) * "]"

  else
    s = repr(x)
  end


  s = replace(s,

    # r"[\,|]"=>"\n]",
    r"x\^?-?\d+" => l -> "{green:$("𝒵"*to_superscript(l[3:end]))}",
    r"x" => l -> "{green:$("𝒵")}",
    # r"z(?<subscr>\_\d+)(?<power>\^-?\d*)"=>s"z\g<power> z\g<subscr>"
  )
  s = replace(s,
    r"\^-?\d+" => l -> "{bright_magenta:$(to_superscript(l[2:end]))}",
    r"z\_\d+" => l -> "{bright_blue:ζ$(to_subscript(l[3:end]))}",
    r"\*" => l -> "",
    r"\d+" => l -> "{code:$l}",
  )

  s = replace(s,
    Regex("\\d{$(2*k),}") => l -> "($(l[1:k]).($(length(l))).$(l[end-k+1:end]))",
  )
  # println("hh")
  StyledStrings.styled(s)
end

function mateval(x::Farey.AbstractMonoidWord, mats::Dict{K,V}; IFUNC=inv2)::V where {K,V}
  di = merge(mats, Dict([inv(l) => IFUNC(r) for (l, r) in pairs(mats)]))

  mapreduce(*, eachgen(x)) do x
    di[x]
  end
end

pretty_printing(big(2)^100) |> printstyled
# printstyled(styled"{red:𝒳}")
# replace("x^-34", r"z"=>l->to_superscript(l[2:end]))

FareyDecomp!(CF(0, [2, 2, 2, 2]))


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
ms = [A, B]
f = CyclicGen("f")
g = CyclicGen("g")
gens = [f, g]
mats = Dict([l => r for (l, r) in zip(gens, ms)])
rWl(Q) = (mateval(farey_word(Q, f, g), mats) |> trace)

pretty_printing(rWl(1 // 2) * rW) |> printstyled;
println("\n");
pretty_printing(rW2) |> printstyled;
println("\n");
Q = @cf0 [2, 6, -2]

pretty_printing(gcd(rWl(1 // 2), rWl(1 // 3))) |> printstyled;
println("\n");
# Φ

MMM = zeros(3, 2)
@info MMM, 0, 0
[l for l in eachindex(MMM)]




pretty_printing(rWl(3 // 4)) |> printstyled;
println("\n");
pretty_printing(rWl(1 // 4)) |> printstyled;
println("\n");





printcomp(x...) = printstyled(IOContext(stdout, :compact => true), x...)

styled"\n{red:| }" * pretty_printing(mod(rW, rW1^2 - 2)) * "\n" |> printcomp
gcdx(3, 7)

rhoW, rhoWt, rhoWr,
qr, qt, Q,
ζ, x, α, β, l,
Kx =
  let (a, b, c) = (2, 3, 7)
    Q = positive_form(CF(0, [2, 2, 2, -2, -2]))
    l = lcm(a, b, c)

    K, ζ = Nemo.cyclotomic_field(2 * l, "z_$(2*l)")
    # K, ζ = QQ, 1
    # Kx, x = K[:x];
    Kx, x = polynomial_ring(K, "x")
    α = Kx(ζ^Int(l // a))
    β = Kx(ζ^Int(l // b))
    γ = Kx(ζ^Int(l // c))
    trW = γ + inv(γ)


    fn, fw = expandhom(Q)
    Wm = i -> (MonoidWord("W$(i)"))
    ws = [Wm(i) for i in eachindex(qnbs(Q))]
    qr, qt, k = qpair(Q)
    qs = qnbs(Q)
    qvec(q) = [numerator(q), denominator(q)]
    vecq(q) = q[1] // q[2]

    # printstyled(styled"{yellow:tr(W)} = $(pretty_printing(trW))\n")
    # printstyled(styled"{yellow:tr(W)²-2} = $(pretty_printing(trW^2-2))\n")

    B = Kx[
      β zero(x);
      x inv(β)
    ]


    A = Kx[
      α one(x);
      zero(x) inv(α)
    ]

    # Ai = A^(a-1)
    # Bi = B^(b-1)

    f = CyclicGen("f")
    g = CyclicGen("g")
    # L1 = [2]
    # L = [L1...,2*c-1,-reverse(L1)...]
    # Q = CF(0, L)|>positive_form
    # Q1 = CF(0,L1)|>positive_form

    printstyled(styled"Q = {code:$Q}\nqr = {green:$qr}\nqt = {green:$qt}\n")
    printstyled(styled"qr⊕$(to_superscript(2*k))qt = {bright_blue:$(vecq(qvec(qt)+2*k*qvec(qr)))}\n")

    printstyled(styled"({red:a},{blue:b},{green:c}) = ({red:$a},{blue:$b},{green:$c})\n")

    printstyled(styled"{bright_magenta:lcm(a,b,c)} = {magenta:$l}\n")

    # printstyled(styled"({red:α},{blue:β},{green:γ}) = ({red:$α},{blue:$β},{green:$γ})\n")
    printstyled(styled"({red:α},{blue:β},{green:γ}) = ({red:$(pretty_printing(α))},{blue:$(pretty_printing(β))},{green:$(pretty_printing(γ))})\n")
    Wr = farey_word(qr, f, g)
    Wt = farey_word(qt, f, g)
    W = farey_word(Q, f, g)
    gens = [f, g]
    ms = [A, B]
    mats = Dict([l => r for (l, r) in zip(gens, ms)])
    rhoWr = (mateval(Wr, mats))
    pretty_printing(rhoWr) * "\n" |> printstyled
    pretty_printing(rhoWr |> trace) * "\n\n" |> printstyled

    rhoWt = (mateval(Wt, mats))
    # println(rhoWt)
    pretty_printing(rhoWt) * "\n" |> printstyled
    pretty_printing(rhoWt |> trace) * "\n\n" |> printstyled

    rhoW = mateval(W, mats) |> trace
    pretty_printing(rhoW) * "\n" |> printstyled
    # rW1 = (mateval(W1^2, mats)|>trace)
    # rW,rW1, ζ, x, α, β, l
    # println.([(k,fw(wl), ql) for (k,(wl, ql)) in enumerate(zip(ws, qs))])

    rhoW, rhoWt, rhoWr,
    qr, qt, Q,
    ζ, x, α, β, l,
    Kx

  end;



q = Q
# FareyDecomp!(q)
p = (rhoW)
S, = residue_ring(Kx, trace(rhoWr));
rhoWr |> trace |> pretty_printing |> printstyled;
print("\n");
(rhoWt) |> trace |> pretty_printing |> printstyled;
print("\n");
# S(p)|>pretty_printing|>printstyled;print("\n")




f, g = MonoidWord("f g")
# collect(eachgen(inv((f*g)^2)))


K, ζl = cyclotomic_field_as_cm_extension(2 * l);
Kx, X = K["x"]
# printstyled(pretty_printing(rW1(X, ζl))*"\n")
# printstyled(pretty_printing(Φ(X, ζl))*"\n")
(pretty_printing(mod(Φ(X, ζl), rW1(X, ζl))) * "\n")


qq = @cf0 [2, 2, 5, 2]
expandhom(qq)
(qnbs(qq)) .|> l -> (l, qpair(l)[1:2]..., qpair(l)[end] * 2)





let qq = @cf0 [2, 2, 5, 2]
  # qq = @cf0 [2, 1,3]
  qq = CF(qq)

  (f, g) = MonoidWord("f g")
  Fw(Q) = farey_word(Q, f, g)
  Ψ(Q) = Fw(Q)[1:denominator(Q)]
  qa, qb, qc = farey_neighbours(qq)
  if qq[end] % 2 == 0

    qu, l = descendQ(qq, qa; l=1)
    k = l ÷ 2
    qr = qa
    qt = qu ⊖ qa
  else
    qu, l = descendQ(qq, qa; l=0)
    k = l ÷ 2
    qr = qa
    qt = qu ⊕ qa

  end
  qr, qt, qq, 2 * k, addK(qt, qr, 2 * k)

  commu(f, g) = f * g * inv(f) * inv(g)

  θa(Q) = f \ Ψ(Q)

  Γ(qq), Γ(qt), Γ(qr), 2 * k

end

qq = @cf0 [2]
qqq = @cf0 [2, 4, -2]
FareyDecomp!(7 // 16)


fww = (g0 * d0)
fwwi = (inv(g0 * d0))
hfww = (d0 * g0)
hfwwi = (inv(d0 * g0))

gww = (g1 * d1)
gwwi = (inv(d1 * g1))
hgww = (d1 * g1)
hgwwi = (inv(d1 * g1))

hom(fww)


let qq = @cf0 [2, 2, 5, 2]
  # qq = @cf0 [2, 1,3]
  qq = CF(qq)

  (f, g) = MonoidWord("f g")
  Fw(Q) = farey_word(Q, f, g)
  Ψ(Q) = Fw(Q)[1:denominator(Q)]
  qa, qb, qc = farey_neighbours(qq)
  if qq[end] % 2 == 0

    qu, l = descendQ(qq, qa; l=1)
    k = l ÷ 2
    qr = qa
    qt = qu ⊖ qa
  else
    qu, l = descendQ(qq, qa; l=0)
    k = l ÷ 2
    qr = qa
    qt = qu ⊕ qa

  end
  qr, qt, qq, 2 * k, addK(qt, qr, 2 * k)

  commu(f, g) = f * g * inv(f) * inv(g)

  θa(Q) = f \ Ψ(Q)

  qnbs(qq) .|> qnbs#farey_neighbours.|>l->(Ψ(l[1]),Ψ(l[2]), Ψ(l[1]⊕l[2]))#.|>l -> [(s.exp) for s in filter(x->(Monoids.id(x)==1), Monoids.gens(l)|>collect)]

end

SFX.SuffixTree("abababc")

farey_word(13 // 92, f, g)
