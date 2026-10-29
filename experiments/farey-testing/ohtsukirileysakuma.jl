function ~(f, g)
  commu(f, g)
end

include("hall_words.jl")

f, g, W1, W2 = FreeGen("f g W1 W2")
z1, z2 = FreeGen("z1 z2")
w1 = ZG.GroupAlgebra(f / f)
cc = FreeGen("c")
zg_map = Dict([
  f => w1 + z1,
  -f => w1 - z1 + z1^2 + cc^10,
  -W1 => w1 - W1 + W1^2 + cc^10,
  -g => w1 - z2 + z2^2 + cc^10,
  -W2 => w1 - W2 + W2^2 + cc^10,
  g => w1 + z2,
  W1 => w1 + W1,
  W2 => w1 + W2,
])


magnus_embed = MonoidReps.MonoidRep(zg_map; identity=w1);
# lw = (W1*W2*(-W1^f))|>magnus_embed
# lw^2|>truncate_zg
# dc = magnus_embed((W1*W2*(-W1^f)))|>truncate_zg
dc
trc(dc, k) = truncate_zg(dc; L=k) - truncate_zg(dc; L=k - 1)
# truncate_zg(ZG.MagnusCoeff(commu(g,f),[g,f]); L=4)

ZG.MagnusCoeff(commu(f, g), [g, f, g])
function compute_magnus(w, n)
  out = []
  zg1 = ZG.GroupAlgebra(one(z1))
  for x in Iterators.product([[f => (zg1 * z1), g => (zg1 * z2)] for i in 1:n]...)
    coeff = ZG.MagnusCoeff(w, [g.first for g in x]) |> (ZG.ϵ)

    coeff == 0 || push!(out, (
      coeff,
      # ZG.MagnusCoeff(w, [g.first for g in x]),
      prod(reverse(last.(x))),
      n
    ))
  end
  # filter(!iszero∘first, out)
  out
end

function magnus_coefficient(w, s)
  ZG.MagnusCoeff(w, [FreeGen("$c") for c in s])
end



magnus_coefficient(f ~ g, "fg")
magnus_coefficient(Fw(2 // 5), "ffgfg")
compute_magnus(Fw(4 // 7), 4)
4 // 7 |> Γ
compute_magnus(commu(f, g), 2)
compute_magnus(commu(-f, g), 7)
compute_magnus(commu(-f, -g), 2)
Δ(w, f, g) = ZG.MagnusCoeff(w, [f, g]) - ZG.MagnusCoeff(w, [g, f])
f, g = FreeGen("f g")
Δ(Fw(1 // 3), f, g)
@memoize function standard_comm(s::String)
  length(s) == 1 && return FreeGen(s)
  length(duval(s)) == 1 || return standard_comm.(duval(s))
  first, last = standard_factor(s)
  commu(standard_comm(first), standard_comm(last))
end

@memoize function standard_comm_tree(s::String)
  length(s) == 1 && return s
  first, last = standard_factor(s)
  [standard_comm_tree(first), standard_comm_tree(last)]
end
standard_comm("ggfg")
ZG.MagnusCoeff(Fw(2 // 5), [f, f, g, f, g])#|>ZG.ϵ
standard_comm_tree("ffgfgfg")
filter(x -> (length(x) == 6), generate_lyndon_words(6, "fg"))
ZG.MagnusCoeff(-wu * Fw(2 // 3), [f, f, f, f, f, g]) |> ZG.ϵ
ZG.MagnusCoeff((f * g * -(f ~ g)) * -Fw(2 // 3), [f, g, g]) |> ZG.ϵ

xf = TD.Replace(Dict(f => w1 + f, -f => w1 - f, g => w1 + g, -g => w1 - g))
wordrep(w) = TD.foldxl(*, xf, w)
wu = (g * f)
standard_comm("gfg")
standard_comm("fgg")
standard_comm_tree("fgg")
compute_magnus(Fw(2 // 3) / wu, 3)
wv = f * -(g) * (g ~ f) * standard_comm("fgg")
wu = standard_comm("fgg") * g * f
wv = g * -f * standard_comm("ffg") * standard_comm("fg")
standard_comm("ffgg")^f
(f ~ g)^2 * standard_comm("ffg")
compute_magnus((g ~ f)^2 * Fw(1 // 4), 3)
compute_magnus(wv \ Fw(1 // 3), 1)
truncate_zg(wordrep(-wu * Fw(2 // 3)))

(f ~ g)^2
((g ~ f) ~ g) * ((g ~ f) ~ f) * ((f ~ (g ~ f)) ~ g)
w1 = ((a ~ b) ~ ((b ~ a)^-f))
w2 = ((a ~ b) ~ (f)) * ((b ~ a) ~ f)
w1 / w2
(f ~ g) * (-f ~ g)

function ∧(u, v)
  u * v - v * u
end

standard_comm("g")

(wf - 1) ∧ (wg - 1) ∧ (wg - 1) * 1 // 2 +
(wg - 1) ∧ (wf - 1) ∧ (wf - 1) * 1 // 2
# commu(-g,-f)|>w->println.(Iterators.flatten([compute_magnus(w, i) for i in 1:4])|>collect);
K = 4
tword = commu(commu(-f, -g), -g)
tword = commu(f, g)^2
tword |>
# w->println.(Iterators.flatten([compute_magnus(w, i) for i in 1:K])|>collect);
# truncate_zg(ZG.MagnusCoeff(commu(f,g),[f,f,g]))
wr = Fw(1 // 3)
Δ(g / f, g, f)
Δ(g / f, f, g)
Δ(g^(f^2), f, g)
dc
zf = magnus_embed(f)
zg = magnus_embed(g)
zf^zg

commu(Γ(1 // 3), Γ(1 // 2))
qa, qb, qc = farey_neighbours(8 // 19)
qa - qb
cayleyΩ(qa), cayleyΩ(qb), cayleyΩ(qc)
commu(cayleyΩ(qa), cayleyΩ(qb)) |> conj_prefix
Ω(qa)
commu(Ω(qa), Ω(qb)).elem |> conj_prefix

q = 3 // 8
f * θa(q)[1:denominator(q)÷2-1]
θa(5 // 12)
te = (ee * oo)^2
(Fw(3 // 8)) |> reverse
(Fw(3 // 8))
k = -1
CF0([2, 2, -2])
CF0([2, 2, -2, -2k])
CF(11 // 28) +
((5) * 17) % 28
17 // 28 |> CF
q = CF0([3, 3, 5])
q = CF0([1, 2, 5])
q = CF0([4, -2])
q = CF0([1, 2, 1, 5])
1 // (4 - 11 // 16)
CF0(q.L |> reverse)
ZG.∂(a \ cayleyΩ(4 // 11) / b)[-a]
cayleyΩ(4 // 11)

truncate_zg(dc)
g ∈ f^2

commu(-cayleyΩ(1 // 3), -cayleyΩ(1 // 2))
commu(-Ω(2 // 5), -Ω(1 // 3))^Ω(1 // 2)
ffh(2 // 5)

Fw(2 // 5)
domega(qa, qb) = (commu(-Ω(max(qa, qb)), -Ω(min(qa, qb)))^Ω(max(qa, qb))).elem(f, g)

Fw(1 // 3) * Fw(1 // 2) * domega(qa, qb)
Q = (24) // 103
Q + 24 // 103^2 |> make_evens

Q + 6 // 17 |> CF
Q + 6 // 13 |> CF
Q - 3 // 13 |> CF
(24 // 103) |> CF
((4 - 4 // 103) * 1 // 17) |> CF
((1 + 1 // 103) // 13) |> CF
2 // 103 |> CF


make_evens(7 // 19)

make_evens.([4k // (17k + 1) for k in 0:7])
12 // 51 |> make_evens
make_evens(Q)
funnel(Q)
qa, qb, qc = farey_neighbours(Q)
domega(qa, qb) |> conj_prefix |> first

Fw(Q) / (
  Fw(min(qa, qb)) * Fw(max(qa, qb)) * domega(qa, qb)
)

Fw(1 // 3) * Fw(1 // 2) * (Fw(1 // 2)^-Fw(0))^Fw(1 // 2)
f * g * commu(-g, commu(g, f))
f / g * commu(g, commu(g, -f))
f / g * commu(g, commu(g, -f))

f * g * commu(-(g * f), g)
f / g * commu(g, f) * commu(commu(f, g), g)
compute_magnus(Fw(1 // 4), 3)
commu(f * g, W1)
commu(g, W1)^-f * commu(f, W1)


f * (g ~ f) * -g * (g ~ (g ~ f)) * (f ~ g)
x, y = FreeGen("a b")
f * -g * (g ~ f) * ((f ~ g) ~ g)^2 * ((g ~ (f ~ g)) ~ (g ~ f))
f * -g * (g ~ f) * ((f ~ g) ~ g)^2 * ((g ~ (f ~ g)) ~ (g ~ f))
f * -g * (g ~ f) * (((f ~ g) ~ g) * (g ~ f))^2 * ((f ~ g))^2
f * -g * (g ~ f) * (((f ~ g) ~ g) * ((f ~ g) ~ g) * ((g ~ (f ~ g)) ~ (g ~ f)))

MagnusCoeff
[g^2, (g ~ f)^2, (g ~ f) ~ g, ((g ~ f) ~ g) ~ (g ~ f)] |> prod
[g^2, (g ~ f)^2, (g ~ f) ~ g, ((g ~ f) ~ g) ~ (g ~ f)] |> prod
[f^2 * g^2, (g ~ f), (g ~ f) ~ g] |> prod
[f^3 * g^3, (g ~ f), (g ~ f) ~ (g), ((g ~ f) ~ (g))^g, (g ~ f), (g ~ f) ~ (g)] |> prod

(f ~ g) * (f ~ g ~ f) * ((f ~ g) ~ f)
((f ~ g) ~ g) * g * f
((f ~ g) ~ f) * (f ~ g)
(f ~ g) * ((f ~ g) ~ f) * ((f ~ (f ~ g)) ~ (g ~ f))

(f / g) * (f ~ g) * ((g ~ f) ~ g)


compute_magnus(((g ~ (f ~ g)) ~ (g ~ f)), 5)
compute_magnus((f * g)^2, 2)
generate_lyndon_words(4, "fg")
ZG.MagnusCoeff(f ~ (f ~ g), [f, f, g]) |> ZG.ϵ
ZG.MagnusCoeff((f ~ g) ~ g, [f, f, g]) |> ZG.ϵ
compute_magnus((f ~ g) ~ g, 3)
compute_magnus(f ~ f ~ g, 3)
compute_magnus((g ~ f)^2 \ g \ f \ Fw(1 // 3), 3)
commu(g, f) * commu(f, commu(f, g)) * commu(commu(commu(f, g), f), commu(f, g))


f * -g * (g ~ f) * ((f ~ g) ~ g)^2 * ((g ~ (f ~ g)) ~ (g ~ f))
f * -g * ((g) ~ (g ~ f))^2 * (g ~ f) * ((g ~ (f ~ g)) ~ (g ~ f))

f / g * (g ~ f) * ((f ~ g) ~ g)^2 * ((g ~ (f ~ g)) ~ (g ~ f))
f / g * (g ~ f) * ((f ~ g) ~ g) * ((g ~ (g ~ f)))

f * g * (-g ~ (g ~ f))
f * g * (g ~ (f ~ g))^(g * (f ~ g))
a^2 * b^2 * (-b ~ -a)^b
a^3 * b^3 * (-b ~ -a^2)^b^2 * (-b ~ -a)^b
a^3 * b^3 * (-b ~ -a)^b^2 * (-b ~ -a)^(a * b^2) * (-b ~ -a)^b

a^3 * b * (-b ~ -a) * -a * (-b ~ -a) * a * b * (-b ~ -a) * b
a^2 * b * (-b ~ -a) * (b ~ a) * a * b^2
a^2 * (-a ~ b) * a * (-a ~ b) * b * (-a ~ b) * b^2

f * (f ~ g)^2 * (-f) * ((g ~ f) ~ g) * ((g ~ (g ~ f)) ~ f) * (f ~ g)



# macro commutator()

# end

# @cf0 [4,-2,4,-2]
# @cf0 [2,2,-2,2,-2,2,2]
# make_evens(5//8)

# (3-5//8)//(8-5//8) |> CF |> make_evens
# k = (2 - 5//8)
# qe = (-3 * k + 2)//
# (-8 * k + 5)
# (qe) |> make_evens
# (qe) |> CF
# (3//8 - 1//64)|>make_evens

function gCF(L)
  1 // reduce((x, y) -> y + 1 // x, reverse(L))
end

CC = CalciumField()
C128 = AcbField(128)

RCC, x = R32[:x]
KCC = fraction_field(RCC)

gCF([2, 2, -2, x])
gCF([2, 2, -2, 2x, -2, 2, 2])
pad = gCF([2, 2, -2, 2x, 2, -2, -2])
invQ(x) = 1 // x
@memoize function epsilon_alg(seq, k, j, ifunc)
  k == -1 && return 0
  k == 0 && return get(seq, j, 0)
  epsilon_alg(seq, k - 2, j + 1, ifunc) + ifunc(epsilon_alg(seq, k - 1, j + 1, ifunc) - epsilon_alg(seq, k - 1, j, ifunc))
end

function epsilon_alg(seq, K; ifunc=invQ)
  [epsilon_alg(seq, 2k, j, ifunc) for (j, k) in Iterators.product(1:K, 1:K)]
end

using Plots
phi_approx = Iterators.cycle([x], 10) |> collect
approx = Iterators.map((y) -> gCF(phi_approx[1:y]) * x, 1:10) |> collect
R32 = ArbField()
xs = BigFloat.((1:10) ./ 10)
ys = epsilon_alg(approx, 3)[3, 3].(xs)
(real.(ys)) .|> BigFloat

plot(xs .|> BigFloat, pad.(xs) .|> BigFloat)
plot(xs, BigFloat.(real.(ys)))

# kB = fraction_field(KB)
# i,j = Hecke.quaternion_algebra(kB, x, w) |> Hecke.gens

# (i - i^2)

# p = x^2 + 1
# [i for i in coefficients(p)]

# KB, w = Kx[:w]

# n0 = w * KB[0 1; 0 0]
# nz = x * KB[0 0; 1 0]
# (n0 * nz)
# (nz * n0)*nz
# n0*(nz * n0)
# (n0 * nz - nz * n0)
# (one(n0) + n0) * (one(n0) + nz) * (one(n0) - n0) * (one(n0) - nz)

# σ3 = KB[1 0; 0 -1]
# σ2 = KB[0 1; -1 0]
# σ1 = KB[0 w; w 0]

# (zf * zg * (-zf)*(-zg)) |> truncate_zg
# f2, g2 = CyclicGen("f2 g2", order=2)
# wf, wg,w1 = ZG.GroupAlgebra.([[f=>one(w)*x],[g=>w],[one(g)=>one(w)]])

# (wf + w1) * (wg + w1) * (-wf + w1) * (-wg + w1)



