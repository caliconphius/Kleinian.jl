include(joinpath(@__DIR__, "../../packages/Farey/working/GroupAlgebra.jl"))
import LinearAlgebra: diagm
import .GroupAlgebras as ZG


function ψf(Q)
  [ϕ^l for (l, ϕ) in zip(CF(Q), ITR.cycle([ld1, lg]))] |> prod
end
ψe(Q) = [ϕ^l for (l, ϕ) in zip(make_evens(CF(Q))[1].L, ITR.cycle([ld1, lg]))] |> prod |> inv
ψe(4 // 7)
ffh(Q) = ψf(Q).rep
ffhi(Q) = (-ψf(Q)).rep
ffhi(1) *
Q = @cf 3 // 13
q1, q2, q3 = farey_neighbours(Q)
Γ(q1)
ZG.∂(Γ(q1), [a, b])
ZG.∂(Γ(q2), [a, b])
ψf(Q).rep |> typeof
function birmanJac(ϕ::MonoidHom)
  gens = [FreeGen(k) for k in keys(ϕ.map)]
  [ZG.∂(ϕ(g), h) for (g, h) in Iterators.product(gens, gens)]
end
ql = 5 // 12
FH = ψf(ql)
lf(FH, v) = ZG.mapt(FH, v)
det2(A) = A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]

fzg(v) = ZG.mapt(FH.rep, v)
fizg(v) = ZG.mapt((-FH).rep, v)
fa, ga = ZG.GroupAlgebra(f), ZG.GroupAlgebra(g)
fz, gz = fa - 1, ga - 1

(fz |> typeof) <: ZG.AbstractGroupRing

D1 = birmanJac(FH.rep)
D2 = birmanJac((-FH).rep)
del1 = fzg.(D2)
det2(del1)
det2(D1)
del1[1, 1]
D1[1, 1]
del1[2, 1]
D1[1, 2]
D1[1, 1] * del1[1, 1]
D1[1, 2] * del1[2, 1]

D1[2, 2]

del1 .|> ZG.ϵ

inv2f(A) = [
  A[2, 2] -A[1, 2];
  -A[2, 1] A[1, 1]
]


ZG.mapt(FH.rep, fz)
D2 = birmanJac((-FH).rep)
D2 = birmanJac(ffhi(3 // 5))
ffh(3 // 5)
ffhi(3 // 5)(christoffel(1 // 2, f, g))
ffh(5 // 12) * ffhi(5 // 12)
ffh(2 // 5)(f)
ZG.mapt(ffh(2 // 5), (fz))
q = @cf0 [2, 2, 2, 2]
fhz = x -> lf(ffh(q), x)
(birmanJac(ffh(q))) * fhz.((birmanJac(ffhi(q))))
(fhz.((birmanJac(ffhi(q))))) .|> l -> ZG.∂(l, [f, g]) .|> ZG.ϵ
(fhz.((birmanJac(ffhi(q))))) .|> ZG.ϵ
(fhz.((birmanJac(ffhi(q))))*[fz, gz])[2] |> l -> ZG.∂(l, [f, g]) .|> ZG.ϵ
fhz.((birmanJac(ffhi(q)))) * [fz, gz]
fhz.((birmanJac(ffhi(q))))[1, 1]
fhz.((birmanJac(ffhi(q))))[2, 2]
((birmanJac(ffhi(q))))[1, 2]
fhz.((birmanJac(ffhi(q))))[2, 1]

fz * gz - gz * fz
wa, wb, wc, wd = ZG.GroupAlgebra.(FreeGen("a b c d"))
a, b, c, d = (FreeGen("a b c d"))
matm = typeof(wa * wb)[
  wa wb;
  wc wd
]
matmi = typeof(wa * wb)[
  wd -wb^a;
  -wc wa
]
matm * matmi
trm = diagm([tr(matm), tr(matm)])
matm * inv2(matm)
(trm * matm - matm * trm)
mu = MonoidHom([f => g * f, g => -g])

lm = ψf(5 // 12) * -lg1
lm.rep.([f, g]) .|> ab
(((lm.rep) | mu).([f, g])) .|> ab
(((lm.rep) | mu).([f, g])) .|> conj_prefix
(⊙)(f, g) = commu(-f, -g)

wu = (f * g)^2 * (f ⊙ g) / (f ⊙ (f ⊙ g)) * ((f ⊙ g) ⊙ g)^2
ZG.MagnusCoeff(-wu * f^2 * g^2, [f, f, g]) |> ZG.ϵ
compute_magnus(-wu * f^2 * g^2, 4) .|> first
Base.empty!(Memoize.memoize_cache(standard_comm))
(standard_comm("ffgfg")) * (standard_comm("fg"))
min_cycle("ffgfg")
standard_comm("fgffg")
(standard_comm("fffg") * standard_comm("ffg"))
standard_comm("fggg")
standard_comm("fhg")
(-prod(standard_comm("hfg"))) * prod(standard_comm("ghf")) * standard_comm("fgh")
standard_comm("fgh")
standard_comm("fhg")
standard_factor("ffgfg")
standard_comm("fffg")^f^2
f ⊙ (f ⊙ g)
(f ⊙ g) ⊙ g

(standard_comm("f") ⊙ standard_comm("gh"))
(standard_comm("f") ⊙ standard_comm("gh"))
standard_comm("fggh")
