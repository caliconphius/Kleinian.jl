using Farey
include("working_setup.jl")
Q = @cf0 [1, 1, 1,2]
Qr = CF(0, reverse(Q.L))
qs = farey_neighbours(Q)
wq = Γ(Q)
wa, wb, wc = [Γ(q) for q in qs]
pa, pb, pc = [Ψ(q) for q in qs]
w = prod([(fh)^l for (l, fh) in zip(Q.L, ITR.cycle([g0, d1]))])
# wh = ψf(Q)
w = prod([(fh)^l for (l, fh) in zip(Q.L, ITR.cycle([d1,g0]))])

ϕ = hom(w)
ϕi = hom(w|>inv)
ϕi
pa|>(ϕi)|>conj_prefix
Ψ(3//5)|>(ϕ^3)|>conj_prefix
Fw(5//8)|>(ϕi)|>conj_prefix

pb|>(ϕ)|>
  conj_prefix     |>
  x->x[2]*F*G*f*g |>
  conj_prefix     |>
  x->x[2]
wb|>(ϕi)
wc|>(ϕi)
wq|>(ϕi)

import .SuffixTrees as sfx
tr = sfx.SuffixTree(pa)

Γ(27//61)
# w = g1*d1|>
    # inv |>
    # x->x

([
  g/f^2*g,
  f*g*f,
  g*f*g,
]).|>x->(x=>(ϕ(x), ϕi(x)))


([
  f,
  f*g,
  f\g,
  f/g,
  g,
]).|>x->(x=>(ϕ(x), ϕi(x)))

(ϕ).([g*f, f/g, f*g, g])

tr = sfx.SuffixTree(pa)
tr = sfx.SuffixTree(pb)
# tr = sfx.SuffixTree("bacabacacabacab")
# tr = sfx.SuffixTree("hey mum how has your day been going... the enigma machine. alksdjaskldjaks")
