using Farey
# Pkg.add("Hecke")
# include("Workdin")
using Hecke: trace, complex_embeddings, radical_extension
using Farey.Monoids
using Nemo
using Nemo: QQ
using StyledStrings, UnicodeFun
include(joinpath(@__DIR__, "../../packages/Farey/src/MonoidHoms.jl"))
# import .modhoms
using .modhoms:MonoidHom
const CF = ContinuedFraction
include("poly_eval.jl")

using .modhoms:MonoidHom
begin
  A = i->CyclicGen("a1"; order = i)
  B = i->CyclicGen("b1"; order = i)
  a, b = MonoidWord("a b")

  hom(i,j) = MonoidHom([
    a => A(i),
    b => B(j)
  ])
  hom(2,5)
  fin_words = [
    1=>(element("a b a b a -b a -b") => (2,3,2)),
    2=>(element("a b -a b a -b") => (2,3,3)),
    3=>(element("a b a -b") => (3,3,2)),
    4=>(element("a b -a -b") => (3,3,2)),
    5=>(element("a b a b 2") => (2,5,2)),
    6=>(element("a b a b a b 4") => (2,5,2)),
    7=>(element("a b a b 2 a b 4") => (2,4,2)),
    8=>(element("a b a b a -b") => (2,4,2)),
    9=>(element("a b a -b") => (2,3,2)),
    10=>(element("a b a b a -b") => (2,3,2)),
    11=>(element("a b a b a b a -b") => (2,3,2)),
    12=>(element("a b a b a -b a b a -b") => (2,3,2)),
    13=>(element("a b a b a b a -b a -b a b a -b a -b") => (2,3,2))
  ]


ws =  last.(fin_words) .|> x->(first(x)|hom(last(x)[1], last(x)[2]), last(x)[3]
)
# wps = last.(fin_words) .|> x->(first(x)^last(x)[3]
# )|hom(last(x)[1], last(x)[2])
end

struct FAut
    word::AbstractMonoidWord
    hom::MonoidHom
    ihom::MonoidHom
end

x, y = MonoidWord("x y")

g0 = MonoidHom([
    x => x,
    y => x*y,
])

g1 = MonoidHom([
    x => x,
    y => y*x,
])

g0i = MonoidHom([
    x => x,
    y => x\y,
])

g1i = MonoidHom([
    x => x,
    y => y/x,
])

d0 = MonoidHom([
    x => y*x,
    y => y,
])

d1 = MonoidHom([
    x => x*y,
    y => y,
])

d0i = MonoidHom([
    x => y\x,
    y => y,
])

d1i = MonoidHom([
    x => x/y,
    y => y,
])

ee =  MonoidHom([
    x => element("$y"),
    y => x,
])

oo =  MonoidHom([
    x => element("$x -1"),
    y => y,
])

Γ(Q) = christoffel(Q, x, y)

il = oo*ee*oo*ee

idhom = oo^2

G,D,G1,D1, E, O = MonoidWord("G D G1 D1 E O")

mapd = Dict([
    G=>g0,
    D=>d0,
    G1=>g1,
    D1=>d1,
    E=>ee,
    O=>oo,
])
idict = Dict([
    g0=>g0i,
    d1=>d1i,
    g1=>g1i,
    d0=>d0i,
    g0i=>g0,
    d1i=>d1,
    g1i=>g1,
    d0i=>d0,
    ee=>ee,
    oo=>oo,
])
# get!(idict, ee, ee)

homev(w) = mateval((w), mapd; IFUNC=l->idict[l])

Q = @cf0 [1,1,3]
# Q = @cf0 [1,1,2]
    Qa, Qb, Qc = farey_neighbours(Q)
    wb = (christoffel(Qb, G, D1))
    wa = (christoffel(Qa, G, D1))
    # wc = christoffel(Qc, cx, cy)
    # w = christoffel(1//3, G1, D)
    cx = (G/G1)
    cy = (D/D1)
    invw = (E*O)^2
    aw = (E*O*E)
    qr, qt = qpair(Q)
    # w = (cx \ w )*
    # -D1
    wc = (christoffel(Qc, G, D1))
    # w = -D1
    w = -G*wc
    fh = homev(-w)
    fhi = homev(w)
    # fh

    # println.([q|fhi for q in [Γ(Qa),Γ(Q),Γ(Qb),Γ(Qc)]])
    # println.([q=>Γ(q)|fhi for q in [Qa,Q,Qb,Qc, qr]])
    # fhl = mateval((w), mapd; IFUNC=l->idict[l])
    w
    # s_seq(Qb)
    wp = (G\D1)^2
    fp= homev(-(wp))
    fpi= homev((wp))
Qpp =  CF(0, [1,1,2,2])
printstyled(styled"{tip:$(Γ(Qpp))}\n")
println.([Rational(q)=>Γ(q)=>(fp*fh)(Γ(q)) for q in [Qa,Q,Qb,qt, qr]])
fhi*fpi


σ = k->(G^k*D1)
τ = k->(G^k)*(D1)*-G
# τ = k->(G1^k)*(D)
ti = k->homev(
τ(k)*
# σ(4)*
one(x)|>inv
)
tt = k-> homev(τ(k))
ss = k-> homev(σ(k))
si = k->homev(
σ(k)*
# σ(4)*
one(x)|>inv
)

Γ(CF(0,[2,2,2]))|ti(2)
Γ(CF(0,[1,2]))|si(1)

Γ(CF(0,[2,1,2]))|si(2)
wcr = Γ(CF(0,[4,7,4])); (wcr|(ti(4)|g1|il|(ti(1)|g1|il)^5)|si(1)|si(4))
Q = @cf0 [4,7,6,2]
qr,qt = qpair(Q)
Qa, Qb, Qc = farey_neighbours(Q)
wcr = Γ(Q); (wcr|(ti(4)|g1|il|(ti(1)|g1|il)^6)|si(7))

wi = inv((τ(4))*(τ(1))^6)
homev(wi)(Γ(Q)),
homev(wi)(Γ(Qb)),
homev(wi)(Γ(Qa)),
homev(wi)(Γ(qt))

Q = @cf0 [5,4,4,5]
qn = @cf0 [5,4]
# qn =
# Q = @cf0 [4,1,1,2,2]
qr,qt = qpair(Q)
Qa, Qb, Qc = farey_neighbours(Q)
wcr = Γ(Q); (wcr|(ti(4)|g1|il|(ti(1)|g1|il)^6)|si(7))


chom(Q) = prod([h^l for (h,l) in zip(ITR.cycle([G, D1]), [Q.L[1:end-1]..., Q[end]])])
wp = (G1\D)
wcr = -chom(Q)
w = wcr
ϕ = homev(-w)
ϕi = homev(w)

[
q.L => ϕi(Γ(q))
for q in [Q, Qa, Qb, Qc, qt, qn]
# ϕi(Γ(Qa)),
# ϕi(Γ(Qb)),
# ϕi(Γ(Qc))
]
# ϕ(Γ(qr)),


Fw = Q->farey_word(Q, x, y)
Wq = Fw(3//8)
# SFX.SuffixTree(Wq)
k = 3
wpj = D1/G1
wp = (wpj^k)*D
hp = homev(wp)
styled"{red:$(hp(Wq))}"|>printstyled
hp

f, g = CyclicGen("f g")
inv2(A::T) where T = A.base_ring[
    A[2,2] -A[1,2];
    -A[2,1] A[1,1];
]


function ρellip(a, b, c)
l = lcm(a,b, c)

K, ζ = Nemo.cyclotomic_field(2*l, "z_$(2*l)")


Kx, x = polynomial_ring(K, "x");

α = Kx(ζ^Int(l/a))
β = Kx(ζ^Int(l/b))
γ = Kx(ζ^Int(l/c))
trW = γ+1/γ
# trW2 = 2 + trW^2
R = Nemo.matrix_space(Kx, 2,2)

# printstyled(styled"{yellow:tr(W)} = $(pretty_printing(trW))\n")
# printstyled(styled"{yellow:tr(W)²-2} = $(pretty_printing(trW^2-2))\n")
A = R([
α BigInt(1);
0 inv(α);
])
B = R([
    β BigInt(0);
    x inv(β);
    ])


# A = Kx[
# 1 BigInt(1);
# 0 1;
# ]


# B = Kx[
#     1 0;
#     x 1;
#     ]




Ai = inv2(A)
Bi = inv2(B)


gens = [f, g]
ms = [A, B]
mats = Dict([l=>r for (l,r) in zip(gens, ms)])



# printstyled((rW)*"\n")

# pretty_printing((trW))|>printstyled;println("\n")

R, Kx, gens, ms, mats, x, ζ, α, β, γ, trW
end;

(a, b, c) = (2,3,7)

R, Kx, gens, ms, mats, x, ζ, alph, bet, gam, trW = ρellip(a,b,c)
l = lcm(a,b,c)
Qbar = algebraic_closure(QQ);
C = CalciumField()
Zq = Nemo.root_of_unity(Qbar, 2*l)
zc = C(Zq);
# cos(C())
# ms = [A, B]
# f = CyclicGen("f")
# g = CyclicGen("g")
# gens = [f, g]
# mats = Dict([l=>r for (l,r) in zip(gens, ms)])
rWl(Q) = (mateval(farey_word(Q,f,g), mats)|>trace)
ρ(w) = (mateval(w, mats)|>trace)
# pretty_printing(rWl(3//5))|>printstyled;println("\n")

# pretty_printing(rWl(1//2))|>printstyled;println("\n")

# pretty_printing(ρ(g))|>printstyled;println("\n")
# Q = @cf0 [2, 6, -2]
# pretty_printing(gcd(rWl(1//2),rWl(1//3)))|>printstyled;println("\n")
# Φ
# complex_embeddings
# em = complex_embeddings(ζ.parent)[1]
# pretty_printing(ρ(f*g*f*g/f))|>printstyled;println("\n")
# CC = em(ζ).parent
# # [polynomial(QQ, coefficients(l), :y)(zc) for l in coefficients(ρ(f*g*f*g/f))]
# # polynomial(CC, [1], :x)
# Nemo.polynomial(CC, [em(l) for l in coefficients(ρ(f*g/f/g))])
# Nemo.polynomial(C, [polynomial(C, coefficients(l))(zc) for l in coefficients(ρ(f*g/f/g))])
# em
# ζ
# to_calc = l->polynomial(C, coefficients(l))(zc)
# PΦ = map_coefficients(AcbField(64)∘to_calc, rWl(4//11)-2)
# pretty_printing(PΦ)|>printstyled;println("\n")
# rts = Nemo.roots(PΦ; target=-2)
P = rWl(4//11)-(2)
dP = derivative(P)
pretty_printing(gcd(P,dP))|>printstyled;println("\n")
pretty_printing(P)|>printstyled;println("\n")
pretty_printing(dP)|>printstyled;println("\n")




p = rWl(2//7).coeffs
roots = [Kx(i) for i in 1:7]
#  monomial_to_newton!(p, roots)
Zx, qx= ZZ["x"];

(a,b,c) = (2,3,5)
pa = cyclotomic(a, qx)
pb = cyclotomic(b, qx)
pc = cyclotomic(c, qx)
# L, Ax = radical_extension(5, QQ(1), "a")
L, bl = Hecke.number_field([pa, pb, pc], "z")
# Ax

b1, b2 = bl
Hecke.basis(L)

# d = Hecke.degree(L)

import Hecke

Zx, qx= Nemo.ZZ["x"];

(a,b,c) = (2,3,5)
pa = cyclotomic(a, qx)
pb = cyclotomic(b, qx)
pc = cyclotomic(c, qx)
# L, Ax = radical_extension(5, QQ(1), "a")
L, bl = Hecke.number_field([pa, pb, pc], "z")
M, c = Hecke.simplified_simple_extension(L)
a = Hecke.gen(M)
# Mmat = Matrix

MM = matrix(QQ, (hcat(Hecke.coordinates.([c(a^i) for i in 0:d-1])...)))


MMpi, Δ = Nemo.pseudo_inv(MM)
MMpi
[MMpi[I] for I in eachindex(MM)]
MMpi * Hecke.basis(L)
# c(a)
[c(a^i) for i in -d÷2:d]


S, X = Nemo.rational_function_field(QQ, :z_0)
RR = matrix_space(S, 2, 2)
M1 = RR([X 1; 0 inv(X)])
pretty_printing(M1^10)|>printstyled;println("\n")


Lx, x = polynomial_ring(Nemo.QQ
, "x");

TY = Lx[
  1 1;1 x
]
# TY.super
TY.base_ring
