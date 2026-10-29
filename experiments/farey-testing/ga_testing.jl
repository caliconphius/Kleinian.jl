using Base: product
include("working_setup.jl")
include(joinpath(@__DIR__, "../../packages/Farey/working/GroupAlgebra.jl"))

import Nemo:tr
using Nemo

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
fg_trace = ρw(Fw(0))|>tr

f, g = FreeGen("f g")
import .GroupAlgebras as ZG
cof = Dict([f*g => u, g/f => ZZ(1)])
cof1 = Dict([f => BigInt(1), g => 1])
l = ZG.GroupAlgebra(cof1)
ll = l*l + l
ZG.GroupAlgebra(f)

Fw(Q) = farey_word(Q, f, g)
Θ(3//7)

il = l * im
il^2

Q = @cf0 [2]
W = Fw(Q)
ZG.∂(W, [f, g])
dWx, dWy = ZG.∂(Fw(Q), [f, g])
ZG.∂(dWy, [f,g])
dWy
ZG.∂(W, g)

using Nemo


fgf =  ZG.GroupAlgebra([f=>w, g=>one(w), one(f)=>v])
GFG =  ZG.GroupAlgebra([-f=>one(w), -g=>w, one(f)=>u])
(fgf*GFG)
# ZG.∂(W,f) - ZG.∂(W,f)*f + ZG.∂(W,g) - ZG.∂(W,g)*g


ZG.∂(f^5)
zf(f) = zero(dWx) + one(f) - f

d1 = [zf(f), zf(g), zf(W)]
d2 = [
  ]
d2 * d1
Q = 0
W = Fw(Q)
f, g, h, c, a, b = FreeGen("f g h c a b")
l, m, n = 3,4,5
r = -f^l * h^n
s = -g^m * h^n
t = W/h

rel = [r, s, t]
gens = [f, g, h]

d2 = [ZG.∂(R, G) for (R,G) in Iterators.product(rel, gens)]
d2[end, end ] = one(g)
d2[1,3]
d1 = [zf(G) for G in gens]
d2 * d1

d3 = [ZG.∂(W,f)*zf(f) ZG.∂(W,g)*zf(g) -zf(c)]

Φ = (d3 * d2)

[x for x in Φ]

Φ[1]
Φ[2]
Φ[3].terms
Fw(Q) = farey_word(Q, f, g)
Rw(Q)  = farey_word(Q, f, g)

function Rw(Q; f=f, g=g)
  typ = Q|>Qtype|>string
  u = typ=="p" ? inv(g) : typ == "v"  ? g : inv(f)
  Θ = rileyw(Q; f=f, g=g)
  f ^-Θ * u
end
f ^-rileyw(1//3)
Rw(3//5)

Q = @cf0 [2,2]
rileyw(Q)
Θ(Q)
import Hecke:structure_constant_algebra
import Hecke
Ktr, xs = polynomial_ring(QQ, [:trf, :trg, :trfg])
trf, trg, trfg = xs

RA = Hecke.zeros_array(Ktr, 4, 4, 4)

IA = identity_matrix(Ktr, 4)

FA = Ktr[
  0 1 0 0; #f*I
  -1 trf 0 0; #f*f = -I + trf f
  0 0 0 1; #f*g
  0 0 -1 trf; #f*fg = -g + trf fg
]

GA = Ktr[
  0 0 1 0; #g*I
  -trf*trg+trfg trg trf -1; #g*f = gfg (trg I - g) = -trf trg I + trg f - fg - trf g + trfg I
  -1 0 trg 0; #g*g = -I + trg g
  -trf 1 trfg 0; #g*fg = -trf I + f + trfg g
]

FGA = Ktr[
  0 0 0 1; #fg*I
  -trg trfg 1 0; #fg*f = -trg I + g + trfg f
  0 -1 0 trg; #fg*g = trg fg - f
  -1 0 0 trfg; #fg*fg = trfg fg - I
]
RA[1,:,:] .= IA
RA[2,:,:] .= FA
RA[3,:,:] .= GA
RA[4,:,:] .= FGA
RA
# P = stack([IA, FA, GA, FGA])

A = structure_constant_algebra(Ktr,RA; one=[1,0,0,0], check = false)
e1, ef, eg, efg = basis(A)

trA(v) = [a*x for (a,x) in zip([2, trf, trg, trfg], v.coeffs)]|>sum



eF = trf*e1 - ef
eG = trg*e1 - eg
eGF = eG * eF

eG * eF - eGF
ef * eg - efg
eF * efg * ef
mons = Dict([f=>ef, g=>eg, -f=>eF, -g=>eG])
eW = monoideval(Fw(1//2), mons)
gs = [one(f), f, g, f*g]


Kz, xs = polynomial_ring(Qb, [:z])
Kx, zxz = polynomial_ring(QQ, :z)
zz = xs[1]
Zx, x = polynomial_ring(ZZ, :x)

# Zx, x = polynomial_ring(QQ,:z)
ms_v = Kz[2 2 zz^2 + 2] |> collect |> x->reshape(x, 3)
l, m, n = (4, 5, 7)
lcm(l,m,n)
Ka, za = Nemo.cyclotomic_real_subfield(2*l)
Kb, zb = Nemo.cyclotomic_real_subfield(2*m)
Kc, zc = Nemo.cyclotomic_real_subfield(2*n)
p1, p2, p3 = minpoly(Kx, za),minpoly(Kx, zb),minpoly(Kx, zc)


L, gs = Nemo.number_field([p1, p2, p3], "a")
uu = Hecke.absolute_primitive_element(L)
Lz, zl = L[:z]
a1, a2, a3 = gs
Hecke.absolute_basis(L)
sL, u = Hecke.simplified_simple_extension(L)
ua = u(gen(sL))

ui = inv(u)
c1,c2,c3 = ui.([a1,a2,a3])


pu = minpoly(Kx, ua)

rtu = roots(Qb, pu)
a0 = CC(rtu[1])
A1 = Kx(Hecke.coordinates(c1))(a0)|>Qb
A2 = Kx(Hecke.coordinates(c2))(a0)|>Qb
A3 = Kx(Hecke.coordinates(c3))(a0)|>Qb
Hecke.splitting_field(p3)
A2
# p1(A1)

Zb, zx = polynomial_ring(Qb, :x)

dd, s, t =gcdx(pu, p1)


pu(p1)
t
s
Hecke.embedding(L, sL)
zz = xs[1]

ms_v = Zb[A1 A2 zx^2+2] |> collect |> x->reshape(x, 3)
Q = @cf0 [1,1,2]
Q = @cf0 [2]
ρw(Fw(Q))|>tr

n1, n2, n3 = gs
gcdx(pu, p2)
p1(ua)
eW = monoideval(Fw(Q), mons);
Fp = evaluate(change_base_ring(Qb, trA(eW)), ms_v)
CC = CalciumField()
# change_base_ring(Qb, Fp-A3)|>roots

eW


cm = Hecke.real_embeddings(L)[1]
fn = Hecke.evaluation_function(cm, 64);
rts = map_coefficients(fn, Fp)|>roots

fg_gens = [one(f), f, g, f*g]
import .GroupAlgebras as ZG
zgW = ZG.GroupAlgebra([
  GEN=>COEFF for (GEN, COEFF) in zip(fg_gens, eW.coeffs)
])

AsL, mapsL = structure_constant_algebra(sL)
# AsL|>Hecke.structure_constant_table

minpoly(QQ[:z][1], Hecke.gen(sL))


eF * eG
trfg*e1 - eg * ef


Qb = algebraic_closure(QQ)
Rx, x = polynomial_ring(Qb, "z")
Kx = fraction_field(Rx)
bb = one(x)#root_of_unity(Qb, 3, 1)
aa = one(x)#oot_of_unity(Qb, 4, 1)
q = @cf0 [1,1,2]
Po = ρ(q)
