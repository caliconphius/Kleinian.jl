
import Hecke:structure_constant_algebra
import Hecke
using Nemo
using Farey
using Farey.Monoids
include(joinpath(@__DIR__, "../../packages/Farey/working/GroupAlgebra.jl"))
import .GroupAlgebras as ZG
Fw(Q) = farey_word(Q, f, g)

function Rw(Q; f=f, g=g)
  typ = Q|>Qtype|>string
  u = typ=="p" ? inv(g) : typ == "v"  ? g : inv(f)
  Θ = rileyw(Q; f=f, g=g)
  f ^-Θ * u
end
Ktr, xs = polynomial_ring(QQ, [:trf, :trg, :trfg])
trf, trg, trfg = xs
f, g = FreeGen("f g")
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

eF = trf*e1 - ef
eG = trg*e1 - eg
eGF = eG * eF
eu = eg*eF
eF*eG

trA(v) = [a*x for (a,x) in zip([2, trf, trg, trfg], v.coeffs)]|>sum
detA(v) = trA(v)*v - v*v
adjA(v) = trA(v) - v
mons = Dict{typeof(f), typeof(ef)}([f=>ef, g=>eg, -f=>eF, -g=>eG])
W = g^-f
W = Fw(1//4)
eQ(Q) = monoideval(Fw(Q), mons)
eW = monoideval(W, mons)
trA(ef*eg*eF*eG)
trA(eW)
eQ(1//4)

zz = 2*ef^2 - trA(ef)^2 + 2
(zz)^2
2*trf * ef - trf^2
(2*ef - trA(ef)) * (2*ef + trA(ef))
detA(eW - e1)
eF = trf*e1 - ef
eG = trg*e1 - eg
eGF = eG * eF
eu = eg*eF

eG * eF - eGF
ef * eg - efg
eF * efg * ef

ef^2
tg = FreeGen("t")
import .GroupAlgebras as ZG
ZG.MagnusCoeff(Fw(1//2), [f,g])
ZG.MagnusCoeff(Fw(1//2), [g,f])
ZG.MagnusCoeff(Fw(1//2), [f,f])
ZG.MagnusCoeff(commu(tg,Fw(1//2)), [g,g])
ZG.MagnusCoeff(commu(tg,Fw(1//2)), [tg,tg])
Abasis = [one(f), f, g, f*g]


w1, wf, wg, wfg =  [ZG.GroupAlgebra([gen*one(f)=>1//1]) for gen in Abasis]
zf, zg = [wf-w1,wg-w1]
zw(w) = w*w1 - w1
w2 = (zw(-Fw(1//2)))

ZG.∂(w2, [f,g])

zg * zf - zw(commu(g,f))*(zf*zg + 1 + zg + zf)
zw(commu(g,f))*(zf*zg+1)+zw(commu(g,f))*(zg+zf)

(ef-e1)^2
θa(Q) = -(f\Fw(Q)|>conj_prefix)[1]


zgW = ZG.GroupAlgebra([
  GEN=>COEFF for (GEN, COEFF) in zip(Abasis, eW.coeffs)
])

eth(Q) = monoideval(θa(Q), mons)
eth(2//5)|>adjA
(eg*ef+ef*eg)
eg - ef * eg * eF
(ef - eg)|>adjA

θa(2//5)

FreeGen("s1")

function braid_rel(i::Int, j::Int)
s1, s2 = FreeGen("s$i s$j")
if abs(i-j) == 1
  return s1 * commu(s2, s1) / s2
end

commu(s1,s2)
end

function wirt(f,g,h)
f * g / f * -h
end

wirt(f,g,tg) * -wirt(-g,tg,f)

braid_rel(1,2)

commu(f,g) * commu(f, tg)^-g
commu(g,tg) * commu(g, f)^-tg
"aba" < "abc"
