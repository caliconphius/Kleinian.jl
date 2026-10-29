include(joinpath(@__DIR__, "../../packages/Farey/working/GroupAlgebra.jl"))
  import .GroupAlgebras as ZG
f, g, h, c, a, b = FreeGen("f g h c a b")
l, m, n = 7,3,2
1//l + 1//m - 1//n
q = 1//2
W = Fw(q)
r = h^n*-f^l
s = h^n * -g^m
t = W/h

# W = Fw(1//2)
# rD = f^l
# sD = g^m
# tD = W^n

rel = [r, s, t]
gens = [f, g, h]
include("fareydecomp.jl")
d2 = [ZG.∂(R, G) for (R,G) in Iterators.product(rel, gens)];
d2
Dxf = w->ZG.∂(w,f)
Dxf(Dxf(W))



𝒜 = matrix_space(ZZ, 3,3)
ZG.ϵ.(d2)|>
  𝒜 |> snf|>
  x->x


zf(f) = zero(ZG.∂(f,f)) + one(f) - f

d1 = [zf(G) for G in gens]
d2 * d1

d3 = [ZG.∂(W,f)*zf(f) ZG.∂(W,g)*zf(g) -zf(c)]

d32 = d3 * d2
d32[2].terms
ZG.∂(Fw(1//3), f)

# Δ(p,q,n) <-- F(p,p,q,q,n) --> Δ(p,q,n)
# Δ(p,p,n) <-- F(p,p,q,q,n) --> Δ(q,q,n)
# Γ ℤ


df(w) = ZG.∂(w, f)
dg(w) = ZG.∂(w, g)
w = Fw(3//5)
commu(u, v) = u * v * -u * -v
-commu(f,commu(g,f)^f)
w0 = commu(f,h)*commu(f,g)^-h
ϕ = MonoidHom([h=>-g^-f])
(w0|ϕ)
commu(f,g)^-h
(f)^2 * (-f*g)^2 * (-g^2)
(f*g)^2 * (-g)^2 /(g*f)

jac(u,v;x=f,y=g) = ZG.∂(u, x) * ZG.∂(v, y) -  ZG.∂(u, y) * ZG.∂(v, x)

jac(f,-f^g)
θa(Q) = -(f\Fw(Q)|>conj_prefix)[1]
θa(1//3)
f * h^-θa(3//8)


a, b, c, d = FreeGen("a b c d")
Ta = MonoidHom([
  a=>d*-a*b
])
Tai = MonoidHom([
  a=>b*-a*d
])

Tb = MonoidHom([
  b=>a*-b*c
])
Tbi = MonoidHom([
  a=>c*-b*a
])

Tc = MonoidHom([
  c=>b*-c*d
])
Tci = MonoidHom([
  c=>d*-c*b
])

fd = a/b
gd = b/c
hd = c/d

decomposeWord(3//8)
# Tbi|Ta^2|Tb
W1, W2, d1, d2 = FreeGen("W1 W2 d1 d2")

dWg = ZG.∂((W1*W2/(W1^f)))
dWg = ZG.∂(((-(W1^-f))*W2*(W1)))
(dWg[W2]*d2 + dWg[W1]*d1 + dWg[f])

dWg = ZG.∂(W1*W2/(W1^f))
(dWg[W2]*d2 + dWg[W1]*d1 + dWg[f])
W1*W2*-f*(zero(f+f) - one(f) + inv(W1))

KZ = ZG.∂(W1^-W2)
KZ = ZG.∂(f*W1*-f*-W1)
(KZ[W1]*d1 + KZ[W2]*d2)
ww = W1*W2/(W1^f)
ww = (W1*W2/W1/W2)
zero(f+f) + one(f) - f
(ZG.∂(ww,[W1,W2,f]).*[d1,d2, one(f)])|>sum


df(Fw(1//2)*W1)
dg(Fw(1//2)*W1)

function cayleyΩ(Q::Number; a=FreeGen("a"), b=FreeGen("b"))
  Q = CF(Q)
  n, d = numerator(Q), denominator(Q)
  K = n+d
  mapreduce(*, 0:K-1) do l
    (l*n)%K > ((l+1)*n%K) && return b
    a
  end
end

cayleyΩ(3//5)

Qhat(3//8)

eoe = (lo^le).rep
eoe
s1, u1 = FreeGen("s u")
hl = MonoidHom([f=>s1, g=>s1\u1])
g1
hli = MonoidHom([s1=>f, u1=>f*g])
Q = 1//3
Q = 2//5
w = (Fw(Q))
(w|hl)

uu = i->u1^(-s1^i)
(uu(1)/uu(2)/uu(0))|hli
(uu(1)*uu(3)/uu(4)/uu(2)/uu(0))|hli
(uu(1)*uu(3)/uu(2)/uu(0))|hli
(uu(1)/uu(0)^2*uu(-1)/uu(0))|hli
(uu(1)/uu(0)^2*uu(-1))|hli

Q = @cf 2//7
Q = @cf0 [2,2,2]

Γ(Q)
W = Fw(Q)
decomposeWord(Q)
δ = Fw(2//5)
W

R, = Nemo.residue_ring(ZZ, 11)
-1/R(4)
Fw(2//5)
@cf0 [3,1,2]
gcdx(4,10)
