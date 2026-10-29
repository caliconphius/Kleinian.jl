using Farey, Farey.Monoids
# using Transducers
includet("MonoidRep.jl")
import .MonoidReps
l, m, n = 2,3,7
# RΔ = KB.ExampleRWS.triangle(l,m,n)
# Δlmn = MON.FPMonoid(RΔ)
fa, fb = FreeGen("a b")
ms = Dict([
  fa=>a,
  fb=>b,
  -fb=>b^(m-1),
  -fa=>a^(l-1),
  ])

# rh(w) = monoideval(w, ms)

function cayleyΩ(Q::Number; a=FreeGen("a"), b=FreeGen("b"))
  Q = CF(Q)
  n, d = numerator(Q), denominator(Q)
  K = n+d
  mapreduce(*, 0:K-1) do l
    (l*n)%K > ((l+1)*n%K) && return b
    a
  end
end
# keys(ms) |> Filter(ispositive) |> collect
θa(Q) = -(f\Fw(Q)|>conj_prefix)[1]

G0 = MonoidReps.MonoidRep(ms; identity=one(a));
# ga,gb = G0(fa),G0(fb)
x, y = FreeGen("f g")
f, g = FreeGen("f g")
# farey_word(2//5,ga,gb)
#
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

Γ(Q) = christoffel(Q, FreeGen("a b")...)

il = oo*ee*oo*ee

idhom = oo^2
G,D,G1,D1, E, O = FreeGen("G D G1 D1 E O")

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

idict[ee]
maphom = merge(mapd,Dict([
  -k=>idict[v] for (k,v) in pairs(mapd)
]))
z1, z2 = FreeGen("z1 z2")
w1 = ZG.GroupAlgebra(z1/z1)
cc = FreeGen("c")
zg_map = Dict([
  f=> w1 + z1,
  -f=> w1 - z1 + z1^2 + cc^10,
  -g=> w1 - z2 + z2^2+ cc^10,
  g=> w1 + z2,
  ])

function truncate_zg(w::ZG.AbstractGroupRing; el = cc, L=4)
	ZG.GroupAlgebra(
	Dict(
	filter(
   	l->(
   	  el∉l.first &&
    		length(l.first) ≤ L
   	),
  	pairs(w.terms)
  	)
	))
end

function truncate_zg(w::MonoidReps.RepElement; el = cc, L=4)
	ZG.GroupAlgebra(
	Dict(
	filter(
   	l->(
   	  el∉l.first &&
    		length(l.first) ≤ L
   	),
  	pairs(w.rep.terms)
  	)
	))
end
chrisMorphs = MonoidReps.MonoidRep(maphom; identity=idhom);
magnus_embed = MonoidReps.MonoidRep(zg_map; identity=w1);
dc = magnus_embed(cayleyΩ(2//5;a=f,b=g))
trc(dc, k) = truncate_zg(dc; L = k) - truncate_zg(dc; L = k-1)
truncate_zg(ZG.MagnusCoeff(commu(g,f),[g,f]))
ZG.MagnusCoeff(commu(f,g),[f,g,f])
ZG.MagnusCoeff(commu(f,g),[g,f,g])
function compute_magnus(w, n)
  out = []
  for x in Iterators.product([[f=>z1,g=>z2] for i in 1:n]...)
    push!(out, (
      ZG.MagnusCoeff(w, [g.first for g in x])|>(ZG.ϵ) ,
      ZG.MagnusCoeff(w, [g.first for g in x]),
      prod(reverse(last.(x))),
      n
    ))
  end
  filter(!iszero∘first, out)
end

Δ(w,f,g) = ZG.MagnusCoeff(w, [f,g]) - ZG.MagnusCoeff(w, [g,f])
Δ(Fw(1//2)^2, f,g)
compute_magnus(Γ(1//3), 1)

function ∧(u,v)
	u * v - v * u
end
wf, wg = ZG.GroupAlgebra.(FreeGen("f g"))
(wf - 1)∧(wg-1)∧(wg-1)*1//2 +
(wg - 1)∧(wf-1)∧(wf-1)*1//2
(w1*z1)∧z2
# commu(-g,-f)|>w->println.(Iterators.flatten([compute_magnus(w, i) for i in 1:4])|>collect);
K = 4
tword = commu(commu(-f,-g),-g)
tword = commu(f,g)^2
tword|>
w->println.(Iterators.flatten([compute_magnus(w, i) for i in 1:K])|>collect);
truncate_zg(ZG.MagnusCoeff(commu(f,g),[f,f,g]))
wr = Fw(1//3)
Δ(g/f, g,f)
Δ(g/f, f,g)
Δ(g^(f^2), f,g)
dc
zf = magnus_embed(f)
zg = magnus_embed(g)
zf ^ zg
zf ^ -zg
q = 3//8
f * θa(q)[1:denominator(q)÷2-1]
θa(5//12)
te = (ee*oo)^2
(Fw(3//8))|>reverse
(Fw(3//8))
k = -1
CF0([2,2,-2])
CF0([2,2,-2,-2k])
CF(11//28) +
((5)*17)%28
17//28|>CF
q = CF0([3,3,5])
q = CF0([1,2,5])
q = CF0([4,-2])
q = CF0([1,2,1,5])
1//(4-11//16)
CF0(q.L|>reverse)
ZG.∂(a\cayleyΩ(4//11)/b)[-a]
cayleyΩ(4//11)
chrisMorphs(G*E*O*D*E)

truncate_zg(dc)
g ∈ f^2

function make_evens(L::Number)
	q = ContinuedFraction(L)
	if !iseven(denominator(q)*numerator(q))
	  q = -(CF0(q.L))

	end
	q = ContinuedFraction(q)
	ll = length(q)
	remains = 0
	s = 1
	Lnew = BigInt[]
	for (k,p) in enumerate(q.L)
	    if remains==1

					p = s*p>0 ? s+p : s-p
					remains = 0
			end
			if isodd(p) && remains==0
	      s*p>0 ? push!(Lnew, s+p) : push!(Lnew, s-p)
				remains = -1
				s = -s
			elseif remains==-1
			  append!(Lnew, vcat([[2s,-2s] for i in 1:(p÷2+1)]...)[1:p-1]...)
				s = isodd(p) ? s : -s
				remains = 1
			else
        s*p>0 ? push!(Lnew, p) : push!(Lnew, -p)

			end
	end
	CF0(Lnew), CF(L)
end

(4 )//(5*11)
4 * 16


cayleyΩ(4//11)
(CF0([2,2,-2,2, -2 ,-2,2,-2,-2]))
4//11 - 87//242
8 * 4 - 11*3
gcdx(11,87)
8 * 11 - 87
4//11 + 1//88
(4 - 8 * 11)//242
(242 * 4 - 87*11)

farey_neighbours(87//242)
@cf0 [-2,2,-2,-2]
[(K,(4//11 - CF0([2,2,-2,2,K,-2,2,-2,-2]) )|> CF,CF(1)* CF0([2,2,-2,2,K,-2,2,-2,-2]) ) for K in -3:3]
qq = 1//2
9*4%11
make_evens(CF0([3,2,2,2]))
make_evens(4//11)
CF0([-2,5])
make_evens(CF0([1,1,4]))
CF0([1,1,4])
p,s = 5,1
vcat([[2s,-2s] for i in 1:(p÷2+1)]...)[1:p-1]
qq
CF0([-2,2,-2,1])
CF0([2,2,4,-2,2,-2,1])
(1-CF0([2,2,3]|>reverse))|>CF

CF0([3,9,3])
CF0([3,8,3])
CF0([4,-2,2,-2,2,-2,2,-2,2,-1,0,-3])
CF0([4,-2,2,-2,2,-2,2,-2,2,-1,0,-3])
CF0([4,-2,2,-2,2,-2,2,-2,1,0,3])

CF0([4,-2,2,-2,2,-2,2,-2,4])
CF0([4,-2,2,-2,2,-2,2,-2,2,-4])
CF0([3,3,2,3])
CF0([4,-1,-2,-2,-3])
CF0([4,-2,1,1,2,3])
CF0([4,-2,2,-4,2,-2])
-CF0([-2,0,-2,-3])

CF0([4,-1,-2,-3,-2])
qq = CF0([4,-2,2,-4,-2])
qq = CF0([4,-2,2,-4,2,-2])
1+CF0(reverse(qq.L))
(24 * (79-56))%79
lg, ld, lg1, ld1, le, lo = [chrisMorphs(X) for X in FreeGen("G D G1 D1 E O")]
commu(f,g) = f / (f^-g)
commu(lg, ld)^lg
commu(ld1, -lg1)
commu(ld, lg)

farey_word(3//5,ld1,-lg)
farey_word(1//4,ld,-lg)
farey_word(1//4,ld1,-lg1)
farey_word(1//5,ld,lg1)

# θa(1//5)

lx = lg/lg1
ly = ld/ld1
lH = ld1/lg1
(le*-lH)
le*lo*le*lo


# A = KB.Alphabet([:f, :g, :t, :F, :G, :T])
# for pair in pairs([:f=>:F, :g=>:G, :t=>:T]|>Dict)
#     KB.setinverse!(A, pair.first, pair.second)
# end

# l,m,n = (2,3,7)
# w = farey_word(5//12,ga, gb)

# f = KB.Word([1])
# F = KB.inv(f, A)
# g = KB.Word([2])
# G = KB.inv(g, A)
# t = KB.Word([3])

# R = KB.RewritingSystem(
#     [(f^l, t),
#     (g^m, t),
#     ((f*g*F*G)^n, t)],
#     KB.LenLex(A),
# )

# taskRkb = Threads.@spawn KB.knuthbendix(R)

# @info taskRkb

# Rkb = taskRkb.result
# Grp = MON.FPMonoid(Rkb)
# Rkb
# f1, g1, t1, F1, G1, T1 = GPC.gens(Grp)
# t1
# (f1*g1)*(F1*G1)*t1
# F1*G1
# f1^4 * g1^1
lw = lg^2*ld1*lg^2
# lw = lg*ld1^2*lg*ld1^2
@cf0 [2,1,2,2,-2,-1,-2,2,2,1,2]

ab_gens = Dict([f=>[1,0],g=>[0,1],-f=>[-1,0],-g=>[0,-1]])
function ab(w::T; ggs=ab_gens) where T<: Farey.AbstractMonoidElement
  sum([ab_gens[l] for l in w],init=[0,0])
end

lw.rep.([f,g]).|>ab.|>l->CF(l[2]//l[1])
(lw).rep(g)



(lw*ld1^2*-lw*ld^2*lw).rep.([f,g]).|>ab.|>l->CF(l[2]//l[1])

lw^2
lw

# W1, W2, d1, d2 = FreeGen("W1 W2 d1 d2")
# dWg = ZG.∂((W1*W2/(W1^f)))
# dWg = ZG.∂(((-(W1^-f))*W2*(W1)))
# (dWg[W2]*d2 + dWg[W1]*d1 + dWg[f])

# dWg = ZG.∂(W1*W2/(W1^f))
# (dWg[W2]*d2 + dWg[W1]*d1 + dWg[f])
# W1*W2*-f*(zero(f+f) - one(f) + inv(W1))

# KZ = ZG.∂(W1^-W2)
# KZ = ZG.∂(f*W1*-f*-W1)
# ww = W1*W2/(W1^f)
# ww = (W1*W2/W1/W2)
# (KZ[W1]*d1 + KZ[W2]*d2)
# zero(f+f) + one(f) - f
# (ZG.∂(ww,[W1,W2,f]).*[d1,d2, one(f)])|>sum

θa(3//8)/f/θa(3//8)
rileyw(3//8)
qq = @cf 4//11
qq = @cf0 [2,1,4,5]
qq = @cf0 [2,1,5]
qr, qt, k = CF(qq)|>qpair

Fw(qr)|>signs_seq
Fw(qt)*g*f|>signs_seq
signs_seq(f*g)
(Fw(qq)^(Fw(qr)^(k)*Fw(qt)))
# -Fw(qq)^(Fw(qr)^k)
(Rw(qq))
θa(qq)
Fw(2//7)
Rw(2//7)
θa(2//7)

CF(2//7)
CF(20//23)
CF(0,[1,3,3])


CF(1//7)
ϕ = MonoidHom([
  f=>g\-f*g,
  g=>g
  ])
farey_word(1//3, f, g)|ϕ
fta = MonoidHom([f=>FreeGen("a"), g=>FreeGen("b")])
g0i(cayleyΩ(3//5;a=f,b=g))
cayleyΩ(3//5;a=f,b=g)
f2 = CyclicGen("f2"; order=2)
g2 = CyclicGen("g")

cyc_mor = MonoidHom([f=>f2, g=>g2])
JR_mor = MonoidHom([f=>g2^-f2, g=>g2])
(Fw(1//5)^2)|cyc_mor
(Fw(2//5)|JR_mor)

Fw(3//8)|> signs_seq
Fw(3//4)^2 |> signs_seq
Γ(3//4)
d1
ψf(Q) = [ϕ^l for (l, ϕ) in zip(CF(Q), ITR.cycle([ld1, lg]))] |> prod
ψe(Q) = [ϕ^l for (l, ϕ) in zip(make_evens(CF(Q))[1].L, ITR.cycle([ld1, lg]))] |> prod |> inv
ψe(CF0([2,3,2]))
ffh(Q) = ψf(Q).rep
ffhi(Q) = (-ψf(Q)).rep
ψf(3//8) *  -ψf(1//2)
ψf(2//5) * -ψf(3//5)
palm(k) = -ψf(CF0([1,k]))
palm(4)
(palm(2)).rep(Fw(1//4))|(inv(palm(2)).rep)
(palm(2)).rep(Fw(1//4))
[(palm(6)).rep(Fw(1//k)) for k in 1:10]
(Fw(4//11))
(palm(2)).rep(Fw(1//3))

Fw(1//2)
