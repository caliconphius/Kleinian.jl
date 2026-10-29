include("working_setup.jl")
include("fareydecomp.jl")
include("sl2algebra.jl")
# include("field_calc_hecke.jl")
include("triangle_rewriting.jl")
include("erics_functions.jl")

f

xxp = f0*hh
f0^xxp
g0p = f0^vv
FareyDecomp!(Q)
funnel(24//103).|>x->(x.qP, qpair(x.qP))
Ωa(Q) = christoffel(Q, f0*vv, hh)

Ωa(1//2)^Ωa(2//3)
(f0\Ωa(1//2), f0\Ωa(2//5), f0*Ωa(1//3))
(f0*Ωa(1//3))^2
commu(a0, b0)

Q = @cf 24//103

p1, p2 = expandhom(Q)
p1*p2
qns = qnbs(Q)

f, g, w, ws... = FreeGen("f g W0 " * join(["W$i " for i in eachindex(qns)]))
(Ωa(2//5))^2
a0^vv
f0 ^ Λ(2//3)
g0p ^ Λ(2//3)
Λ(2//3)
a0^zz
b0/a0*b0
f0*g0p

a0^yy
decomposeWord(4//7)
Ω(1//2)^(Ω(1//2)Ω(1//1))
-Ω(1//2)^f0
Ω(1//2)^Ω(1//1)
Ω(1//3)^Ω(2//5)
Ω(1//3)^Ω(1//2)

decomposeWord(CF0([2,3,3,2]))
Ω(10//23)^Ω(3//7)
@cf 7//16
ffh(1//3)


Q = @cf0 [1,4,2,3]
FareyDecomp!(Q)
Qr, Qt = qpair(Q)

ffhi(Q)(Γ(Qt))
ffhi(Q)(Γ(Qr))
ffhi(Q)(Γ(Q))
Γ(Qt)
Γ(Qr)

b2 = CyclicGen("b0";order= 2)
b
b2
pl = MonoidHom([
    fb => b2,
])

Γ(1//2)
Λ(1//2)
xx^zz

lg1

Mn(n) = (ψf(CF0([n])))
(Mn(2))

m3 = -lx*ψf(1//2)

(-m3).rep(farey_word(1//3,a,b))
wQ = farey_word(1//3,a,b)
(-Mn(3)).rep(wQ)

λ(Q) = a\Γ(Q)/b
λ(3//8)
λ(4//7)
λ(7//12)
λ(5//14)
a0^zz
θa(4//11)[1:5]
θa(7//12)
θa(5//18)
Λ(4//11)
@cf0 [1,1,2,3]

1//2 < 7//12 < 10//17 < 3//5
(a0*b0)^zz
Q = 3//5
θa(Q) = -(f\Fw(Q)|>conj_prefix)[1]
λw(Q) = θa(Q)[1:(denominator(Q)-1)÷2]
pλw(Q) = monoideval(λw(Q), Dict([f=>xx*yy*zz,-f => zz*yy*xx,g=>zz*xx*yy,-g => yy*xx*zz]))
cλw(Q) = monoideval(λ(Q), Dict([a=>xx*yy,-a => yy*xx,b=>zz*yy,-b => yy*zz]))
λw(1//10)
λw(2//5)
λw(3//5)
rev_conj(3//5)
rev_cong

cλw(5//18)
θa(3//8)
ww = xx^hh
ww^-pλw(3//8)
ww^-pλw(3//16)

(xx^zz)^-pλw(5//18)
θa(5//18)
((a0*b0*a0)^2)^zz




rev_conj(5//7)
f\Fw(5//7)|>conj_prefix

-lg1

(lg1).rep(Fw(1//5))

uw(i) = FreeGen("u$i")
s1 = g/f
umap = MonoidHom([
  uw(i)=>s1^(-f^i) for i in 0:30
  ])

function torus_word(n::Integer)
  up = one(uw(0))
  un = one(uw(0))

  for i in 0:n-1
    if isodd(i)
      up = up * uw(i)
    else
      un = uw(i)\un
    end
	end

	up * un
end

uw(3)^(-uw(1))|umap

torus_word(4)|umap
decomposeWord(2//7)
expandhom(2//7)

-(g*commu(f,g)/f)
commu(f,g*f*g)

f * -f^-g
