using Farey
# loaded = true
if loaded
include("working/kassel_morphs.jl")
include("working/autF2.jl")
include("working/PillowGroup.jl")
using .PillowGroup
loaded = false
end


function (w::Farey.FreeGroupElement)(subs...)
    FGr = w.parent
    SS = Dict([
        [x=>s for (x,s) in zip(GPC.gens(FGr), subs)]...
        [inv(x)=>inv(s) for (x,s) in zip(GPC.gens(FGr), subs)]...
    ])
    prod([SS[FGr(x)] for x in w.word])
end

include("pi1o_tests.jl")
PillowGroup.Klein4(1)
const CF = ContinuedFraction
function Qtype(Q)
    n=numerator(Q);d=denominator(Q)
    if d%2==0
        return PillowGroup.Klein4(1)
    elseif n%2==0
        return PillowGroup.Klein4(2)
    else
        return PillowGroup.Klein4(3)
    end
end

hK = PillowGroup.Klein4(1)
vK = PillowGroup.Klein4(2)
pK = PillowGroup.Klein4(3)

xx = f0 * hh
yy = pp
zz = vv
a0 = xx*yy
c0 = xx*zz


rQ = @cf0 [2,2,3,3]
ContinuedFraction(5//8)
s_seq(rQ)[1:numerator(rQ)]
# rQ = @cf0 [2,1,3,3]
rtyp = Qtype(rQ)
function Qhat(rQ; α=2)
rtyp = Qtype(rQ)
rhat =  rtyp==hK ?  CF(α*rQ) :
        rtyp==vK ?  CF(rQ/α) :
                error("r∉H∪V")
    rhat
end

rhat = Qhat(rQ)
t0Q, tQ, t1Q = farey_neighbours(rhat)

sQ, cQ =  tQ⊕t1Q==rhat ?    (t0Q, t1Q) :
                            (t1Q, t0Q)

s0Q, c0Q = rtyp==hK ?  (CF(sQ/2),CF(cQ/2)) :
         rtyp==vK ?  (CF(sQ*2),CF(cQ*2)) :
                     error("r∉H∪V")

Ω(rQ)
Ω(rhat)
Ω(c0Q)

function S_Palin(Q, a, c)
    Q = CF(Q)

    S1, S2, S3, S4 = a, one(a), one(a), c
    for (k,l) in enumerate(Q.L)
        S1, S2, S3, S4 = k%2==1 ?   ((S1*S2)^(l)*S3,S4, (S1*S2)^(l-1)*S3,S4) :
                                    (S1,S2*(S3*S4)^(l), S1,S2*(S3*S4)^(l-1))
    end
    (S1,S2,S3,S4)
end

S_Palin(CF(3//7),a,b)
S_Palin(CF(3//5),a,b)
CF(0,[2,2])|>Γ
CF(0,[2,1,2,2])|>Γ
CF(0,[2,1,2,2])|>Γ|>G0i
CF(0,[2,1,2])|>Γ|>G0i
CF(0,[2,1,3])|>Γ|>G0i
Q0, Q1, Q2 = farey_neighbours(Q)
Q3, _ = descendQ(Q1, Q0)
Q4, _ = descendQ(Q0,Q3)
Q5, _ = descendQ(Q3,Q4)

G0i
QQ = CF(0,[2,2,2])
Γ(QQ)
QQ = CF(0,[3])
Γ(QQ)|>G0i

f0\Ω(Q)

function sw_descent(Q)
    Q = CF(Q)
    ITR.flatmap(enumerate(Q.L)) do (k,l)
        Lk = Q.L[1:(k-1)]
        [CF(0,[Lk..., l]), CF(0,[Lk...,l+1])]
    end|>collect

end
Q
sw = sw_descent(rQ).|>x->(x,Γ(x))


ψf(Q) = [ϕ^l for (l, ϕ) in zip(CF(Q), ITR.cycle([d1i, g0]))] |> prod
ψfi(Q) = [ϕ^l for (l, ϕ) in zip(CF(Q), ITR.cycle([d1i, g0]))] |> prod |> inv

ψfi(3//10)|>hom|>x->D1i*x
ψf(3//10)|>hom|>x->x*D1


QQ = 30//103
QQ = 1//4
ffi = ψfi(QQ)|>hom
ff = ψf(QQ)|>hom

farey(1//3)|>ffi



sw.|>(x)->ff(x[2])
sw.|>(x)->ff(x[2])

sw.|>first
((ff*G1).(sw.|>last))
((ff).(sw.|>last))
CF(3//10)
riley(3//10)[1]




function aab_dec(rhat)
t0Q, tQ, t1Q = farey_neighbours(rhat)

sQ, cQ =  tQ⊕t1Q==rhat ?    (t0Q, t1Q) :
                            (t1Q, t0Q)
sQ,cQ
end
Q = @cf0 [2,5,2]
aab_dec(3//10)
aab_dec(3//4)

Q = @cf0 [2,2,3]
Qb, Qa = aab_dec(Q)
Qs = sw_descent(Q)
aab_dec(Qa)
Q
Wr(Qa)*Ω(Qb)
Wr(Qa)|>inv
(Wr(Qa)^-b0)^-
(Ω(1//3).elem)(f0,h0)

(Wr(Qa)*b0)^Ω(1//3)
Ω(Qa)
Ω(0//1)
Ω(1//2)
Ω(1//2)^3
Ω(3//5)*b0*c0
Ω(2//5)^2

f0\Ω(2//5)
Ω(0//1)^3
@cf 4//11
@cf0 [2,2,3]

(ψf(3//8)|>hom)(Γ(1//2))

sqbak(w) = (w^(one(f)⋊w.hom)).elem


Ω(1//3)*Ω(2//5)
Ω(1//3)^2*Ω(2//5)
Ω(4//11)/Ω(1//2)
(ψf(4//7)|>hom)(Γ(1//2))
(ψfi(3//5)|>hom)((Ω(1//2)^2).elem)
(ψfi(7//17)|>hom)(sqbak(Ω(1//5))|>inv)
(ψfi(7//17)|>hom)(Ω(1//5).elem)
(ψf(4//7)|>hom)(Γ(4//7))
Γ(4//7)(a0,c0)
Γ(1//2)(a0,c0)


s_seq(4//7)
(ψf(4//7)|>hom)(Γ(1//2))
(ψf(4//7)|>hom)(Γ(3//5))
Ω(3//7)
Ω(1//3)*Ω(1//2)
Ω(0//1)*Ω(1//2)*Ω(1//2)
Ω(2//5)


w = Ω(0//1)*Ω(1//2)*Ω(1//2)*Ω(0//1)
w = Ω(1//3)*Ω(1//2)*Ω(1//2)*Ω(1//3)
Ω(1//2)*Ω(1//2)
Ω(3//7)|>sqr
Ω(2//5)^2
w.elem



Q = @cf0 [2,1,2,2]
sw_descent(Q)
Qb, Qa = aab_dec(Q)
Qa⊕Qb
s_seq(Qa)
s_seq(Q)[1:numerator(Qa)*3+numerator(Qb)]
s_seq(Q)[numerator(Qa)+1:end]



fin =  w -> w * (w|zz)
(Γ(Qb)*Γ(Qa)^2*Γ(Qb))
Γ(Qa)^2

Q = @cf 4//11
Q = @cf 7//17
Q = @cf0 [2,1,1,3,4]
Q = @cf0 [2,1,1,2,3]
# Q = @cf 3//8
Qs = sw_descent(Q)
Qb, Qa = aab_dec(Q)
Qab = Qa⊕Qb
Qs[8]
Qb, Qa = aab_dec(Q)
(Ω(Q)/f0) ^ Ω(Qa)
((Ω(Q)/f0) ^ Ω(Qb)) ^ inv(c0*a0)
zz^a0
(Ω(3//5)/f0)^(Ω(2//3))
(Ω(4//7)/f0)^(Ω(2//3)/f0)
(Ω(3//8)/f0)^(Ω(2//5))
Γ(3//5)|>G0i |> x->x(a0,c0)
Γ(2//5)|>(ψf(3//8)|>hom)
Γ(3//5)|>(ψf(3//5)|>hom)




Q1 = CF(0,[1,1,2])
Q2 = CF(0,[1,1,2,4])
Q3 = CF(0,[1,1])

ψfi(Q3)|>hom
g1/d1|>hom

fh = (ψf(Q2)|>hom)
fhi = (ψfi(Q2)|>hom)
fh1 = (ψf(Q1)|>hom)
fh1i = (ψfi(Q1)|>hom)
Q3, l = descendQ(Q2,Q1)
Q4, l3 = descendQ(Q1,Q3)
Γ(Q1)|>fh
Γ(Q1)|>fh1
Γ(Q2)|>fh
Γ(Q2)|>fh1
Γ(Q3)|>fh
Γ(Q3)|>fh1
Ω(Q3)^l3 * Ω(Q4)
Ω(Q3)*Ω(Q1)^l
Ω(Q2)

(G1i^2*D0)
(D1*G0)
(fhi((a*b)*a))
(fhi((a*b)^2*a))
(a*b)^2

c1 = Ω(Q1)
a1 = Ω(Q2)
n1 = Ω(Q3)
aab_dec(Q2)
Q3 > Q1
s_seq(Q2)


fAdd = (q1,q2,n) -> foldl(⊕, ITR.cycle([q2], n), init=q1)

fAdd(Q3,Q1,l)
Q2


Q = CF(0,[2,1,2])
Q0 = CF(0,[1,1,2,4])
Q1 = CF(0,Q0[2:end-1])
Q2, l = descendQ(Q0,Q1)
Q3, l3 = descendQ(Q1,Q2)

g\(Ψ(Q0)/g)^Ψ(Q1)/g|>conj_prefix
Ψ(Q0)
(Ψ(Q1)|>Qtype(Q1)|>inv)/g

ff = q->ψf(q)|>hom
ffi = q->ψfi(q)|>hom

w = farey(Q0)
(ff(Q0))(w)
D1*ffi(1//2)
ffi(Q0)

Ω(Q)
farey(1//3)|>(G0i*D0i*G0i*D1*G0i)
Ψ(1//3)|>(G0i*D0i*G0i*D1*G0i)
Ψ(1//2)|>(G0i*D0i*G0i*D1*G0i)

(f^g)|>G0*D1i*G0*D0*G0
(f^2*g)|>G0*D1i*G0*D0*G0
(f^2*g)|>G0*D1i*G0*D0*G0
(f^3*g^3)|>G0*D1i*G0*D0*G0
(f/g)|>G0*D0*G1*D1i*G0i
(g)|>G0*D1i*G0*D0
(f)|>G0*D1i*G0*D0

g|>G0*D1i*G0*D0*G0
(G0i*D1*G0i*D0*G1i)(Wr(1//3).elem)
(G0i*D1*G0i*D0*G1i)(Wr(1//3).elem)
Wr(1//4).elem|>G0i|>D0i|>G0i|>D1|>G0i


g\Ψ(Q0)|>conj_prefix
