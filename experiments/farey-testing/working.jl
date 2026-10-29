import Farey

CF = Farey.ContinuedFraction
Farey.farey_word(1//2,f,g)
F2 = Farey.FreeGroup((;:f=>:F, :g=>:G))
f, g = Farey.GPC.gens(F2)
F, G = Farey.GPC.gens(F2) .|> inv
include("../Farey/working/PillowGroup.jl")
using .PillowGroup
include("../Farey/working/autF2.jl")
f0, h0, hh, vv, pp = Farey.GPC.gens(pi1OGroup())
a0 = f0*vv
b0 = hh
c0 = a0*b0
xx = f0 * hh
function (w::Farey.FreeGroupElement)(subs...)
    FGr = w.parent
    SS = Dict([
        [x=>s for (x,s) in zip(GPC.gens(FGr), subs)]...
        [inv(x)=>inv(s) for (x,s) in zip(GPC.gens(FGr), subs)]...
    ])
    prod([SS[FGr(x)] for x in w.word])
end

function sw_descent(Q)
    Q = CF(Q)
    Farey.ITR.flatmap(enumerate(Q.L)) do (k,l)
        Lk = Q.L[1:(k-1)]
        [CF(0,[Lk..., l]), CF(0,[Lk...,l+1])]
    end|>collect

end

Ω(Q) = Farey.christoffel(Q, a, b)(a0,b0)
Γ(Q) = Farey.christoffel(Q, a, b)
Wr(Q) = farey_word(Q, a, b)(f0,h0)
Fw(Q) = farey_word(Q, a, b)(f0,h0).elem
Ψ(Q) = Ω(Q).elem
Ψa(Q) = (one(f)⋊Ω(Q).hom * Ω(Q)).elem
Ω(1//2)
Ψa(3//5)

Ω(Q⊕Qb)/Ω(Qb)
Ω(Q⊕Qb)/Ω(Qb)
Ω(Qa)*Ω(Qb)
Γ(Q)
Ψ(Q⊕Qb)|>ab
Ψ(Q)
hom(w::Farey.FreeGroupElement)::Hom = w==one(w) ? idh : [mapauts[AutGens(x)] for x in w.word] |> prod
hom

idh = Farey.Hom(F2,F2,(a=>a,b=>b))
ψf(Q) = [l>0 ? ϕ^l : one(d1) for (l, ϕ) in zip(CF(Q), ITR.cycle([D1g, G0g]))] |> x->prod(x, init=one(d1))
fp(0//2)

Q = @cf0 [2,1,3]
Qa, Qb,Qc = Farey.farey_neighbours(Q)

fp(Q) = ψf(Q)|>hom
fpi(Q) = inv(ψf(Q))|>hom
ff1 = ψf(Q)|>inv|>hom

Farey.ITR.map(sw_descent(Q)) do x
    x=>(G0^2*fp(Q)*G0^2)(Ψ(x))
end|>collect
D1*G0i*D0i*G0i
Φ = O*D1*G0*D1i*O
Φi = O*D1*G0i*D1i*O

(Φ)(Fw(2//5))
(Φ*G1i)(f*g*f)/
(Φ*G1i)(g*f*g)
(fpi(Qc))(Ψ(Q))
(fpi(Qc))(Ψa(Q))
(fpi(Q))(Fw(Q))
Fw(Q)

(fp(Qa)^2*fpi(Q))(Ψa(Q))|g|f|f
(fp(Qa)^2*fpi(Q))(Fw(Q))

fpi(Q)
fp(Q)(Fw(2//5))
