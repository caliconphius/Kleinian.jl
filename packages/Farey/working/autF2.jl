include("kassel_morphs.jl");


AutGens = FreeGroup(:D0, :D1, :G0, :G1, :E, :O)
d0, d1, g0, g1, e, o = AutGens|>GPC.gens
d0i, d1i, g0i, g1i, ei, oi = AutGens|>GPC.gens.|>inv
mapauts = Dict([
    AutGens(1)=>D0, AutGens(2)=>D1,
    AutGens(3)=>G0, AutGens(4)=>G1,
    AutGens(5)=>E, AutGens(6)=>O,
    AutGens(1)|>inv=>D0i, AutGens(2)|>inv=>D1i,
    AutGens(3)|>inv=>G0i, AutGens(4)|>inv=>G1i,
    AutGens(5)|>inv=>E, AutGens(6)|>inv=>O,

])



riley(Q, a, b) = conj_prefix(a\farey_word(Q,a,b))
riley(Q) = riley(Q, a, b)


Q = @cf0 [2,1,2]
farey_neighbours(Q).|>s_seq
s_seq(Q)
riley(Q)|>prod
farey(Q)

A*Γ(Q)*B



w = riley(2//5, d0, g1i)
hom(w) = [mapauts[AutGens(x)] for x in w.word] |> prod

Q = 3//5

auta = d0^-1; autb = g1
w = christoffel(Q, auta, autb)
ff1 = w|>hom
ff2 = w|>inv|>hom

comm = a * A^B

ff1(a)

(ff2)(comm^a*comm)
(ff1)(farey(1//4))
ff2((ff1)(a)^b)
(ff1)(a)^b
# ab(g) = mapreduce(
#     x->
#     Int(x)==1 ? [1,  0] :
#     Int(x)==2 ? [0,  1] :
#     Int(x)==3 ? [-1, 0] :
#                 [0, -1]
#     , + , g.word
# )
riley(3//8, f, g)
riley(3//8, g, f)

G1 * D1i

ϕ = (G1i ∘ D1)
ψ = (D1i ∘ G1)
ϕ0 = (G0∘D0i)
ψ0 = (D0∘G0i)
n = 4



ϕupn(n) = G1i↑n ∘ D1
ψupn(n) =  D1i∘G1↑n

(D0i^2 * G0^2)
(D0i^2 * G0^2)
D1i^n * G1

function ab(ϕ::Farey.Hom)
    Gr = ϕ.dom
    [ϕ(x)=>ab(ϕ(x)) for x in Gr|>GPC.gens]
end

G1*(D0^2*G0)
G1*(D0*G0)
(G0i*D0i^2)


(O^3 * E*O)

eoe = E*O*E

D0
G0i



G0



Q = 3//5

auta = g1; autb = d1^-1
w = christoffel(Q, auta, autb)

wrd(n) = g1 * d0^n * g0
N = 1

wrd(0) |> hom
ψ = wrd(N) |> inv |> hom

(d1)|>hom


σ(m) = Farey.Hom(F2, F2, (
    a => a^(m+1) * b,
    b => a^(m) * b,
))

σi(m) = Farey.Hom(F2, F2, (
    a => a/b,
    b => (a/b)^m \ b,
))

σ1 = σ(1)
Σ1 = σi(1)

Q = @cf0 [1, 1, 2]
Q1 = @cf0 [1,1,1, 2]

σ1(Γ(Q))
s_seq(Q1)

farey_neighbours(Q1)
w, w1, w2, w3 = (Q1, farey_neighbours(Q1)...).|> Γ

w2 |> Σ1
w1 |> Σ1
Q
# farey(Q) |>

D1*G0
(D1*G0) * (D1*G0)

s1 = d1 * g0
s2 = d1 * g0^2

w = s1^2


w * g0 |> hom
Φ = w |> inv |> hom

Γ(3//5) |> Φ
Γ(2//3) |> Φ

Psia(Q) = farey(Q)[1:denominator(Q)]
Psia(1//2)

vhom = E
hhom = eoe*O
phom = vhom*hhom


