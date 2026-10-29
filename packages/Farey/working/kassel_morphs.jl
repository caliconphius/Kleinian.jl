using Farey
↑(ϕ, n) = foldl(|>, ITR.take(ITR.cycle([ϕ]), n))

F2 = FreeGroup(:a, :c)
id = one(F2)
a, b = GPC.gens(F2)

A, B = GPC.gens(F2) .|> inv
farey(q) = farey_word(q, a, b)
Γ(p, q) = christoffel(p, q, a, b)
Γ(Q) = christoffel(Q, a, b)
E = Farey.Hom(F2, F2, (
    a => b,
    b => a,
))


D0 = Farey.Hom(F2, F2, (
    a => b * a,
    b => b,
))

D1 = Farey.Hom(F2, F2, (
    a => a * b,
    b => b,
))

G0 = E ∘ D0 ∘ E

G1 = E ∘ D1 ∘ E



D0i = Farey.Hom(F2, F2, (
    a => b \ a,
    b => b,
))

D1i = Farey.Hom(F2, F2, (
    a => a / b,
    b => b,
))

O = Farey.Hom(F2, F2, (
    a => inv(a),
    b => b,
))

O ∘ D0 ∘ O

F2fg = FreeGroup(:f, :g)
f, g = GPC.gens(F2fg)
F, G = GPC.gens(F2fg) .|> inv

G0i = E ∘ D0i ∘ E
G1i = E ∘ D1i ∘ E
((G1i ∘ D1)↑6)(a)
((G1i ∘ D1)↑6)(g)
ϕ = (G1i ∘ D1)
# ϕ = (G0∘D0)
ψ = (D1i ∘ G1)
# ψ = (D1i∘G1)


# w = (ϕ↑3)(f)
# ϕ(w) * ψ(w)


ψ↑2

# u, v = (ϕ↑(6-k)).([f, G])
# u, v = (ϕ↑(k)).([f,g])
k = 3
u = (ϕ↑(k))(b)
u = (ϕ↑(k))(a)
v = (ψ↑(k))(B)
v = (ψ↑(k))(b)
u * v

Q = @cf0 [2, 2]
farey(Q)
ϕ(b)

(ψ↑(k))(v)
(ϕ↑(k))(inv(v)) * (ψ↑(k))(u * v)

# PHI IS PERIODIC OF PERIOD 6 (in outF2)!!
# 1\    /2
#   \  /
#    \/3
#    |
#    |
#   /\4
#  /  \
# /5   \6


bkl = G0 ∘ D0i ∘ G1
bkli = G1i ∘ D0 ∘ G0i
G0 ∘ G1i
D0i ∘ G0 ∘ G1i ∘ D0
bkl ∘ G0 ∘ bkli

Delta = G0 ∘ D0i ∘ G0 ∘ G1 ∘ D0i ∘ G0
Delta

bkl ∘ Delta ∘ bkli

w0 = farey(2 // 3)
w1 = farey(1 // 4)
w = farey(1 // 2)
w = Γ(2 // 3)
((G1i ∘ D1i↑1)).([w0, w1])
((D1↑3 ∘ G1)).([w0, w1])
((D0↑3 ∘ G0i))

D0 ∘ O

reduce(|>, [G0, D1i↑3])
joined = reduce(∘, [
    D0i↑2,
    G0,
    D1i↑2,
    G1
])

w0 = farey(3 // 5)

θ = (E |> O)↑2

joined = foldl(|>, ITR.cycle([
        O, G0, D0, G1, O
    ], 1))

joinedi = reduce(|>, [
    G1, D1, G0
])


joined(a * b * a * B * A * B)
joined1 = reduce(∘, [
    O, D1, G1i↑2, D0, O
])


u = (a) |> joined1
v = (joinedi(b))

joined.([a, b])
joined(a * B * a * B * a)
joined(a * B * a * B * a)
joined(A * b * A * b * A * b * A)
k = 6

joined1(a)
((D0↑4 ∘ G0) ∘ D1i)
((D0↑4 ∘ G0) ∘ E ∘ (D1↑3 ∘ G1))

G1i ∘ D1
ϕ↑3


(G1i ∘ D1i ∘ G1i ∘ D1i↑2)
(G1i ∘ D1i ∘ G1i ∘ D1i↑2)
(D1↑2 ∘ G1 ∘ D1)(w)
