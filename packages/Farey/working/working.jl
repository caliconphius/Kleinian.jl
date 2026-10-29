using Farey
using LinearAlgebra
# computing in FreeGroups/Monoids!

F2ab = FreeGroup(:a, :b)
F2 = FreeGroup(:f, :g)
F2n = FreeGroup(:u, :v)
id = one(F2)
idab = one(F2ab)
f, g = GPC.gens(F2)
F, G = GPC.gens(F2) .|> inv
a, b = GPC.gens(F2ab)
A, B = GPC.gens(F2ab) .|> inv
farey(q) = farey_word(q, f, g)
Γ(p, q) = christoffel(p, q, a, b)
Γ(Q) = christoffel(Q, a, b)
# riley(Q) = 
w = A*farey(3//5)



s_seq(3//5)
u = christoffel(3//5,a,b)
u = christoffel(1//10,a,b)
palindrome(f,g,10)
palindrome(f,g,9)


Q =(3//5)

Q = @cf 31//278
Q = @cf0 [1,1,1,1,1,1]
farey_word(Q, f, g)
s_seq(Q)


Q = @cf 3//8
W3_8 = farey(Q)

A3 = palindrome(f,g,3)
A2 = palindrome(f,g,2)
An(n) = palindrome(f,g,n)
n = 4
nielsen = Farey.Hom(F2, F2, (
    f=>g^f,
    g=>f^-1
)) 

ϕ = Farey.Hom(F2, F2, (
    f=>f*g,
    g=>F
)) 

ψ = Farey.Hom(F2, F2, (
    f=> G,
    g =>g*f 
))

(ϕ∘ψ)

Q
h = Farey.Hom(F2, F2, (
    f=> f,
    g => g*f 
))
farey(Q)
h(farey(Q))<<4

switch = Farey.Hom(F2, F2, (f=>g, g=>f))

Q = @cf 5//7
Q0 = @cf0 [3, 2]
Γ(Q)

christoffel(Q, f, g)
farey(Q)
F2ab(1)
((christoffel(Q, F2n(1), F2n(2))|>ϕ) << 7 ) 
Q |> farey_neighbours
s_seq(Q)'

MM(F2n(1))

NN(f*g)

k_iter = 8
n2 = accumulate(∘, repeat([nielsen], k_iter))|>collect


nielsen_i = Farey.Hom(F2, F2, (
    f=>g * f, 
    g=>f
)) 

n_iter = n2[end]
[hom.([f,g]) .<< k_iter for hom in n2]

nielsen_i.(nielsen.([f,g]))



ϕ = Farey.Hom(F2, F2, (
    f=>An(n+1),#*(A(n)^-1),
    g=>An(n)
)) 

ff = ϕ(f)
gg = ϕ(g)

gg^-1 * ff
newword(f,g) = g^-1 * f

newword(ff,gg) 



((gg*ff)^-1 * gg) 

# ψ = Farey.Hom(F2ab, F2ab,(
#     a=>a
# ))


ϕ.image
ϕ(f*g*f^(g^2))
ww = f*g*f^(g^2) 
Γ(3//5)
christoffel(3, 5, a, b^a)<<5
christoffel(3, 8, a, b)
christoffel(3, 8, a, b^a)<<10

Q = @cf0 [1, 2, 4, -2, -1]
Q = positive_form(Q)
q = shrink_cf(Q)

farey(Q)
Γ(q)

(q1, q2, q3) = farey_neighbours(Q)

Γ(q1), Γ(q2)

W = farey(1//2)
ϕ(W)<<5

ww = ϕ(f^inv(g)) 
ϕ(f^inv(g)) ^ A(4)
ϕ(g)
ϕ(f)

(f*g)


Γ(3, -4)
