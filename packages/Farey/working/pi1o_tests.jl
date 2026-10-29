

f0, h0, hh, vv, pp = GPC.gens(pi1OGroup())
a0 = f0*vv
b0 = hh
c0 = a0*b0
xx = f0 * hh


Ω(Q) = christoffel(Q, a, b)(a0,b0)
Wr(Q) = farey_word(Q, a, b)(f0,h0)
q = 3//5
(Ω(q))
(Ω(q)) 
vv*inv(a0)*vv

# (vv | (f0*h0)) * f0
Ω(2//3)



a0^2 * vv / a0^2

(f0*pp)

christoffel(q, a0, b0)
π1 = pi1OGroup()

SS = Dict([
    0x1=>a0,
    0x2=>a0*b0,
    0x3=>a0^-1,
    0x4=>b0/a0,
])

p1o(w) = prod([SS[x] for x in w.word])
Ψ(Q) = Ω(Q).elem
Ψ(1//2)
Q = @cf0 [2,1,2]
Ω(Q)/f0

(a0^3*b0)

Ω(1//2)^3 * h0

Wr(1//2)|>inv
Ω(1//2) 

xx 
(h0 * hh)
hh
xx*vv*pp

Ω(3//8)
Ω(1//3)^2 * Ω(1//2)
riley(Q, a, b) = conj_prefix(a\farey_word(Q,a,b))
riley(Q, a, b) = conj_prefix(a\farey_word(Q,a,b))
riley(Q) = riley(Q, a, b)[1]
riley(3//5)
prod([π1O(x) for x in riley(3//5)[1].word])
Riab(Q) = riley(Q)(a0,c0)
Ri(Q) = riley(Q)(f0,h0)
Λ(Q) = f0\(Γ(Q))(a0,b0)
Λ(3//8)
Λ(3//5)

farey(3//5)(a0,c0)
farey(3//8)(a0,b0)

n = 10
Q = @cf0 [2,2,3, 4]
Q1, Q2, Q0 = farey_neighbours(Q)
rw = 
Q1
Q0
Q
Q2>Q1
descendQ(Q, Q1; l=0) = (foldl(⊖, ITR.cycle([Q1], Q[end]+l), init=Q), Q[end])

Q1⊖

Q
Q1
Q2, l = descendQ(Q, Q1)
Q3, l = descendQ(Q1, Q2)

Γ(Q1)
Γ(Q2)^3 * Γ(Q3)

Γ(Q)
s_seq(Q)
s_seq(Q1)
a0^3
a0*b0

Γ(2//5) |>  G0i
Γ(3//8) |>  G0i
Γ(3//8) |>  G0i


Ω(1//2)


a0^2
a0^-2
(c0)^2
(c0)^-2
(c0^a0)^2
(c0^a0)^-2
(a0^c0)^2
(a0^c0)^-2
sqr = x->x^2
a*b^2|>p1o |> sqr
G0i(b) 





a|>p1o|>x->x^2
# |>x->(f0\x)^2

D0i^2 * G1i
G0i^4 * D1i


Ω(Q0)*Ω(Q1)
Ω(Q1)*Ω(Q2)
Γ(Q2)*Γ(Q1)
Γ(Q1)*Γ(Q2)
Γ(Q)
Γ(Q1)^2*Γ(Q0)
Ω(Q2)^Ω(Q1)
Ω(Q0)^2
Ω(Q) 

Ω(Q1)
f0\Ω(Q1)

w = Γ(3//8)|>G0i
w(a0,c0)
(a0, c0) |> s->f0\w(s...)
Q
rw = Ri(Q1) 
rw = Ri(Q0) 

Ω(Q)
f0*Λ(Q1)

b0^-f0
b0^-h0

(a*b)(a0,c0)


Q = @cf0 [2,1,3]
# Q = @cf0 [2,4,3,3]
begin
    Q0, Q1, Q2 = farey_neighbours(Q)
    Q3, l = descendQ(Q1, Q0)
    Q4, l = descendQ(Q0,Q3)
    Q5, l = descendQ(Q3,Q4)
    Γ(Q3)*Γ(Q0)^Q[end]==Γ(Q)
    Γ(Q0)^Q[end] * Γ(Q3)==Γ(Q)

    Γ(Q4) * Γ(Q3)

    ψ = [ϕ^l for (l, ϕ) in zip(Q, ITR.cycle([D1i, G0i]))] |> prod
    ψi = [ϕ^l for (l, ϕ) in zip(Q, ITR.cycle([D1, G0]))] |>reverse|>  prod
    ψi      |>
        ab .|>
        (y)->ContinuedFraction(last(y)[2]//last(y)[1])

    ψ       |>
        ab .|>
        (y)->ContinuedFraction(last(y)[2]//last(y)[1])

    wf = (ψ)
    wfi = (ψi)
end

Q
Q1
Q0
Q2
Q3
Q4
Q5

Γ(Q⊕Q0) |> wf 
Γ(Q) |> wf
Γ(Q0) |> wf
Γ(Q1) |> wf 




wf

Q⊕Q0,Q, Q0
Γ(Q)|>wf
Γ(Q0)|>wf
Γ(Q1)|>wf
Γ(Q2)|>wf
Γ(Q3)|>wf
Γ(Q4)|>wf
Γ(Q5)|>wf

a1 = wf(a)(a0,c0)
c1 = wf(b)(a0,c0)
wfi(a)(a1,c1)
wfi(b)(a1,c1)
wfi(a*b)(a1,c1)|>sqr



(a1^2*c1*a1*c1)^hh |> inv |> sqr
(a1^2*c1)^hh |> inv |> sqr
(a1*c1)^hh |> inv |> sqr
c1^hh |> inv |> sqr
(a1\c1)^hh |> inv |> sqr


Γ(Q⊕Q0)(a1,c1)
Γ(Q)(a1,c1)
u1 = Γ(Q0)(a1,c1)

Γ(Q0)(a1,c1)|>sqr
Γ(Q)(a1,c1)|>sqr
Γ(Q3)
Γ(Q1)
Γ(Q1)(a1,c1)
Γ(Q2)(a1,c1)
v1 = Γ(Q3)(a1,c1)

Ω(Q)
Ω(Q1)
Ω(Q0)*Ω(Q2) 
Ω(Q0)^2*Ω(Q2) 
Ω(Q0)^2


(Ω(1//3)/f0)^(a0^2)
(Ω(2//3)/f0)^(a0*c0)

(Ω(1//4)/f0)^(a0^2)
(Ω(3//4)/f0)^(a0*c0)

(Ω(1//6)/f0)^(a0^3)
(Ω(5//6)/f0)^(a0*c0^2)

(Ω(1//5)/f0)^(a0^3)
(Ω(3//5)/f0)^(a0*c0*a0)
(Ω(4//5)/f0)^(a0*c0*a0)

(Ω(1//7)/f0)^(a0^4)
(Ω(2//7)/f0)^(a0^3*c0)
(Ω(3//7)/f0)^(a0^2*c0^2)
(Ω(4//7)/f0)^((a0*c0)^2)
(Ω(5//7)/f0)^(a0*c0^3)
(Ω(6//7)/f0)^(a0*c0^3)

(Ω(1//8)/f0)^(a0^3*c0)
(Ω(3//8)/f0)^(a0^2*c0*c0)
(Ω(5//8)/f0)^((a0*c0)^2)
(Ω(7//8)/f0)^((a0*c0^3))

pp*f0*vv

a0^2

@cf 5//8

Ω(Q)/f0 
Ω(Q0)/f0 

Q0⊕Q1
Q0⊕Q2⊕Q0
Q0
Q2
Q


# Q0,Q,Q⊕Q0
u1^2*v1

