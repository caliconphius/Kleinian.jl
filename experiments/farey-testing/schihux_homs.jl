include("working_setup.jl")
hmphw = o * e / g0 * d1 * e * o
hmphi = hmphw|>inv|>hom
hmph = hmphw|>hom

pw = g0 * hmphw * g1
pf = hom(pw)
pfi = hom(inv(pw))
(D0*pf*D1i)
pfi(a*b)#*a*B^2)
pf(a*b)#*a*B^2)
(G1 * hmph * G0)(b*a)
(G1 * hmph * G0)(a*B)
((G1 * hmph * G0))(b*a*b*a*b)
((G0i * hmphi * G1i))(a*b*a*b*a)

rileyw(3//5) * f
g * rileyw(3//5)

Rw(Q) = rileyw(Q) * (
  denominator(Q) % 2 == 0 ? G :
  numerator(Q)   % 2 == 0 ? F :
                            f
) / (g*rileyw(Q))



fww = (g0*d0)
fwwi = (inv(g0*d0))
hfww = (d0*g0)
hfwwi = (inv(d0*g0))

gww = (g1*d1)
gwwi = (inv(d1*g1))
hgww = (d1*g1)
hgwwi = (inv(d1*g1))


fw = fww|>hom
fwi = fwwi|>hom
hfw = hfww|>hom
hfwi = hfwwi|>hom
gw = gww|>hom
gwi = gwwi|>hom
hgw = hgww|>hom
hgwi = hgwwi|>hom

  phif(Q, A, B) = [ϕ^l for (l, ϕ) in zip(CF(Q).L, ITR.cycle([A, B]))]|> prod 
  phifr(Q, A, B) = [ϕ^l for (l, ϕ) in zip(CF(Q).L|>reverse, ITR.cycle([A, B]))]|> prod 


let X = fww; Y=gww; Q = @cf0 [2,1,2]
  X, Xi, Xp, Xv =([X, inv(X), reverse(X), inv(reverse(X))])
  Y, Yi, Yp, Yv =([Y, inv(Y), reverse(Y), inv(reverse(Y))])

  c1, c2, c3, c4 = f*g, g*f, f*G, G*f
  ts = [
    Fw(Q),signs_dec(Fw(Q))...
    ]
  ts = [f, g
  # , ((Fw(Q)))
  , (Fw(3//5))
  ]


  o2 = 
  (g1)*
  o
  ws = [
    
    
    (g0/d0)/(d0/g0)|>
    inv|>
    x->x
    # ((g0)|d0)
    # phif(3//8, (g1), (d0))|>inv
    # phifr(3//8, (g0),(d0))|>inv
    # phifr(1//2, inv(g0),inv(d0))
    # phif(3//8, inv(g1), (d0))|>inv
    # phifr(3//5, inv(g0), inv(d1))|>inv
    # phifr(Q, inv(d0), (g1))|>inv
    
    # (phif(Q, inv(g0),(d0)))|>inv
    # phif(3//5, inv(d0), g1)|>inv
    # e*phif(3//5, (g1), inv(d1))|>inv
    # phif(3//5, (g0), (d0))|>inv
  ]
println("----")
  [
  """
  $t = \n$(["($s)\n" for s in signs_dec(t)]|>join)
  ↓  $w\n
  $(hom(w)(t)) = \n$(["($(hom(w)(s)))\n" for s in signs_dec((t))]|>join)
  \n$(hom(w)(t)) = \n$(["($((s)))\n" for s in signs_dec(hom(w)(t))]|>join)
  """
  for (t,w) in ITR.product(ts, ws)].|>

  # begin
  #   hs = hom.(ws)

  #   prt =
  #   # " $("⎸"^(((length(ts)-1)*(length(ws)))))$("╲"^(length(ts)))\n"*
  #   join([
  #     "$("⎸"^(((length(ts))*(length(ws)-kh))))$(" "^(k-1))$("╲"^(length(ts)-k))"*
  #     join([
  #       # "$("⎸"^((length(ts)-1)*(length(ws)-kh+1)))  $("╲"^(length(ts)-k+1))",
  #       "╲",
  #       "╲",
  #       " $t",
  #       " ↧ $n",
  #       " $(x(t))",
  #       " = $(join(signs_seq(x(t)).|>l->"S$(l)","*"))"
  #     ],"\n$("⎸"^((length(ts))*(length(ws)-kh)))$(" "^(k))$("⎸"^(length(ts)-k))")
    
  #   for (kh, (n, x)) in enumerate(zip(ws, hs)) for (k, t) in enumerate(ts)], "\n")
  #   """
  #   u ∈ {$(join([repr(t) for t in ts], ", "))}
  #   ϕ ∈ {$(join([repr(t) for t in ws], ", "))}
  #   --------
  #   """*prt*
  #   """\n  """
  # end.|>
  println
# "
#     $(join(["$(length(x)): $x ↤ $(hom(inv(w))(x))\n" for x in ((t)|>signs_dec)]))⤓\n$((t))\n↓ $w\n$(hom(w)(t))\n⇣\n$(join(["$(length(x)): $x \t↤ $(hom(inv(w))(x))\n" for x in (hom(w)(t)|>signs_dec)]))\n"
  ws[1]|>inv|>hom,
  ws[1]|>hom
end
(f\Fw(3//5))

D0*G0i*D0i
D1i*G0

signs_seq(Γ(3//5))

ff(-1)
ψf(-5//2)

# ((hgw*fw)^((a*b*a*B*A*B)))(a)
  # ((hgw*fw)^(a*b/a/b))(b)

  # ((hgw*fw)^((a*b*a)))(A*b)
  # ((hgw*fw)^((a*b*a)))((B*a)^3)
  # ((hgw*fw)^((a*b*a)))(a)

  # (hgw*fw)(((B*a)))
  # (hgw*fw)(B*a)
  # (fwi*hgwi)
  # a|>(hgw*fw)^(a*b*a)
  # hfw*hgwi
  # ((hgw*hfw))(a*b*A*B)
  # ((gwi))(a*b*A*B)
  # gw
  # fwi(a^2*b)
  # fwi(a^3*b)
  # fwi(a*b*a*B*A*B)
  # (gwi*fwi)(a*b*a*b*A*B*A*B)
  # (gwi*fwi)(a*b*a*b*a*B*A*B*A*B)

  # (b*a)^2 |> hfwi
  # (a*b)^2 |> hgwi

  # ((fwi)(Fw(3//5)))^(F*g)|>signs_seq
  # Ψ(3//5)
  # ((hfw)(Fw(2//7)))

  # (Rw(3//8))
  # Fw(3//8)^rev_conj(3//8)

  # Rw(2//7)
  # Fw(2//7)|rev_conj(2//7)


  # Rw(3//11)
  # (Fw(3//11)^rev_conj(3//11))


  # Ψ(3//8)
  # Ψa(3//8)
  # Rw(3//8)

  # Fw(Q)|>signs_seq

  # Ψ(3//5)
  # # import .SuffixTrees:SuffixTree
  # # SuffixTree(ϕi(Ψ(Q)))
  # # SuffixTree((Ψ(Q)))


  # let Q = CF(0, [1,1,2])
  # Ψ(Q)^rev_conj(Q)
  # # ψf(Q)|>hom
  # w = prod([(fh)^l for (l, fh) in zip(Q.L, ITR.cycle([d0*inv(g1), inv(d0)]))])
  # HH = w|>hom
  # HHi = w|>inv|>hom
  # HHi(Ψ(3//5)),HHi(Ψa(3//5))
  # # (HH(Ψ(3//5))HH(Ψa(3//5)))^(g*F*G)
  # end

  # G1i * D0

  # O * D0i * O
  # D1i * G0
  # hmph = O * E * G0i * D1 * E * O
  # hmphw = o * e * inv(g0)*d1 * e * o
  # hmphi = hmphw|>inv|>hom
  # hmph = hmphw|>hom

  # (D0i*hmphi)(Ψa(3//5))
  # Ψa(3//5) |> (D0 * hmph)
  # (D0i*hmphi)(Ψa(3//5))
  # hmph(Ψa(1//4))
  # ((G1i*hmph)(Fw(3//5)))#* Fw(1//2)
  # ((hmph^4)(Fw(5//8)))^(f*g)

  # (D0 * hmph)
  # (G1 * hmph)

  # (D1i^2*hmphi*G1i*(hmphi*G1i)^2)(Fw(3//5))

  # Fw(1//2)\
  # (((hmph^3))(Ψa(3//5)))*
  # Fw(1//2)

  # (hmph^3)(a*b*a*B*A)
  # (hmph^3)(a*b*a*b*a)
  # (hmph^3)(a*b*a)


  # homhom = (g0/(d1^2))*g0/(d0^2)
  # fpw = homhom|>inv|>hom


  # fpw(rileyw(3//5;f=g,g=f))
  # fpw(Ψa(3//5))
  # fpw(Ψa(3//5))
  # hmph
  # ((G1*hmph)^3)
  # ((G1*hmph))
  # ((G1*hmph))(a*b*a*b*a*b)
  # ((G1*hmph))(b*a)
  # ((hmphi*G1i))
  # ((a*b*a)/(b*a*b))

  # (hmph*G1)^2

  # ((G1i*hmphi)^2)(
  #   (f*g*f*g)*
  #   (F*G*F)*
  #   (g*f*g)*
  #   (F*G*F*G)*
  #   (g*f*g)*
  #   (F*G*F)
  # )
  # (
  # (hmph*G1i)^2*
  # (hmph*D1i)
  # )(
  #   (f*g*f*g)*
  #   (F*G*F)*
  #   (g*f*g)*
  # (F*G*F*G)  
  # )

  # (
  # # (hmph*G1i)^2*
  # (hmph*G0i)*
  # (hmph*D1i)*
  # (hmph*G0i^2)
  # )
  # (
  # Fw(3//5)
  # )


  # ((hmphi*G1i)^2)(b*a*b*a)
  # ((hmphi*G1^2))(Ψ(3//5))
  # ((hmphi*G1i)^3)(Ψa(3//5))
  # (hmphi*D0i)(Ψa(3//5))
  # (hmph)(a*b*a*B*A*B)

  # fw*gw
  # (hw*gw)(a)

  # ((a*b) * (A*B)^2 )|> (fw*fw1i*fwi*fw1)
  # (fw*fw1i*fwi*fw1)(a*b*a*b*A*B*A*B)
  # (fw*fw1i*fwi*fw1)(B*A*B)
  # (fw*fw1i*fwi*fw1)(A*(B*A)^2)
  # (fw*fw1i*fwi*fw1)(b^A)
  # (fw1i*fw*fw1*fwi)(a^B)
  # (fw1i*fw*fw1*fwi)(a^B)
  # (fw1i*fw*fw1*fwi)(a)

(D0*O)
((O*E)^3*D1*O*E)
((O*E)*D0*O*E)
((O*E)*D1*O*E)*D0

(O*E)^2
G0*D0*G1i
# ()(f*g*F*G*f)
dgh = (d1*g1i)|>hom
dgh0 = (d0*g0i)|>hom
dgh0i = inv(d0*g0i)|>hom
dghi = inv(d1*g1i)|>hom
@cf0 [2,2,3]
out = let k = 1; qq = @cf0 [1,1,2,2]
println(CF(qq))
println(farey_neighbours(CF(qq)))
# println(farey_neighbours(3//7))
qa, qb, qc = farey_neighbours(qq)
(qa⊕qc⊕qa)|>println
ts = [

  (Fw(qq))\
  inv(Fw(qa)^F)*
  (Fw(qc))*
  (Fw(qa)),
  # (Fw(qq))*F,
  # (rev_conj(qq)),
  # (Fw(qa)),
  # (Fw(qc))*F*
  # inv(Fw(qa))*f,
  # inv(Fw(qa)),
  # (F*Ψ(qq)),

  # inv(Fw(qa)),
  # (Fw(qc))^f,
  # (Fw(qa))^f,
  # (Fw(qq)),
  # (Fw(qc)),
  # (Fw(qc)),

  # Fw(qa)*Fw(qc)*inv(Fw(qa)^f),

  # ((E)(rileyw(qq))),
  # Ψa(qq)*g/f/g
  # f,
  # g,
  # f*g,
  # (Γ(qq)),#Fw(3//7),
  # (Fw(qq)),#Fw(3//7),
  # (Fw(1//2))*((Fw(1//2))*g/(Fw(1//2))),
  # ((Fw(3//8))^rev_conj(3//8))|Fw(1//3),
  # (rev_conj(qq)),
  # ((Fw(qa))*(Fw(qc)))*F*inv(Fw(qa))*f,
  # (Ψa(qq)*f),
  # (Fw(qq)),
  # f*inv(Fw(qa)),
  # F*(Fw(qc)),
  # (Fw(qa)),
  # one(Fw(1//1)),
]
hf1, hf1i = ((dgh^k),(dghi^k))
hf0, hf0i = ((dgh0^k),(dgh0i^k))
println.(["""\n
---------
$(t)\n=\n
$(signs_seq(t))
$(signs_seq(t)|>ITR.flatten|>collect|>x->(sum(x), length(x)))
---------
""" for t in ts])
println()
println("qa = $(qa),qc = $(qc), qa<qc ? $(qa<qc), $(qa⊕qc⊕qa)")
# println.([hf1(t) for t in ts])
# println.([signs_seq(hf1i(t)) for t in ts])
[hf1(t) for t in ts]
# [signs_seq(hf1(t)) for t in ts]
end;
Q = @cf 7//16
qa,qb,qc = farey_neighbours(Q)
qa<qb
Γ(qa)*Γ(qb)*Γ(qa)*Γ(qb)

if qa <= qc
frw = Fw(qa)*Fw(qc)*F*inv(Fw(qa))*f
else
frw = f*inv(Fw(qa))*F*Fw(qc)*(Fw(qa))
frw1 = f*inv(Fw(qa))*F*Fw(qc)

end

frw
Fw(qa)