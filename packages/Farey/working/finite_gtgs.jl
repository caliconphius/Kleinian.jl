module GTGS

begin
using Farey
  
  using Farey.Monoids
  include("../src/MonoidHoms.jl")
  import .modhoms
  const CF = ContinuedFraction

  function qpair(qq)
    qq = CF(qq) 
    qr = CF(0,qq[2:end-1])
    qt = qq[end] % 2 == 0 ? CF(0,qq[2:end-2]) : CF(0,[qq[2:end-1]..., 1])
    k = qq[end]÷2

    qr, qt, k  
  end

  function qnbs(Q::T) where T<:Number
    Q = CF(Q)
    qs = CF.([1//(Q[2])])
    for (k,l) in enumerate(Q.L[2:end])
      append!(qs, [positive_form(CF(0, [Q[2:k+1]..., x])) for x in 1:l])
    end
    qs
  end

  function farey_expr(qq::Number)
    qq = CF(qq)

    idxr = 2+sum(qq[3:end-1]); idxt = qq[end] % 2 == 1 ? idxr+1 : 1+sum(qq[3:end-2]) 
    k = qq[end]÷2

    w1, w2, f = FreeGen("W$idxr W$idxt f")

    # length(qq) % 2 == 0 ? w2 * inv(w1^f)^k*w1^k : w2 * w1^k*inv(w1^-f)^k
    length(qq) % 2 == 1 ? w1^k * w2 * inv(w1^f)^k : inv(w1^-f)^k *w2 * w1^k
  end
end
begin
  A = i->CyclicGen("x"; order = i)
  B = i->CyclicGen("y"; order = i)
  a, b = FreeGen("a b")

  hom(i,j) = modhoms.MonoidHom([
    a => A(i),
    b => B(j)
  ])
  hom(2,5)
  fin_words = [
    1=>(element("a b a b a -b a -b") => (2,3,2)),
    2=>(element("a b -a b a -b") => (2,3,3)),
    3=>(element("a b a -b") => (3,3,2)),
    4=>(element("a b -a -b") => (3,3,2)),
    5=>(element("a b a b 2") => (2,5,2)),
    6=>(element("a b a b a b 4") => (2,5,2)),
    7=>(element("a b a b 2 a b 4") => (2,4,2)),
    8=>(element("a b a b a -b") => (2,4,2)),
    9=>(element("a b a -b") => (2,3,2)),
    10=>(element("a b a b a -b") => (2,3,2)),
    11=>(element("a b a b a b a -b") => (2,3,2)),
    12=>(element("a b a b a -b a b a -b") => (2,3,2)),
    13=>(element("a b a b a b a -b a -b a b a -b a -b") => (2,3,2))
  ]

 
ws =  last.(fin_words) .|> x->(first(x)|hom(last(x)[1], last(x)[2]), last(x)[3]
)
wps = last.(fin_words) .|> x->(first(x)^last(x)[3]
)|hom(last(x)[1], last(x)[2])
end
k = 2
k = 1
# for k in eachindex(fin_words)
p,q,n = last.(last.(fin_words[k]))
x, y = hom(p,q).([a,b])
X, Y = FreeGen("X$p Y$q")
uu = (wps[k])^-y
gcp(uu, -uu)
first.(ws)[k]^-y
ffeq = (w1, w2) -> (map(x->(x[1]=>x[2]), zip(gens(w1), eachgen(w2))))


one(uu)
(y^x)^y
inv(x*y)
xyho = 
modhoms.MonoidHom([
    x => x^y,
    y => y^-(x*y)
])#|
modhoms.MonoidHom([
    x => X,
    y => Y
])

g1i = modhoms.MonoidHom([
  x=>a\y,
  y=>y^(a\y),
])

g1 = modhoms.MonoidHom([
  a => x*y,
  y => x*y*x,
])

  xyho(ws[k][1])
  (ws[k][1])
  println("$((ws[k][1])^n|xyho) -> $((wps[k]|xyho|g1i))  ($p, $q, $n)")
# end
((ws[k][1]|xyho))

Γ = Q->christoffel(Q, a, b)
reverse(Γ(3//5))
# g1i|g1

include("../src/Palindromes.jl")
import .Palindromes as pal

Q = @cf0 [1,2,2]
pal.pal_suff(Γ(Q)/b)
Ty = typeof(Γ(Q)*x)
Ty([1])
pal.pal_conj(Γ(Q))
# pal.Pal()
function Pal(w::T; ids=(T(1), T(2))) where T<: AbstractMonoidWord
      
    idg, idd = ids
    G0 = 
    out = one(w)
    for x in eachgen(w)
        out = x==idg ? G0(out)*Gr(g) : D0(out)*Gr(g)
    end
    out
end
end