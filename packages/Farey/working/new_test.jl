include("../src/MonoidHoms.jl")
import .modhoms

using Farey
using Farey.Monoids
include("src/MonoidHoms.jl")
import .modhoms
const CF = ContinuedFraction

using .modhoms:MonoidHom
f1, g1 = FreeGen("f0 g1")
f, g = FreeGen("f g")

g2 = CyclicGen(2; order = 3)


f, g, h = FreeGen("f g h")
a, b, c = FreeGen("a b c")
x,y,z = CyclicGen("x y z"; order = 2)
X,Y,Z = FreeGen("X Y Z")
p1 = MonoidHom([
f=>a*z,
g=>z*a,
])

expand!(p1(g*f)|>inv)
(expand!(g*f|>inv))|p1


expand!(p1(inv(f)))
expand!(p1(inv(g)))
expand!(p1(f*(g)))
expand!(p1(inv(g*f)))

# p1 = MonoidHom([
#   f=>X*Y,
#   g=>Z*Y,
# ])

# p1 = MonoidHom([
#   X=>f*h,
#   g=>Z*Y,
# ])
Γ(p, q) = christoffel(p, q, a, b)
Γ(Q) = christoffel(Q, a, b)

Q=@cf0 [1,2,2]
Q=@cf0 [1,1,2]
Γ(Q)

f2 = CyclicGen("f"; order = 3)
g3 = CyclicGen("g"; order = 3)
g3^-1
p3 = MonoidHom([
f=>f2,
g=>f2*g3,
])
inv(a)
p1 = MonoidHom([
f=>a*z,
g=>-(z*a),
])
w=-(a*z)
w=(-a)
dump(w)
www = farey_word(3//5,f,g)


Q = @cf0 [2, 2, 2]
(s_seq(Q))|>sum
r, s = numerator(Q), denominator(Q)
www = farey_word(Q,f,g)[1:s]
# aaa, bbb = www|>conj_prefix
www = farey_word(Q,f,g)[s+1:end]|p1
www|>conj_prefix
ww = farey_word(3//5,f2,g3)
SFX.SuffixTree(ww)
ww|p3
# p1 = MonoidHom([
#   f=>a*z,
#   g=>z*a,
# ])

dhnx = MonoidHom([
Y => inv(Z)*inv(Y)*inv(X),
])
dhny = MonoidHom([
X => Z/X/Y,
])

dhnz = MonoidHom([
Z=>inv(Y)/Z*X
])

p2 = MonoidHom([
x=> x,
y=> c*z*a,
z=> z,
])
[1,1]
-(a*z)
ww = (farey_word(1//2,x, g))
ww|(dhny)
f2 * inv(g3)
-f

inv(f2)
(-farey_word(12//29,f1,g1))
q1 = @cf0 [2, 3, 2, 4,3]
q2 = @cf0 [2, 2, 2]


gens = [FreeGen("W$i") for i in 1:10]
# @cf 4//11


qq = @cf0 [2, 2]
# qq = @cf0 [2, 1,3]
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


MonoidWord([1,2])
element("f g 2")
q1 = CF(0,[2,2,1,3,5])
qns = qnbs(q1)
collect(enumerate(qns))
qpair.(filter(x->length(x)>=3, qns))

idxr = sum(q1[1:end-1]); idxt = q1[end] % 2 == 1 ? idxr+1 : sum(q1[1:end-2])
function farey_expr(qq::Number)
qq = CF(qq)

idxr = 2+sum(qq[3:end-1]); idxt = qq[end] % 2 == 1 ? idxr+1 : 1+sum(qq[3:end-2])
k = qq[end]÷2

w1, w2, f = FreeGen("W$idxr W$idxt f")

# length(qq) % 2 == 0 ? w2 * inv(w1^f)^k*w1^k : w2 * w1^k*inv(w1^-f)^k
length(qq) % 2 == 1 ? w1^k * w2 * inv(w1^f)^k : inv(w1^-f)^k *w2 * w1^k
end

# qs = Q_path(q1)
# qts = [qpair(qq)[2] for qq in filter(x->0<x<1, qs)]
# qts1 = [qpair(qq)[2] for qq in filter(x->0<x<1, qts)]
function expandhom(q1)
qns = qnbs(q1)

f,g,  w, ws... = FreeGen("f g W "*join(["W$i " for i in eachindex(qns)]))
t0, t1 = FreeGen("t$(q1[2]) t$(q1[2]+1)")
u, v = q1[2] % 2 == 0 ? (-f, -g) : (-g, -f)
MonoidHom([[
    ws[1]=>f*t0*u/t0,
    ws[2]=>f*t1*v/t1,
]..., w=>farey_word(0, f, g)]),
MonoidHom([
    ws[i] => farey_expr(q) for (i, q) in Iterators.filter(x->length(x[2])>2, enumerate(qns))
])
end
q1 = CF(0,[2,5,5])
fn, fw = expandhom(q1)
fn, fw^3


Wm = i-> (element("W$(i)"))
ww17 = (fw^4)(element("W17"))
ww(15)
K = 5
fw^2
th1 = -gcp(-ww(K), -ww(K-1))
th2 = gcp(ww(K), ww(K-1))
length(ww(10))
ww(10)
#  gp = MonoidHom([
#     f1=>element("a b -1"),
#     g1=>element("b c -1"),
#   ])
fw^2
fn

f1\((ww(10)|(fw^5))) |> conj_prefix

farey_word(2//5, f1, g1)

comm(ww(2)^2, f1/ww(1))
fw * fn
fw^2

ww(4)

w6 = ww(6)|fw|fw|fn

cm = element("f t4 f -1 t4 -1")^2 * element("f")

conj_prefix(f\w6)

collect(enumerate(qnbs(q1)))
