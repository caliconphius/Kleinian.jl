function nregtree(k, depth)
    gr = Gjl.SimpleDiGraphFromIterator(
        [Gjl.Edge(((i)÷(k))+2=>i+3) for i in 0:((k)^(depth+1)-1)÷(k-1)-2]
        )
    # Gjl.add_vertex!(gr)
    # Gjl.add_edge!(gr, Gjl.Edge(1=>2))
    gr
end
gra = nregtree(2,4)
Gjl.outneighbors(gra, 10)
nbs = i->Gjl.outneighbors(gra, i).-2*(i-1)
[nbs(i) for i in 1:20]
struct VertexData{T}
    left::T
    center::T
    right::T
    op
end

Lshift(v::VertexData{T}) where T = VertexData{T}(v.left, v.op(v.left, v.center), v.center, v.op)
Rshift(v::VertexData{T}) where T = VertexData{T}(v.center, v.op(v.center, v.right), v.right, v.op)

@present SchWeightedGraph <: SchGraph begin
  T::AttrType
  weight::Attr(V,T)
end

@acset_type WeightedGraph(SchWeightedGraph, index=[:src,:tgt, :weight])<:AbstractGraph


otimes = Farey.:⊕
initV = VertexData{ContinuedFraction}(CF(0//1), CF(1//1), CF(1//0), otimes)|>Lshift
initVfg = VertexData{Farey.FreeGroupElement}(f, f*g, g, *)|>Lshift

function farey_tree(depth, initV, Lshift, Rshift)
    gra = nregtree(2,depth)
    
    
    gr = @acset WeightedGraph{typeof(initV)} begin
        V = 0
        E = 0
        src = []
        tgt = []
        weight = []
    end

    add_vertices!(gr, Gjl.nv(gra))
    set_subpart!(gr, 2, :weight, initV)

    for e in Gjl.edges(gra)
        e.src==1 && continue
        add_edge!(gr, e.src, e.dst)
        shift = e.dst%2 ==0 ? Rshift : Lshift
        newV = shift(gr[:weight][e.src])
        set_subpart!(gr, e.dst, :weight, newV)
    end

    gr


end



dpt = 5
grfg = farey_tree(dpt, initVfg, Lshift, Rshift)
gr = farey_tree(dpt, initV, Lshift, Rshift)
nbs = outneighbors(gr,10)|>collect
gr[:weight] .|> x->isnothing(x) ? nothing : x.center|>Fw
θ = g1\d1^3|>inv|>hom
θi = g1\d1^3|>hom
Φi
gr[:weight][2:end] .|> x->(isnothing(x) ? nothing : fp(x.left)(f))
gr[:weight][2:end] .|> x->Ω.(farey_neighbours(x.center))
gr[:weight][2:end] .|> x->(Ω(x.left),Ω(x.right))
gr[:weight][2:end] .|> x->x.center=>fp(x.center)(Γ(x.center))

fpi(2//7)

tw = k-> g0\d1*g0^k
sw = k-> d1*g0^k

fp(Q)(Γ(Q))


si(k) = sw(k)|>inv|>hom
s0(k) = sw(k)|>hom
ti(k) = tw(k)|>inv|>hom
t0(k) = tw(k)|>hom

a0^2*c0
a0*c0
chr(k) = d1*g0^k|>hom
ti(4)(Γ(CF(0,[4,1,2])))
s0(2)
ϕ = t0(1)*s0(1)*ti(2)
t0(3)*si(1)*ti(1)
(ϕ)(Ψ(1//3))
(ϕ)(Ψa(1//3))
s0(4)(Ψ(1//4))
si(4)(Γ(CF(0,[4,1,2])))
(Γ(CF(1,[2])))

mapp(Q) = fpi(Q)



gr[:weight][2:end] .|> x->Rational.([x.left, x.right])=>Fw(x.center)|>l->(mapp(x.right))(l)

gr[:weight] .|> x->isnothing(x) ? nothing : x.center


Λ(Q) = f0\Ω(Q)
χ(Q) = Ω(Q)/f0
Q = @cf0 [2,1,1]
Qa, Qb, Qc = farey_neighbours(Q)
(Λ(Qa)*χ(Qb))^2
Λ(Qc)*Λ(Qa)
function aab_dec(rhat)
t0Q, tQ, t1Q = farey_neighbours(rhat)

sQ, cQ =  otimes(tQ,t1Q)==rhat ?    (t0Q, t1Q) :
                            (t1Q, t0Q)
sQ,cQ    
end
# Q

Qs = sw_descent(Q)
Qb⊕Qa⊕Qa
Qa<Qb

Ω(Qb)*Wr(Qa)
Wr(Qa)*Ω(Qb)/Ω(Q)
Λ(Q)|(Λ(Qs[1])*Ω(Qs[3])*Ω(Qs[5]))

Λ(Qs[end])|(Λ(Qs[2])*Ω(Qs[4]))

χ(Q)
(Λ(Qa)*f0*Λ(Qa))
(Λ(Qa)*h0*Λ(Qa))
Ω(Qc)
descendQ(Q, Q1; l=0) = (foldl(Farey.:⊖, Farey.ITR.cycle([Q1], Q[end]+l), init=Q), Q[end])

Q = @cf0 [2,2,4,2]
Qb, Qa = aab_dec(Q)
Qb,l = descendQ(Q, Qa)

Q<Qa ? Γ(Q)/Γ(Qb)^l|>fp(Qa)*G0 : Γ(Qb)^l\Γ(Q)|>fp(Qa)*G0
Q<Qa ? Γ(Q)/Γ(Qa)^l : Γ(Qa)^l\Γ(Q)

import Iterators as ITR

function rileyw(Q; f=f, g=g)
    p, q = denominator(Q), numerator(Q)
    ϵi = i-> (-1)^Int(floor(i*q//p))
    prod(
        [u^ϵi(i) for (u, i) in zip(ITR.cycle([f,g]), 1:(p-1))]
    ,init=one(f))
end


rileyw(1//2; f=f, g=g)*f
f*rileyw(1//3; f=g, g=f)
rileyw(2//3; f=f, g=g)*f
QQ = @cf 5//7
# QQ = @cf0 [2,2,1,3]
u = denominator(QQ)%2==0    ? f :
    numerator(QQ)%2==0      ? G :
    g
long_n(QQ)
V = rev_conj(QQ)|>signs_seq

f*rileyw(QQ;f=g,g=f)
(rileyw(QQ;f=g,g=f)*u)^-1

f*rileyw(QQ;f=g,g=f)/(rileyw(QQ;f=g,g=f)*u)


function long_n(Q)
    p, q = denominator(Q), numerator(Q)
    l = gcdx(-q,2*p)[2]
    l = (l%(2*p))
    l < 0 ? l+2*p : l

end

function rev_conj(Q)
    p, q = denominator(Q), numerator(Q)
    l = long_n(Q)
    ϵi = i-> (-1)^Int(floor(i*q//p))
    prod(
        [u^ϵi(i) for (u, i) in zip(ITR.cycle([f,g]), 0:(l))]
    )
end

long_n(3//5)
invQ(Q) = CF(0, [collect(Q[2:end-2])..., Q[end-1]])
#  Q⊖(sw_descent(Q)[end]⊖Q)#⊖(sw_descent(Q)[end]⊖Q)

long_n(3//7)
rileyw(1//2)


rw(Q) = rileyw(Q)*f/rileyw(Q)/g
rw(1//3)
begin


Q = @cf0 [2,3,3,4]
Qs = sw_descent(Q)
Qn = Qs[end]⊖Q
long_n(Q)
Q1 = Q⊖Qn⊖Qn⊖Qn
rev_conj(Q)*Fw(Qn),Q,Qn
end
Qn
Q
QL = @cf0 [2, 3,3]
long_n(Q)
signs_seq(Ψ(QL)/g)
signs_seq(rev_conj(Q))
Ψ(QL)\rev_conj(Q)
begin
Q = @cf0 [3,6,3,3]
long_n(Q)
Qs = sw_descent(Q)
Qn = Qs[end]⊖Q
Q1 = (Qs[end])
Ψ(Q1)\rev_conj(Q),Q, Q1
end



begin
Q = @cf0 [3,2,3,3]
long_n(Q)
Qs = sw_descent(Q)
Qn = Qs[end]⊖Q
Q1 = (Qs[end])
rev_conj(Q)*(g*rileyw(Q1)),Q, Q1
end


begin
Q = @cf0 [3,2,4,3]
long_n(Q)
Qs = sw_descent(Q)
Qn = Qs[end]⊖Q
Q1=Qs[end]
# Q1 = @cf [0,3,2,2,2]
rev_conj(Q)*(g*rileyw(Q1)),Q,Q1
end

# begin

[a⊖b for (a,b) in zip(Q_abv[1:end-1], Q_abv[2:end])]
long_n.(Q_abv)

Q = @cf0 [3,5,4,3]
long_n(Q)
Qa,Qb,Qc = farey_neighbours(Q)
Qs = sw_descent(Q)
Qn = Qs[end]⊖Q
Q1 = (Q⊕Qa⊕Qa)
Q*long_n(Q) |> x->x - floor(x) |> CF
Q⊕Qa⊖Q1
@cf0 [3,5]
rev_conj(Qa)|>signs_seq
Qa
Qd = @cf0 [3,5]
Qc = Qa⊖Qd
Ψ(Qc)
(Qa⊕Qc⊕Qc)|>long_n
addK(Qa,Qb,k) = foldl(⊕, ITR.cycle([Qb],k), init=Qa)
Qk = addK(Qa,Qc,10)
Qa
descendQ(Qk,Qc)
descendQ(Qc,Q2)
addK(Q2, Qc, 11)
Qk⊕Q1
Q1⊕Q2
Qa, Qb, Qc = farey_neighbours(CF(0,[2,3,7,3]))
begin
Q = addK(Qb, Qc, 4)
Q_abv = [CF(0,Q[2:i]) for i in 2:length(Q)]
Q_abv.|>x->(rev_conj(x)|>l->[l,Fw(x)].|>signs_seq)
[a⊕b for (a,b) in zip(Q_abv[1:end-1], Q_abv[2:end])].|>Ψ.|>signs_seq
end
long_n.(Q_abv)
Q_abv
rw
Q_abv.|>q->Rational(q)=>rileyw.(q)|>l->signs_seq(l)
Q_abv.|>q->Rational(q)=>Fw(q)|>l->signs_seq(l)=>l[1]
revs = Q_abv.|>q->Rational(q)=>rev_conj.(q)|>l->signs_seq(l)
Q_abv.|>q->Rational(q)=>rev_conj.(q)\Fw.(q)#|>l->signs_seq(l)

q=CF(0,[2,2,2,2])
q1 = CF(0,q[2:end-1])
q2, q3,q4,q5 = ((q1⊖q),q1⊕q,q1⊕q⊕q1,q1⊖q⊖q1)



A1, A2, A3,A4 = let q=CF(0,[2,2,2,2])
q1 = CF(0,q[2:end-1])
q2, q3,q4,q5 = ((q1⊖q),q1⊕q,q1⊕q⊕q1,q1⊖q⊖q1)
[println(x) for x in (q, q1, q2,q3,q4,q5)]
A1, A2,A3,A4 = [rev_conj(q),(rev_conj(q)\Fw.(q)), rev_conj(q4), rev_conj(q4)\Fw(q4)]

println(A1|>signs_seq|>l->(length(l),sum(l)))
println(A2|>signs_seq|>l->(length(l),sum(l)))
println(Fw(q)|>signs_seq|>l->(length(l),sum(l)))
println(Ψ(q4)|>signs_seq|>l->(length(l),sum(l)))
A1,A2, A3,A4
end;
Fw(5//13)|>signs_seq|>l->(length(l),sum(l))
Ψ(q4)|>signs_seq|>l->(length(l),sum(l))
rev_conj(q)\Ψ(q4)|>signs_seq|>l->(length(l),sum(l))
Ψ(q4)|>signs_seq|>l->(length(l),sum(l))
rev_conj(q)|>signs_seq|>l->(length(l),sum(l))
Fw(11//27)|>signs_seq|>l->(length(l),sum(l))
(rev_conj(q)|>signs_seq)-((Ψ(q4)/g)|>signs_seq)
@cf 11//27
q4
(((A1*(A3))/((A3)*A1)|>signs_seq).!=(A3*A4|>signs_seq))
[A1,A3, A4,A3*A4,A1*A3/(A3*A1) ].|>x->(signs_seq(x)|>l->(length(l),sum(l)))

rev_conj(15//22)
revs[4][2]|>length
revs[4][2]|>signs_seq

rev_conj(2//5)
rw(2//5)
Fw(1//2)

rev_conj(Qk)|>signs_seq
rev_conj(Q1)|>signs_seq
Ψ(Q2)|>signs_seq
Q1 = Qc
Q2,l = descendQ(Q1,Qa)
Ψ(Q1)|>signs_seq
rev_conj(Qk)|>signs_seq
rev_conj(Qk)\Fw(Qk)|>signs_seq
Qk = addK(Qc,Qa,3)

rev_conj(Qa⊕Qc⊕Qc)|>signs_seq
Fw(Qa⊕Qc⊕Qc)|>signs_seq
Qa, Qc,Qa⊕Qc
# end




(g*rileyw(Q1))

Qn = Qs[end]⊖Q
Qp = Qn⊖Qs[end-2]
rev_conj(Q)\(Fw(Q)), Q

reverse(w::Farey.AbstractGroupElement) = prod([inv(w.parent(x)) for x in w.word])|>inv
reverse(rw(Qn))
(f*reverse(rileyw(Qn; f=g, g=f)))


long_n(Q)
s_seq(Q)
rev_conj(Q)|>x->(length(x))/2
rev_conj(Q)

rev_conj(Q)|>conj_prefix
(rw(Q)*rev_conj(Q))

f*rileyw(3//5; f=g, g=f)
rev_conj(Q)
f|rileyw(Q)|rev_conj(Q)
g*rileyw(3//5)

rileyw(Q)|>inv
rileyw(Q)
rev_conj(Q)\Fw(Q)
function signs_seq(w)
    letters = [w.parent(x) for x in w.word]
    genrs = Farey.GPC.gens(w.parent)

    out = foldl(letters, init=(0, [], true)) do x, y
        if ((y∈genrs)==x[3])
            return (x[1]+1, x[2], x[3])
        else
            return (1, [x[2]..., x[1]], !x[3])
        end
    end

    [out[2]..., out[1]]
end
Q,Q1

signs_seq(rileyw(Q))
signs_seq(Ψ(Q))
signs_seq(rev_conj(Q))
signs_seq(rileyw(Q1))
f\rev_conj(Q)*g*rileyw(Q1)

Q = @cf0 [ 2, 3, 5]
begin
 Q = @cf0 [ 4, 2,3,4]
long_n(Q)
Qa,Qb,Qc = farey_neighbours(Q)
Q⊖Qb
Q*long_n(Q) |> x->x - floor(x) |> CF
rev_conj(Q)*(Fw(Qa))
end
Q,Qa

begin
Qa = @cf0 [ 2, 3, 4, 2]
Q1
@cf0 [2, 3, 5]
end
rev_conj(Q)
(Fw(Q1)/f)|>inv



