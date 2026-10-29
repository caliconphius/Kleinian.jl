# includet("../src/monoid.jl")
using Farey


import Farey.Monoids as mond


using Farey.Monoids


FreeGen("-a -1")

ja = element("a1 -2 b -a")

gens(ja)|>collect
s = "a1 b -a"
gens(ja)|>collect
ms = eachmatch(r"(?<gen>-?(?:[A-Z]|[a-z])\d*)(?:[\s\^](?<power>\d+|-\d+))?[\s\*]*", s)
collect(ms)

matches = map(ms) do m
    exp =  m[:power]
    sgn = m[:minus]
    pow = (isnothing(exp) || isempty(exp)) ? 1 : parse(Int, exp)

    inv = xor(!(isnothing(sgn)||isempty(sgn)), pow<0)


    FreeGen(sgn*m[:gen]; p=abs(pow))
end


prod(matches)

f, g, a = FreeGen("f g -a1")
inv(f)
ja*ja

ff = MonoidHom([
    f=>element("f g"),
    FreeGen("g")=>element("f -1"),
]; strict=true)

w = ff(f*g);
one(w)
w.word
w
(f*g)^g^a

fi = MonoidHom([
    g=>element("g f"),
    f=>element("g -1"),
]; strict=true)
fi * ff
E = MonoidHom([f=>g, g=>f])
fi*E
ff * E
ff*fi
k = 11
u = (ff^k)(f)
v = (fi^k)(g)
E(v)
(u) * (v)
k=1
ff(f^k)/fi(g^k)
ff(f^k)
f / f^-g

(f*g)|fi
f*(g|fi)
((farey_word(4//7, f, g)))^(f*(-f)^-g)

ge = FreeGen(0)
ge = FreeGen(10)

element("a b a1 -2")
ge1 = FreeGen(1; p=-3)
gei = FreeGen("x")
ge1 * gei
f, g = FreeGen("f g")
expand!(inv(farey_word(3//5, f,g)))
expand!(inv(farey_word(4//11, f,g)))

w = farey_word(3//8, f,g)
tr = SFX.SuffixTree(w)

x = ge1 * ge

MonoidWord([1,1,2,10])
w = farey_word(3//5, f,g)
element("f -1")

ww = w^2
length(ww)
# [g for g in eachgen(ww[1:27])]
w = farey_word(3//10, f, g)

f\w|>conj_prefix
#
    ge * one(ge)

MonoidWord(ge)

    a, b = FreeGen("a b")
    x = christoffel(3//5, a, b)

    ss = a^2 * b
    ss.exp
    gei
    a, b = CyclicGen("a b")
    x = farey_word(3//5, a, b)
    [x for x in (gens(inv(x)))]
    [x for x in (gens((x)))]
    inv(x)*x
    # ITR.take(ITR.drop(eachgen(y),2), 3)|>collect
    length(x)
    inv(ge)
    (ge^2)/ge1
    ge^-ge1

    dump(-ge)
    dump(ge*-ge)

    ge * -ge
    gcp(ge1^2, ge1^3)
    ge1
    x * one(ge)
    Fw(Q) = farey_word(Q, f, g)
    gcp(Fw(1//3), Fw(3//8))


    w = inv(x)

    expand!(w)

    ms = eachmatch(r"(?<gen>[A-Z]|[a-z]\d*)(?:[\s\^](?<power>\d*))?[\s\*]*", "a b a")
    [x[:power] for x in ms]


function save_dot_svg(fw; cmd="dot", out="working/output.svg")
    run(`cmd /c "C:\\Program Files\\Graphviz\\bin\\$cmd.exe" -Tsvg $fw -o $out`)
end

using Farey
begin
tr = SFX.SuffixTree("We_follow_Bridson's_proof")
# tr = SFX.SuffixTree("We_may_assume_that_G")
SFX.dotfile!(tr, "working/tree"; direction="LR",  to_string=x->x)

save_dot_svg("working/tree.dot")
end

# Monoids.MonoidHom(Dict([1=>w]))
f^-1 * g

f|ff|E|fi
ff|E|fi
g|ff|E|fi
# SFX.text(tr)
    # [x for x in gens(x)
