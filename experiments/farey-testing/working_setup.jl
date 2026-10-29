using Farey
using Farey.Monoids
const CF = Farey.ContinuedFraction
CF0(L) = CF(0, L)
F2 = Farey.FreeGroup((; :f => :F, :g => :G))
F2ab = Farey.FreeGroup((; :a => :A, :b => :B))
f, g = Farey.GPC.gens(F2)
a1, b1 = Farey.GPC.gens(F2ab)
A, B = Farey.GPC.gens(F2ab) .|> inv
F, G = Farey.GPC.gens(F2) .|> inv
using StyledStrings
using UnicodeFun
include(joinpath(@__DIR__, "../../packages/Farey/working/PillowGroup.jl"))
include(joinpath(@__DIR__, "../../packages/Farey/working/autF2.jl"))
# include("SuffixTrees.jl")

begin # Word Building Functions
  function (w::Farey.FreeGroupElement)(subs::Vararg{T,N}; Upper::Bool=false)::T where {T<:AbstractString,N}
    subsi = Upper ? uppercase.(subs) : (subs .|> x -> x * "⁻")
    FGr = w.parent
    SS = Dict([
      [(x => s) for (x, s) in (zip(GPC.gens(FGr), subs))]...,
      [(inv(x) => s) for (x, s) in (zip(GPC.gens(FGr), subsi))]...
    ])
    prod([SS[FGr(x)] for x in w.word], init=one(subs[1]))
  end

  function (w::Farey.FreeGroupElement)(subs::Vararg{T,N}) where {T<:Farey.AbstractElement,N}
    FGr = w.parent
    SS = Dict([
      [(x => s) for (x, s) in (zip(GPC.gens(FGr), subs))]...,
      [(inv(x) => inv(s)) for (x, s) in (zip(GPC.gens(FGr), subs))]...
    ])
    prod([SS[FGr(x)] for x in w.word], init=one(T))
  end


  function (w::MonoidWord{S})(mats::T; IFUNC=inv2) where {S<:AbstractMonoidGen,T<:AbstractDict}
    # subsi = invFn ? inv.(subs) :  inv2.(subs)

    di = merge(mats, Dict([inv(l) => IFUNC(r) for (l, r) in pairs(mats)]))

    mapreduce(*, ITR.drop(eachgen(w), 1); init=(di[w[1]])) do x
      di[x]
    end
  end



  ↑(f, n) = foldl(∘, ITR.cycle([f], n), init=l -> l)

  function Qtype(Q)
    n = numerator(Q)
    d = denominator(Q)
    if d % 2 == 0
      return PillowGroup.Klein4(1)
    elseif n % 2 == 0
      return PillowGroup.Klein4(2)
    else
      return PillowGroup.Klein4(3)
    end
  end

  function Qhat(rQ; α=2)
    rtyp = Qtype(rQ)
    rhat = rtyp == hK ? CF(α * rQ) :
           rtyp == vK ? CF(rQ // α) :
           CF((1+rQ)/α)
    rhat
  end

  f0, h0, hh, vv, pp = GPC.gens(PillowGroup.pi1OGroup())

  #  a, b = Farey.GPC.gens(F2ab)
  f, g = FreeGen("f g")
  a, b = FreeGen("a b")
  Ω(Q) = christoffel(Q, a0, b0)
  Fw(Q) = ((Ω(Q)^2).elem)(f, g)
  Ψ(Q) = Ω(Q).elem
  Wr(Q) = farey_word(Q, a, b)(f0, h0)
  hK = PillowGroup.Klein4(1)
  vK = PillowGroup.Klein4(2)
  pK = PillowGroup.Klein4(3)

  xx = f0 * hh
  yy = pp
  zz = vv
  a0 = xx * yy
  b0 = zz * yy
  c0 = xx * zz

  Ψa(Q) = Ψ(Q) \ Fw(Q)

  descendQ(Q, Q1; l=-1) = (foldl(⊖, ITR.cycle([Q1], Q[end] + l), init=Q), Q[end])
  Λ(Q) = f0 \ Ω(Q)
  χ(Q) = Ω(Q) / f0
  function aab_dec(rhat)
    t0Q, tQ, t1Q = farey_neighbours(rhat)

    sQ, cQ = ⊕(tQ, t1Q) == rhat ? (t1Q, t0Q) :
             (t0Q, t1Q)
    sQ, cQ
  end

  function sw_descent(Q)
    Q = CF(Q)
    ITR.flatmap(enumerate(Q.L)) do (k, l)
      Lk = Q.L[1:(k-1)]
      [CF(0, [Lk..., l]), CF(0, [Lk..., l + 1])]
    end |> collect

  end

  function rileyw(Q; f=f, g=g)
    p, q = denominator(Q), numerator(Q)
    ϵi = i -> (-1)^Int(floor(i * q // p))
    prod(
      [u^ϵi(i) for (u, i) in zip(ITR.cycle([f, g]), 1:(p-1))], init=one(f))
  end


  function long_n(Q)
    p, q = denominator(Q), numerator(Q)
    l = gcdx(-q, 2 * p)[2]
    l = (l % (2 * p))
    l < 0 ? l + 2 * p : l

  end

  function rev_conj(Q)
    p, q = denominator(Q), numerator(Q)
    l = long_n(Q)
    ϵi = i -> (-1)^Int(floor(i * q // p))
    prod(
      [u^ϵi(i) for (u, i) in zip(ITR.cycle([f, g]), 0:(l))]
    )
  end

  function signs_dec(w; POWERS=false)
    letters = [w.parent(x) for x in w.word]
    genrs = Farey.GPC.gens(w.parent)
    s_seq = []
    y = letters[1]
    out = foldl(
      letters[:end],
      init=(;
        w=y,
        pos=(y ∈ genrs),
        z=(y)
      )
    ) do x, y

      if y == (x.z)
        push!(s_seq, x.w)
        return (;
          w=y,
          pos=(x.pos),
          z=y
        )
        # end

      elseif ((y ∈ genrs) == x.pos)
        return (;
          w=x.w * y,
          pos=x.pos,
          z=y
        )
      else
        push!(s_seq, x.w)
        return (;
          w=(y),
          pos=!x.pos,
          z=y
        )
      end

    end
    push!(s_seq, out.w)
    s_seq

  end

  function signs_seq(w::Farey.AbstractMonoidElement)
    ss = Int[]
    out = foldl(w, init=(; sgn=0, pos=true)) do x, y
      if ispositive(y)==x.pos
        return (; sgn=x.sgn+1, pos=x.pos)
      else
        push!(ss, x.sgn)
        return (; sgn=1, pos=!x.pos)
      end
    end
    push!(ss, out.sgn)
    ss
  end

  function signs_seq(w)
    w == one(w) && return [0]
    letters = [w.parent(x) for x in w.word]
    genrs = Farey.GPC.gens(w.parent)
    s_seq = Int[]
    out = foldl(letters, init=(; sgn=0, pos=true, z=inv(letters[1]))) do x, y

      if y == x.z
        push!(s_seq, x.sgn)
        push!(s_seq, 0)
        return (; sgn=1,
          pos=x.pos,
          z=y
        )

      elseif ((y ∈ genrs) == x.pos)
        return (;
          sgn=x.sgn + 1,
          pos=x.pos,
          z=y
        )
      else
        push!(s_seq, x.sgn)
        return (;
          sgn=1,
          pos=!x.pos,
          z=y
        )
      end

    end
    push!(s_seq, out.sgn)
    s_seq

  end

  ab_gens = Dict([f=>[1,0],g=>[0,1],-f=>[-1,0],-g=>[0,-1]])
  function ab(w::T; ggs=ab_gens) where T<: Farey.AbstractMonoidElement
    sum([ab_gens[l] for l in w],init=[0,0])
  end

  Base.reverse(x::Farey.FreeGroupElement)::Farey.FreeGroupElement = x.parent.([l for l in x.word]) |> reverse |> s -> prod(s, init=one(x))

  # ψf(Q) = [ϕ^l for (l, ϕ) in zip(CF(Q), ITR.cycle([d1, g0]))] |> prod

  # ffh(Q) = ψf(Q) |> hom
  # ffhi(Q) = ψf(Q) |> inv |> hom

  addK(Qa, Qb, k) = foldl(⊕, ITR.cycle([Qb], k), init=Qa)
  Q_path(Q) = [CF(0, Q[2:i]) for i in 2:length(Q)]

end

begin # SL2Z functions

  mobA(Q) = CF(Q / (Q + 1))
  moba(Q) = CF(-Q / (Q - 1))
  mobB(Q) = CF(1 - Q)
  mobd(Q) = CF(1 - 2 * Q)
  mobC = mobA ∘ mobB
  mobc = mobB ∘ moba
  inv2(A) = [
    A[2, 2] -A[1, 2];
    -A[2, 1] A[1, 1]
  ]
  mob(A, x) = CF(
    (A[1, 1] * numerator(x) + A[1, 2] * denominator(x)) // (A[2, 1] * numerator(x) + A[2, 2] * denominator(x))
  )
  Ac(n)::Matrix{BigInt} = [
    n 1;
    -n^2+3n-1 3-n
  ]
  Bc(n)::Matrix{BigInt} = Ac(n) * Ac(n + 1)
  # function matrEval(w::Farey.AbstractElement, subs::Vector{Matrix{S}})::Matrix{S} where S
  #     Dict{UInt, Matrix{S}}([
  #         [x=>]
  #     ])
  # end
  # ITR.flatmap(x->[x, inv(x)], GPC.gens(F2))|>collect

end

begin # palindrome functions
  ispalindrome(w) = (w == reverse(w))
  # function conj_prefix(x)
  #     xi = inv(x)
  #     Gr = x.parent

  #     prefix = zip(x.word, xi.word) |>
  #         itr->ITR.takewhile((x)->x[1]==x[2], itr) |>
  #         l->mapreduce(x->Gr(x[1]), *, l, init=one(Gr))

  #     prefix, x^prefix

  # end
  function pal_suff(w)
    w == one(w) && return w
    for i in 1:(length(w)-1)
      wsuff = w[i:length(w)]
      ispalindrome(wsuff) && return wsuff
    end
  end

  function Pal(w)
    out = one(w)
    Gr = w.parent
    for g in reverse(w).word
      out = g == 1 ? G0(out) * Gr(g) : D0(out) * Gr(g)
    end
    out
  end


  letters(ww::Farey.AbstractElement) = ww[1:length(ww)]

  function lcp(w1, w2)
    mapreduce((a) -> a[1], *, ITR.takewhile(a -> a[1] == a[2], zip(letters.((w1, w2))...)), init=one(w1))
  end

  function lcs(w1, w2)
    reverse(mapreduce((a) -> a[1], *, ITR.takewhile(a -> a[1] == a[2], zip(letters.(reverse.((w1, w2)))...)), init=one(w1)))
  end

  function pref_array(w1, w2)
    [lcp(w1, w2[length(w2)-i:length(w2)]) for i in 1:length(w2)-1]
  end

  function suff_array(w1, w2)
    [lcs(w1, w2[1:i]) for i in 2:length(w2)]
  end

  pal_conj(w) = (w / pal_suff(w)) * reverse(w)




  # rileyw(Qp)|>signs_seq
end

begin #decomposing farey words
  descendQ(Q, Q1; l=0) = (foldl(⊖, ITR.cycle([Q1], Q[end] + l), init=Q), Q[end])

  function FareyDecomp!(qq)
    f, g = FreeGen("f g")

    qq = CF(qq)
    qa, qb, qc = farey_neighbours(qq)
    if qq[end] % 2 == 0

      qu, l = descendQ(qq, qa; l=1)
      k = l ÷ 2
      qr = qa
      qt = qu ⊖ qa
    else
      qu, l = descendQ(qq, qa; l=0)
      k = l ÷ 2
      qr = qa
      qt = qu ⊕ qa

    end
    pow(k) = to_superscript("$k")
    npow(k) = to_superscript("-$k")
    # qr = addK(qv, qu, l1); qt = addK(qv, qu, l1+1);
    printstyled(styled"{yellow:$(addK(qt, qr, k*2))}   {code:qq = $(qq)}\n{tip:k = $k}; {bright_blue:(k-$(l%2))/2 = $(k÷2)}; {magenta:l = $l}\n\n")
    printstyled(styled"{red:qr = $(qr)}   {green:qt = $(qt)}\n{bright_blue:qc = $(qc)} {yellow:qu = $(qu)} {magenta:qa = $(qa)}\n\n")

    fwr = farey_word(qr,f,g)
    fwt = farey_word(qt,f,g)
    nt = 2 * numerator(qt)
    dt = 2 * denominator(qt)
    nr = 2 * numerator(qr)
    dr = 2 * denominator(qr)
    if qr < qt
      printstyled("" *
                  (styled"{bright_magenta:$(fwr)}" *
                   styled"{bright_blue:$(fwr)}")^(k ÷ 2) *
                  styled"{red:$(fwr)}"^(k % 2) *
                  styled"{green:$(fwt)}" *
                  styled"{yellow:$(inv(f))}" *
                  styled"{bright_red:$(inv(fwr))}"^(k % 2) *
                  (styled"{magenta:$(inv(fwr))}" *
                   styled"{code:$(inv(fwr))}")^(k ÷ 2) *
                  styled"{yellow:f}\n"
      )
      printstyled(styled"{emphasis:$(farey_word(qq,f,g))}\n")
      printstyled(styled"\n{red:Wr[$(Rational(qq))]} \t= \t...{emphasis:$(farey_word(qq,f,g)[k*dr+1:k*dr+dt+1])}...\n")
      printstyled(styled"{red:[$(k*dr+1)...$(k*dr+dt+1)]}\n")
      printstyled(styled"{bright_blue:Wr[$(Rational(qt))]/f} \t= \t   {magenta:$(fwt)}|{yellow:$(F)}\n")

      printstyled(styled"\n{red:Wr[$(Rational(qq))]} \t= \t{emphasis:$(farey_word(qq,f,g)[1:dr])}...\n")
      printstyled(styled"{red:[1...$(dr+1)]}\n")
      printstyled(styled"{bright_blue:Wr[$(Rational(qr))]} \t= \t{magenta:$(fwr)}\n")
      printstyled(styled"\n{code:Wr[$(Rational(qq))]} \t= {bright_red:Wr}{bright_blue:$(pow(k))}{bright_red:[$(Rational(qr))]}.{green:Wr[$(Rational(qt))]}.{yellow:$(inv(f))}.{bright_red:Wr}{bright_blue:$(npow(k))}{bright_red:[$(Rational(qr))]}.{yellow:$(f)}\n")

      printstyled(styled"\n \t\t= ({bright_red:$(fwr)})$(pow(k)){green:$(fwt)}({yellow:$(inv(f))}{bright_red:$(inv(fwr)*f)})$(pow(k))\n")
    else

      printstyled("" *
                  styled"{yellow:$(f)}" *
                  styled"{bright_red:$(inv(fwr))}"^(k % 2) *
                  (styled"{magenta:$(inv(fwr))}" *
                   styled"{code:$(inv(fwr))}")^(k ÷ 2) *
                  styled"{yellow:$(inv(f))}" *
                  styled"{green:$(fwt)}" *
                  (styled"{bright_magenta:$(fwr)}" *
                   styled"{bright_blue:$(fwr)}")^(k ÷ 2) *
                  styled"{red:$(fwr)}\n"^(k % 2)
                  * "\n"
      )
      printstyled(" ↓\n" *
                  styled"{yellow:$(f)}" *
                  styled"{bright_red:$(inv(fwr))}"^(k % 2) *
                  (styled"{magenta:$(inv(fwr))}" *
                   styled"{code:$(inv(fwr))}")^(k ÷ 2) *
                  styled"{green:$(fwt[2:end])}" *
                  (styled"{bright_magenta:$(fwr)}" *
                   styled"{bright_blue:$(fwr)}")^(k ÷ 2) *
                  styled"{red:$(fwr)}"^(k % 2)
                  * "\n"
      )

      printstyled(styled"{emphasis:$(farey_word(qq,f,g))}\n")


      printstyled(styled"\n{red:Wr[$(Rational(qq))]} \t= \t  ...{emphasis:$(farey_word(qq,f,g)[k*dr+2:k*dr+dt])}...\n")
      printstyled(styled"{red:[$(k*dr+2)...$(k*dr+dt)]}\n")
      printstyled(styled"{bright_blue:$(F).Wr[$(Rational(qt))]} \t= \t{yellow:$(F)}|{magenta:$(fwt)}\n")

      printstyled(styled"\n{red:Wr[$(Rational(qq))]} \t= \t{emphasis:$(farey_word(qq,f,g)[1:dr+1])}...\n")
      printstyled(styled"{red:[1...$(dr+1)]}\n")
      printstyled(styled"{bright_blue:f/Wr[$(Rational(qr))]} \t= \t{yellow:$(f)}|{magenta:$(inv(fwr))}\n")
      printstyled(styled"\n{code:Wr[$(Rational(qq))]} \t= {yellow:$(f)}{bright_red:Wr}{bright_blue:$(npow(k))}{bright_red:[$(Rational(qr))]}{yellow:$(inv(f))}{green:Wr[$(Rational(qt))]}{bright_red:Wr}{bright_blue:$(pow(k))}{bright_red:[$(Rational(qr))]}\n")

      printstyled(styled"\n \t\t= ({yellow:$((f))}{bright_red:$(inv(fwr))})$(pow(k)){green:$((-f)*fwt)}({bright_red:$(fwr)})$(pow(k))\n")
    end
  end
end


function qpair(qq)
  qq = CF(qq)
  qr = CF(0, qq[2:end-1])
  qt = qq[end] % 2 == 0 ? CF(0, qq[2:end-2]) : CF(0, [qq[2:end-1]..., 1])
  k = qq[end] ÷ 2

  qr, qt, k
end

function qnbs(Q::T) where T<:Number
  Q = CF(Q)
  qs = CF.([1 // (Q[2])])
  for (k, l) in enumerate(Q.L[2:end])
    append!(qs, [positive_form(CF(0, [Q[2:k+1]..., x])) for x in 1:l])
  end
  qs
end

function commu(f,g)
f/(f^-g)
end


function farey_expr(qq::Number)
  qq = CF(qq)

  idxr = 1 + sum(qq[3:end-1])
  idxt =  qq[end] % 2 == 1 ?  idxr + 1 :
          length(qq)>3 ?     1 + sum(qq[3:end-2]) :
          0

  k = qq[end] ÷ 2

  w1, w2, f = FreeGen("W$idxr W$idxt f")

  # length(qq) % 2 == 0 ? w2 * inv(w1^f)^k*w1^k : w2 * w1^k*inv(w1^-f)^k
  length(qq) % 2 == 0 ? w1^k * w2 * inv(w1^f)^k : inv(w1^-f)^k * w2 * w1^k
end
function expandhom(q1)
  q1 = CF(q1)
  qns = qnbs(q1)

  f, g, w, ws... = FreeGen("f g W0 " * join(["W$i " for i in eachindex(qns)]))
  t0, t1 = FreeGen("t$(q1[2]) t$(q1[2]+1)")
  u, v = q1[2] % 2 == 0 ? (-f, -g) : (-g, -f)
  MonoidHom([[
      ws[1] => f * t0 * u / t0,
      ws[2] => f * t1 * v / t1,
    ]..., w => farey_word(0, f, g)]),
  MonoidHom([
    ws[i] => farey_expr(q) for (i, q) in Iterators.filter(x -> length(x[2]) > 2, enumerate(qns))
  ])
end


function pretty_printing(x::T; k=5) where T
  if T <: MatElem || T<:AbstractMatrix
    s = "[\n" * join(["$(Tuple(l)):\t" * repr(x[l]) * "\n" for l in eachindex(x)]) * "]"

  else
    s = repr(x)
  end


  s = replace(s,

    # r"[\,|]"=>"\n]",
    # r"(x|z)\^?-?\d+" => l -> "{magenta:$("𝒵"*to_superscript(l[3:end]))}",
    r"x|z" => l -> "{magenta:$("𝒵")}",
    # r"z(?<subscr>\_\d+)(?<power>\^-?\d*)"=>s"z\g<power> z\g<subscr>"
  )
  s = replace(s,
    r"\^-?\d+" => l -> "{bright_magenta:$(to_superscript(l[2:end]))}",
    r"z\_\d+" => l -> "{bright_blue:ζ$(to_subscript(l[3:end]))}",
    r"\*" => l -> "",
    r"\d+" => l -> "{code:$l}",
  )

  s = replace(s,
    Regex("\\d{$(2*k),}") => l -> "($(l[1:k]).($(length(l))).$(l[end-k+1:end]))"
  )
  # println("hh")
  StyledStrings.styled(s)

end

function mateval(x::Farey.AbstractMonoidElement, mats::Dict{K,V}; IFUNC=inv2)::V where {K,V}
  di = merge(mats, Dict([inv(l) => IFUNC(r) for (l, r) in pairs(mats)]))

  mapreduce(*, eachgen(x)) do x
    di[x]
  end
end


function monoideval(x::Farey.AbstractMonoidElement, di::Dict{K,V}; init=one(first(values(di))))::V where {K,V}
  mapreduce(*, eachgen(x); init=init) do x
    di[x]
  end
end
