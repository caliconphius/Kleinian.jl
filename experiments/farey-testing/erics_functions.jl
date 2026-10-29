using Memoize
using Farey
using Farey: Monoids
CC = algebraic_closure(QQ)
# using Oscar
#
Kx, x = polynomial_ring(QQ, :z)
using Memoize

function dQ(Q; x=x)
  (-1)^numerator(Q) * x^denominator(Q)
end
dQ(3 // 5)

function aab_dec(rhat)
  t0Q, tQ, t1Q = farey_neighbours(rhat)

  sQ, cQ = ⊕(tQ, t1Q) == rhat ? (t1Q, t0Q) :
           (t0Q, t1Q)
  sQ, cQ
end

aab_dec(1 // 1)


qnbs(3 // 5) .|> qpair

qpair(3 // 10)

aab_dec


function funnel(Q::T) where T<:Number
  Q = CF(Q)
  qs = SBEdge[SBEdge(0, 1 // i) for i in 1:(Q[2])]
  qhub = CF(1 // Q[2])
  for (k, l) in enumerate([Q.L[2:end-1]..., Q[end] - 1])
    spokes = [positive_form(CF(0, [Q[2:k+1]..., x])) for x in 1:l]
    append!(qs, [SBEdge(qhub, s) for s in spokes])
    qhub = spokes[end]
  end
  qs
end


@memoize function Qpoly(Q::Number; x=x)::QQPolyRingElem
  Q == 1 // 0 && return zero(x)
  0 <= Q <= 1 || error("$Q ∉ ℚ₀")
  Q in [0 // 1, 1 // 1] && return one(x)

  qa, qb = aab_dec(Q)
  -dQ(qa; x=x) * Qpoly(qb; x=x) + Qpoly(qa; x=x) * Qpoly(qa ⊕ qb; x=x)
end

@memoize function Npoly(Q::Number; x=x)::QQPolyRingElem
  # @info "computed $Q"
  Q in [1 // 0, 0 // 1, 1 // 1] && return one(x)

  qa, qb = aab_dec(Q)
  -dQ(qa; x=x) * Npoly(qb; x=x) + Qpoly(qa; x=x) * Npoly(qa ⊕ qb; x=x)
end
Qpoly(2 // 5)
Qpoly(3 // 5)
Qpoly(3 // 10)

V(Q) = Npoly(Q) // Qpoly(Q)
V(9 // 31)

Qx = fraction_field(Kx)

Wmat(Q) = Qx[
  1 0;
  -Qpoly(Q)^2//dQ(Q) 1
]

Umat(Q) = Qx[
  1 -V(Q);
  0 1
]
QMx = Umat(1 // 4) |> parent
inv2x(A) = QMx(inv2(A)) * 1 // Nemo.det(A)
Wmat(1 // 4) |> inv2x

Tmat(Q) = inv2x(Umat(Q)) * inv2x(Wmat(Q)) * (Umat(Q))
Smat(Q) = inv2x(Tmat(Q)) * inv2x(Umat(0))
Tmat(1 // 2)
Q = @cf 3 // 8
Tmat(Q)
q1, q2, q3 = farey_neighbours(Q)

struct SBEdge
  qL::ContinuedFraction
  qR::ContinuedFraction
  qP::ContinuedFraction
  SBEdge(q1, q2) = new(sort([CF(q1), CF(q2)])..., q1 ⊕ q2)
end



function SBEdge(q)
  q1, q2, _ = farey_neighbours(q)
  SBEdge(q1, q2)
end

function highlow(s::SBEdge)
  sort([s.qL, s.qR], by=denominator)
end

function lr(s::SBEdge)
  [s.qL, s.qR]
end

Q = @cf 24 // 103
highlow.(funnel(Q))
lr.(funnel(Q))

SBEdge(3 // 8)
Z(e::SBEdge) = (V(e.qR) - V(e.qP)) // (V(e.qL) - V(e.qP))
Z(e::SBEdge) = -(2 - tr(Smat(e.qL))) // (2 - tr(Smat(e.qR)))

Z(q::Number) = Z(SBEdge(q))
Z(q1::Number, q2::Number) = Z(SBEdge(q1, q2))
Z(3 // 8)



function discriminantQ(Q::Number; x=x)
  Qpoly(Q; x=x)^2 - 4 * dQ(Q; x=x)
end

tr(Smat(1 // 5))
discriminantQ(1 // 5)
discriminantQ(3 // 8)

Xmat(Q) = Qx[
  0 1;
  -dQ(Q) Qpoly(Q)
]

function Ymat(Q)
  e = SBEdge(Q)
  b0 = e.qR
  b1 = e.qR ⊕ e.qP
  Qx[
    -dQ(Q)*Npoly(b0) Npoly(b1);
    -dQ(Q)*Qpoly(b0) Qpoly(b1)
  ]
end

Rmat(Q) = Ymat(Q) * Xmat(Q) * inv2x(Ymat(Q))
(Rmat(2 // 3)^2 + Smat(2 // 3))
Tmat(2 // 3)
Umat(0) |> inv2x


TEST_S(Q) = Qx[
  dQ(Q)-Qpoly(Q)*Npoly(Q) Npoly(Q)^2+dQ(Q)-Qpoly(Q)*Npoly(Q);
  -Qpoly(Q)^2 -Qpoly(Q)^2+dQ(Q)+Qpoly(Q)*Npoly(Q)
]


q = 5 // 12
(Smat(q) |> tr) - 2 +
Qpoly(q)^2 // dQ(q)

# Rmat(q)

[s(-1 // x) for s in (Smat(3 // 5))]
[s(-1 // x) for s in (Rmat(3 // 5))]
Tmat(3 // 5)
q0 = 0
q1 = 2 // 7
monsUW = Dict([f => inv2x(Umat(q0)), g => Wmat(q0), -f => Umat(q0), -g => inv2x(Wmat(q0))])
monsUW0 = Dict([f => inv2x(Umat(q1)), g => Wmat(q1), -f => Umat(q1), -g => inv2x(Wmat(q1))])
eW = monoideval(Fw(1 // 2), monsUW)
ev(w; monsUW=monsUW) = monoideval(w, monsUW)
evM(w; monsUW=monsUW) = [s(-1 // x) for s in (ev(w; monsUW=monsUW))]
tr(evM(f * g))
(evM(farey_word(q1, -f, -g); monsUW=monsUW0))
(evM(Fw(1 // 2)))

# ev(Θ(3//5); monsUW=mons)

sw = funnel(5 // 12)[end]
qL, qR, qP = sw.qL, sw.qR, sw.qP
Smat(qL)
ev(Fw(qL))
ev(f), ev(g)


Tmat(qP)
Xm = parent(x)[
  1 1;
  0 1
]
Ym = parent(x)[
  1 0;
  x 1
]

Xm * inv(Ym)
Xm * (Ym) * inv(Ym * (Xm))

function make_evens(L::Number)
  perpFlag = false
  q = ContinuedFraction(L)
  if !iseven(denominator(q) * numerator(q))
    q = (CF(q.leading, [q.L[1:end-1]..., q.L[end] - 1]))
    perpFlag = true
  end
  q = ContinuedFraction(q)
  ll = length(q)
  remains = 0
  s = 1
  Lnew = BigInt[]
  for (k, p) in enumerate(q.L)
    if remains == 1

      p = s * p > 0 ? s + p : s - p
      remains = 0
    end
    if isodd(p) && remains == 0
      s * p > 0 ? push!(Lnew, s + p) : push!(Lnew, s - p)
      remains = -1
      s = -s
    elseif remains == -1
      append!(Lnew, vcat([[2s, -2s] for i in 1:(p÷2+1)]...)[1:p-1]...)
      s = isodd(p) ? s : -s
      remains = 1
    else
      s * p > 0 ? push!(Lnew, p) : push!(Lnew, -p)

    end
  end
  perpFlag ? CF0([Lnew..., sign(Lnew[end])]) : CF0(Lnew)
end


function cf_matrix(q; z=x, evens=true)

  Q = evens ? ContinuedFraction(q) |> make_evens |> first : ContinuedFraction(q)
  ms = [parent(z)[
    1 l;
    0 1
  ] for l in Q.L]

  mapreduce(
    M -> iseven(M[1]) ? M[2] : transpose(M[2]),
    *,
    enumerate(ms); init=one(ms[1])
  ) * parent(z)[
    1 z;
    0 1
  ]
end
funnel(1 // 4)
@memoize function geometric_root(Q)
  Q = Rational(Q)
  CC = algebraic_closure(QQ)
  rlist = roots(CC, Qpoly(Q))
  for fe in funnel(Q)[2:end-1]
    rlist = filter(x -> imag(Z(fe)(x)) > 0, rlist)
  end
  rlist
end
# geometric_root(24 // 103)

# roots(CalciumField(), Qpoly(2 // 5))

make_evens(3 // 5)

cf_matrix(2 // 11)
cf_matrix(2 // 5)
cf_matrix(5 // 24)
cf_matrix(11 // 27)
q = 8 // 19 |> CF
q1 = CF0(reverse(CF(q).L |> collect))
T = cf_matrix(3 // 5)
(cf_matrix(-q))

CC = Nemo.algebraic_closure(QQ)
eigenspaces(cf_matrix(24 // 103; z=CC(0))) |> keys |> first |> minimal_polynomial
make_evens(7 // 19)
make_evens(7 // 17)
Zx, xz = ZZ[:x]
zz = Kx([0, 1])
cos_minpoly(15, xz)
T, (u, v, w) = QQ[:u][:v][:w]
Q = 8 // 19
cf_matrix(Q; z=CC(0)) |> eigenspaces |> values .|> x -> minimal_polynomial.(x)
cf_matrix(Q; z=CC(0)) |> eigenspaces |> keys .|> x -> minimal_polynomial.(x)

l, m, n = (2, 3, 7)
lmn = lcm(l, m, n)



# L, (zl, zm, zn) = number_field([cyclotomic(2l, xx),cyclotomic(2m, xx),cyclotomic(2n, xx)], [:zl, :zm, :zn])
# Lz, zlx = L[:z]


# evM(Fw(3//5))
# Wmat(0)
# zq = Lz|>gen
# W0 = Lz[
#   zm 0 ;
#   zq 1/zm
#   ]

# U0 = Lz[
#   zl -1 ;
#   0 1/zl
#   ]

# LM = parent(U0)

# U0i = inv2(U0)|>LM
# W0i = inv2(W0)|>LM

# monsUW0 = Dict([f=>U0, g=>W0, -f=>U0i, -g=>W0i])

# ev(Fw(3//5); monsUW=monsUW0)|>tr

# sL, u = Hecke.simplified_simple_extension(L)
# ua = u(gen(sL))
# sLz, sz = sL[:z]

# ui = inv(u);
# szl,szm,szn = ui.([zl, zm, zn]);
# szq = gen(sLz)
# sW0 = sLz[
#   szm 0 ;
#   szq 1/szm
#   ];

# sU0 = sLz[
#   szl -1 ;
#   0 1/szl
#   ];




# sLM = parent(sU0)
# sU0i = inv2(sU0)|>sLM;
# sW0i = inv2(sW0)|>sLM;
# smonsUW0 = Dict([f=>sU0, g=>sW0, -f=>sU0i, -g=>sW0i]);

# qb = @cf 24//103
# Φ = ev(Fw(qb); monsUW=smonsUW0)|>tr;
# Φ|>coefficients
# scms = Hecke.complex_embeddings(sL)
# scms
# pfc = map_coefficients(z->scms[3](z,4*128), Φ-szn)
# rtf = pfc|>roots
# z0 = rtf[6]
# pfc(z0)

# Ct = parent(z0)
# Crt = parent(real(z0))
# for (k,z0) in enumerate(rtf)|>collect
# vvsr = ([round(imag(Z(e)(z0))|>BigFloat; sigdigits=4)|>Float64 for e in funnel(qb)])
# vvs = ([(imag(Z(e)(z0))) for e in funnel(qb)])

# printstyled(
#   styled"\n{yellow:$k}:{code:$(vvsr)}{red:$(all((vvs.>-1e-4)))}\n"
# )
# end


