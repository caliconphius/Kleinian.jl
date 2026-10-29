Qb = algebraic_closure(QQ)
Rx, x = polynomial_ring(QQ, "z")

Kx = fraction_field(Rx)
bb = one(x)#root_of_unity(Qb, 3, 1)
aa = one(x)#oot_of_unity(Qb, 4, 1)
q = @cf0 [1,1,2]
CC = CalciumField()


B = Kx[
  bb BigInt(0);
  x^2 1/bb
]


A = Kx[
  aa BigInt(1);
  0 1/aa
]

J = Kx[
  Qb(im)*x 0;
  0 -Qb(im)*1//x;
  ]

Θ(Q) = -conj_prefix(f\farey_word(Q, f, g))[1]
Fw(Q) = farey_word(Q, f, g)

ms = [A, B]
Θ(1//3)
ZZx, zx = polynomial_ring(ZZ,:zx)

function ρw(W; gens=gens, ms=ms, Fk = Qb)
  mats = Dict([l => r for (l, r) in zip(gens, ms)])
  rW = (mateval(W, mats)) |>
    # z->map_coefficients(BigInt, numerator(z))
    z->z
  Zx, x = polynomial_ring(Fk,:z)
  FF = z->map_coefficients(BigInt, numerator(z))
  # FF = z->z
  rW = [el|>FF|>l->Zx(l.coeffs) for el in rW]
  # pretty_printing(rW) |> printstyled
  # println("\n")
  return rW
end
# J *
import Nemo:tr

Kz, xs = polynomial_ring(Qb, [:z])
Kx, x = polynomial_ring(QQ, :z)
zz = xs[1]
Zx, x = polynomial_ring(ZZ, :x)

# Zx, x = polynomial_ring(QQ,:z)
ms_v = Kz[2 2 zz^2 + 2] |> collect |> x->reshape(x, 3)
l, m, n = (4, 5, 7)
lcm(l,m,n)
Ka, za = Nemo.cyclotomic_real_subfield(2*l)
Kb, zb = Nemo.cyclotomic_real_subfield(2*m)
Kc, zc = Nemo.cyclotomic_real_subfield(2*n)
p1, p2, p3 = minpoly(Kx, za),minpoly(Kx, zb),minpoly(Kx, zc)


L, gs = Nemo.number_field([p1, p2, p3], "a")
uu = Hecke.absolute_primitive_element(L)
Lz, zl = L[:z]
a1, a2, a3 = gs
Hecke.absolute_basis(L)
sL, u = Hecke.simplified_simple_extension(L)
ua = u(gen(sL))

ui = inv(u)
c1,c2,c3 = ui.([a1,a2,a3])


pu = minpoly(Kx, ua)

rtu = roots(Qb, pu)
a0 = CC(rtu[1])
A1 = Kx(Hecke.coordinates(c1))(a0)|>Qb
A2 = Kx(Hecke.coordinates(c2))(a0)|>Qb
A3 = Kx(Hecke.coordinates(c3))(a0)|>Qb
# Hecke.splitting_field(p3)
A2
# p1(A1)

Zb, zx = polynomial_ring(Qb, :x)

dd, s, t =gcdx(pu, p1)


pu(p1)
t
s
zz = xs[1]

ms_v = Zb[A1 A2 zx] |> collect |> x->reshape(x, 3)
Q = @cf0 [1,1,2]
Q = @cf0 [2]

n1, n2, n3 = gs
gcdx(pu, p2)
p1(ua)

function zalpha(Q)
  rtyp = Qtype(Q)
  rhat =
    rtyp == hK ?
    g  :
    rtyp == vK ?
    -f :
    f
  rhat
end

Rw(Q) = zalpha(Q)^-rileyw(Q) / g

Rw(3//8)
Fw(3//8)

Rw(2//7)
Fw(2//7)
