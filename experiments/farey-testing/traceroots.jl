# al, am, an = Int(lmn//l), Int(lmn//m), Int(lmn//n)
# ul = 2 * chebyshev_t(al, 1//2 * zz)(zl)
# um = 2 * chebyshev_t(am, 1//2 * zz)(zl)
# un = 2 * chebyshev_t(an, 1//2 * zz)(zl)
l,m,n = (3,4,7)
lmn = lcm(l,m,n)
L, (zl, zm, zn) = number_field([cos_minpoly(2l, xx),cos_minpoly(2m, xx),cos_minpoly(2n, xx)], [:zl, :zm, :zn])

Lz, zlx = L[:z]
ms_t = Lz[zl zm zlx] |> collect |> x->reshape(x, 3)
# ms_t = Lz[2 2 zlx] |> collect |> x->reshape(x, 3)
bigpoly = Threads.@spawn ev(Fw(3//5); monsUW = mons)|>trA;
bigpoly
Fqtr = bigpoly.result
Flmn = evaluate(Fqtr, ms_t)#|>coefficients|>collect
Flmn
cms = Hecke.complex_embeddings(L)
map_coefficients(cms[1], Flmn-zn)|>roots
Flmn(1//2)

sL, u = Hecke.simplified_simple_extension(L)
ua = u(gen(sL))

ui = inv(u)
xcl,xcm,xcn = ui.([zl, zm, zn])



pu = minpoly(Kx, ua)
rtu = roots(Qb, pu)
a0 = CC(rtu[1])
A1 = Kx(Hecke.coordinates(xcl))(a0)|>Qb
A2 = Kx(Hecke.coordinates(xcm))(a0)|>Qb
A3 = Kx(Hecke.coordinates(xcn))(a0)|>Qb
Qbx, xb = Qb[:x]
ms_cc = Qbx[A1 A2 xb] |> collect |> x->reshape(x, 3)
Flmn = evaluate(change_base_ring(Qb,Fqtr), ms_cc)-A3#|>coefficients|>collect
rtslmn = Threads.@spawn Flmn|>roots
rtslmn
ttx  = rtslmn.result[1]|>minimal_polynomial
rtsttx = Threads.@spawn ttx|>x->roots(Qb, x)
rtslmn.result
rtsttx.result
ttx
Flmn
A3|>minimal_polynomial
