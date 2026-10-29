
include("sl2algebra.jl")
include("field_calc_hecke.jl")
include("erics_functions.jl")


l, m, n = (3, 5, 2)
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

Zb, zx = polynomial_ring(Qb, :x)
Kx, x = polynomial_ring(QQ, :z)
Lz, zl = L[:z]
a1, a2, a3 = gs

Q = @cf0 [1,1,2]
ms_vb = Zb[A1 A2 zx] |> collect |> x->reshape(x, 3)
ms_v = Kx[2 2 x] |> collect |> x->reshape(x, 3)
ms_v = Lz[a1 a2 zl] |> collect |> x->reshape(x, 3)
eW = monoideval(Fw(Q), mons);
Fp = evaluate(change_base_ring(QQ, trA(eW)), ms_v)
Fp = evaluate(change_base_ring(Qb, trA(eW))-A3, ms_vb)

Cb = complex_field()

# rt_sols = Threads.@spawn roots(Qb, Fp)
rt_sols
fproot = rt_sols.result
# Threads.@spawn
z0 = fproot[1]

psol = z0|>minimal_polynomial
# fcs = factor(Qb, change_base_ring(QQ, psol))
fcs.fac
Fm = change_base_ring(Qb, psol)
d, s, t = gcdx(Fm, Fp)
Fm//Fp

# factor(Cb, change_base_ring(Qb, (Fm//Fp)|>numerator))
pp = coefficients(psol)
pp|>collect
CC.(fproot)


# save_object( "($l,$m,$n,$(Rational(Q)|>x->"$(numerator(x))-$(denominator(x))")).jld2", rt_sols)
z0|>dump
CC(z0)
minimal_polynomial(z0)
roots(Cb, Fm)

# Q = @cf0 [2,2,4,-2,-2,2]
Cb = complex_field()
#
Q = 3//5
eW = monoideval(Fw(Q), mons);

CC = CalciumField()
Acb = AcbField(128)
Arb = ArbField(128)
Fp = evaluate(change_coefficient_ring(Qb, trA(eW)), ms_v)
FpAcb = change_coefficient_ring(Acb,Fp)
function big_calc(;K=Acb, Fp = FpAcb)
  r = roots(Fp-A3)
  @info "finished: $r"
  r
end
resultQb2 = Threads.@spawn big_calc(;K=Qb, Fp=Fp)
resultC = Threads.@spawn big_calc(;K=Acb, Fp=FpAcb)
z00 = resultQb2.result[1]
Pf = minimal_polynomial(z00)
gcdx(change_base_ring(Qb, Pf), Fp)

resultQb2.result
resultC
# z00 = resultC.result|>first
# z01 = resultC.result[2]



ZZx(cfs)(z00)

# result.result
# using JLD2
# save("big_trace_obj.jld2", Dict("eW"=>eW))
# JLD2.save("fourtitwo.jld2", Dict(["$k"=>r for (k,r) in enumerate(resultC.result)]));

f1 = resultC.result.fac|>keys|>collect|>first
p00 = z00|>minimal_polynomial|>x->change_coefficient_ring(ZZ, x)
# p00|>x->factor(Cb, x)

trA(eW)
