module SakumaWeeks

import Farey
import Nemo
using Kleinian, Poincare

export PolynomialSystem, SBEdge, funnel, interior_edges, aab_dec, dQ, Qpoly,
       Npoly, V, Z, discriminantQ, geometric_roots, geometric_root,
       Umat, Wmat, Tmat, Smat, boundary_isometry, developed_domain

const Slope = Rational{BigInt}
const INF = big(1)//big(0)
slope(q::Union{Integer,Rational,Farey.ContinuedFraction}) = Slope(Rational(q))
function unit_slope(q)
    r = slope(q)
    (r == INF || 0 <= r <= 1) || throw(ArgumentError("slope must lie in [0,1] or be infinity"))
    r
end
mediant(a,b) = slope(Farey.:⊕(a,b))

"A polynomial ring and caches scoped to this calculation, not to global REPL state."
struct PolynomialSystem{P,F}
    x::P
    field::F
    qcache::Dict{Slope,P}
    ncache::Dict{Slope,P}
end
function PolynomialSystem()
    R,x = Nemo.polynomial_ring(Nemo.QQ, :x)
    PolynomialSystem(x,Nemo.fraction_field(R),
        Dict(big(0)//1=>one(x),big(1)//1=>one(x),INF=>zero(x)),
        Dict(big(0)//1=>one(x),big(1)//1=>one(x),INF=>one(x)))
end

# Adapted from farey-testing/erics_functions.jl: r = b ⊕ a ⊕ a.
function aab_dec(q)
    r = unit_slope(q)
    0 < r < 1 || throw(ArgumentError("decomposition requires 0 < slope < 1"))
    t0,t,t1 = slope.(Farey.farey_neighbours(r))
    mediant(t,t1) == r ? (t1,t0) : (t0,t1)
end
dQ(s::PolynomialSystem,q) = (r=unit_slope(q); (-1)^numerator(r)*s.x^denominator(r))
function Qpoly(s::PolynomialSystem,q)
    r = unit_slope(q)
    get!(s.qcache,r) do
        a,b = aab_dec(r)
        -dQ(s,a)*Qpoly(s,b)+Qpoly(s,a)*Qpoly(s,mediant(a,b))
    end
end
function Npoly(s::PolynomialSystem,q)
    r = unit_slope(q)
    get!(s.ncache,r) do
        a,b = aab_dec(r)
        -dQ(s,a)*Npoly(s,b)+Qpoly(s,a)*Npoly(s,mediant(a,b))
    end
end
V(s::PolynomialSystem,q) = s.field(Npoly(s,q))/s.field(Qpoly(s,q))
discriminantQ(s::PolynomialSystem,q) = Qpoly(s,q)^2-4*dQ(s,q)

"A validated Stern–Brocot edge, ordered left to right."
struct SBEdge
    qL::Slope
    qR::Slope
    qP::Slope
    function SBEdge(a,b)
        l,r = sort([unit_slope(a),unit_slope(b)])
        r <= 1 || throw(ArgumentError("edge endpoints must be finite"))
        abs(numerator(l)*denominator(r)-numerator(r)*denominator(l)) == 1 ||
            throw(ArgumentError("endpoints must be Farey neighbours"))
        new(l,r,mediant(l,r))
    end
end
opposite(e::SBEdge) = (numerator(e.qR)-numerator(e.qL))//(denominator(e.qR)-denominator(e.qL))
nu(e::SBEdge) = denominator(e.qR) < denominator(e.qL) ? 1 : -1
highlow(e::SBEdge) = sort([e.qL,e.qR];by=denominator)

"Top-to-bottom funnel edges, including the two boundary edges."
function funnel(q)
    r = unit_slope(q)
    0 < r < 1 || throw(ArgumentError("funnel requires 0 < slope < 1"))
    l,h = big(0)//1,big(1)//1
    es = SBEdge[]
    while true
        e = SBEdge(l,h)
        push!(es,e)
        e.qP == r && return es
        if r < e.qP
            h = e.qP
        else
            l = e.qP
        end
    end
end
interior_edges(q) = funnel(q)[2:end-1]

# Lemma 5.5 / Definition 4.17; avoids the notebook's overwritten Z methods.
Z(s::PolynomialSystem,e::SBEdge) = -s.field(dQ(s,e.qR))*s.field(Qpoly(s,e.qL))^2 /
    (s.field(dQ(s,e.qL))*s.field(Qpoly(s,e.qR))^2)

"Exact algebraic roots satisfying positivity for every interior shape."
function geometric_roots(s::PolynomialSystem,q)
    r = unit_slope(q)
    (0 < r < 1//2 && numerator(r) != 1) ||
        throw(ArgumentError("use a hyperbolic two-bridge slope in (0,1/2), excluding 1/n"))
    es = interior_edges(r)
    C = Nemo.algebraic_closure(Nemo.QQ)
    filter(Nemo.roots(C,Qpoly(s,r))) do z
        all(es) do e
            f = Z(s,e)
            den = Nemo.evaluate(Nemo.denominator(f),z)
            !iszero(den) && imag(Nemo.evaluate(Nemo.numerator(f),z)/den) > 0
        end
    end
end
function geometric_root(s::PolynomialSystem,q)
    rs = geometric_roots(s,q)
    length(rs) == 1 || error("expected one positive-shape root; found $(length(rs))")
    only(rs)
end
polyeval(p,z) = Nemo.evaluate(p,z)
function polyeval(p,z::Union{AbstractFloat,Complex{<:AbstractFloat}})
    result = zero(z)
    for k in Nemo.degree(p):-1:0
        result = result*z + typeof(z)(Rational{BigInt}(Nemo.coeff(p,k)))
    end
    result
end
function value(f,z)
    den = polyeval(Nemo.denominator(f),z)
    iszero(den) && throw(ArgumentError("rational function has a pole at the chosen parameter"))
    polyeval(Nemo.numerator(f),z)/den
end

function Umat(s::PolynomialSystem,q)
    F = s.field
    [F(1) -V(s,q); F(0) F(1)]
end
function Wmat(s::PolynomialSystem,q)
    F = s.field
    [F(1) F(0); -F(Qpoly(s,q))^2/F(dQ(s,q)) F(1)]
end
inv2(A) = [A[2,2] -A[1,2]; -A[2,1] A[1,1]] ./ (A[1,1]*A[2,2]-A[1,2]*A[2,1])
Tmat(s::PolynomialSystem,q) = inv2(Umat(s,q))*inv2(Wmat(s,q))*Umat(s,q)
Smat(s::PolynomialSystem,q) = inv2(Tmat(s,q))*inv2(Umat(s,0))

"Convert a paper boundary-action matrix to Kleinian's conjugate spinor convention."
boundary_isometry(A) = MobiusIsometry(conj.(A))
complexfloat(z) = ComplexF64(Float64(real(z)),Float64(imag(z)))

"""
    developed_domain(system, slope; parameter=geometric_root(system,slope))

Develop the tetrahedron pairs from Definitions 5.2 and 5.8. Exact algebraic
coordinates and shape parameters are retained; geometry is converted to Float64
for visualization. This is not yet a face-paired TetrahedralComplex certificate.
"""
function developed_domain(s::PolynomialSystem,q;parameter=geometric_root(s,q))
    r = unit_slope(q)
    (0 < r < 1//2 && r != 1//3) || throw(ArgumentError("domain requires 0 < slope < 1/2, excluding 1/3"))
    es = interior_edges(r)
    isempty(es) && throw(ArgumentError("no interior tetrahedron pairs"))
    shapes = [value(Z(s,e),parameter) for e in es]
    all(z -> imag(z)>0,shapes) || throw(ArgumentError("parameter gives a flat or inverted layer"))
    tetrahedra = HyperbolicTetrahedron[]
    triples = []
    for e in es
        hi,lo = highlow(e)
        v(q) = value(V(s,q),parameter)
        epsilon = nu(e)==1 ? -1 : 0
        pair = ((epsilon+v(opposite(e)), -1+v(e.qL), v(e.qR)),
                (v(e.qP),v(hi),v(lo)))
        for triple in pair
            push!(triples,triple)
            # Preserve the source vertex order; idealshape may be negative.
            push!(tetrahedra,HyperbolicTetrahedron([
                complexidealpoint(nothing),complexidealpoint.(complexfloat.(triple))...]))
        end
    end
    (; slope=r,parameter,edges=es,shapes,triples,tetrahedra,
       complete=iszero(polyeval(Qpoly(s,r),parameter)))
end
end
