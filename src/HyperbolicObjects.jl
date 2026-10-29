# Value-based field operations: no assumption that an element subtypes Number,
# or that its parent field can be recovered from its Julia type.
_two(x) = one(x) + one(x)
_fielddiv(x,y) = x/y
_fielddiv(x::Integer,y::Integer) = big(x)//big(y)
_divide(a::Multivector,s) = mapcoefficients(x -> _fielddiv(x,s),a)
_sqrt_preserving(x) = sqrt(x)
function _sqrt_preserving(x::Integer)
    r = isqrt(x)
    r*r == x ? r : sqrt(x)
end
function _sqrt_preserving(x::Rational)
    n,d = isqrt(numerator(x)),isqrt(denominator(x))
    n*n == numerator(x) && d*d == denominator(x) ? n//d : sqrt(x)
end
function _sqrt_preserving(x::Complex)
    if iszero(imag(x))
        a = real(x)
        return a < zero(a) ? complex(zero(a),_sqrt_preserving(-a)) : _sqrt_preserving(a)
    end
    sqrt(x)
end

# Roundoff allowance is local to validation of built-in floating coefficients.
# Exact and foreign coefficient backends use their exact predicates.
_negligible(x) = iszero(x)
_negligible(x::AbstractFloat) = abs(x) <= 64eps(typeof(x))
function _geometric_input(v::Multivector,r::Integer)
    all(isreal,coefficients(v)) || throw(ArgumentError("geometry needs real-valued coefficients"))
    all(coefficient_isfinite,coefficients(v)) || throw(ArgumentError("geometry needs finite coefficients"))
    v == grade(v,r) || throw(ArgumentError("expected a homogeneous grade-$r multivector"))
    w = mapcoefficients(x -> _predicate_value(real(x)),v)
    if eltype(w) <: AbstractFloat && !iszero(w)
        w = w/maximum(abs,coefficients(w))
    end
    w
end
function _square_scalar(v::Multivector)
    vv = v*v
    all(_negligible,coefficients(vv)[2:16]) || throw(ArgumentError("blade square must be scalar"))
    real(scalarpart(vv))
end

"""A future null ray, stored as an unnormalized grade-one vector."""
struct IdealPoint{T}
    vector::Multivector{T}
    function IdealPoint(v::Multivector{T}) where {T}
        w = _geometric_input(v,1)
        real(coefficient(w,1)) > zero(real(coefficient(w,1))) ||
            throw(ArgumentError("ideal point must be future directed"))
        _negligible(_square_scalar(w)) || throw(ArgumentError("ideal point must be null"))
        new{T}(v)
    end
end
IdealPoint(t,x,y,z) = IdealPoint(clvector(t,x,y,z))
IdealPoint(c::Union{Tuple,AbstractVector}) = IdealPoint(clvector(c))
representative(p::IdealPoint) = p.vector
coordinates(p::IdealPoint) = _pointcoordinates(p.vector)
function _projective_key(p::IdealPoint)
    c = map(x -> _equality_value(real(x)),coordinates(p))
    ntuple(i -> c[i+1]/c[1],3)
end
Base.:(==)(p::IdealPoint,q::IdealPoint) = _projective_key(p) == _projective_key(q)
Base.isequal(p::IdealPoint,q::IdealPoint) = isequal(_projective_key(p),_projective_key(q))
Base.hash(p::IdealPoint,h::UInt) = hash(_projective_key(p),hash(:IdealPoint,h))
function Base.isapprox(p::IdealPoint,q::IdealPoint;kwargs...)
    a,b = coordinates(p),coordinates(q)
    isapprox(a/a[1],b/b[1];kwargs...)
end

"""
    Geodesic(L)
    Geodesic(p, q)

An oriented hyperbolic geodesic, represented by a nonzero simple timelike
bivector L (L² > 0). For two finite/ideal points, L=p∧q; reversing their
order reverses the axis orientation. Positive scale does not change orientation.
"""
struct Geodesic{T}
    blade::Multivector{T}
    function Geodesic(L::Multivector{T}) where {T}
        w = _geometric_input(L,2)
        all(_negligible,coefficients(w∧w)) || throw(ArgumentError("geodesic bivector must be simple"))
        _square_scalar(w) > zero(real(scalarpart(w))) || throw(ArgumentError("geodesic must be timelike and nonzero"))
        new{T}(L)
    end
end
Geodesic(p::Union{HyperbolicPoint,IdealPoint},q::Union{HyperbolicPoint,IdealPoint}) =
    Geodesic(representative(p)∧representative(q))
representative(g::Geodesic) = g.blade
Base.reverse(g::Geodesic) = Geodesic(-g.blade)

"""An oriented plane specified by an unnormalized spacelike normal vector."""
struct HyperbolicPlane{T}
    normal::Multivector{T}
    function HyperbolicPlane(n::Multivector{T}) where {T}
        w = _geometric_input(n,1)
        _square_scalar(w) > zero(real(scalarpart(w))) || throw(ArgumentError("plane normal must be spacelike"))
        new{T}(n)
    end
end
representative(p::HyperbolicPlane) = p.normal

"""
    Horosphere(h)

A horosphere defined by the future null vector h, with its scale retained.
For a unit point x its equation is -x⋅h=1. For a projective point v,
incidence is (v⋅h)²=-v². Rescaling h changes the horosphere.
"""
struct Horosphere{T}
    normal::Multivector{T}
    function Horosphere(h::Multivector{T}) where {T}
        IdealPoint(h) # Validate without discarding scale.
        new{T}(h)
    end
end
Horosphere(c::IdealPoint) = Horosphere(representative(c))
representative(h::Horosphere) = h.normal
idealcenter(h::Horosphere) = IdealPoint(h.normal)
Base.:(==)(h::Horosphere,k::Horosphere) = h.normal == k.normal
Base.isequal(h::Horosphere,k::Horosphere) = isequal(h.normal,k.normal)
Base.hash(h::Horosphere,s::UInt) = hash(h.normal,hash(:Horosphere,s))
Base.isapprox(h::Horosphere,k::Horosphere;kwargs...) = isapprox(h.normal,k.normal;kwargs...)

function _blade_key(v::Multivector;oriented=false)
    c = map(x -> _equality_value(real(x)),coefficients(v))
    pivot = first(filter(x -> !iszero(x),c))
    sign = oriented && pivot < zero(pivot) ? -1 : 1
    (sign,Tuple(map(x -> _fielddiv(x,pivot),c)))
end
for G in (:Geodesic,:HyperbolicPlane)
    @eval begin
        Base.:(==)(a::$G,b::$G) = _blade_key(representative(a);oriented=true) == _blade_key(representative(b);oriented=true)
        Base.isequal(a::$G,b::$G) = isequal(_blade_key(representative(a);oriented=true),_blade_key(representative(b);oriented=true))
        Base.hash(a::$G,h::UInt) = hash(_blade_key(representative(a);oriented=true),hash($(QuoteNode(G)),h))
    end
end
"""Compare geodesic axes without regard to orientation."""
sameaxis(a::Geodesic,b::Geodesic) = _blade_key(a.blade) == _blade_key(b.blade)

_incident_zero(x,atol) = iszero(x) || (!iszero(atol) && abs(x) <= atol)
"""Incidence in homogeneous coordinates; atol defaults to zero (exact incidence)."""
function isincident(p::Union{HyperbolicPoint,IdealPoint},g::Geodesic;atol=0)
    all(x -> _incident_zero(x,atol),coefficients(representative(p)∧g.blade))
end
function isincident(p::Union{HyperbolicPoint,IdealPoint},g::HyperbolicPlane;atol=0)
    a,b = map(_predicate_value,coordinates(p)),map(_predicate_value,_pointcoordinates(g.normal))
    _incident_zero(_lorentz(a,b),atol)
end
function isincident(p::HyperbolicPoint,h::Horosphere;atol=0)
    a,b = map(_predicate_value,coordinates(p)),map(_predicate_value,_pointcoordinates(h.normal))
    _incident_zero(_lorentz(a,b)^2+_lorentz(a,a),atol)
end
Base.in(p::Union{HyperbolicPoint,IdealPoint},g::Union{Geodesic,HyperbolicPlane}) = isincident(p,g)
Base.in(p::HyperbolicPoint,h::Horosphere) = isincident(p,h)

for G in (:IdealPoint,:Geodesic,:HyperbolicPlane,:Horosphere)
    @eval begin
        mapcoefficients(f,g::$G) = $G(mapcoefficients(f,representative(g)))
        function Base.show(io::IO,g::$G)
            print(io,$(string(G)),"(")
            show(io,representative(g))
            print(io,")")
        end
    end
end
