# Hyperbolic 3-space consists of future timelike rays in R^(3,1).
_pointcoordinates(v::Multivector) = SVector(ntuple(i -> coefficient(v, 1 << (i-1)), 4))
_lorentz(a, b) = -a[1]*b[1] + a[2]*b[2] + a[3]*b[3] + a[4]*b[4]

# Predicate arithmetic must not overflow on fixed-width integers/rationals.
_predicate_value(x) = x
_predicate_value(x::Integer) = big(x)
_predicate_value(x::Rational) = big(x)
function _predicate_coordinates(v::SVector{4,T}) where {T}
    T <: AbstractFloat ? v / maximum(abs, v) : map(_predicate_value, v)
end

"""
    HyperbolicPoint(v::Multivector)
    HyperbolicPoint(t, x, y, z)
    HyperbolicPoint(coordinates)

A point of H³ represented by a future timelike grade-one vector in Cl(3,1).
The metric is (-,+,+,+), so t > 0 and -t²+x²+y²+z² < 0.
The supplied representative and its real coefficient type are preserved:
the constructor does not normalize. Positive multiples represent the same point.
Null, spacelike, past-directed, nonfinite, and mixed-grade inputs are rejected.
"""
struct HyperbolicPoint{T}
    vector::Multivector{T}
    function HyperbolicPoint{T}(v::Multivector{T}) where {T}
        all(isreal, coefficients(v)) || throw(ArgumentError("point coefficients must be real-valued"))
        all(coefficient_isfinite, coefficients(v)) || throw(ArgumentError("point coefficients must be finite"))
        v == grade(v, 1) || throw(ArgumentError("a point must be a grade-one vector"))
        c = map(real, _pointcoordinates(v))
        c[1] > zero(c[1]) || throw(ArgumentError("a point must be future directed (t > 0)"))
        w = _predicate_coordinates(c)
        _lorentz(w,w) < zero(w[1]) || throw(ArgumentError("a point must be strictly timelike"))
        new{T}(v)
    end
end

HyperbolicPoint(v::Multivector{T}) where {T} = HyperbolicPoint{T}(v)
HyperbolicPoint(t,x,y,z) = HyperbolicPoint(clvector(t,x,y,z))
HyperbolicPoint(c::Union{Tuple,AbstractVector}) = HyperbolicPoint(clvector(c))

"""Return the original, unnormalized grade-one multivector representing `p`."""
representative(p::HyperbolicPoint) = p.vector
"""Return the original homogeneous coordinates `(t,x,y,z)` as an SVector."""
coordinates(p::HyperbolicPoint) = _pointcoordinates(p.vector)

# Exact projective comparison for built-in integer, rational, and floating types.
# Converting finite floats to rationals preserves their represented binary value.
_equality_value(x) = x
_equality_value(x::Union{Integer,Rational,AbstractFloat}) = Rational{BigInt}(x)
_equality_value(x::Complex) = complex(_equality_value(real(x)), _equality_value(imag(x)))
function _projective_key(p::HyperbolicPoint)
    c = map(x -> _equality_value(real(x)), coordinates(p))
    ntuple(i -> c[i+1]/c[1], 3)
end
Base.:(==)(p::HyperbolicPoint,q::HyperbolicPoint) = _projective_key(p) == _projective_key(q)
Base.isequal(p::HyperbolicPoint,q::HyperbolicPoint) = isequal(_projective_key(p),_projective_key(q))
Base.hash(p::HyperbolicPoint,h::UInt) = hash(_projective_key(p),hash(:HyperbolicPoint,h))
function Base.isapprox(p::HyperbolicPoint,q::HyperbolicPoint; kwargs...)
    a, b = coordinates(p), coordinates(q)
    isapprox(a/a[1], b/b[1]; kwargs...)
end

"""Return a new representative of `p` scaled by a finite, strictly positive real `s`."""
function rescale(p::HyperbolicPoint,s)
    isreal(s) && coefficient_isfinite(s) && real(s) > zero(real(s)) || throw(ArgumentError("projective scale must be finite and positive"))
    HyperbolicPoint(s*representative(p))
end
Base.:*(s::Real,p::HyperbolicPoint) = rescale(p,s)
Base.:*(p::HyperbolicPoint,s::Real) = rescale(p,s)

_analytic_value(x) = x
_analytic_value(x::Union{Integer,Rational}) = float(x)
function _distance_coordinates(p::HyperbolicPoint)
    c = map(x -> _analytic_value(real(x)), coordinates(p))
    all(coefficient_isfinite,c) || throw(DomainError(c,"coordinate conversion overflowed; use BigFloat"))
    c = c/maximum(abs,c)
    _lorentz(c,c) < zero(c[1]) ||
        throw(DomainError(c,"timelike norm lost at working precision; use BigFloat"))
    c
end

"""
    unitvector(p::HyperbolicPoint)

Compute a future unit timelike multivector representing `p`, with square -1.
Floating arithmetic is used for the square root. The stored representative is unchanged.
"""
function unitvector(p::HyperbolicPoint)
    c = _distance_coordinates(p)
    clvector(c/sqrt(-_lorentz(c,c)))
end

"""
    distance(p::HyperbolicPoint, q::HyperbolicPoint)

Hyperbolic distance of projective points:
acosh(-⟨v,w⟩ / (sqrt(-⟨v,v⟩)*sqrt(-⟨w,w⟩))).
No stored representative is normalized. Working copies are rescaled for numerical
safety; tiny roundoff below the acosh domain boundary is clamped to 1.
"""
function distance(p::HyperbolicPoint,q::HyperbolicPoint)
    p == q && return zero(_analytic_value(real(coordinates(p)[1])))
    a,b = _distance_coordinates(p),_distance_coordinates(q)
    p == q && return zero(promote(a[1],b[1])[1])
    c = (-_lorentz(a,b)/sqrt(-_lorentz(a,a)))/sqrt(-_lorentz(b,b))
    if c < one(c)
        c isa AbstractFloat && one(c)-c <= 16eps(typeof(c)) ||
            throw(DomainError(c,"distance argument below 1; use higher precision"))
        c = one(c)
    end
    acosh(c)
end

mapcoefficients(f,p::HyperbolicPoint) = HyperbolicPoint(mapcoefficients(f,representative(p)))

function Base.show(io::IO,p::HyperbolicPoint)
    print(io,"HyperbolicPoint(")
    show(io,Tuple(coordinates(p)))
    print(io,")")
end
