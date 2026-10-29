"""An orientation-preserving hyperbolic isometry. g*h acts by g(h(x))."""
abstract type AbstractIsometry end

# Remove only irrelevant projective scale in floating working copies. Exact
# and foreign coefficient backends retain their original arithmetic.
_working_vector(v::Multivector) = v
function _working_vector(v::Multivector{T}) where {T<:AbstractFloat}
    s=maximum(abs,coefficients(v))
    iszero(s) ? v : v/s
end
function _working_matrix(A::SMatrix{2,2,T,4}) where {T}
    if T <: AbstractFloat || T <: Complex{<:AbstractFloat}
        s=maximum(abs,A)
        return iszero(s) ? A : A/s
    end
    A
end
_projective_mul(a::Multivector,b::Multivector) = _working_vector(a)*_working_vector(b)
_projective_mul2(A,B) = _mul2(_working_matrix(A),_working_matrix(B))

"""A half-turn about a geometric algebra geodesic axis."""
struct HalfTurn{G<:Geodesic} <: AbstractIsometry
    axis::G
end
HalfTurn(L::Multivector) = HalfTurn(Geodesic(L))
versor(h::HalfTurn) = dual(representative(h.axis))
Base.inv(h::HalfTurn) = h

"""
    VersorIsometry(R)

A real even Clifford versor with positive scalar reverse norm. R is stored
without normalization; nonzero real rescalings give the same isometry.
"""
struct VersorIsometry{T} <: AbstractIsometry
    value::Multivector{T}
    function VersorIsometry(R::Multivector{T}) where {T}
        all(isreal,coefficients(R)) && all(coefficient_isfinite,coefficients(R)) ||
            throw(ArgumentError("versor coefficients must be finite and real-valued"))
        R == gradeinvolution(R) || throw(ArgumentError("isometry versor must be even"))
        w = mapcoefficients(x -> _predicate_value(real(x)),R)
        if eltype(w) <: AbstractFloat && !iszero(w)
            w = w/maximum(abs,coefficients(w))
        end
        s = w*reverse(w)
        all(_negligible,coefficients(s)[2:16]) || throw(ArgumentError("versor reverse norm must be scalar"))
        real(scalarpart(s)) > zero(real(scalarpart(s))) || throw(ArgumentError("versor reverse norm must be positive"))
        new{T}(R)
    end
end
versor(g::VersorIsometry) = g.value
Base.inv(g::VersorIsometry) = VersorIsometry(reverse(g.value))

_det2(A) = A[1,1]*A[2,2]-A[1,2]*A[2,1]
_adjugate2(A) = SMatrix{2,2}(A[2,2],-A[2,1],-A[1,2],A[1,1])
# Explicit 2x2 kernels use only value arithmetic, including for foreign fields.
function _mul2(A,B)
    SMatrix{2,2}(A[1,1]*B[1,1]+A[1,2]*B[2,1],
                 A[2,1]*B[1,1]+A[2,2]*B[2,1],
                 A[1,1]*B[1,2]+A[1,2]*B[2,2],
                 A[2,1]*B[1,2]+A[2,2]*B[2,2])
end
_dagger2(A) = SMatrix{2,2}(conj(A[1,1]),conj(A[1,2]),conj(A[2,1]),conj(A[2,2]))

"""
    MobiusIsometry(A)

An invertible 2x2 complex matrix up to nonzero complex scale (PGL₂(C),
representing PSL₂(C)). No determinant normalization is performed.
The coefficient backend supplies conjugation, real/imaginary parts, and arithmetic.
"""
struct MobiusIsometry{T} <: AbstractIsometry
    value::SMatrix{2,2,T,4}
    function MobiusIsometry(A::AbstractMatrix)
        size(A) == (2,2) || throw(DimensionMismatch("expected a 2x2 matrix"))
        B = SMatrix{2,2}(A)
        all(coefficient_isfinite,B) || throw(ArgumentError("matrix coefficients must be finite"))
        iszero(_det2(_working_matrix(B))) && throw(ArgumentError("isometry matrix must be invertible"))
        new{eltype(B)}(B)
    end
end
matrix(g::MobiusIsometry;imaginary_unit=nothing) = g.value
Base.inv(g::MobiusIsometry) = MobiusIsometry(_adjugate2(g.value))

"""Ordered factors; composition retains its construction instead of refactorizing."""
struct IsometryProduct{F<:Tuple} <: AbstractIsometry
    factors::F
    function IsometryProduct(factors::Tuple)
        !isempty(factors) && all(g -> g isa AbstractIsometry,factors) ||
            throw(ArgumentError("a product requires at least one isometry factor"))
        new{typeof(factors)}(factors)
    end
end
_factors(g::AbstractIsometry) = (g,)
_factors(g::IsometryProduct) = g.factors
Base.:*(g::AbstractIsometry,h::AbstractIsometry) = IsometryProduct((_factors(g)...,_factors(h)...))
Base.:∘(g::AbstractIsometry,h::AbstractIsometry) = g*h
Base.inv(g::IsometryProduct) = IsometryProduct(map(inv,reverse(g.factors)))
function versor(g::IsometryProduct)
    reduce(_projective_mul,map(versor,g.factors))
end
function matrix(g::IsometryProduct;imaginary_unit=nothing)
    reduce(_projective_mul2,map(f -> matrix(f;imaginary_unit),g.factors))
end

"""
    Translation(axis, length)
    Translation(axis; multiplier=q)

Signed translation along an oriented GA geodesic. The multiplier form uses
q>0 with signed length log(q) and avoids evaluating a logarithm or exponential.
Parameters are retained and evaluated only by versor/matrix/action.
"""
struct Translation{G<:Geodesic,S} <: AbstractIsometry
    axis::G
    parameter::S
    by_multiplier::Bool
    function Translation(axis::Geodesic,s,by_multiplier::Bool=false)
        isreal(s) && coefficient_isfinite(s) || throw(ArgumentError("translation parameter must be finite and real-valued"))
        by_multiplier && !(real(s)>zero(real(s))) && throw(ArgumentError("multiplier must be positive"))
        new{typeof(axis),typeof(s)}(axis,s,by_multiplier)
    end
end
Translation(axis::Geodesic;multiplier) = Translation(axis,multiplier,true)
function versor(g::Translation)
    B = representative(g.axis)
    b = _sqrt_preserving(real(scalarpart(B*B)))
    if g.by_multiplier
        q = g.parameter
        return scalar((q+one(q))*b)+(q-one(q))*B
    end
    t = _fielddiv(g.parameter,_two(g.parameter))
    scalar(cosh(t)) + _fielddiv(sinh(t),b)*B
end
Base.inv(g::Translation) = g.by_multiplier ?
    Translation(reverse(g.axis);multiplier=g.parameter) : Translation(g.axis,-g.parameter)

"""
    Elliptic(axis, order; numerator=1)

Rotation about an oriented GA axis by 2π*numerator/order. The integers are
retained exactly. They must be coprime for the specified order to be the order.
Evaluation uses the coefficient backend's sinpi/cospi (and sqrt for the axis).
"""
struct Elliptic{G<:Geodesic} <: AbstractIsometry
    axis::G
    order::Int
    numerator::Int
    function Elliptic(axis::Geodesic,n::Integer;numerator::Integer=1)
        n>0 && gcd(n,numerator)==1 || throw(ArgumentError("order must be positive and numerator coprime to it"))
        new{typeof(axis)}(axis,Int(n),Int(mod(numerator,n)))
    end
end
function versor(g::Elliptic)
    B = representative(g.axis)
    u = one(scalarpart(B))
    angle = _fielddiv(u*g.numerator,u*g.order)
    b = _sqrt_preserving(real(scalarpart(B*B)))
    scalar(cospi(angle)) + _fielddiv(sinpi(angle),b)*dual(B)
end
Base.inv(g::Elliptic) = Elliptic(g.axis,g.order;numerator=-g.numerator)

"""
    Parabolic(horosphere, u)

Null rotation with N=h∧u and R=1+N/2, where h is the horosphere's future
null vector, u⋅h=0, and N≠0. u modulo multiples of h is the displacement
parameter. Its sign follows this rotor convention. The supplied horosphere
is preserved; replacing h by s*h requires u/s for the same transformation.
"""
struct Parabolic{H<:Horosphere,T} <: AbstractIsometry
    horosphere::H
    direction::Multivector{T}
    function Parabolic(h::Horosphere,u::Multivector{T}) where {T}
        w = _geometric_input(u,1)
        hv = _geometric_input(representative(h),1)
        _negligible(_lorentz(_pointcoordinates(w),_pointcoordinates(hv))) ||
            throw(ArgumentError("parabolic direction must be orthogonal to the null centre"))
        iszero(representative(h)∧u) && throw(ArgumentError("parabolic direction must not be parallel to its centre"))
        new{typeof(h),T}(h,u)
    end
end
function versor(g::Parabolic)
    N = representative(g.horosphere)∧g.direction
    one(N)+_divide(N,_two(scalarpart(N)))
end
Base.inv(g::Parabolic) = Parabolic(g.horosphere,-g.direction)

# Cl⁺(3,1) ≅ Mat₂(C): e0ei maps to the Pauli matrix σi, I maps to i.
function matrix(g::AbstractIsometry;imaginary_unit=nothing)
    R = _working_vector(versor(g))
    i = isnothing(imaginary_unit) ? coefficient_imaginary_unit(scalarpart(R)) : imaginary_unit
    iszero(i*i+one(i)) || throw(ArgumentError("imaginary_unit must square to -1"))
    c0 = coefficient(R,0)+i*coefficient(R,15)
    c1 = coefficient(R,3)+i*coefficient(R,12)
    c2 = coefficient(R,5)-i*coefficient(R,10)
    c3 = coefficient(R,9)+i*coefficient(R,6)
    SMatrix{2,2}(c0+c3,c1+i*c2,c1-i*c2,c0-c3)
end

"""Compute an SL₂(C) lift explicitly; the backend's sqrt chooses a sign/branch."""
function sl2matrix(g::AbstractIsometry;imaginary_unit=nothing)
    A = _working_matrix(matrix(g;imaginary_unit))
    determinant = _det2(A)
    d = isreal(determinant) && real(determinant)<zero(real(determinant)) ?
        coefficient_imaginary_unit(determinant)*_sqrt_preserving(-real(determinant)) :
        _sqrt_preserving(determinant)
    map(x -> _fielddiv(x,d),A)
end
function versor(g::MobiusIsometry)
    A = sl2matrix(g)
    a,b,c,d = A[1,1],A[1,2],A[2,1],A[2,2]
    half(x) = _fielddiv(x,_two(x))
    z = zero(real(a))
    coeff = (half(real(a)+real(d)),z,z,half(real(b)+real(c)),z,
             half(imag(c)-imag(b)),half(imag(a)-imag(d)),z,z,
             half(real(a)-real(d)),half(real(c)-real(b)),z,
             half(imag(b)+imag(c)),z,z,half(imag(a)+imag(d)))
    Multivector(coeff)
end

function _linear_action(g::AbstractIsometry,v::Multivector;normalized=false)
    R = _working_vector(versor(g))
    source = normalized ? v : _working_vector(v)
    w = grade(R*source*reverse(R),1)
    normalized ? _divide(w,real(scalarpart(R*reverse(R)))) : w
end
function _linear_action(g::MobiusIsometry,v::Multivector;normalized=false)
    # Construct H in a compatible complex field via the coefficient protocol.
    A = _working_matrix(g.value)
    t,x,y,z = _pointcoordinates(normalized ? v : _working_vector(v))
    H = SMatrix{2,2}(t+z,coefficient_complex(x,y),coefficient_complex(x,-y),t-z)
    W = _mul2(_mul2(A,H),_dagger2(A))
    result = clvector(real(W[1,1])+real(W[2,2]),real(W[1,2])+real(W[2,1]),
                      imag(W[2,1])-imag(W[1,2]),real(W[1,1])-real(W[2,2]))
    # result is twice the Hermitian-coordinate vector; projective objects
    # can retain that scale, but a horosphere must divide by 2*|det A|.
    normalized || return result
    d = _det2(A)
    magnitude = isreal(d) ? abs(real(d)) : _sqrt_preserving(real(d*conj(d)))
    _divide(result,_two(magnitude)*magnitude)
end
(g::AbstractIsometry)(p::HyperbolicPoint) = HyperbolicPoint(_linear_action(g,representative(p)))
(g::AbstractIsometry)(p::IdealPoint) = IdealPoint(_linear_action(g,representative(p)))
(g::AbstractIsometry)(p::HyperbolicPlane) = HyperbolicPlane(_linear_action(g,representative(p)))
(g::AbstractIsometry)(h::Horosphere) = Horosphere(_linear_action(g,representative(h);normalized=true))
function (g::AbstractIsometry)(L::Geodesic)
    R=_working_vector(versor(g))
    Geodesic(grade(R*_working_vector(L.blade)*reverse(R),2))
end
function (g::MobiusIsometry)(L::Geodesic)
    sample = scalarpart(L.blade)
    e = ntuple(i -> _linear_action(g,basisvector(i-1,one(sample))),4)
    result = zero(e[1]∧e[2])
    for i in 0:3, j in i+1:3
        result += coefficient(L.blade,(1<<i)|(1<<j))*(e[i+1]∧e[j+1])
    end
    Geodesic(result)
end
for G in (:HyperbolicPoint,:IdealPoint,:Geodesic,:HyperbolicPlane,:Horosphere)
    @eval function (g::IsometryProduct)(x::$G)
        foldr((f,y) -> f(y),g.factors;init=x)
    end
end

Base.copy(g::AbstractIsometry) = g
Base.one(g::AbstractIsometry) = VersorIsometry(one(versor(g)))
Base.one(g::MobiusIsometry) = MobiusIsometry(SMatrix{2,2}(one(g.value[1]),zero(g.value[1]),zero(g.value[1]),one(g.value[1])))
Base.one(g::IsometryProduct) = one(first(g.factors))

# Keep powers as construction data; do not allocate a word of length n.
struct IsometryPower{G<:AbstractIsometry,N<:Integer} <: AbstractIsometry
    base::G
    exponent::N
end
function Base.:^(g::AbstractIsometry,n::Integer)
    iszero(n) && return one(g)
    n == 1 && return g
    n < 0 ? IsometryPower(inv(g),-big(n)) : IsometryPower(g,n)
end
function _power_value(base,n,multiply,result)
    while n>0
        isodd(n) && (result=multiply(result,base))
        n >>= 1
        n>0 && (base=multiply(base,base))
    end
    result
end
function versor(g::IsometryPower)
    R=versor(g.base)
    _power_value(R,g.exponent,_projective_mul,one(R))
end
function matrix(g::IsometryPower;imaginary_unit=nothing)
    A=matrix(g.base;imaginary_unit)
    u,z=one(A[1]),zero(A[1])
    _power_value(A,g.exponent,_projective_mul2,SMatrix{2,2}(u,z,z,u))
end
Base.inv(g::IsometryPower) = inv(g.base)^g.exponent
Base.one(g::IsometryPower) = one(g.base)
_uses_matrix(::AbstractIsometry) = false
_uses_matrix(::MobiusIsometry) = true
_uses_matrix(g::IsometryProduct) = any(_uses_matrix,g.factors)
_uses_matrix(g::IsometryPower) = _uses_matrix(g.base)
for G in (:HyperbolicPoint,:IdealPoint,:Geodesic,:HyperbolicPlane,:Horosphere)
    @eval function (g::IsometryPower)(x::$G)
        evaluated = _uses_matrix(g.base) ? MobiusIsometry(matrix(g)) : VersorIsometry(versor(g))
        evaluated(x)
    end
end

function _matrix_key(g::AbstractIsometry)
    c = map(_equality_value,matrix(g))
    pivot = first(filter(x -> !iszero(x),c))
    Tuple(map(x -> _fielddiv(x,pivot),c))
end
Base.:(==)(g::AbstractIsometry,h::AbstractIsometry) = _matrix_key(g)==_matrix_key(h)
Base.isequal(g::AbstractIsometry,h::AbstractIsometry) = isequal(_matrix_key(g),_matrix_key(h))
Base.hash(g::AbstractIsometry,h::UInt) = hash(_matrix_key(g),hash(:HyperbolicIsometry,h))
function Base.isapprox(g::AbstractIsometry,h::AbstractIsometry;kwargs...)
    A,B = matrix(g),matrix(h)
    # Normalize by the largest entry of B rather than an arbitrarily tiny pivot.
    k = argmax(abs.(B))
    isapprox(map(x -> x*B[k],A),map(x -> x*A[k],B);kwargs...)
end
function isidentity(g::AbstractIsometry;atol=0)
    A = matrix(g)
    all(x -> _incident_zero(x,atol),(A[1,2],A[2,1],A[1,1]-A[2,2]))
end

mapcoefficients(f,g::HalfTurn) = HalfTurn(mapcoefficients(f,g.axis))
mapcoefficients(f,g::VersorIsometry) = VersorIsometry(mapcoefficients(f,g.value))
mapcoefficients(f,g::MobiusIsometry) = MobiusIsometry(map(f,g.value))
mapcoefficients(f,g::Translation) = Translation(mapcoefficients(f,g.axis),f(g.parameter),g.by_multiplier)
mapcoefficients(f,g::Elliptic) = Elliptic(mapcoefficients(f,g.axis),g.order;numerator=g.numerator)
mapcoefficients(f,g::Parabolic) = Parabolic(mapcoefficients(f,g.horosphere),mapcoefficients(f,g.direction))
mapcoefficients(f,g::IsometryProduct) = IsometryProduct(map(x -> mapcoefficients(f,x),g.factors))
mapcoefficients(f,g::IsometryPower) = mapcoefficients(f,g.base)^g.exponent
