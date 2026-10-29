"""
    Multivector(coefficients)
    Multivector{T}(coefficients)

An immutable Cl(3,1) multivector backed by `SVector{16,T}`. Coefficients
are in ascending bit-mask order (mask 0 is the scalar; masks 1,2,4,8 are
e0,e1,e2,e3). `T` is unconstrained; arithmetic assumes commuting coefficients.
"""
struct Multivector{T}
    coeffs::SVector{16,T}

    Multivector{T}(xs::SVector{16,T}) where {T} = new{T}(xs)
end

Multivector(xs::SVector{16,T}) where {T} = Multivector{T}(xs)
Multivector(xs::Union{Tuple,AbstractVector}) = Multivector(SVector{16}(xs))
Multivector{T}(xs::Union{Tuple,AbstractVector}) where {T} =
    Multivector{T}(SVector{16,T}(xs))

"""Return the immutable, fixed-size coefficient vector (ordinary 1-based indexing)."""
coefficients(a::Multivector) = a.coeffs
"""Map coefficients explicitly, preserving their backend's arithmetic and parent."""
mapcoefficients(f, a::Multivector) = Multivector(map(f, a.coeffs))
Base.eltype(::Type{Multivector{T}}) where {T} = T
Base.eltype(::Multivector{T}) where {T} = T

function _checkmask(mask::Integer)
    0 <= mask <= 15 || throw(ArgumentError("blade mask must be in 0:15"))
    Int(mask)
end

"""`coefficient(a, mask)` returns the coefficient for a blade mask in `0:15`."""
coefficient(a::Multivector, mask::Integer) = a.coeffs[_checkmask(mask) + 1]
scalarpart(a::Multivector) = a.coeffs[1]

"""Construct a single basis blade with the given coefficient (default `1`)."""
function basisblade(mask::Integer, value=1)
    m = _checkmask(mask)
    Multivector(SVector{16}(ntuple(i -> i == m + 1 ? value : zero(value), Val(16))))
end

"""Construct the scalar multivector corresponding to `value`."""
scalar(value) = basisblade(0, value)

"""Construct `value*ei`, where `i` is in `0:3` and `e0^2 = -1`."""
function basisvector(i::Integer, value=1)
    0 <= i <= 3 || throw(ArgumentError("basis vector index must be in 0:3"))
    basisblade(1 << Int(i), value)
end

"""Construct the grade-one multivector with coordinates `(t,x,y,z)`."""
function clvector(t, x, y, z)
    v = promote(t, x, y, z)
    Multivector(SVector{16}(ntuple(Val(16)) do i
        m = i - 1
        ispow2(m) ? v[trailing_zeros(m) + 1] : zero(v[1])
    end))
end
function clvector(v::Union{Tuple,AbstractVector})
    length(v) == 4 || throw(DimensionMismatch("expected four vector coordinates"))
    clvector(v...)
end

Base.zero(::Type{Multivector{T}}) where {T} = scalar(zero(T))
Base.one(::Type{Multivector{T}}) where {T} = scalar(one(T))
Base.zero(a::Multivector) = scalar(zero(scalarpart(a)))
Base.one(a::Multivector) = scalar(one(scalarpart(a)))
Base.copy(a::Multivector) = a
Base.iszero(a::Multivector) = all(iszero, a.coeffs)
Base.isone(a::Multivector) = isone(scalarpart(a)) && all(iszero, a.coeffs[2:16])
Base.:(==)(a::Multivector, b::Multivector) = a.coeffs == b.coeffs
Base.isequal(a::Multivector, b::Multivector) = isequal(a.coeffs, b.coeffs)
Base.hash(a::Multivector, h::UInt) = hash(a.coeffs, hash(:Cl31Multivector, h))
Base.isapprox(a::Multivector, b::Multivector; kwargs...) =
    isapprox(a.coeffs, b.coeffs; kwargs...)

Base.:+(a::Multivector, b::Multivector) = Multivector(a.coeffs + b.coeffs)
Base.:-(a::Multivector, b::Multivector) = Multivector(a.coeffs - b.coeffs)
Base.:-(a::Multivector) = Multivector(-a.coeffs)
Base.:+(a::Multivector) = a
Base.:*(s, a::Multivector) = Multivector(map(x -> s*x, a.coeffs))
Base.:*(a::Multivector, s) = Multivector(map(x -> x*s, a.coeffs))
Base.:/(a::Multivector, s) = Multivector(map(x -> x/s, a.coeffs))
Base.:+(a::Multivector, s::Number) = a + scalar(s)
Base.:+(s::Number, a::Multivector) = scalar(s) + a
Base.:-(a::Multivector, s::Number) = a - scalar(s)
Base.:-(s::Number, a::Multivector) = scalar(s) - a

# Moving each generator of b left across the higher generators in a gives
# the permutation sign; the shared e0 contributes the negative metric sign.
function _bladesign(a::Int, b::Int)
    swaps = sum(count_ones(a >> i) for i in 1:4 if !iszero(b & (1 << (i-1))); init=0)
    isodd(swaps + (!iszero(a & b & 1))) ? -1 : 1
end
const _PRODUCT_SIGNS = ntuple(i -> ntuple(j -> _bladesign(i-1, j-1), 16), 16)

# For output mask k and left mask a, exactly one right mask b=a xor k
# contributes. Immutable accumulators also work for non-isbits coefficients.
@inline _keep(::Val{:geometric}, a, b) = true
@inline _keep(::Val{:wedge}, a, b) = iszero(a & b)
@inline _keep(::Val{:leftcontract}, a, b) =
    count_ones(a) <= count_ones(b) && count_ones(xor(a,b)) == count_ones(b)-count_ones(a)

function _product(a::Multivector, b::Multivector, kind)
    Multivector(SVector{16}(ntuple(Val(16)) do i
        k = i - 1
        value = zero(a.coeffs[1] * b.coeffs[1])
        for amask in 0:15
            bmask = xor(amask, k)
            if _keep(kind, amask, bmask)
                term = a.coeffs[amask+1] * b.coeffs[bmask+1]
                value = _PRODUCT_SIGNS[amask+1][bmask+1] == 1 ? value + term : value - term
            end
        end
        value
    end))
end

"""Geometric product in Cl(3,1); also available as `a*b`."""
geometricproduct(a::Multivector, b::Multivector) = _product(a, b, Val(:geometric))
Base.:*(a::Multivector, b::Multivector) = geometricproduct(a, b)

"""
    dot(a::Multivector, b::Multivector)
    a ⋅ b

Symmetrized geometric product `(a*b + b*a)/2`, returning a multivector.
For vectors this is their scalar Lorentz product. For general multivectors
the result may contain other grades. Coefficients are not conjugated;
this operation is symmetric, including for complex coefficients.
Division follows the coefficient type's ordinary `/` rules.
"""
dot(a::Multivector, b::Multivector) = (a*b + b*a)/2

"""Outer product, extended bilinearly across grades; also available as `a ∧ b`."""
wedge(a::Multivector, b::Multivector) = _product(a, b, Val(:wedge))
const ∧ = wedge

"""
    leftcontract(a, b)

For homogeneous grades r,s, retain grade s-r of a*b when r <= s;
otherwise return zero. Extend bilinearly to mixed grades.
"""
leftcontract(a::Multivector, b::Multivector) = _product(a, b, Val(:leftcontract))

"""Project onto grade `r` in `0:4`, retaining multivector storage."""
function grade(a::Multivector, r::Integer)
    0 <= r <= 4 || throw(ArgumentError("grade must be in 0:4"))
    Multivector(map((m,x) -> count_ones(m) == r ? x : zero(x), SVector{16}(0:15), a.coeffs))
end

function _gradesigns(a::Multivector, sign)
    Multivector(map((m,x) -> sign(count_ones(m)) ? -x : x, SVector{16}(0:15), a.coeffs))
end

"""Reverse basis-generator order, with grade-r sign (-1)^(r(r-1)/2)."""
Base.reverse(a::Multivector) = _gradesigns(a, r -> isodd(r*(r-1)÷2))
"""Negate odd grades."""
gradeinvolution(a::Multivector) = _gradesigns(a, isodd)
"""Compose reversion and grade involution; coefficient conjugation is not applied."""
cliffordconjugate(a::Multivector) = gradeinvolution(reverse(a))

"""Right dual `a*I^-1`, with `I=e0*e1*e2*e3` and `I^-1=-I`."""
dual(a::Multivector) = a * basisblade(15, -one(scalarpart(a)))
"""Inverse right dual `a*I`; `undual(dual(a)) == a`."""
undual(a::Multivector) = a * basisblade(15, one(scalarpart(a)))

"""
    checked_inverse(a)
    inv(a)

Return reverse(a)/s only if both a*reverse(a) and reverse(a)*a are exactly
the same nonzero scalar s. Covers non-null vectors and suitable versors.
This is a restricted inverse: other invertible multivectors are rejected.
No tolerance silently discards nonscalar residuals in floating-point inputs.
"""
function checked_inverse(a::Multivector)
    r = reverse(a)
    ar, ra = a*r, r*a
    all(iszero, coefficients(ar)[2:16]) && all(iszero, coefficients(ra)[2:16]) &&
        scalarpart(ar) == scalarpart(ra) ||
        throw(ArgumentError("reverse products must be the same scalar; general inverse is not implemented"))
    s = scalarpart(ar)
    iszero(s) && throw(DomainError(a, "zero reverse norm: no inverse by this formula"))
    Multivector(map(x -> x/s, coefficients(r)))
end
Base.inv(a::Multivector) = checked_inverse(a)
Base.:^(a::Multivector, n::Integer) = n < 0 ?
    Base.power_by_squaring(inv(a), -big(n)) : Base.power_by_squaring(a, n)

# Grade order followed by lexicographic basis-index order; storage is unchanged.
const _DISPLAY_MASKS = (0, 1, 2, 4, 8, 3, 5, 9, 6, 10, 12, 7, 11, 13, 14, 15)

# Unordered coefficient rings keep their signs inside the parentheses.
_negative_coefficient(c) = false
_negative_coefficient(c::Real) = (c < zero(c)) === true

function Base.show(io::IO, a::Multivector)
    firstterm = true
    for mask in _DISPLAY_MASKS
        c = coefficient(a, mask)
        iszero(c) === true && continue
        negative = _negative_coefficient(c)
        print(io, firstterm ? (negative ? "-" : "") : (negative ? " - " : " + "))
        print(io, '(')
        # Widen fixed-width integers before negation, including typemin(Int).
        magnitude = negative ? (c isa Integer && isbitstype(typeof(c)) ? -big(c) : -c) : c
        show(io, magnitude)
        print(io, ')')
        iszero(mask) || print(io, '*')
        for i in 0:3
            iszero(mask & (1 << i)) || print(io, "e$i")
        end
        firstterm = false
    end
    firstterm && print(io, '0')
end

Base.show(io::IO, ::MIME"text/plain", a::Multivector) = show(io, a)
