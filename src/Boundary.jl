# Projective working copies preserve exact coefficient fields and put numerical
# residuals on the t=1 chart. No stored vertex is normalized or converted.
function _chartvector(p::Union{HyperbolicPoint,IdealPoint})
    c = map(x -> _predicate_value(real(x)),coordinates(p))
    map(x -> _fielddiv(x,c[1]),c)
end

# Homogeneous complex boundary coordinate z=x+iy. This is the conjugate of
# the spinor convention used by matrix(g); see docs/isometries.md.
function _boundarypair(p::IdealPoint)
    t,x,y,z = _chartvector(p)
    abs(t+z) > abs(t-z) ? (t+z,coefficient_complex(x,-y)) : (coefficient_complex(x,y),t-z)
end
_detpair(a,b) = a[1]*b[2]-a[2]*b[1]
"""
    crossratio(p, q, r, s)

Complex ideal-point cross-ratio in the boundary coordinate x+iy, with
crossratio(infinity,0,1,z)=z. Homogeneous evaluation includes infinity.
"""
function crossratio(p::IdealPoint,q::IdealPoint,r::IdealPoint,s::IdealPoint)
    a,b,c,d = map(_boundarypair,(p,q,r,s))
    denominator = _detpair(a,d)*_detpair(b,c)
    iszero(denominator) && throw(ArgumentError("cross-ratio denominator is zero"))
    _fielddiv(_detpair(a,c)*_detpair(b,d),denominator)
end

"""
    complexcoordinate(p::IdealPoint)

Boundary coordinate `(x+im*y)/(t-z)` in the upper half-space chart.
Returns `nothing` at infinity; finite values retain the coefficient backend.
Uses a second homogeneous chart near the north pole to avoid cancellation.
This is the conjugate of the spinor coordinate for `matrix(g)`: the action
on this coordinate uses `conj.(matrix(g))` in the usual Möbius formula.
"""
function complexcoordinate(p::IdealPoint)
    a,b = _boundarypair(p)
    iszero(b) ? nothing : _fielddiv(a,b)
end

"""
    complexidealpoint(z)
    complexidealpoint(nothing)

Lift a finite complex boundary coordinate to the existing null-ray
`IdealPoint` representation. `nothing` denotes the single point at infinity.
Exact inputs retain exact arithmetic; floating inputs are scaled before
squaring to avoid overflow for large complex coordinates.
"""
complexidealpoint(::Nothing) = IdealPoint(1,0,0,1)
function complexidealpoint(z)
    coefficient_isfinite(z) || throw(ArgumentError("use nothing for boundary infinity; coordinates must be finite"))
    x,y = real(z),imag(z)
    if x isa AbstractFloat
        s = max(one(x),abs(x),abs(y))
        u,v,w = x/s,y/s,inv(s)
        r2 = u*u+v*v
        return IdealPoint(r2+w*w,2u*w,2v*w,r2-w*w)
    end
    r2 = x*x+y*y
    IdealPoint(r2+one(r2),2x,2y,r2-one(r2))
end

