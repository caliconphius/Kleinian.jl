# Chart convention: e3 is vertical in upper half-space; (0,0,1) is infinity.
function _checkmodel(model)
    model in (:klein, :ball, :halfspace) ||
        throw(ArgumentError("model must be :klein, :ball, or :halfspace"))
    model
end

"""
    modelcoordinates(p; model=:ball)

Map a finite or ideal point to three coordinates in `:klein`, `:ball`
(Poincaré), or `:halfspace` (upper half-space, third coordinate vertical).
Ideal points have unit radius in the balls and height zero in half-space.
The ideal point `(1,0,0,1)` is excluded from the half-space chart and throws
`DomainError`. Computation uses floating arithmetic without changing `p`.
"""
function modelcoordinates(p::Union{HyperbolicPoint,IdealPoint}; model=:ball)
    _checkmodel(model)
    c = map(x -> _analytic_value(real(x)), coordinates(p))
    all(coefficient_isfinite,c) || throw(DomainError(c,"coordinate conversion overflowed"))
    c = c / maximum(abs,c)
    v = SVector(c[2],c[3],c[4])
    model === :klein && return v/c[1]
    norm2 = -_lorentz(c,c)
    p isa HyperbolicPoint && norm2 <= 0 &&
        throw(DomainError(p,"timelike norm lost at working precision; use BigFloat"))
    s = p isa IdealPoint ? zero(c[1]) : sqrt(norm2)
    model === :ball && return v/(c[1]+s)
    d = c[1]-c[4]
    d > zero(d) || throw(DomainError(p,"point is at half-space infinity or chart precision was lost"))
    SVector(c[2]/d,c[3]/d,s/d)
end

function _modelinput(c,model)
    _checkmodel(model)
    length(c) == 3 || throw(ArgumentError("expected three model coordinates"))
    v = SVector{3}(map(x -> _analytic_value(real(x)),c))
    all(isreal,c) && all(coefficient_isfinite,v) ||
        throw(ArgumentError("model coordinates must be finite and real"))
    v
end

"""
    modelpoint(c; model=:ball)

Construct a `HyperbolicPoint` from three finite model coordinates.
Ball/Klein inputs must be strictly inside the unit ball; half-space height
must be positive. The returned representative need not be unit normalized.
"""
function modelpoint(c; model=:ball)
    v = _modelinput(c,model)
    r2 = sum(abs2,v)
    if model === :halfspace
        v[3] > 0 || throw(ArgumentError("half-space height must be positive"))
        return HyperbolicPoint((r2+1)/2,v[1],v[2],(r2-1)/2)
    end
    r2 < 1 || throw(ArgumentError("point must be inside the unit ball"))
    model === :klein && return HyperbolicPoint(one(r2),v...)
    HyperbolicPoint(1+r2,2v[1],2v[2],2v[3])
end

"""
    modelidealpoint(c; model=:ball)

Construct an ideal point from unit-sphere coordinates in `:klein`/`:ball`,
or `(x,y,0)` in `:halfspace`. Use `IdealPoint(1,0,0,1)` for chart infinity.
Sphere inputs within floating roundoff of unit length are normalized.
"""
function modelidealpoint(c; model=:ball)
    v = _modelinput(c,model)
    if model === :halfspace
        iszero(v[3]) || throw(ArgumentError("ideal half-space height must be zero"))
        r2 = v[1]^2+v[2]^2
        return IdealPoint((r2+1)/2,v[1],v[2],(r2-1)/2)
    end
    r = sqrt(sum(abs2,v))
    isapprox(r,one(r);rtol=64eps(typeof(r)),atol=0) ||
        throw(ArgumentError("ideal ball coordinates must have unit length"))
    IdealPoint(one(r),(v/r)...)
end
