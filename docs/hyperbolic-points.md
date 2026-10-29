# Projective points of hyperbolic 3-space

Implementation: `src/HyperbolicGeometry.jl`.

`HyperbolicPoint{T}` wraps a `Multivector{T}` of grade one, with real
coefficients, positive time coordinate, and negative Lorentz square. The
ordered basis is `(e0,e1,e2,e3)`, so coordinates are `(t,x,y,z)` and
the metric is `diag(-1,1,1,1)`. The zero vector, null or spacelike vectors,
past-directed vectors, nonfinite coefficients, and mixed grades are rejected.

Points represent future timelike rays: `v` and `s*v` describe the same point
for any positive real `s`. Construction stores the supplied vector unchanged,
including its coefficient type and scale. No unit-length condition is imposed.

```julia
using Kleinian

p = HyperbolicPoint(clvector(5,3,0,0))
o = HyperbolicPoint(1,0,0,0)

coordinates(p)                 # SVector(5,3,0,0)
representative(p)              # original grade-one multivector
p == 7p                       # true: same projective point
distance(o,p) ≈ log(2)         # true
distance(3o,7p) ≈ distance(o,p)
u = unitvector(p)             # computes a unit representative on demand
scalarpart(u*u) ≈ -1          # true
coordinates(p)                 # still SVector(5,3,0,0)
```

The constructor also accepts a tuple or vector of four coordinates. `rescale(p,s)`,
`s*p`, and `p*s` return a new point and require finite positive `s`.
Point addition is not defined. Use representatives for geometric algebra
calculations, then validate any resulting point with `HyperbolicPoint`.

## Equality and distance

`==` compares projective points rather than their stored representatives.
Hashing and `isequal` use the same projective key, making scaled equivalents
interchangeable in dictionaries and sets. For built-in integer, rational,
and floating coefficient types, exact rational coordinate ratios avoid overflow
and make equality exact for the values actually stored. Rounded floating
rescaling can change a ray slightly; use `isapprox` (`≈`) in that case.
Approximate equality compares representatives in the time-one chart, without
changing storage or constructing unit timelike vectors.

For arbitrary future timelike representatives `v,w`, distance is

```text
acosh(-⟨v,w⟩ / (sqrt(-⟨v,v⟩) * sqrt(-⟨w,w⟩)))
```

This follows by applying the unit-hyperboloid distance formula to rescaled
representatives; see [Brice Loustau, Hyperbolic geometry](https://brice.loustau.eu/ressources/book.pdf),
the chapter on the hyperboloid model. The square roots are used only during
the calculation; stored vectors remain unchanged. `unitvector(p)` is the
explicit API for requesting a unit representative.

## Arithmetic and numerical limits

Validation uses widened integer/rational arithmetic to avoid overflow in
causal predicates. Floating-point working coordinates are rescaled before
squaring to avoid overflow or underflow caused solely by projective scale.
These working copies do not alter the stored representative.

Distance and explicit normalization use floating arithmetic. Float32 and
BigFloat working types are preserved. Floating-point classification near the
light cone and distance near coincidence remain precision-sensitive; use
BigFloat when necessary. A distance acosh argument below 1 by at most 16
machine epsilons is clamped to 1; larger violations raise DomainError.
Converting an exact timelike vector to insufficient precision can also raise
DomainError. No tolerance admits null or spacelike points during construction.

Point coefficients no longer require `T <: Real`: each coefficient must be
real-valued according to its backend's `isreal`, with real-part ordering and
finiteness predicates available. Value-based arithmetic preserves coefficient
parents. Foreign backends are not automatically converted to Float64 for
distance or unit-vector calculations; they must supply the required analytic
operations, or the user can explicitly convert with `mapcoefficients`.

See [geometry and isometries](isometries.md) for ideal points, planes,
geodesics, horospheres, and the optional Nemo algebraic-number adapter.
