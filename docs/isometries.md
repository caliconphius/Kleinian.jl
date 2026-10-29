# Geometric algebra objects and hyperbolic isometries

The geometry is specified in Cl(3,1), with basis `(e0,e1,e2,e3)` and metric
`(-,+,+,+)`. Implementations are in `src/HyperbolicObjects.jl` and
`src/Isometries.jl`; coefficient hooks are in `src/Coefficients.jl`.

## Geometric objects

| Type | Stored data | Scale and orientation |
| --- | --- | --- |
| `HyperbolicPoint(v)` | Future timelike vector | Positive multiples are equal |
| `IdealPoint(h)` | Future nonzero null vector | Positive multiples are equal |
| `Geodesic(L)` | Nonzero simple timelike bivector, L²>0 | Positive multiples are equal; sign reverses orientation |
| `HyperbolicPlane(n)` | Spacelike normal vector | Positive multiples are equal; sign reverses the plane's orientation |
| `Horosphere(h)` | Future nonzero null vector | **Scale is retained**; rescaling changes the horosphere |

`Geodesic(p,q)` constructs `representative(p) ∧ representative(q)` from two
distinct finite or ideal points. `reverse(L)` reverses its orientation;
`sameaxis(L,M)` compares axes without orientation. `representative(obj)`
returns the stored multivector. No object normalizes its representative.

`isincident(p,obj)` and `p in obj` use homogeneous equations. A point v is on
a geodesic L when v∧L=0, and on a plane n when v⋅n=0. The horosphere defined
by h satisfies -x⋅h=1 on the unit hyperboloid. For our projective points this
becomes (v⋅h)²=-v²; this incidence test needs no square root. Both vectors are
future directed, which fixes the sign lost by squaring.
`idealcenter(H)` returns the null ray corresponding to the horosphere's centre.

Incidence is exact by default, with an optional `atol` on the algebraic
residual. This tolerance is not projectively invariant: rescaling inputs
rescales the residual. Constructors validate exact backend values exactly.
For built-in floating coefficients, nullity, simplicity, and scalar reverse
norm checks allow residuals up to 64 machine epsilons after scaling the
largest input coefficient to magnitude one. Timelike/spacelike signs remain
strict. Near-degenerate inputs require adequate precision.

## Constructing isometries

```julia
using Kleinian
e0,e1,e2,e3 = (basisvector(i,big(1)//1) for i in 0:3)
axis = Geodesic(e0 ∧ e3)

h = HalfTurn(axis)
t = Translation(axis; multiplier=2//1) # signed length log(2)
r = Elliptic(axis,7)                    # angle 2π/7
H = Horosphere(e0+e3)
p = Parabolic(H,e1)

g = h*t*p                             # p acts first, then t, then h
x = HyperbolicPoint(5,1,2,3)
y = g(x)
inv(g)(y) == x
```

Write `2*e1`, not `2e1`: the latter is Julia's scientific-notation literal 20.0.

### Half-turns and compositions

A half-turn about L uses `R=dual(L)` and acts on projective points by
`v ↦ R*v*reverse(R)`. The input geodesic need not be normalized. A composition
`g*h` or `g ∘ h` means `g(h(x))`. Products retain their factors and act through
those factors. Integer powers retain the base and exponent and evaluate by
repeated squaring; they do not allocate a word of length n.

Products of two half-turns include screw motions about an axis, as demonstrated
in `examples/hyperbolic_isometries.jl`. The implementation accepts a pair of
half-turns directly; it does not automatically refactor arbitrary isometries
back into two half-turns.

### Translations and elliptics

`Translation(axis,length)` stores a signed hyperbolic length ℓ. Evaluation
uses `cosh(ℓ/2) + sinh(ℓ/2)*L/sqrt(L²)`.

`Translation(axis;multiplier=q)` stores q>0 with signed length log(q).
Its projective versor is `(q+1)*sqrt(L²) + (q-1)*L`. This avoids logarithms,
exponentials, and sqrt(q); it can remain exact with rational data when sqrt(L²)
is rational. Reversing the oriented axis reverses the translation direction.

`Elliptic(axis,n;numerator=k)` retains integers n,k with n>0 and gcd(n,k)=1.
Its angle is 2πk/n and its versor is
`cospi(k/n) + sinpi(k/n)*dual(L)/sqrt(L²)`. For the positive z-axis
`L=e0∧e3`, positive rotation takes the positive x-direction toward y.
The coefficient backend supplies the square root and trigonometric operations.
The order is exact construction data; evaluating with floating coefficients
does not make the computed nth power exactly the identity.

### Parabolics and horospheres

`Parabolic(H,u)` uses the null vector h stored by H and a nonparallel vector u
with h⋅u=0. It stores these inputs, then evaluates `N=h∧u`, `R=1+N/2`.
Since N²=0, this is a finite exact expression for exp(N/2). Adding a multiple
of h to u leaves the transformation unchanged. Rescaling h by s requires
rescaling u by 1/s for the same transformation. The direction parameter is
in the quotient h-perpendicular/span(h), with sign fixed by this rotor formula.

An isometry can act on points, ideal points, geodesics, planes, and horospheres.
**Horospheres need an actual Lorentz action, not just a projective action**:
for a GA versor R their null vector is divided by the positive scalar
`scalarpart(R*reverse(R))` after the sandwich. Omitting this factor would move
the horosphere even when the underlying transformation is the identity.

## Matrix representation

`matrix(g)` returns an `SMatrix{2,2,T,4}` representing the same isometry up to
nonzero complex scale. `MobiusIsometry(A)` wraps an invertible matrix without
normalizing its determinant. This PGL₂(C) representation describes the same
orientation-preserving isometries as PSL₂(C). `sl2matrix(g)` explicitly requests
a determinant-one lift, using the coefficient backend's square root branch.

The bridge uses the Hermitian matrix

```text
H(t,x,y,z) = [ t+z   x-iy ]
            [ x+iy  t-z  ]
H ↦ A H A†.
```

For Clifford conversion, e0ei maps to the Pauli matrix σi and I=e0e1e2e3
maps to i. These choices fix the boundary-coordinate and rotation conventions.
The boundary spinor coordinate is ζ=(x-iy)/(t-z) where defined; its action is
ζ ↦ (aζ+b)/(cζ+d).

The default action for GA constructors uses GA sandwiches. An explicitly
constructed `MobiusIsometry` uses Hermitian matrix arithmetic. Products act
factor by factor; a power containing matrix factors is evaluated as a matrix.
`versor(M::MobiusIsometry)` requests a normalized spin lift, which can require
a square root and a coefficient extension. Matrix actions themselves do not
require that conversion. A matrix action on a horosphere divides by |det A|
to remove the matrix representative's scale, which may require a square root.

`inv` uses reversion for GA versors and an adjugate for projective matrices.
`==`, `isequal`, and hashing compare projective matrix representatives.
`isapprox` is available for numerical coefficients. `isidentity(g;atol=0)`
tests a scalar matrix; its optional absolute tolerance depends on matrix scale.

## Coefficients and precision

Geometry and isometry types have no Real or Number subtype restriction on
their coefficient parameters. Supported backends must supply appropriate
value-based arithmetic (`zero(x)`, `one(x)`, +, -, *, /), realness and ordering
for real values, conjugation, equality, and the analytic operations actually
requested. An abstract number field needs a specified compatible embedding
and conjugation to act as a complex coefficient field.

The three extensible hooks are:

- `coefficient_isfinite(x)`: defaults to `isfinite(x)`.
- `coefficient_imaginary_unit(x)`: an i in a compatible complex domain.
- `coefficient_complex(x,y)`: constructs x+iy using that domain.

General coefficient multiplication uses value identities rather than
`one(T)` or `zero(T)`, since a type alone may not determine a field parent.
Explicit 2x2 multiplication also avoids array operations that assume Number.
When manipulating a returned matrix over a foreign field, use coefficient-wise
`map` or supported broadcasting instead of assuming scalar-array multiplication
is supplied by StaticArrays.

Use BigInt or Rational{BigInt} for growing exact computations: machine integers
still have their ordinary overflow limits. Native square root and analytic
operations can introduce approximations where their number type cannot express
the result. Rational perfect-square roots and exact integer divisions needed
by geometric constructors are preserved explicitly.

Construction parameters and factors are retained. To choose numerical precision:

```julia
setprecision(BigFloat,512) do
    numerical = mapcoefficients(BigFloat,Elliptic(axis,7))
    A = matrix(numerical)
end
```

`mapcoefficients` converts the source geometry and parameters; versors are then
recomputed. It cannot recover information already lost in Float64 inputs.
Mixed Float64 and exact big-number inputs may produce BigFloat outputs that
still carry only Float64 input accuracy; choose comparison tolerances accordingly.
Construct decimal parameters from strings inside the desired precision context.
For exact algebraic finite-order rotations, use the optional Nemo integration.

There is no general claim that half-turn/GA evaluation loses fewer digits than
2x2 matrices. Conditioning, scaling, operation count, and evaluation order all
matter. Exact coefficients eliminate rounding while supported operations stay
in that exact domain; BigFloat controls precision, not error certification.
Ball arithmetic and uncertain predicates need dedicated certification handling
and are not claimed supported by the present exact/numerical validation policy.
Floating working copies of projective versors, matrices, and point vectors
are rescaled before action/composition to reduce overflow from irrelevant scale.
Stored representatives are unchanged. Horosphere vectors retain their scale
through the correctly normalized Lorentz action. These rescalings do not cure
ill-conditioning or recover precision already lost near the ideal boundary.

## Optional Nemo integration

Nemo is a weak dependency. Loading it alongside Kleinian activates
`ext/KleinianNemoExt.jl`; ordinary use of Kleinian does not install or load Nemo.
The adapter targets Nemo 0.56's exact `QQBarFieldElem`, providing finiteness and
an imaginary unit in its algebraic closure. Other Nemo field types are not
covered by this adapter unless they already satisfy the documented protocol.

See `examples/nemo_isometries.jl` for an exact order-seven rotation about a
non-unit GA axis. The optional `test/nemo_integration.jl` suite runs in an
environment containing both packages. The regular test target exercises
built-in types and a parent-dependent field type outside Number, with no
type-based zero or one methods.

General translations specified by transcendental lengths may leave the
algebraic numbers. Use an exact multiplier when available, or explicitly map
to an analytic numerical backend; exact algebraic arithmetic is not a symbolic
representation of arbitrary transcendental functions.

## Validation record

On 2026-10-28, the full package test target passed 5,226 assertions on Julia
1.13.3. These include 114 isometry checks, 23 geometric-object checks, and
12 checks using a parent-dependent coefficient type outside Number. The
optional integration suite passed 58 checks with Nemo 0.56.2, including exact
finite-order rotations and matrix/GA agreement with algebraic coefficients.
Both `examples/hyperbolic_isometries.jl` and `examples/nemo_isometries.jl` ran
successfully. The minimum declared Julia version has not been tested here.
