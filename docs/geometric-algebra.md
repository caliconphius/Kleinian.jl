# Geometric algebra kernel

This first implementation fixes Cl(3,1) with ordered basis `(e0,e1,e2,e3)`
and metric `diag(-1,1,1,1)`. **e0 is timelike**. Coordinates passed to
`clvector` are `(t,x,y,z)`. Basis indices run from 0 to 3.

## Representation for review

```julia
struct Multivector{T}
    coeffs::SVector{16,T}
end
```

The actual struct also has a typed constructor. StaticArrays supplies immutable
fixed-size storage; `T` has no subtype bound. Vectors, bivectors, and other grades
currently use this same representation, with zero coefficients outside their
grades. Separate compact grade types and a signature parameter are future API
choices, not required to use this kernel.

Coefficient order is ascending blade mask, rather than grouped by grade:

| Mask | Blade | Mask | Blade |
| --- | --- | --- | --- |
| 0 | 1 | 8 | e3 |
| 1 | e0 | 9 | e0e3 |
| 2 | e1 | 10 | e1e3 |
| 3 | e0e1 | 11 | e0e1e3 |
| 4 | e2 | 12 | e2e3 |
| 5 | e0e2 | 13 | e0e2e3 |
| 6 | e1e2 | 14 | e1e2e3 |
| 7 | e0e1e2 | 15 | e0e1e2e3 |

`coefficients(a)` returns the `SVector`, with Julia's ordinary **1-based** array
indices. `coefficient(a, mask)` uses blade masks **0:15**. `basisvector(i)`
uses basis indices **0:3**. These are separate indexing conventions.

## Display

Multivectors print as sums such as `(3) + (2)*e0 - (1)*e0e1`. Coefficients
are enclosed in parentheses, including unit coefficients and scalars. Zero
terms are omitted; the zero multivector prints as `0`. Terms are grouped by
grade, then ordered lexicographically by basis indices within each grade.
Negative real coefficients use an external minus sign; complex and unordered
coefficients retain their own representation inside parentheses.
This applies to both ordinary printing and REPL display.

## Usage

```julia
using Kleinian

e0 = basisvector(0)
e1 = basisvector(1)
v = clvector(5//3, 4//3, 0//1, 0//1)
B = e0 ∧ e1

v*v                          # scalar(-1)
grade(v + B, 2)              # B
leftcontract(v, e0)          # scalar(-5//3)
reverse(B)                   # -B
undual(dual(B))              # B
inv(v)*v                     # scalar(1)
```

See [the runnable example](../examples/cl31_basics.jl) for a reflection checked
against a coordinate formula. Geometry wrappers and point validation belong to
the next milestone; `clvector` constructs any vector, without normalization.

## Operations and conventions

| Operation | Contract |
| --- | --- |
| `a*b`, `geometricproduct(a,b)` | Clifford product |
| `a ⋅ b`, `dot(a,b)` | Symmetrized geometric product `(a*b+b*a)/2` |
| `a ∧ b`, `wedge(a,b)` | Outer product, extended bilinearly |
| `leftcontract(a,b)` | On grades r,s: grade s-r of a*b for r≤s; zero otherwise |
| `grade(a,r)` | Grade projection for r in 0:4, returning a Multivector |
| `reverse(a)` | Reversion: grade sign (-1)^(r(r-1)/2) |
| `gradeinvolution(a)` | Grade sign (-1)^r |
| `cliffordconjugate(a)` | Grade involution composed with reversion |
| `dual(a)` | Right multiplication by -I, where I=e0e1e2e3 and I²=-1 |
| `undual(a)` | Right multiplication by I; inverse of dual |
| `scalarpart(a)` | The coefficient of mask 0 |
| `checked_inverse(a)`, `inv(a)` | Reverse-norm formula with exact algebraic checks |

`dual(dual(a)) == -a`. Clifford conjugation does not complex-conjugate
coefficients. The dot product is symmetric and does not conjugate coefficients.
For vectors it returns their scalar Lorentz product as a multivector; for
mixed grades it may return nonscalar grades. It extends `LinearAlgebra.dot`
for multivectors, with `⋅` as the same function. Division by 2 follows the
coefficient type's `/` rules; rational inputs preserve exact arithmetic.
Integer powers use repeated squaring; negative powers require the restricted inverse.

## Generic coefficients and inversion limits

Addition, subtraction, and products assume a commutative coefficient ring with
`zero`, `+`, `-`, and `*`; constructors such as `basisvector` use an explicitly
supplied coefficient or default to Int. Type-based `zero` and `one` additionally
need those operations for `T`. Mixed coefficient types must provide Julia
promotion/arithmetic rules. The kernel does not convert coefficients to Float64.

Tests exercise integers, rationals, BigInt, Float32, Float64, BigFloat, complex
numbers, and a custom coefficient ring which does not subtype Number. This is
not a claim of compatibility with every symbolic or exact-arithmetic package.
Symbolic equality/zero tests may require simplification outside this kernel.

Scalar multiplication and division accept coefficient types outside Number.
For scalar addition with a non-Number coefficient, use an explicit scalar
multivector, e.g. `scalar(c)+a`. `mapcoefficients(f,a)` provides explicit
coefficient conversion. Value-based identities preserve coefficient parents.

Inversion requires coefficient division and decidable scalar/zero checks.
It verifies that both reverse products are the same nonzero scalar before
dividing. Zero reverse norm raises DomainError; nonscalar reverse products
raise ArgumentError. This rejects some invertible general multivectors and
floating-point versors with nonscalar rounding residuals. It deliberately
does not guess a tolerance or implement a general multivector inverse.
Integer division follows Julia's `/` rules; use rational coefficients when
an exact inverse is required.

## Implementation and verification

### Multiplication design review (2026-10-29)

Decision: retain the current table-based multiplication implementation. We
considered expressing multiplication as a matrix acting on the coefficient
vector, but agreed to leave the implementation unchanged.

Each basis blade `E_a` is indexed by a four-bit mask identifying its generators.
The defining relations are `e0² = -1`, `ei² = 1` for `i = 1,2,3`, and
`ei*ej = -ej*ei` for distinct generators. Therefore

```text
E_a * E_b = sign(a,b) * E_xor(a,b).
```

XOR removes shared generators. `_bladesign` counts the generator swaps needed
to restore ascending order and adds a negative metric sign when both blades
contain `e0`. `_PRODUCT_SIGNS` precomputes these signs for all 256 blade pairs.
For example, `(e0*e1)*(e0*e2) = e1*e2`: the generator swap and `e0² = -1`
contribute two cancelling minus signs.

For each output mask `k`, `_product` sums over left masks `a`. The corresponding
right mask is uniquely determined by `b = xor(a,k)`. In zero-based blade indices,

```text
(A*B)[k] = sum(sign(a, xor(a,k)) * A[a] * B[xor(a,k)] for a in 0:15).
```

Actual coefficient-array accesses add one to these masks for Julia indexing.
The shared accumulator uses `_keep` to select the operation:

- Geometric product: retain every pair.
- Wedge product: retain pairs with no shared generators (`a & b == 0`).
- Left contraction: retain the grade `s-r` part for input grades `r <= s`.

`dot` is implemented separately as `(A*B + B*A)/2`. It is the scalar Lorentz
product for vectors, but can contain nonscalar grades for general multivectors.
It does not conjugate coefficients. This existing convention is unchanged by
the review.

The alternative matrix formulation fixes the left operand and constructs a
16-by-16 matrix `L(A)` such that `coefficients(A*B) = L(A)*coefficients(B)`.
Its column `j` contains the coefficients of `A*E_j`; explicitly,

```text
L(A)[k,j] = sign(xor(k,j), j) * A[xor(k,j)].
```

This is again zero-based mathematical indexing. Wedge and contraction can use
the same construction with excluded entries set to zero. The matrix approach
makes application concise, but moves rather than eliminates blade-index and
sign logic. Adopting it would require checking generic coefficient support,
performance, and floating-point differences from changed summation order.
No matrix implementation or performance comparison was made in this review.

### Existing implementation and checks

The product uses a precomputed blade-sign table and 16 contributions per output
coefficient. Results are constructed directly as immutable SVectors, avoiding
mutable scratch storage that would restrict non-isbits coefficients such as
BigFloat. Product arithmetic never branches on whether a coefficient is zero.
This is a small auditable kernel; performance claims require measurement.

Tests compare all 256 blade products against a separate ordered-generator
reduction algorithm and check all 4,096 basis triples for associativity.
Additional tests cover grade conventions, mixed multivectors, involutions,
duality, coefficient types, null vectors, zero divisors, and reflections.

Validation on 2026-10-28: `Pkg.test()` passed the scaffold assertion and all
5,010 kernel assertions; `examples/cl31_basics.jl` completed successfully.
A warmed Float64 geometric product stored into a preallocated `Ref` allocated
0 bytes in a local `@allocated` check. This checks one operation and scalar
type, not general performance or the minimum supported Julia version.
