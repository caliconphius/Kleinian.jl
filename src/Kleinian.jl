"""
    Kleinian

Tools for studying Kleinian groups and hyperbolic geometry.

The geometric algebra kernel uses the metric (-,+,+,+). See
`docs/geometric-algebra.md` for storage and operation conventions.
"""
module Kleinian

using StaticArrays: SVector, SMatrix
import LinearAlgebra: dot, ⋅

# All files share the Kleinian namespace. Load foundations before their users.
include("Coefficients.jl")
include("GeometricAlgebra.jl")
include("HyperbolicGeometry.jl")
include("HyperbolicObjects.jl")
include("Isometries.jl")
include("Models.jl")
include("Boundary.jl")
# Plotting interfaces live here; Makie methods live in ext/KleinianMakieExt.jl.
include("Visualization.jl")

export AbstractIsometry, basisblade, basisvector, checked_inverse, cliffordconjugate,
       clvector, coefficient, coefficient_complex, coefficient_imaginary_unit, coefficient_isfinite,
       coefficients, complexcoordinate, complexidealpoint, coordinates, crossratio,
       distance, dot, dual, Elliptic, Geodesic,
       geometricproduct, grade, gradeinvolution, HalfTurn, Horosphere,
       HyperbolicPlane, hyperbolicplot, hyperbolicplot!, HyperbolicPoint, idealcenter,
       IdealPoint, isidentity, isincident, IsometryProduct, leftcontract,
       mapcoefficients, matrix, MobiusIsometry, modelboundary!, modelcoordinates,
       modelidealpoint, modelpoint, Multivector, Parabolic, plotgeometry,
       representative, rescale, sameaxis, scalar, scalarpart,
       sl2matrix, Translation, undual, unitvector, versor,
       VersorIsometry, wedge, ∧, ⋅

end
