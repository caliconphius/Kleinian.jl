module Poincare
using Kleinian
using StaticArrays: SVector, SMatrix
import Kleinian: mapcoefficients
# Coefficient/predicate hooks shared with the geometry foundation.
import Kleinian: _fielddiv, _lorentz, _pointcoordinates, _predicate_value,
                 _sqrt_preserving, _incident_zero, _chartvector
include("Polyhedra.jl")
include("Gluings.jl")
include("LimitSetHulls.jl")
export GeodesicSegment, HyperbolicPolygon, AbstractHyperbolicPolyhedron, HyperbolicPolyhedron, HyperbolicTetrahedron,
       vertices, faceindices, face, faceplane, edges,
       interiorpoint, interiorangles, polygonarea, edgelength, grammatrix,
       dihedralcos, dihedralangle, idealshape, edgeparameter, FacePairing,
       TetrahedralComplex, EdgeCycle, boundaryfaces, edgecycles, cycleword,
       verify_edge_gluings, LimitSetHull, limitsethull, hullbendingangles, hulledgesgeometry,
       pleatedges
end
