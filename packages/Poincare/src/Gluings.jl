"""
    FacePairing(source, target, isometry; label="g", internal=false)

A map from face `source=(cell,face)` to `target=(cell,face)`. Its inverse is
implicit. Tetrahedral face i is opposite vertex i. `internal=true` marks an
identity gluing across a shared face in the displayed domain.
"""
struct FacePairing{G<:AbstractIsometry}
    source::Tuple{Int,Int}
    target::Tuple{Int,Int}
    isometry::G
    label::String
    internal::Bool
end
FacePairing(source::Tuple{Int,Int},target::Tuple{Int,Int},g::AbstractIsometry;label="g",internal=false) =
    FacePairing(source,target,g,String(label),Bool(internal))

function _identity_atol(g,atol)
    A = matrix(g)
    s = maximum(abs,A)
    isidentity(MobiusIsometry(map(x -> _fielddiv(x,s),A));atol)
end

"""
    TetrahedralComplex(cells, pairings; atol=nothing)

A collection of tetrahedra with every face paired exactly once. The
constructor computes vertex correspondences from the isometries and checks
face matching and the interior-to-exterior side condition. It does not infer
global disjointness or certify a fundamental-domain theorem.
"""
struct TetrahedralComplex{C<:Tuple,P<:Tuple,V<:Tuple}
    cells::C
    pairings::P
    vertexmaps::V
end
function TetrahedralComplex(cells,pairings;atol=nothing)
    cs,ps = Tuple(cells),Tuple(pairings)
    !isempty(cs) && all(c -> c isa HyperbolicTetrahedron,cs) ||
        throw(ArgumentError("cells must be tetrahedra"))
    !isempty(ps) && all(p -> p isa FacePairing,ps) || throw(ArgumentError("face pairings are required"))
    tol = _poly_atol(Tuple(p for c in cs for p in vertices(c)),atol)
    used = Set{Tuple{Int,Int}}(); maps = Tuple[]
    for p in ps
        for f in (p.source,p.target)
            1 <= f[1] <= length(cs) && 1 <= f[2] <= 4 || throw(ArgumentError("invalid cell/face index"))
            f in used && throw(ArgumentError("a face has been paired more than once"))
            push!(used,f)
        end
        a,b = cs[p.source[1]],cs[p.target[1]]
        af,bf = faceindices(a)[p.source[2]],faceindices(b)[p.target[2]]
        correspondence = Pair{Int,Int}[]
        for i in af
            q = p.isometry(vertices(a)[i])
            matches = filter(j -> _samevertex(q,vertices(b)[j],tol),bf)
            length(matches)==1 || throw(ArgumentError("pairing $(p.label) does not match its face vertices"))
            push!(correspondence,i=>only(matches))
        end
        length(unique(last.(correspondence)))==3 || throw(ArgumentError("face vertex map is not bijective"))
        _planeside(faceplane(b,p.target[2]),p.isometry(interiorpoint(a))) > tol ||
            throw(ArgumentError("pairing $(p.label) does not carry the interior to the exterior"))
        p.internal && !_identity_atol(p.isometry,tol) && throw(ArgumentError("internal shared-face gluing must be the identity"))
        push!(maps,Tuple(correspondence))
    end
    length(used)==4length(cs) || throw(ArgumentError("every tetrahedral face must be paired"))
    TetrahedralComplex(cs,ps,Tuple(maps))
end

"""Exterior (cell,face) incidences, omitting identity-glued internal faces."""
boundaryfaces(d::TetrahedralComplex) = [f for p in d.pairings if !p.internal for f in (p.source,p.target)]
contains(d::TetrahedralComplex,p::PolyhedralPoint;kwargs...) = any(c -> contains(c,p;kwargs...),d.cells)
Base.in(p::PolyhedralPoint,d::TetrahedralComplex) = contains(d,p)

function _face_links(d::TetrahedralComplex)
    links = Dict{Tuple{Int,Int},Any}()
    for (i,p) in enumerate(d.pairings)
        mapping = Dict(d.vertexmaps[i])
        links[p.source] = (target=p.target,mapping=mapping,isometry=p.isometry,crossing=(i,1))
        links[p.target] = (target=p.source,mapping=Dict(v=>k for (k,v) in mapping),isometry=inv(p.isometry),crossing=(i,-1))
    end
    links
end

"""A circuit of tetrahedral edge sectors, with its ordered face crossings and full holonomy."""
struct EdgeCycle{H<:AbstractIsometry,A,S}
    sectors::Vector{Tuple{Int,Tuple{Int,Int}}}
    crossings::Vector{Tuple{Int,Int}}
    holonomy::H
    angle::A
    shape_product::S
    reverses_edge::Bool
end

"""
    edgecycles(complex; T=Float64)

Walk across paired faces, switching to the other face at each edge. Derive
the quotient edge classes, their total dihedral angle, full isometry product,
and (for ideal cells) complex shape product. Each tetrahedral edge sector is
counted once, including both sides of an internal face.
"""
function edgecycles(d::TetrahedralComplex;T=Float64)
    links = _face_links(d)
    visited = Set{Tuple{Int,Tuple{Int,Int}}}()
    result = EdgeCycle[]
    for ci in eachindex(d.cells),edge in edges(d.cells[ci])
        (ci,edge) in visited && continue
        outgoing = first(i for i in 1:4 if !(i in edge))
        start = (ci,edge,outgoing)
        cell,ordered,f = ci,edge,outgoing
        sectors = Tuple{Int,Tuple{Int,Int}}[]; crossings = Tuple{Int,Int}[]
        h = one(first(d.pairings).isometry)
        angle = zero(T); shape = nothing; allideal = true
        while true
            sector = (cell,minmax(ordered...))
            sector in visited && throw(ArgumentError("edge link is not a closed circular circuit"))
            push!(visited,sector); push!(sectors,sector)
            tetrahedron = d.cells[cell]
            angle += dihedralangle(tetrahedron,sector[2];T)
            if all(p -> p isa IdealPoint,vertices(tetrahedron))
                z = edgeparameter(tetrahedron,sector[2])
                shape = isnothing(shape) ? z : shape*z
            else
                allideal = false
            end
            link = links[(cell,f)]
            push!(crossings,link.crossing)
            h = link.isometry*h
            ordered = (link.mapping[ordered[1]],link.mapping[ordered[2]])
            cell,incoming = link.target
            f = only(i for i in 1:4 if !(i in ordered) && i != incoming)
            (cell,minmax(ordered...),f)==start && break
            length(sectors) <= 6length(d.cells) || throw(ArgumentError("edge circuit failed to close"))
        end
        push!(result,EdgeCycle(sectors,crossings,h,angle,allideal ? shape : nothing,ordered==reverse(edge)))
    end
    result
end

"""The cycle's matrix-composition word (rightmost factor acts first), omitting internal identities."""
function cycleword(d::TetrahedralComplex,c::EdgeCycle)
    factors = String[]
    for (i,s) in reverse(c.crossings)
        p = d.pairings[i]
        p.internal && continue
        push!(factors,s==1 ? p.label : p.label*"⁻¹")
    end
    isempty(factors) ? "1" : join(factors," ")
end

"""
    verify_edge_gluings(complex; orders=nothing, atol=nothing, angle_atol=1e-9)

Check every edge circuit: angle sum 2π/m, individually fixed endpoints,
holonomy of exact order m, and (for ideal cells) shape_product^m = 1. Shape
and matrix equations use exact arithmetic when atol=0; angles are numerical
and independently check the winding number. `orders` supplies one positive
integer per cycle, or is inferred from the angle sums. A failed inference is
reported with order zero. Returns `(valid, edge_cycles)` with detailed checks.
This verifies local edge gluings, not global tiling, cusp completeness, or
faithfulness of a group presentation.
"""
function verify_edge_gluings(d::TetrahedralComplex;orders=nothing,atol=nothing,angle_atol=1e-9)
    isfinite(angle_atol) && angle_atol >= 0 || throw(ArgumentError("angle_atol must be nonnegative and finite"))
    tol = _poly_atol(Tuple(p for c in d.cells for p in vertices(c)),atol)
    cycles = edgecycles(d)
    if !isnothing(orders)
        length(orders)==length(cycles) && all(m -> m isa Integer && m>0,orders) ||
            throw(ArgumentError("orders must specify one positive integer per edge cycle"))
    end
    reports = map(enumerate(cycles)) do (i,c)
        candidate = max(1,round(Int,2pi/c.angle))
        m = isnothing(orders) ? (abs(c.angle-2pi/candidate)<=angle_atol ? candidate : 0) : orders[i]
        residual = m>0 ? c.angle-2pi/m : Inf
        angle_ok = abs(residual)<=angle_atol
        cell,e = first(c.sectors)
        endpoints_ok = !c.reverses_edge && all(j -> _samevertex(c.holonomy(vertices(d.cells[cell])[j]),vertices(d.cells[cell])[j],tol),e)
        power_ok = m>0 && _identity_atol(c.holonomy^m,tol)
        primitive = power_ok && all(k -> m%k != 0 || !_identity_atol(c.holonomy^k,tol),1:m-1)
        shape_ok = isnothing(c.shape_product) ? nothing : m>0 && _incident_zero(c.shape_product^m-one(c.shape_product),tol)
        shape_residual = isnothing(c.shape_product) || m==0 ? nothing : abs(c.shape_product^m-one(c.shape_product))
        valid = angle_ok && endpoints_ok && power_ok && primitive && shape_ok !== false
        (cycle=c,word=cycleword(d,c),order=m,angle=c.angle,angle_residual=residual,
         endpoints_fixed=endpoints_ok,holonomy_power_identity=power_ok,
         holonomy_order_correct=primitive,shape_equation=shape_ok,shape_residual=shape_residual,valid=valid)
    end
    (valid=all(r -> r.valid,reports),edge_cycles=reports)
end

function (g::AbstractIsometry)(d::TetrahedralComplex)
    ps = map(p -> FacePairing(p.source,p.target,g*p.isometry*inv(g);label=p.label,internal=p.internal),d.pairings)
    TetrahedralComplex(map(g,d.cells),ps)
end
function mapcoefficients(f,d::TetrahedralComplex)
    ps = map(p -> FacePairing(p.source,p.target,mapcoefficients(f,p.isometry);label=p.label,internal=p.internal),d.pairings)
    TetrahedralComplex(map(c -> mapcoefficients(f,c),d.cells),ps)
end
