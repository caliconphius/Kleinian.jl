using Farey
using Graphs, TikzGraphs, TikzPictures
const FG = Farey.FreeGroup
const FGe = Farey.FreeGroupElement

struct FGroupItr end

function BFS(Gr::Farey.FreeGroup; depth = 1)
    gs = MON.elements(Gr.monoid,depth)[1]
Gr.(gs)
end

BFS(F2; depth=3)

[x for x in ITR.take(F2.alphabet,10)]

Base.eltype(::Type{FG})  = FGe
Base.iterate(Gr::FG) = one(Gr), 1

function Base.iterate(c::ContinuedFraction, state::Int)
    return state < c.length ? (c.L[state],state+1) : nothing
end

Base.IteratorSize(::Type{ContinuedFraction{S}}) where S<:Integer= Base.HasLength()

B = binary_tree(3)

t = B|>plot;
save(SVG("wheel10.svg"), t)