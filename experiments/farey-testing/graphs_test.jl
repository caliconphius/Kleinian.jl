using Catlab.Graphs, Catlab.Graphics
using Catlab.Theories, Catlab.
CategoricalAlgebra
using DataStructures
SchLabeledGraph.generators
import Graphs as Gjl

# gr = @acset WeightedGraph{Farey.FreeGroupElement} begin
#   V = 4
#   E = 5
#   src = [1,1,1,2,3]
#   tgt = [2,3,4,4,4]
#   weight = [f, g, F, G]
# end


add_edge!(gr, add_vertex!(gr),add_vertex!(gr))
k = 4; depth=2

[(((i)÷(k))+2=>i+3) for i in 0:((k)^(depth+1)-1)÷(k-1)-2]
function nregtree(k, depth)
    gr = Gjl.SimpleDiGraphFromIterator(
        [Gjl.Edge(((i)÷(k))+2=>i+3) for i in 0:((k)^(depth+1)-1)÷(k-1)-2]
        )
    # Gjl.add_vertex!(gr)
    Gjl.add_edge!(gr, Gjl.Edge(1=>2))
    gr
end
gra = nregtree(2,4)
Gjl.outneighbors(gra, 10)
nbs = i->Gjl.outneighbors(gra, i).-2*(i-1)
[nbs(i) for i in 1:20]
struct VertexData{T}
    left::T
    center::T
    right::T
    graph::Gjl.DiGraph
    op
end



# iterate(v::VertexData) = 

Gjl.prufer_encode(gra)
Gjl.d

gra = Gjl.prufer_decode([3,3,3,3])
edg = [x for x in Gjl.edges(gra)]

subpart(gr, 1, :tgt)
subpart(gr, 1, :src)


(D1*G0^2).((f,f*g,g,f/g))
w = (g1^d1)
w = (g0^d0)

((w)|>hom).((f,f*g,g,(g*f)^-1))
(w|>inv|>hom).((f,f*g*f,g*f*g,g,(g*f)^-1))
(w|>inv|>hom).((
    f,
    g*f*g*f*g,
    g*f*g,
    f*g,
    g,
    )
)
(w|>hom).((
    f,
    g,
    f*g*f,
    g*f*g,
    (f*g),
    g*f,
    )
)
D1.((f,f*g,g))