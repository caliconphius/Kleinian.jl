using AbstractTrees
using DataStructures: Trie, subtrie, partial_path
import DataStructures as datas
include("working_setup.jl")
using Farey
using Farey.Monoids

a = MonoidWord("a")

w = element("a 10 b")
# tt = Trie([w, w^2])
# keys(tt)

# keys(subtrie(tt, a^9))
# qs = qnbs(7//19)
# ws = [Γ(q) for q in qs]
# ws = [Fw(q) for q in qs]


# tt = Trie(ws, qs)
# keys(tt)
# st = subtrie(tt, ws[3])
# keys(st).|> x->prod(x, init=one(ws[1]))

# keys(st)
