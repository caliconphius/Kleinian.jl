using Farey
include("working_setup.jl")
qs = qnbs(7 // 31)
ws = [Fw(q) for q in qs]
tr = SFX.SuffixTree(ws[3]^2)
Fw(3 // 5)

SFX.atr.printnode(stdout, tr[2])
SFX.atr.print_tree(tr[0]; maxdepth=10)
# tr = sftr.SuffixTree("abaaaba")
# tr = sftr.SuffixTree("abacabaabaca")

# function crochemore(tre::sftr.SuffixTree{S}) where S
#   croche = Int[]
#   ix0 = 1; ix1 = 1
#   node = 0
#   while ix1 <= tre.lw
#     nxt = sftr.next(tre, node, ix1)
#     if (nxt > 0 || tre[nxt].position) && push!(croche, ix1)
#   end
# end

# tr
# sftr.nodepath(tr, sftr.walktree(tr, Int.([x for x in "baa"])))
# tr
# tr = sftr.SuffixTree("abaaaaba")
# tr = sftr.SuffixTree("abaaba")
# croch = let tr = tr
# croch = []

# ix0 = 1
# node = 0
# txt = sftr.text(tr)

# for ix1 in eachindex(txt)[2:end]
#   nxt = sftr.next(tr, node, ix1-1)
#   if (nxt>0) && sftr.nodepath(tr, nxt) <= ix0
#       node = nxt
#   else
#     push!(croch, ix1)
#     ix0 = ix1
#     node = 0
#   end
# end
# croch
# end


import .SuffixTrees as sftr

sftr.nodepath(tr, 7)
sftr.subword(tr, 4)
txt = sftr.text(tr)
node = sftr.walktree(tr, txt[3:4])
(tr[node].position)
(tr[node])
sftr.printpath(tr, 1) |> println

tr.tree
