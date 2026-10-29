include("working_setup.jl")
import .SuffixTrees as SFX

tr = SFX.SuffixTree(Fw(3//5))
tra = SFX.SuffixTree("abacabaabaca")

function crochemore(tr)
    node = 0
    cro = []
    N = 0
    edge = 0
    pos = 1
    len = 0
    while true
        if len == 0 
            edge = t
            nxt = SFX.next(tr, pos; idx = node)
            len = SFX.edgelength(tr, nxt)
        end

        if tr[nxt].position > N
            push!(cro, k)  
            
            node = 0
            pos = 0 
        else
            pos = pos + SFX.edgelength(tr, nxt)
            len += 
            node = nxt
            nxt = SFX.next(tr, pos; idx = node)

        end

        pos <= tr.lw || return cro
    end
    error("aaa")
end
# SFX.printpath(tr, 11)|>print
SFX.text(tr)

filter(x->x.suffix, tr.tree)
ndk = SFX.walktree(tr, SFX.text(tr)[4:end])

# SFX.printpath(tr, ndk)|>print
SFX.nodepath(tr,4)

SFX.text(tr, 11)

idx = 7
SFX.subword(tra, idx),SFX.finish(tra, idx), tra.word[SFX.nodepath(tra, idx):SFX.finish(tra, idx)]

# crochemore(tra)

SFX.walktree(tra, SFX.text(tra)[1:2])

# SFX.printpath(tra, 3)|> print

SFX.finish(tra, 1)

SFX.πw(tra, [98])

# tra[12]
SFX.finish(tra,4)
SFX.edgelength(tra,4)
SFX.next(tra, [98, 97, 99, 97])|>k->(k, SFX.subword(tra, k[1]))
# SFX.next(tra, [98, 97])|>k->(k, SFX.subword(tra, k[1]))


SFX.text(tra)
tra