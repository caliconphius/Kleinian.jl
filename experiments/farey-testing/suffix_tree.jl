module SuffixTreeCore

export SuffixTree
using Farey
using FixedSizeArrays
const NextNode = FixedSizeArray{Int}
const ALPHABET_SIZE = 256

mutable struct SuffixNode
    start_idx::Int
    end_idx::Union{Int, Nothing}
    slink_idx::Union{Int, Nothing}
    next::NextNode
    label::String
    issuffix::Bool
    function SuffixNode(start; finish=nothing, slink=0, next=NextNode(undef,ALPHABET_SIZE), label="", suffix=false)
        next.=0
        new(start, finish, slink, next, label, suffix)
    end
end


struct SuffixTree{T} 
    text::T
    tree::Vector{SuffixNode}
end


function edge_length(nd::SuffixNode, idx)
    finish = isnothing(nd.end_idx) ? idx+1 : nd.end_idx
    finish - nd.start_idx
end

function Base.push!(tree::SuffixTree, nd::SuffixNode)
    push!(tree.tree, nd)
    # print("""
    
    # Added $nd at idx $(length(tree.tree))
    
    # """)
end

function Base.getindex(tree::SuffixTree, idx::Integer)::SuffixNode
    tree.tree[idx+1]
end

function text(tr::SuffixTree)
    [[Int(x) for x in tr.text]..., Int('#')]
end

# function extendtree!(
#     tree::SuffixTree, 
#     c::Integer, 
#     active_node::Integer, 
#     active_edge::Integer, 
#     active_len::Integer, 
#     remainder::Integer, 
#     idx::Integer, 
#     needSL::Integer
#     )

#     println.([
#     c, 
#     active_node, 
#     active_edge, 
#     active_len, 
#     active_len, 
#     idx, 
#     # needSL
#     ])

#     remainder=1
#     txt = text(tree)
#     while remainder>0 
#         if (active_len==0) active_edge=idx end
#         if tree[active_node].next[txt[active_edge]]==0
#             leaf = length(tree.tree)+1
#             push!(tree, SuffixNode(idx))
#             tree[active_node].next[txt[active_edge]] = leaf
#             if (needSL > 0) tree[needSL].slink_idx = active_node end
#             needSL = active_node

#         else
#             nxt = tree[active_node].next[txt[active_edge]]
#             if (active_len >= edge_length(tree[nxt], idx))
#                 @show edge_length(tree[nxt], idx)
#                 active_edge += edge_length(tree[nxt], idx)
#                 active_len -= edge_length(tree[nxt], idx)
#                 active_node = nxt
#                 continue
#             elseif (txt[tree[nxt].start_idx + active_len]==c)
#                 active_len+=1
#                 if (needSL > 0) tree[needSL].slink_idx = active_node end
#                 needSL = active_node
#                 break                
#             end
#             split = length(tree.tree)+1
#             push!(
#                 tree, 
#                 SuffixNode(
#                     tree[nxt].start_idx;
#                      finish=tree[nxt].start_idx+active_len
#                     )
#                 )
#             tree[active_node].next[txt[active_edge]] = split

#             leaf = length(tree.tree)+1
#             push!(tree, SuffixNode(idx))
#             tree[split].next[c] == leaf
#             tree[nxt].start_idx += active_len
            

#             tree[split].next[txt[tree[nxt].start_idx]] = nxt 
#             if (needSL > 0) tree[needSL].slink_idx = split end
#             needSL = split

#         end

#         remainder-=1
#         if active_node==1 && active_len>0
#             active_len -=1
#             active_edge = idx - remainder + 1
#         else
#             active_node = tree[active_node].slink_idx > 0 ? tree[active_node].slink_idx : 1
#         end
#     end
#     # tree, active_node, active_edge, active_len, remainder, idx, needSL
# end

function BuildSuffixTree(word::String) 
    needSL = 0; last_added = 1; idx = -1; remainder = 0
    active_node = 0; active_edge = 0; active_len = 0
    root = SuffixNode(0; finish=0, label="root")
    treel = [root]

    tree = SuffixTree(word, treel)
    txt = text(tree)
    for c in txt[1:end]
    # tr,
    edge = active_edge == 0 ? '\x0' : Char(txt[active_edge])
    [(;
        nd=tree[active_node].label, 
        edix=active_edge,
        ed=edge,
        len=active_len,
        rem=remainder,
        
        # c=Char(c), ix=idx,slink=needSL,
        ), 
    # active_len,
    # remainder, 
    # idx, 
    ].|>println
    remainder += 1
    idx += 1
    # needSL = 
        while remainder>0 

            if (active_len==0) active_edge=idx+1 end
            edge = txt[active_edge]

            if tree[active_node].next[txt[active_edge]]==0
                leaf = length(tree.tree)
                push!(tree, SuffixNode(idx+1;label = "<$leaf>"))
                tree[active_node].next[txt[active_edge]] = leaf

                if (needSL > 0) 
                    tree[needSL].slink_idx = active_node 
                end
                needSL = active_node

            else
                nxt = tree[active_node].next[txt[active_edge]]

                if (active_len >= edge_length(tree[nxt], idx))
                    active_edge += edge_length(tree[nxt], idx)
                    active_len -= edge_length(tree[nxt], idx)
                    active_node = nxt
                    continue
                
                end
                if (txt[tree[nxt].start_idx + active_len]==c)
                    active_len+=1
                    if (needSL > 0) 
                        tree[needSL].slink_idx = active_node 
                    end
                    needSL = active_node
                    break                
                end
                split = length(tree.tree)
                split_start = tree[nxt].start_idx
                split_end = tree[nxt].start_idx + active_len
                push!(
                    tree, 
                    SuffixNode(
                        split_start;
                        finish=split_end,
                        label = "($split)"
                    )
                    )
                tree[active_node].next[txt[active_edge]] = split

                leaf = length(tree.tree)
                push!(tree, SuffixNode(idx+1; label="<$leaf>"))
                tree[split].next[c] = leaf
                tree[nxt].start_idx += active_len
                

                tree[split].next[txt[tree[nxt].start_idx]] = nxt 
                
                if (needSL > 0) 
                    tree[needSL].slink_idx = split 
                end
                needSL = split

            end

            remainder-=1
            if active_node==0 && active_len>0
                active_len -=1
                active_edge = idx - remainder +1
            else
                active_node = tree[active_node].slink_idx > 0 ? tree[active_node].slink_idx : 0
            end
            # println()
            # println.(tree.tree)
            # println()

        end
    end

    tree
end

function children(nd::SuffixNode)
    [(k,nd.next[k]) for k in findall(nd.next.!=0)]
end

function Base.show(io::IO, nd::SuffixNode)
    print(io, """$(nd.label)($(nd.start_idx), $(nd.end_idx), $(nd.slink_idx), $([Char(k)=>nd.next[k] for k in findall(nd.next.!=0)]))""")
end

function Base.show(io::IO, tr::SuffixTree)
    root = tr[0]
    stack = [(1,x...) for x in children(root)]
    out =  " $(root.label) +\n"
    depth = 1
    path = [" "*prod(["-" for i in 2:length(root.label)])]
    while !isempty(stack)
        (parent, cr, ndix) = pop!(stack)
        nd = tr[ndix]
        chdn = children(nd)
        depth+=1
        
        if isempty(chdn)
            # out *= "\n"*"\t"^(depth)
            # out *= "\n|\n"
            depth-=1
            out*=prod([path[1:parent]..., "> $(Char(cr))|$(nd.label)"])*"\n"
            # print(io, "\n"*"\t"^depth*"")
        else
            out*=prod([path[1:parent]..., "> $(Char(cr))|$(nd.label) +"])*"\n"

            push!(path, "  | "*prod(["-" for i in 2:length(nd.label)]))
            append!(stack, [(depth, x...) for x in chdn])
        end

# "*prod(Char.(txt[split_start:split_end-1]))*", 
        
    end
    print(io, out)
    

end

function tree_print(tr::SuffixTree)
    root = tr[0]
    stack = [tr[x] for x in children(root)]
    out =  "\n+--$(root.label)\t"
    depth = 0
    while !isempty(stack)
        nd = pop!(stack)
        out*="+-$(nd.label)"
        chdn = children(nd)
        depth+=1
        
        if isempty(chdn)
            out *= "\n"*"\t"^(depth)
            out *= isempty(stack) ? "" : "|\n"*"\t"^(depth)
            depth-=1
            # print(io, "\n"*"\t"^depth*"")
        else
            append!(stack, [tr[x] for x in chdn])
        end
# "*prod(Char.(txt[split_start:split_end-1]))*", 
        
    end
    print(out)
    

end

end


