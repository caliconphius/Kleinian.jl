shape = [5,4]

using Farey

function uptableaux(tab)
    nrows = length(tab)
    n = length(vcat(tab...))
    ups = [
    [tab..., [n+1]]
    ]
    previous_length = length(tab[1])+1
    for (i, row) in enumerate(tab)
        newtab = tab|>copy
        if length(row) < previous_length
            newtab[i] = [row..., n+1]
            push!(ups, newtab)
        end
        previous_length = length(row)
    end
    ups
end

tab = [[1,2,3],[4,5]]
tab = [[1,2],[3,4], [5,6]]

uptableaux(tab)



nrows = length(tab)
n = length(vcat(tab...))
ups = [
    [tab..., [n+1]]
]




# newtab = tab
# newtab[1,:] = newrow
# ups = cat(ups, newtab)
for i in 1:nrows
    # if length(tab[i+1,:]) < length(tab[i,:])
    #     newrow = cat(tab[i+1,:],n+1)
    #     newtab = tab
    #     newtab[i+1,:] = newrow
    #     ups = cat(ups,newtab)
    # end

    newtab = tab|>copy
    newtab[i] = [tab[i]..., n+1]
    push!(ups, newtab)
end
ups

1:10                        |> 
    x->map(l->[l, l+1], x)  |> 
    x->filter(l->l[1]>3, x) 


