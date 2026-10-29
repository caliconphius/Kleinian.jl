# Adapted from suffixes_test.jl: repeated words and an ordinary text tree.
using Farey
using Farey.Monoids

function main()
    f, g = FreeGen("f g")
    word = (f * g)^3 * f
    tree = SFX.SuffixTree(word)
    @show word length(tree.tree) length(SFX.suffixes(tree))
    query = Int[SFX.genId(letter) for letter in eachgen(f * g * f)]
    node, _ = SFX.walktree(tree, query)
    @assert node != 0
    @assert SFX.walktree(tree, Int[Int(id(FreeGen("z")))])[1] == 0
    @show SFX.subword(tree, tree[node])

    text_tree = SFX.SuffixTree("abacabaabaca")
    @assert SFX.walktree(text_tree, Int[Int(c) for c in "abaca"])[1] != 0
    @show length(text_tree.tree)
    # Export to a temporary directory; Graphviz is not required to write DOT.
    mktempdir() do directory
        target = joinpath(directory, "text_tree")
        SFX.dotfile!(text_tree, target; to_string=identity)
        @assert isfile(target * ".dot")
        println("DOT export succeeded (temporary file).")
    end
    # The current symbol encoding uses a fixed alphabet and limited generator IDs.
    # Custom-word DOT export currently encounters a missing `word` import;
    # use the string tree for export until that adapter is repaired.
end

main()
