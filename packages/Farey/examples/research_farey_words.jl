# Adapted from working_setup.jl: Christoffel/Farey words and signed Riley words.
using Farey
using Farey.Monoids

function riley_word(slope, f, g)
    p, q = denominator(slope), numerator(slope)
    signs(i) = (-1)^Int(floor(i * q // p))
    prod((letter^signs(i) for (letter, i) in
          zip(Iterators.cycle((f, g)), 1:(p-1))); init=one(f))
end

function main()
    f, g = FreeGen("f g")
    slope = 3 // 5
    christoffel_word = christoffel(slope, f, g)
    farey = farey_word(slope, f, g)
    sequence = s_seq(slope)
    @show slope christoffel_word sequence farey riley_word(slope, f, g)
    # For this positive slope, the Christoffel word has q copies of f and p of g.
    letters = collect(eachgen(christoffel_word))
    @assert count(==(f), letters) == denominator(slope)
    @assert count(==(g), letters) == numerator(slope)
    @assert sum(sequence) == 2 * denominator(slope)
    @assert christoffel(0, f, g) == f
    @assert s_seq(0) == [2]
    @show palindrome(f, g, 5) shrink_cf(slope)
    # The library's S-sequence/Farey path is restricted to slopes in [0,1].
end

main()
