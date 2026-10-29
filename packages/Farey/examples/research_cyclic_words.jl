# Inspired by triangle_rewriting.jl and prove_farey_formula.jl: finite-order letters.
using Farey
using Farey.Monoids

function main()
    involution = CyclicGen("b0"; order=2)
    rotation = CyclicGen("c"; order=3)
    @show involution rotation inv(rotation)
    @assert isequal(involution^2, one(involution))
    @assert isequal(rotation^3, one(rotation))
    @assert isequal(inv(rotation), rotation^2)
    # Use structural isequal here: cyclic-generator == currently calls a
    # missing _exponent method. This is recorded as a cleanup target.
    # This is a word in cyclic generators, not a general triangle-group solver:
    # additional relations between different generators are not imposed here.
    word = involution * rotation
    @show word collect(eachgen(word))
    @assert isempty(collect(eachgen(word * inv(word))))
    # Negative integer powers and mixed-order arithmetic remain cleanup review targets.
end

main()
