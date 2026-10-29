# Adapted from hall_words.jl, BirmanJacobian.jl and working_setup.jl: partial maps.
using Farey
using Farey.Monoids

function main()
    f, g, h, k = FreeGen("f g h k")
    substitution = MonoidHom([h => -f, k => -g])
    word = h * f * k
    @show word substitution(word) (word | substitution)
    @assert collect(eachgen(substitution(word))) == [inv(g)]
    # Compare letters: a one-letter MonoidWord and a FreeGen do not currently
    # share a reliable cross-representation == implementation.
    @assert substitution(g) == g # Unmapped letters remain unchanged by default.

    # A Nielsen-style substitution and its explicit inverse, represented as words.
    forward = MonoidHom([f => f * g])
    backward = MonoidHom([f => f / g])
    original = f * g * f
    @assert backward(forward(original)) == original
    @show forward(original) forward^2
    @assert (forward^2)(original) == forward(forward(original))
    # Use explicit nested calls here. General composition and identity powers
    # need separate migration checks before treating them as supported defaults.
end

main()
