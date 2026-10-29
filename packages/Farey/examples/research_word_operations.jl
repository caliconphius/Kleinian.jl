# Adapted from working_setup.jl and hall_words.jl: words, conjugation and reversal.
using Farey
using Farey.Monoids

function main()
    f, g = FreeGen("f g")
    word = f^2 * g * inv(f)
    @show word collect(gens(word)) collect(eachgen(word)) word[2:3]
    # gens returns power blocks; eachgen returns individual signed letters.
    @assert length(collect(eachgen(word))) == length(word)
    @assert word * inv(word) == one(word)
    @assert word^g == inv(g) * word * g
    @assert word / g == word * inv(g)
    @assert g \ word == inv(g) * word
    @show word^g word^(-g) reverse(word)

    palindrome_word = f * g * f
    @assert palindrome_word == reverse(palindrome_word)
    # conj_prefix returns the removed outside prefix and the reduced core:
    # word = inv(prefix) * core * prefix in this convention.
    conjugated = inv(g) * (f * g) * g
    conjugator, core = conj_prefix(conjugated)
    @show conjugator core
    @assert conjugated^inv(conjugator) == core
    @show element("f^2 g -f")
    # Uppercase letters name distinct generators; inverse letters use inv or '-'.
end

main()
