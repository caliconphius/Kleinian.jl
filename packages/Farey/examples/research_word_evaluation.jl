# Adapted from working_setup.jl: mateval, monoideval and abelianization.
# This evaluator is an example helper, not an exported Farey API.
using Farey
using Farey.Monoids

function evaluate_letters(word, images; init)
    mapreduce(letter -> images[letter], *, eachgen(word); init=init)
end

# Exact inverse for the determinant-one 2x2 matrices used below.
inverse_sl2(matrix) = [matrix[2, 2] -matrix[1, 2]; -matrix[2, 1] matrix[1, 1]]

function main()
    f, g = FreeGen("f g")
    f_matrix = BigInt[1 1; 0 1]
    g_matrix = BigInt[1 0; 1 1]
    images = Dict(f => f_matrix, g => g_matrix,
                  inv(f) => inverse_sl2(f_matrix), inv(g) => inverse_sl2(g_matrix))
    identity_matrix = BigInt[1 0; 0 1]
    word = farey_word(2 // 5, f, g)
    result = evaluate_letters(word, images; init=identity_matrix)
    @show word result
    @assert result[1, 1] * result[2, 2] - result[1, 2] * result[2, 1] == 1
    @assert evaluate_letters(one(word), images; init=identity_matrix) == identity_matrix
    @assert evaluate_letters(inv(word), images; init=identity_matrix) == inverse_sl2(result)

    # Abelianization forgets order and sums signed generator counts.
    exponent_vectors = Dict(f => [1, 0], g => [0, 1],
                            inv(f) => [-1, 0], inv(g) => [0, -1])
    abelianization = mapreduce(letter -> exponent_vectors[letter], +,
                              eachgen(word); init=[0, 0])
    @show abelianization
end

main()
