# Adapted from working_setup.jl: Q_path, sw_descent, qpair and rational neighbours.
# Run from the repository root: julia --project=. examples/research_continued_fractions.jl
using Farey

function main()
    slope = ContinuedFraction(24 // 103)
    path = [CF(0, slope[2:i]) for i in 2:length(slope)]
    # Prefixes describe the continued-fraction route used in research decompositions.
    @show slope Rational.(path)
    left, right, difference = farey_neighbours(slope)
    @show left right difference
    @assert Rational(left ⊕ right) == Rational(slope)
    @assert abs(numerator(Rational(left)) * denominator(Rational(right)) -
                numerator(Rational(right)) * denominator(Rational(left))) == 1
    @show positive_form(ContinuedFraction(0, [2, 2, -1]))
end

main()
