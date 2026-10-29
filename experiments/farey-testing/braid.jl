# BraidGroup.jl – Garside normal form computation
# This script demonstrates how to compute braid group elements using
# the Garside normal form. The Garside element (Δ) is the least common
# multiple of all standard generators and serves as the "atomic" element
# for normal form decompositions.

# Install and load the GAP package (provides BraidGroup with Garside normal forms)
using GAP
BraidGroup = GAP.BraidGroup

# -----------------------------------------------------------
# Example 1: Build a braid directly in Garside normal form
# -----------------------------------------------------------
n = 3  # number of strands
Δ = BraidGroup.GarsideElement(n)   # the Garside element
a = BraidGroup.Generator(1)        # σ₁
b = BraidGroup.Generator(2)        # σ₂

# Normal form: a^2 · b · a (i.e., σ₁² σ₂ σ₁)
nf = [a, a, b, a]                  # list of Garside atoms in order

# Construct the braid and display it
braid = BraidGroup.NormalForm(nf)
println("Braid in Garside normal form:")
println(braid)

# -----------------------------------------------------------
# Example 2: Compute the normal form from a random braid word
# -----------------------------------------------------------
# Define a function to convert a word in the standard generators
# σ₁, …, σ_{n−1} to its Garside normal form.
function garside_normal_form(word::Vector{BraidGroup.BraidWord})
  # Ensure each entry is a word (not already a BraidWord object)
  return BraidGroup.NormalForm(word)
end

# Example word: (σ₂ σ₁)⁴ Δ⁻¹ (σ₁ σ₂)²
w = [BraidGroup.Generator(2), BraidGroup.Generator(1),
  BraidGroup.Generator(2), BraidGroup.Generator(1)]

# For n = 3, the Garside element is Δ = (σ₁σ₂)(σ₂σ₁)
Δ_word = [BraidGroup.Generator(1), BraidGroup.Generator(2),
  BraidGroup.Generator(2), BraidGroup.Generator(1)]

# Compose the word
total = word ∪ Δ_word ∪ [BraidGroup.Generator(1), BraidGroup.Generator(2)]
nf = garside_normal_form(total)

println("\nNormal form of a constructed braid:")
println(nf)
println("Length (number of simple factors): ", length(nf))
