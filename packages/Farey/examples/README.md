# Representative research examples

Run each script from the repository root, for example:

```sh
julia --project=. examples/research_farey_words.jl
```

Each `research_*.jl` file is self-contained and uses the current Monoids
implementation. The scripts adapt selected ideas from the sibling
`farey-testing` checkout without depending on its globals, external includes,
Nemo/Hecke, plotting packages or presentation setup. Assertions check selected
mathematical properties; these examples do not replace a full regression suite.

| Script | Library area / demonstrated workflow | Research source |
| --- | --- | --- |
| [research_continued_fractions.jl](research_continued_fractions.jl) | Rational conversion, CF normalization, coefficient-prefix paths, Farey neighbours and mediants | `working_setup.jl`: Q_path, sw_descent and neighbour-based decompositions |
| [research_word_operations.jl](research_word_operations.jl) | Generator blocks versus letters, reduction, indexing, conjugation, divisions, reversal, common conjugating prefix and parsing | `working_setup.jl`, `hall_words.jl` |
| [research_cyclic_words.jl](research_cyclic_words.jl) | Finite-order generators and mixed cyclic words | `triangle_rewriting.jl`, `prove_farey_formula.jl` |
| [research_homomorphisms.jl](research_homomorphisms.jl) | Partial substitution, inverse images, word application and iterated maps | `hall_words.jl`, `BirmanJacobian.jl`, `working_setup.jl` |
| [research_farey_words.jl](research_farey_words.jl) | Christoffel counts, S-sequences, Farey words, alternating palindromes and a local Riley-word helper | `working_setup.jl` |
| [research_word_evaluation.jl](research_word_evaluation.jl) | Evaluating words in exact determinant-one matrices, target identity and abelianization | `working_setup.jl`: mateval, monoideval and ab |
| [research_suffix_trees.jl](research_suffix_trees.jl) | Trees of repeated words and text, pattern lookup, subword inspection and text DOT export | `suffixes_test.jl` |

`continued_fractions_ex.jl` and `free_group_ex.jl` are earlier examples.
The latter documents the legacy FreeGroup/Hom backend scheduled for removal;
the new scripts do not require that backend directly.

## Understanding the research workflow

Continued fractions organize rational slopes and neighbour/decomposition paths.
Christoffel words turn those slopes into ordered generator words. S-sequences
record blocks used to construct Farey words with alternating inverse letters.
Custom generators and words then supply reduction, powers, conjugation,
letter inspection and symbolic names. Partial homomorphisms rewrite these
symbolic expressions; generic evaluation sends the same letters into matrices
or another target algebra. Suffix trees inspect repeated subwords. Reversal
and palindromic/conjugating-prefix calculations provide additional word analysis.

The researcher-facing algebra thus has three distinct layers: word arithmetic,
substitution into other words, and evaluation into a different algebra.
The GroupInterface cleanup should preserve the shared arithmetic contracts;
Monoids supplies word storage and substitutions. Matrix evaluation helpers in
these examples remain local adaptations, not new library APIs.

## Current limitations exposed while preparing examples

- Cyclic-generator `==` calls a missing `_exponent` method; the cyclic example
  uses `isequal` for structural comparisons. Negative cyclic powers still need review.
- A reduced one-letter MonoidWord and a FreeGen may compare unequal even when
  their letter sequences agree. The substitution example compares letters.
- Custom-word DOT export encounters an undefined `word` binding in SuffixTrees.
  The suffix example exports its text tree with `to_string=identity` instead.
- `conj_prefix` returns the inverse of the element used to conjugate into the
  core. The word example states and checks that orientation explicitly.

The unloaded `src/Palindromes.jl` experiments and external PillowGroup,
automorphism and group-algebra implementations are not standalone supported
parts of this checkout. Their research uses informed these adaptations, but
cannot be reproduced fully without the missing implementations.
