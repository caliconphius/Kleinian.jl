# Farey research experiments

Moved on 2026-10-29 from `/home/tcsiz/Code/Math-Code/farey-testing`.
The original scripts and local environment are retained. References to the old
Farey checkout (including Windows paths) now point to `../../packages/Farey`.
The neighbouring `windowfarey_testing` directory was not imported.

These are exploratory REPL scripts, not a runnable test suite or a stable API.
Several depend on variables established by earlier interactive work, repeat
module definitions, use older Farey interfaces, or contain unfinished expressions.
Inspect and run selected definitions rather than including every file.

## Work developed here

| Area | Main files | What the code explores |
| --- | --- | --- |
| Farey words and continued fractions | `farey_equations.jl`, `fareydecomp.jl`, `erics_functions.jl`, `farey_tree.jl`, `tree_builder.jl`, `prove_farey_formula.jl` | Neighbour pairs, recursive word decompositions, continued-fraction substitutions, Stern–Brocot paths and tree visualizations. `prove_farey_formula.jl` is a working notebook of identities, not a checked proof. |
| Trace and Riley polynomials | `sl2algebra.jl`, `farey_polynomials.jl`, `trace.jl`, `poly_eval.jl`, `traceroots.jl` | Evaluating symbolic words in two-generator matrix representations, polynomial traces and roots. `sl2algebra.jl` constructs a four-dimensional algebra on `1,f,g,fg` with parameters `tr(f),tr(g),tr(fg)` and Cayley–Hamilton multiplication rules. |
| Trace fields and algebraic representations | `field_calc_hecke.jl`, `trace_field_working.jl`, `so31.jl`, `nemoGA.jl` | Hecke/Nemo calculations with cyclotomic fields, primitive elements, embeddings, polynomial roots, and matrix/Lorentz representations. These are candidate calculations, not a general discreteness test. |
| Automorphisms, Fox derivatives and group rings | `BirmanJacobian.jl`, `cohomology_working.jl`, `ga_testing.jl`, `standard_comms.jl`, `farey_commutators.jl` | Word substitutions, Jacobian matrices built from group-ring derivatives, augmentation and commutator calculations. `birmanJac` applies `GroupAlgebras.∂` to images of generators. |
| Relators and triangle groups | `magnus_onerelator.jl`, `ohtsukirileysakuma.jl`, `triangle_rewriting.jl`, `trianglegroupstudy.jl`, `sakuma_weeks.jl` | Magnus expansions, continued-fraction manipulations, finite-order generator substitutions, and experimental relations motivated by two-bridge and triangle groups. |
| String algorithms for words | `SuffixTrees.jl`, `suffix_tree.jl`, `sfx_tree_1.jl`, `crochemore.jl`, `tries.jl`, `suffixes_test.jl`, `transducerstest.jl` | Alternative suffix-tree/trie implementations, repeated substring searches, factorization and streaming transformations. These are overlapping prototypes rather than one chosen implementation. |
| Shared scratch setup | `working_setup.jl`, `setup_scripts.jl`, `working.jl`, `workingmaaaaa.jl` | REPL imports, helpers and exploratory calculations; the setup scripts draw on Farey's retained `working/` sources. |

Some files are placeholders: `TraceAlgebra.jl` is commented-out scaffolding;
`reidemeister_schreier.jl` and `backupsufiiixtree.jl` are empty.

## Environments and known limitations

For the maintained word-to-isometry example, start in the repository root:

```sh
julia experiments/setup.jl
julia --project=experiments experiments/farey_riley.jl
```

The archive has a separate environment with heavier research dependencies
(Oscar, Hecke, Transducers and others). To prepare it, optionally run:

```sh
julia experiments/farey-testing/setup.jl
julia --project=experiments/farey-testing
```

This setup selects the local Farey package; it does not execute the notebooks.
The archived manifest was retained with its Farey path updated. It may need
resolution on a different Julia version. Individual scripts also import packages
not declared in the original project, so this environment is not yet a complete
reproducibility specification.

Known migration work includes old `.modhoms` imports in `farey_equations.jl`,
`poly_eval.jl` and `sakuma_weeks.jl`: pointing their includes at the new location
does not make those old module interfaces compatible. `graphs_test.jl` contains
an unfinished import. Several files use REPL state before defining it.
The research scripts have been inspected, but have not been executed as a suite.

## Possible extraction into the libraries

- **Farey:** tested continued-fraction/word decomposition algorithms and reusable
  word operations; choose one string algorithm only after profiling a real use case.
- **Riley:** trace-polynomial evaluation and parameter/root exploration, using
  Farey words through the optional adapter. Slice or cusp claims need independent
  mathematical criteria beyond finding a polynomial root.
- **Poincare:** geometric domain constructions and verified face pairings derived
  from representations, once the scratch calculations have explicit inputs and checks.
- **Kleinian:** only reusable geometry and representation primitives.

Before extraction, turn a selected calculation into a function with explicit
inputs, establish its conventions, and add checks against small known examples.
