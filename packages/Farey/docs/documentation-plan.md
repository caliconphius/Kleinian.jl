# Documentation and cleanup plan

Use [the codebase map](codebase-map.md) as the initial inventory. Keep cleanup
changes small enough that mathematical behavior can be checked independently.

The agreed first cleanup is consolidation onto Monoids and removal of the
outdated FreeGroup implementation, KnuthBendix and GroupsCore. The
[cleanup tracker](cleanup-tracker.md) takes precedence for that migration's
scope, sequence and validation; the phases below remain the broader roadmap.

## 1. Establish intended scope and a behavioral baseline

- Identify the thesis scripts/notebooks and mathematical definitions this
  package must continue to support, including notation and word conventions.
- Decide the intended status of custom monoids, suffix trees and the unloaded
  palindrome experiments: supported API, research extension or retired work.
- Add focused tests for reduction/inversion, homomorphism application and
  composition, Christoffel/Farey words, S-sequences and suffix reconstruction.
  Use independently known results or mathematical identities rather than
  reproducing the algorithms in tests.
- Run both example scripts and record discrepancies before rewriting them.

Done when the supported features and mathematical conventions are recorded
and their key behavior has reproducible checks. The existing four tests are
the starting baseline, not sufficient protection for broad cleanup.

## 2. Write the user documentation

- Expand the README with purpose, supported scope, setup from this checkout,
  a minimal working example, and links to topic guides.
- Write guides for continued fractions; free groups and homomorphisms;
  custom monoid/cyclic words; Christoffel/Farey words and S-sequences;
  and suffix trees if retained as supported functionality.
- State rational and generator conventions, domains, return types, operator
  meanings, mutation behavior and known limitations alongside examples.
- Add docstrings to supported exported functions and constructors; explain
  dependency aliases and distinguish public helpers from implementation details.
- Connect algorithms to thesis definitions/results where references are
  available. Check example outputs rather than copying unchecked comments.

Start with Markdown guides and Julia docstrings. Introduce a documentation
generator later if a published API reference is wanted.

Done when a new reader can install the package and reproduce a representative
calculation for each supported topic without reading implementation code.

## 3. Remove confirmed unused material

- Review the evidence table in the map against external research code.
- Remove unused imports/dependencies first; move development tools into a
  separate environment if still needed.
- Resolve the unloaded palindrome module and placeholder types separately.
- Review commented algorithms, inactive settings and the nested environment.
  Preserve mathematical rationale in documentation where useful.
- Check package loading and meaningful tests after each change; verify a
  clean environment and the minimum supported Julia version for dependency changes.

Done when each deletion has a stated reason and validation, and retained
experimental code has an explicit status. Lack of a local caller alone is
not sufficient justification to remove externally accessible functionality.

## 4. Finish maintenance documentation

- Document the loading structure, both word representations, and how to add
  tests and examples; keep the codebase map synchronized with changes.
- Ensure CI exercises supported Julia versions and the examples that become
  documentation. Add documentation checks if a generator is introduced.
- Record any API changes and migration advice in a useful changelog.

The next concrete implementation step is to validate the two existing examples
and add mathematical regression tests for the main word-construction path.
That establishes reliable examples for the first documentation guides and
protects the subsequent cleanup.
