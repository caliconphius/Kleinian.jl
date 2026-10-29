# Repository and development workflow

This is one Git repository containing four separately loadable Julia packages.
Keep changes that affect several packages in one commit so their interfaces
remain synchronized. There are no nested repositories or Git submodules.

| Package | Location | Responsibility |
| --- | --- | --- |
| Kleinian | repository root | Geometric algebra and basic hyperbolic geometry |
| Poincare | `packages/Poincare` | Polyhedra, face gluings, convex hulls and bending |
| Riley | `packages/Riley` | Limit-set sampling and parameter families |
| Farey | `packages/Farey` | Symbolic words and combinatorics |

Poincare and Riley depend on Kleinian. Farey is independent; Riley's Farey
adapter is optional. `experiments/` combines these packages and contains
research work that is not part of the maintained package APIs.

## Git workflow

`main` is the shared baseline. Use short-lived branches such as
`codex/geometric-algebra-readability` for changes, then review and merge them.
The initial commit is a development snapshot, not a release certification.
Generated manifests, local Julia caches, coverage output and editor state are
ignored. Documentation figures remain tracked deliberately.

Farey was imported as a working-tree snapshot, including local modifications.
Its original history remains in the original checkout and at
<https://github.com/caliconphius/farey-words>. See
[the import record](../packages/Farey/IMPORT.md) for the source commit.
The imported license remains under `packages/Farey/LICENSE`; it does not assign
that license to the other packages. The nested historical Farey workflow files
are retained as imported material; only root `.github/workflows/` files run here.

## Setup after cloning

From the repository root:

```sh
julia --startup-file=no experiments/setup.jl
julia --project=experiments experiments/check.jl
```

The setup script resolves the four packages from the checkout's actual location;
no account-specific absolute paths or published package releases are needed.
Julia 1.10 users should run the setup script rather than relying on `[sources]`.

## Automated checks

The same runner works locally and in GitHub Actions:

```sh
julia --startup-file=no scripts/test.jl core
julia --startup-file=no scripts/test.jl plotting
julia --startup-file=no scripts/test.jl exact
```

Each invocation creates a temporary environment and resolves all four local
packages together. Each test entry point runs in a separate Julia process.
`core` runs the four package suites and the Farey adapter integration tests.
`plotting` adds CairoMakie for headless rendering and viewer checks. `exact`
adds Nemo for algebraic geometry and trefoil checks. These checks do not execute
all research notebooks or certify unfinished experiments.

CI runs core checks on Julia 1.10 and current stable Julia; plotting and exact
checks run on current stable Julia. Remote CI results are available only after
publishing the repository to GitHub.

## Publishing

The local repository can be published as one GitHub repository. Choose the
owner, repository name and visibility before creating it, then configure the
remote and push `main`. GitHub access is unnecessary for local development.
Individual package repositories and independent releases can be introduced
later if needed; do not initialize nested repositories in this checkout.
