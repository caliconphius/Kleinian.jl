# Farey

<!-- [![Build Status](https://github.com/caliconphius/Farey.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/caliconphius/Farey.jl/actions/workflows/CI.yml?query=branch%3Amaster) -->


### Setup Instructions 
Follow the instructions at 
"https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project"


### Functionality
Currently there are 3 main sets of functions.

`ContinuedFractions.jl` contains code for working with continued fractions

`FreeGroup.jl` and `Homomorphisms.jl` has code for working with words in the free group and homomorphisms between free groups

`FareyWords.jl` has code for computing Farey and Christoffel words

See `examples/` for how to use some of these functions

See the [representative research examples](examples/README.md) for standalone
scripts using the preferred Monoids implementation across the library's main areas.

### Development documentation

- [Cleanup tracker](docs/cleanup-tracker.md): migration to Monoids and removal of KnuthBendix/GroupsCore, with dependency inventories and ordered tasks.
- [Codebase map](docs/codebase-map.md): modules, loading structure, API entry points, validation baseline and cleanup candidates.
- [Documentation and cleanup plan](docs/documentation-plan.md): proposed sequence and completion criteria.
