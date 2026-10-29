# Optional research environment: installs heavier algebra dependencies such as Oscar.
using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(path=normpath(joinpath(@__DIR__, "../../packages/Farey")))
Pkg.instantiate()
