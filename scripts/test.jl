# Resolve all local packages together, then test each in a fresh Julia process.
# Usage: julia --startup-file=no scripts/test.jl [core|plotting|exact]
using Pkg

suite = isempty(ARGS) ? "core" : only(ARGS)
suite in ("core", "plotting", "exact") || error("Expected core, plotting, or exact")
root = normpath(joinpath(@__DIR__, ".."))
packages = [root, [joinpath(root, "packages", name) for name in
                   ("Poincare", "Riley", "Farey")]...]
suites = Dict(
    "core" => ["test/runtests.jl", "packages/Poincare/test/runtests.jl",
               "packages/Riley/test/runtests.jl", "packages/Farey/test/runtests.jl",
               "test/package_integration.jl"],
    "plotting" => ["test/makie.jl"],
    "exact" => ["test/nemo_integration.jl", "test/trefoil_exact.jl"],
)

mktempdir() do environment
    Pkg.activate(environment)
    Pkg.develop([PackageSpec(path=path) for path in packages])
    dependencies = [PackageSpec(name=name) for name in ("Test", "Random", "StaticArrays")]
    if suite == "plotting"
        append!(dependencies, [PackageSpec(name="CairoMakie", version="0.15"),
                               PackageSpec(name="GeometryBasics", version="0.5")])
    elseif suite == "exact"
        push!(dependencies, PackageSpec(name="Nemo", version="0.56"))
    end
    Pkg.add(dependencies)
    Pkg.instantiate()
    for file in suites[suite]
        println("\nTesting ", file)
        run(`$(Base.julia_cmd()) --startup-file=no --project=$environment $(joinpath(root, file))`)
    end
end
