# Separate processes also verify clean package loading and avoid test namespace collisions.
root = normpath(joinpath(@__DIR__,".."))
for file in ("test/runtests.jl", "packages/Poincare/test/runtests.jl",
             "packages/Riley/test/runtests.jl", "packages/Farey/test/runtests.jl",
             "test/package_integration.jl")
    run(`$(Base.julia_cmd()) --startup-file=no --project=$(@__DIR__) $(joinpath(root,file))`)
end
