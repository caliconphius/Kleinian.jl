using Pkg
root = normpath(joinpath(@__DIR__,".."))
Pkg.activate(@__DIR__)
Pkg.develop([PackageSpec(path=p) for p in
    (root,joinpath(root,"packages","Poincare"),joinpath(root,"packages","Riley"),joinpath(root,"packages","Farey"))])
Pkg.instantiate()
