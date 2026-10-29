function setup()
  include("working_setup.jl")
  print("working_setup finished")

  include("sl2algebra.jl")
  print("sl2algebra fin")
  include("triangle_rewriting.jl")
  print("Δ finished")
  include("erics_functions.jl")
  print("erics_functions finished")
end
task2 = Threads.@spawn setup()
include(joinpath(@__DIR__, "../../packages/Farey/working/GroupAlgebra.jl"))
import .GroupAlgebras as ZG