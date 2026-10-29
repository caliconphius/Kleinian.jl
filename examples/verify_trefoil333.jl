include("trefoil333_data.jl")
data = trefoil333()
report = verify_edge_gluings(data.domain;atol=1e-10)
@assert report.valid
for e in report.edge_cycles
    println(e.word," | angle/π = ",e.angle/pi," | order = ",e.order,
            " | shape product = ",e.cycle.shape_product)
end
@assert all(g -> isidentity(g;atol=1e-10),data.relations)
println("PASS: face matching, opposite sides, four edge cycles, shape equations and presentation relations.")
