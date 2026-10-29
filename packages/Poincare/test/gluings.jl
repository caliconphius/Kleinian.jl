include(joinpath(@__DIR__,"..","..","..","examples","trefoil333_data.jl"))

@testset "Trefoil (3,3,3) edge gluings" begin
    data = trefoil333()
    d = data.domain
    a,b,c = data.pairings
    f,g = data.generators
    @test length(boundaryfaces(d))==6
    @test all(t -> interiorpoint(t) in d,d.cells)
    @test idealshape(d.cells[1])≈conj(idealshape(d.cells[2]))
    @test sort([dihedralangle(d.cells[1],e)/pi for e in edges(d.cells[1])])≈[1/6,1/6,1/6,1/6,2/3,2/3]
    @test sum(dihedralangle(t,(2,3)) for t in d.cells)≈4pi/3
    r = verify_edge_gluings(d;atol=1e-10)
    @test r.valid
    @test length(r.edge_cycles)==4
    @test sort([e.order for e in r.edge_cycles])==[1,3,3,3]
    @test sort([length(e.cycle.sectors) for e in r.edge_cycles])==[1,1,4,6]
    @test sum(length(e.cycle.sectors) for e in r.edge_cycles)==12
    @test all(e -> e.endpoints_fixed && e.holonomy_order_correct && e.shape_equation,r.edge_cycles)
    omega = complex(-.5,sqrt(3)/2)
    @test all(e -> isapprox(e.cycle.shape_product,e.order==1 ? 1 : omega),r.edge_cycles)
    @test all(h -> isidentity(h;atol=1e-10),data.relations)
    @test isapprox(b,f*g*f*inv(g)*inv(f);atol=1e-10)
    @test isapprox(c,f*g;atol=1e-10)
    @test isidentity(f^3;atol=1e-10)
    @test isidentity(g^3;atol=1e-10)
    @test isidentity((f*g*f*inv(g)*inv(f)*inv(g))^3;atol=1e-10)
    # Deliberately incorrect cone orders cannot pass on angle sums or holonomy.
    wrong = verify_edge_gluings(d;orders=ones(Int,4),atol=1e-10)
    @test !wrong.valid
    @test count(e -> !e.valid,wrong.edge_cycles)==3
    @test_throws ArgumentError verify_edge_gluings(d;orders=[3,3])
    @test_throws ArgumentError verify_edge_gluings(d;angle_atol=-1)
    @test_throws ArgumentError TetrahedralComplex(d.cells,d.pairings[1:3])
    @test_throws ArgumentError TetrahedralComplex(d.cells,(d.pairings...,d.pairings[1]))
    bad = FacePairing((2,2),(2,3),one(a);label="incorrect")
    @test_throws ArgumentError TetrahedralComplex(d.cells,(bad,d.pairings[2:end]...))
    # Geometry and gluing certificates survive conjugation and recentering.
    moved = MobiusIsometry([2.0 -1.0;0.0 1.0])(d)
    @test verify_edge_gluings(moved;atol=1e-9).valid
    @test modelcoordinates(vertices(moved.cells[1])[2];model=:ball)≈[-1,0,0]
end
