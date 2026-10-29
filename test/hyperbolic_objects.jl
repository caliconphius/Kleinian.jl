@testset "Geometric algebra objects" begin
    e0,e1,e2,e3 = (basisvector(i,1//1) for i in 0:3)
    o = HyperbolicPoint(e0)
    a,b = IdealPoint(e0-e3),IdealPoint(e0+e3)
    L = Geodesic(a,b)
    @test L == Geodesic(e0∧e3)
    @test sameaxis(L,reverse(L))
    @test L != reverse(L)
    @test o in L && a in L && b in L
    @test !(HyperbolicPoint(2,1,0,0) in L)
    @test IdealPoint(2,0,0,2) == b
    @test hash(b) == hash(IdealPoint(2,0,0,2))
    @test Geodesic(o,HyperbolicPoint(5,0,0,3)) == L
    plane = HyperbolicPlane(e1)
    @test o in plane && b in plane
    @test HyperbolicPlane(2*e1) == plane
    @test HyperbolicPlane(-e1) != plane
    H = Horosphere(e0+e3)
    @test o in H && HyperbolicPoint(10,0,0,0) in H
    @test !(o in Horosphere(2*(e0+e3)))
    @test idealcenter(H) == b
    @test H != Horosphere(2*(e0+e3))
    @test representative(H) == e0+e3
    @test_throws ArgumentError IdealPoint(e0)
    @test_throws ArgumentError IdealPoint(-e0-e3)
    @test_throws ArgumentError Geodesic(e1∧e2)
    @test_throws ArgumentError Geodesic((e0∧e1)+(e2∧e3))
    @test_throws ArgumentError Geodesic(o,o)
    @test_throws ArgumentError HyperbolicPlane(e0)
    @test_throws ArgumentError Horosphere(e0)
end
