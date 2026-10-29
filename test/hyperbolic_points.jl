@testset "Projective hyperbolic points" begin
    v = clvector(5,3,0,0)
    p = HyperbolicPoint(v)
    o = HyperbolicPoint(1,0,0,0)
    @test representative(p) === v
    @test coordinates(p) == SVector(5,3,0,0)
    @test p isa HyperbolicPoint{Int}
    @test HyperbolicPoint([5,3,0,0]) == p
    @test HyperbolicPoint{Int}(v) == p
    @test HyperbolicPoint((5//1,3//1,0//1,0//1)) == p
    @test p == 7p == p*2 == rescale(p,1//2)
    @test p == HyperbolicPoint(5.0,3.0,0.0,0.0)
    @test isequal(p,7p)
    @test hash(p) == hash(7p) == hash(HyperbolicPoint(5.0,3.0,0.0,0.0))
    @test length(Set([p,2p,7p])) == 1
    @test p != o
    @test p ≈ HyperbolicPoint(5.0,3.0+1e-12,0.0,0.0)
    @test distance(o,p) ≈ log(2)
    @test distance(o,p) ≈ distance(17o,23p)
    @test distance(p,o) ≈ distance(o,p)
    @test distance(p,7p) == 0
    @test scalarpart(unitvector(p)^2) ≈ -1
    @test HyperbolicPoint(unitvector(p)) ≈ p
    @test representative(p) === v # Neither distance nor unitvector mutates p.
    @test scalarpart(v*v) == -16
    @test distance(HyperbolicPoint(1f0,0f0,0f0,0f0),HyperbolicPoint(5f0,3f0,0f0,0f0)) isa Float32
    @test distance(HyperbolicPoint(big"1.0",0,0,0),HyperbolicPoint(big"5.0",3,0,0)) isa BigFloat
    @test HyperbolicPoint(typemax(Int),1,0,0) isa HyperbolicPoint{Int}
    for s in (1e-200,1e200)
        @test distance(HyperbolicPoint(s,0.0,0.0,0.0),HyperbolicPoint(5s,3s,0.0,0.0)) ≈ log(2)
    end
    @test HyperbolicPoint(big(10)^400,0,0,0) == o
    for c in ((0,0,0,0),(-5,3,0,0),(1,1,0,0),(1,2,0,0),
              (Inf,0,0,0),(1,NaN,0,0))
        @test_throws ArgumentError HyperbolicPoint(c)
    end
    @test_throws ArgumentError HyperbolicPoint(v+scalar(1))
    @test_throws ArgumentError HyperbolicPoint(v+basisblade(3))
    @test_throws ArgumentError HyperbolicPoint(clvector(5+1im,3,0,0))
    for s in (0,-1,Inf,NaN)
        @test_throws ArgumentError rescale(p,s)
    end
end
