@testset "Hyperbolic model charts" begin
    origin = HyperbolicPoint(1,0,0,0)
    @test modelcoordinates(origin;model=:klein) == [0,0,0]
    @test modelcoordinates(origin;model=:ball) == [0,0,0]
    @test modelcoordinates(origin;model=:halfspace) == [0,0,1]
    for c in ((2,1,0,0),(5,1,2,3),(13,0,0,-12)), m in (:klein,:ball,:halfspace)
        p = HyperbolicPoint(c)
        v = modelcoordinates(p;model=m)
        @test isapprox(modelpoint(v;model=m),p;atol=1e-13)
        @test modelcoordinates(rescale(p,17);model=m) ≈ v
    end
    for c in ((1,1,0,0),(1,0,0,-1),(5,3,4,0)), m in (:klein,:ball,:halfspace)
        p = IdealPoint(c)
        v = modelcoordinates(p;model=m)
        @test isapprox(modelidealpoint(v;model=m),p;atol=1e-13)
        @test m == :halfspace ? v[3] == 0 : sum(abs2,v) ≈ 1
    end
    @test modelcoordinates(IdealPoint(1,0,0,1);model=:ball) == [0,0,1]
    @test_throws DomainError modelcoordinates(IdealPoint(1,0,0,1);model=:halfspace)
    @test_throws ArgumentError modelcoordinates(origin;model=:typo)
    @test_throws ArgumentError modelpoint((1,0,0);model=:ball)
    @test_throws ArgumentError modelpoint((2,0,0);model=:klein)
    @test_throws ArgumentError modelpoint((0,0,0);model=:halfspace)
    @test_throws ArgumentError modelpoint((0,0,-1);model=:halfspace)
    @test_throws ArgumentError modelpoint((0,0,Inf);model=:halfspace)
    @test_throws ArgumentError modelpoint((0,0,1im))
    @test_throws ArgumentError modelpoint((0,0))
    @test_throws ArgumentError modelidealpoint((0,0,0))
    @test_throws ArgumentError modelidealpoint((0,0,1);model=:halfspace)
    p = HyperbolicPoint(big"5.0",big"1.0",big"2.0",big"3.0")
    for m in (:klein,:ball,:halfspace)
        v = modelcoordinates(p;model=m)
        @test eltype(v) == BigFloat
        @test isapprox(modelpoint(v;model=m),p)
    end
    # Independent half-space distance formula, fixing the vertical convention.
    a,b = (0.2,-0.3,0.7), (-0.5,0.4,1.8)
    expected = acosh(1+sum(abs2,a.-b)/(2a[3]*b[3]))
    @test distance(modelpoint(a;model=:halfspace),modelpoint(b;model=:halfspace)) ≈ expected
end
