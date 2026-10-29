@testset "Hyperbolic isometries" begin
    e0,e1,e2,e3 = (basisvector(i,big(1)//1) for i in 0:3)
    o = HyperbolicPoint(e0)
    p = HyperbolicPoint(clvector(big(5)//1,1,2,3))
    L = Geodesic(e0∧e3)
    h = HalfTurn(L)
    H = Horosphere(e0+e3)
    P = Parabolic(H,e1)
    T = Translation(L;multiplier=2//1)
    @test h(p) == HyperbolicPoint(clvector(5,-1,-2,3))
    @test h(h(p)) == p
    @test h(o) == o
    @test h(IdealPoint(e0+e3)) == IdealPoint(e0+e3)
    @test isidentity(h*h)
    @test inv(h) === h
    @test isidentity(h^0) && isidentity(h^-2)
    @test h == HalfTurn(reverse(L))
    @test h == VersorIsometry(3versor(h))
    @test h == MobiusIsometry(matrix(h))
    @test hash(h) == hash(MobiusIsometry(matrix(h)))
    @test T(o) == HyperbolicPoint(5,0,0,3)
    @test distance(o,T(o)) ≈ log(2)
    @test inv(T)(T(p)) == p
    @test T(IdealPoint(e0+e3)) == IdealPoint(e0+e3)
    @test T(H) == Horosphere(2*(e0+e3))
    @test T(H) != H
    @test P(o) == HyperbolicPoint(3,2,0,1)
    @test P(H) == H
    @test P(o) in H
    @test P(IdealPoint(e0+e3)) == IdealPoint(e0+e3)
    @test inv(P)(P(p)) == p
    @test iszero((representative(H)∧e1)^2)
    @test scalarpart(versor(P)*reverse(versor(P))) == 1
    @test Parabolic(H,e1+3*(e0+e3)) == P
    @test Parabolic(Horosphere(2*(e0+e3)),e1/2) == P
    @test (T*P)(p) == T(P(p))
    @test (P*T)(p) != (T*P)(p)
    @test (T∘P)(p) == T(P(p))
    @test inv(T*P)((T*P)(p)) == p
    @test (T*P).factors == (T,P)
    @test (P^3)(p) == P(P(P(p)))
    for g in (h,P,T,T*P,inv(T*P))
        M = MobiusIsometry(matrix(g))
        @test M(p) == g(p)
        @test M(H) == g(H)
        @test M(L) == g(L)
        @test M(HyperbolicPlane(e1)) == g(HyperbolicPlane(e1))
        @test inv(M)(M(p)) == p
        @test MobiusIsometry((2+3im)*matrix(g))(p) == g(p)
        @test MobiusIsometry(3matrix(g))(H) == g(H)
        @test VersorIsometry(versor(M))(p) ≈ g(p)
        @test distance(g(o),g(p)) ≈ distance(o,p)
    end
    # Independent skew-axis example: two half-turns yield translation + twist.
    h1 = HalfTurn(Geodesic(e0∧e1))
    h2 = HalfTurn(Geodesic((5*e0+3*e3)∧(e1+e2)))
    G = h2*h1
    A = matrix(G)
    @test iszero(A[1,2]) && iszero(A[2,1])
    @test A[1,1]/A[2,2] == -4im
    @test distance(o,G(o)) ≈ log(4)
    @test G(p) == MobiusIsometry(A)(p)
    # Finite order and axis orientation, with numerical trig evaluated on demand.
    E = Elliptic(mapcoefficients(Float64,L),4)
    @test isidentity(E^4;atol=1e-13)
    @test !isidentity(E^2;atol=1e-13)
    @test E(HyperbolicPoint(2.0,1.0,0.0,0.0)) ≈ HyperbolicPoint(2.0,0.0,1.0,0.0)
    @test inv(E)(E(p)) ≈ p
    # Mixed exact/Float64 arithmetic can promote to BigFloat without gaining
    # accuracy beyond that of the Float64 length parameter.
    @test isapprox(Translation(mapcoefficients(Float64,L),log(2.0))(o),T(o);rtol=1e-14)
    @test_throws ArgumentError Translation(L;multiplier=0)
    @test_throws ArgumentError Elliptic(L,4;numerator=2)
    @test_throws ArgumentError Parabolic(H,e0)
    @test_throws ArgumentError Parabolic(H,e0+e3)
    @test_throws ArgumentError VersorIsometry(e0)
    @test_throws ArgumentError VersorIsometry(e0*e1)
    @test_throws ArgumentError MobiusIsometry([1 2;2 4])
    @test_throws DimensionMismatch MobiusIsometry([1 2 3])
    @test_throws ArgumentError IsometryProduct(())
    Mnegative = MobiusIsometry([1 0;0 -1])
    @test VersorIsometry(versor(Mnegative))(p) == h(p)
    @test (Mnegative^2)(p) == p
    @test (T^100)(o) == MobiusIsometry(matrix(T^100))(o)
    @test length(Test.detect_ambiguities(Kleinian)) == 0
    # Reevaluate construction parameters at the chosen precision.
    setprecision(BigFloat,256) do
        Eb = mapcoefficients(BigFloat,Elliptic(L,7))
        @test eltype(versor(Eb)) == BigFloat
        @test isidentity(Eb^7;atol=BigFloat(2)^(-230))
        @test precision(real(matrix(Eb)[1])) == 256
    end
    @test eltype(matrix(h)) == Complex{Rational{BigInt}}
    @test eltype(matrix(P)) == Complex{Rational{BigInt}}
    pf=mapcoefficients(Float64,p)
    Hf=mapcoefficients(Float64,H)
    hf=mapcoefficients(Float64,h)
    for scale in (1e-200,1e200)
        hs=HalfTurn(Geodesic(scale*representative(hf.axis)))
        Ms=MobiusIsometry(scale*matrix(hf))
        @test hs(pf) ≈ hf(pf)
        @test Ms(pf) ≈ hf(pf)
        @test hs(Hf) ≈ hf(Hf)
        @test Ms(Hf) ≈ hf(Hf)
        @test representative(hs.axis) == scale*representative(hf.axis)
    end
end
