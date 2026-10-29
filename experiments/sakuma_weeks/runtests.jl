using Test, Kleinian, Poincare
import Nemo, Farey
include("SakumaWeeks.jl")
using .SakumaWeeks

@testset "Farey rational regression" begin
    @test Rational(Farey.ContinuedFraction(big(2)//5)) == 2//5
    @test Rational(Farey.ContinuedFraction(big(0)//1)) == 0
end
@testset "Chesebro polynomial and funnel calculations" begin
    s = PolynomialSystem(); x = s.x
    @test Qpoly(s,1//2) == 1
    @test Qpoly(s,1//3) == 1-x
    @test Qpoly(s,1//4) == 1-2x
    @test Qpoly(s,2//5) == 1-x+x^2
    @test Qpoly(s,3//8) == 1-2x+2x^2
    @test Npoly(s,1//2) == 1-x
    @test Npoly(s,1//3) == 1-2x
    @test Z(s,only(interior_edges(2//5))) == s.field(x)
    @test [(e.qL,e.qR) for e in funnel(2//5)] == [(0,1),(0,1//2),(1//3,1//2)]
    @test length(interior_edges(3//8)) == 2
    @test_throws ArgumentError SBEdge(0,2//5)
    @test_throws ArgumentError funnel(0)
    @test_throws ArgumentError Qpoly(s,-1//2)
    @test_throws ArgumentError geometric_root(s,1//4)
    @test_throws ArgumentError developed_domain(s,2//5;parameter=0.5+0im)
    for q in (2//5,3//8,5//12,10//33,24//103)
        for e in interior_edges(q)
            @test Z(s,e) == (V(s,e.qR)-V(s,e.qP))/(V(s,e.qL)-V(s,e.qP))
        end
    end
    for n in 2:20, p in 1:n-1
        gcd(p,n)==1 || continue
        q=p//n; a,b=aab_dec(q)
        @test SakumaWeeks.mediant(a,SakumaWeeks.mediant(a,b)) == q
        @test Nemo.coeff(Qpoly(s,q),0)==1
        @test Nemo.coeff(Npoly(s,q),0)==1
    end
end
@testset "Developed geometry and matrix convention" begin
    s = PolynomialSystem()
    for (q,nt) in ((2//5,2),(3//8,4),(10//33,12))
        d = developed_domain(s,q)
        @test d.complete
        @test length(d.tetrahedra)==nt
        @test all(z -> imag(z)>0,d.shapes)
        for (i,t) in enumerate(d.tetrahedra)
            shape = SakumaWeeks.complexfloat(d.shapes[cld(i,2)])
            # Every tetrahedron must have the prescribed shape on some edge;
            # order in Definition 5.2 is not an oriented vertex ordering.
            @test any(e -> isapprox(edgeparameter(t,e),shape;atol=1e-8),edges(t))
        end
    end
    @test SakumaWeeks.complexfloat(geometric_root(s,2//5)) ≈ cis(pi/3)
    @test SakumaWeeks.complexfloat(geometric_root(s,3//8)) ≈ (1+im)/2
    @test real(SakumaWeeks.complexfloat(geometric_root(s,10//33))) ≈ .451069913 atol=1e-8
    for q in (0,1//3,2//5)
        @test Smat(s,q)[1,1]+Smat(s,q)[2,2] == 2-s.field(Qpoly(s,q))^2/s.field(dQ(s,q))
    end
    z=cis(pi/3)
    for M in (Umat(s,0),Wmat(s,0),Smat(s,0))
        A = map(f -> SakumaWeeks.value(f,z),M)
        g = boundary_isometry(A)
        w=.2+.4im
        @test complexcoordinate(g(complexidealpoint(w))) ≈ (A[1,1]*w+A[1,2])/(A[2,1]*w+A[2,2])
        t = first(developed_domain(s,2//5).tetrahedra)
        @test idealshape(g(t)) ≈ idealshape(t)
    end
end
