# Optional integration suite: run in an environment containing Kleinian and Nemo.
# Nemo is a weak dependency, not required by the ordinary test target.
using Test
using Kleinian
import Nemo

@testset "Nemo exact algebraic geometry" begin
    K=Nemo.algebraic_closure(Nemo.QQ)
    e0,e1,e2,e3=(basisvector(i,K(1)) for i in 0:3)
    @test !(eltype(e0)<:Real)
    o=HyperbolicPoint(e0)
    p=HyperbolicPoint(clvector(K(5),K(1),K(2),K(3)))
    L=Geodesic(e0∧e3)
    H=Horosphere(e0+e3)
    h=HalfTurn(L)
    P=Parabolic(H,e1)
    T=Translation(L;multiplier=K(2))
    E=Elliptic(L,7)
    @test isidentity(E^7)
    @test !isidentity(E^1) && !isidentity(E^2)
    @test P(H)==H
    @test T(H)==Horosphere(K(2)*(e0+e3))
    @test scalarpart(representative(p)^2)==K(-11)
    @test representative(p)===p.vector
    @test h(h(p))==p
    @test inv(T)(T(p))==p
    @test P(o) in H
    @test eltype(versor(P))==typeof(K(1))
    for g in (h,P,T,E,T*P*E,Elliptic(Geodesic(e0∧(e1+e2)),5))
        A=Kleinian.matrix(g)
        @test eltype(A)==typeof(K(1))
        M=MobiusIsometry(A)
        @test M(p)==g(p)
        @test M(H)==g(H)
        @test M(L)==g(L)
        @test inv(M)(M(p))==p
        @test VersorIsometry(versor(M))(p)==g(p)
        @test MobiusIsometry(map(x -> K(2+3im)*x,A))(H)==g(H)
    end
    @test isidentity(HalfTurn(Geodesic(e0∧e1))^2)
    @test HyperbolicPoint(clvector(K(10),K(2),K(4),K(6)))==p
    M=MobiusIsometry(K.([1+im 2-im; 3im 1-im]))
    @test VersorIsometry(versor(M))(p)==M(p)
    @test inv(M)(M(p))==p
    @test M(H)==VersorIsometry(versor(M))(H)
end
