using Test, Kleinian, Poincare
import Nemo
include(joinpath(@__DIR__,"..","examples","trefoil333_data.jl"))

@testset "Exact algebraic trefoil edge equations" begin
    K = Nemo.algebraic_closure(Nemo.QQ)
    data = trefoil333(root3=sqrt(K(3)))
    report = verify_edge_gluings(data.domain;atol=0)
    @test report.valid
    @test sort([e.order for e in report.edge_cycles])==[1,3,3,3]
    omega = (-K(1)+sqrt(K(-3)))/2
    for e in report.edge_cycles
        @test e.endpoints_fixed
        @test e.holonomy_order_correct
        @test e.shape_equation
        @test iszero(e.shape_residual)
        @test e.cycle.shape_product==(e.order==1 ? K(1) : omega)
        @test isidentity(e.cycle.holonomy^e.order)
    end
    @test all(isidentity,data.relations)
    a,b,c = data.pairings
    f,g = data.generators
    @test b==f*g*f*inv(g)*inv(f)
    @test c==f*g
    @test isidentity(f^3) && isidentity(g^3)
    @test isidentity((f*g*f*inv(g)*inv(f)*inv(g))^3)
    @test !isidentity(a) && !isidentity(a^2)
    @test !isidentity(b) && !isidentity(b^2)
    @test !isidentity(b*inv(a)) && !isidentity((b*inv(a))^2)
    println("Exact shape products: ",[e.cycle.shape_product for e in report.edge_cycles])
    println("Exact holonomy orders: ",[e.order for e in report.edge_cycles])
end
