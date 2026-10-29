@testset "Symmetric dot product" begin
    e0, e1, e2, e3 = (basisvector(i, 1//1) for i in 0:3)
    @test e0 ⋅ e0 == scalar(-1)
    @test e1 ⋅ e1 == scalar(1)
    @test e0 ⋅ e1 == scalar(0)
    @test dot(e1,e2) == e1 ⋅ e2
    v, w = clvector(1//1,2//1,3//1,4//1), clvector(5//1,6//1,7//1,8//1)
    @test v ⋅ w == scalar(-5+12+21+32)
    @test v*w == (v ⋅ w) + (v ∧ w)
    a, b = scalar(2//1)+e0+e1∧e2, e3+e0∧e1
    @test a ⋅ b == (a*b+b*a)/2
    @test a ⋅ b == b ⋅ a
    @test (a+v) ⋅ b == a ⋅ b + v ⋅ b
    @test (2a) ⋅ b == 2*(a ⋅ b)
    @test scalar(2//1) ⋅ (e1∧e2) == 2*(e1∧e2)
    @test e0 ⋅ (e0∧e1) == scalar(0)
    @test basisvector(1,1+2im) ⋅ e1 == scalar(1+2im)
    @test eltype(basisvector(0, big"1.5") ⋅ e0) == BigFloat
end
