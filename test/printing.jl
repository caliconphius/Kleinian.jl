@testset "Multivector display" begin
    a = scalar(3) + basisvector(0,2) - basisblade(3) + basisblade(7,4)
    expected = "(3) + (2)*e0 - (1)*e0e1 + (4)*e0e1e2"
    @test repr(a) == expected
    @test sprint(show, MIME"text/plain"(), a) == expected
    @test sprint(print, a) == expected
    @test repr(scalar(0)) == "0"
    @test repr(basisvector(0)) == "(1)*e0"
    @test repr(basisvector(3,-1)) == "-(1)*e3"
    @test repr(scalar(-2)+basisblade(15)) == "-(2) + (1)*e0e1e2e3"
    @test repr(basisblade(5, -2//3)) == "-(2//3)*e0e2"
    @test repr(basisvector(1, 1+2im)) == "(1 + 2im)*e1"
    @test repr(basisvector(2, typemin(Int))) == "-($(big(typemax(Int))+1))*e2"
    @test repr(basisvector(1, RingCoefficient(2))) == "(RingCoefficient(2))*e1"
    @test repr(basisblade(9)+basisblade(6)+basisvector(3)) ==
        "(1)*e3 + (1)*e0e3 + (1)*e1e2"
end
