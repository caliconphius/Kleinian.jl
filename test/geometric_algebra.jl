# Independent reference: concatenate ordered generator lists and reduce
# adjacent inversions and equal pairs, rather than using bit-mask swaps.
function reference_product(a, b)
    word = vcat([i for i in 1:4 if !iszero(a & (1 << (i-1)))],
                [i for i in 1:4 if !iszero(b & (1 << (i-1)))])
    sign = 1
    while true
        i = findfirst(j -> word[j] >= word[j+1], 1:length(word)-1)
        isnothing(i) && break
        if word[i] == word[i+1]
            sign *= word[i] == 1 ? -1 : 1
            deleteat!(word, i:i+1)
        else
            word[i], word[i+1] = word[i+1], word[i]
            sign = -sign
        end
    end
    mask = sum(i -> 1 << (i-1), word; init=0)
    basisblade(mask, sign)
end

# A coefficient ring that deliberately does not subtype Number.
struct RingCoefficient
    value::Int
end
Base.zero(::Type{RingCoefficient}) = RingCoefficient(0)
Base.zero(::RingCoefficient) = RingCoefficient(0)
Base.one(::Type{RingCoefficient}) = RingCoefficient(1)
Base.one(::RingCoefficient) = RingCoefficient(1)
Base.:+(a::RingCoefficient, b::RingCoefficient) = RingCoefficient(a.value+b.value)
Base.:-(a::RingCoefficient, b::RingCoefficient) = RingCoefficient(a.value-b.value)
Base.:-(a::RingCoefficient) = RingCoefficient(-a.value)
Base.:*(a::RingCoefficient, b::RingCoefficient) = RingCoefficient(a.value*b.value)
Base.iszero(a::RingCoefficient) = iszero(a.value)
Base.:(==)(a::RingCoefficient, b::RingCoefficient) = a.value == b.value

@testset "Cl(3,1) kernel" begin
    blades = [basisblade(m) for m in 0:15]
    e = [basisvector(i) for i in 0:3]
    @testset "Construction and coefficient types" begin
        @test coefficients(blades[1]) isa SVector{16,Int}
        @test Multivector(Tuple(1:16)) == Multivector(collect(1:16))
        @test Multivector{BigFloat}(collect(1:16)) isa Multivector{BigFloat}
        @test coefficient(basisblade(3, 7), 3) == 7
        @test clvector((1,2,3,4)) == e[1]+2e[2]+3e[3]+4e[4]
        @test eltype(clvector(1,2.0,3,4)) == Float64
        @test zero(e[1]) == scalar(0)
        @test one(e[1]) == scalar(1)
        @test isone(one(e[1]))
        @test iszero(zero(e[1]))
        @test_throws ArgumentError basisblade(16)
        @test_throws ArgumentError coefficient(e[1], -1)
        @test_throws ArgumentError basisvector(-1)
        @test_throws ArgumentError grade(e[1], 5)
        @test_throws DimensionMismatch clvector([1,2])
        @test_throws DimensionMismatch Multivector([1,2])
        for T in (Int, Rational{Int}, BigInt, Float32, Float64, BigFloat, ComplexF64, RingCoefficient)
            v = clvector(T(1),T(2),T(3),T(4))
            @test v*v == scalar(T(28))
            @test dual(undual(v)) == v
        end
        @test eltype(e[1] + basisvector(1, 1//2)) == Rational{Int}
        @test eltype(e[1] * basisvector(1, big"1.2")) == BigFloat
        @test hash(e[1]) == hash(basisvector(0,1.0))
        @test isequal(e[1], basisvector(0,1.0))
    end
    @testset "Products and associativity" begin
        for a in 0:15, b in 0:15
            @test blades[a+1]*blades[b+1] == reference_product(a,b)
        end
        for a in blades, b in blades, c in blades
            @test (a*b)*c == a*(b*c)
        end
        for i in 1:4, j in 1:4
            @test e[i]*e[j]+e[j]*e[i] == scalar(i == j ? (i == 1 ? -2 : 2) : 0)
        end
        a, b, c = Multivector(Tuple(-7:8)), Multivector(Tuple(1:16)), clvector(1,3,-2,4)
        @test a*(b+c) == a*b+a*c
        @test (a+b)*c == a*c+b*c
        @test a*one(a) == a == one(a)*a
        @test a^0 == one(a)
        @test a^3 == a*a*a
        @test a+2 == 2+a
        @test 2-a == -(a-2)
        @test 2*a == a*2
        @test (a/2)*2 == a
        @test geometricproduct(a,b) == a*b
    end
    @testset "Grades, exterior product, and contraction" begin
        a = Multivector(Tuple(1:16))
        @test sum(grade(a,r) for r in 0:4) == a
        for x in e, y in e
            @test x ∧ y == wedge(x,y) == -wedge(y,x)
            @test wedge(x,y) == (x*y-y*x)/2
            @test leftcontract(x,y) == grade(x*y,0)
        end
        for x in blades, y in blades
            r = count_ones(findfirst(!iszero, coefficients(x))-1)
            s = count_ones(findfirst(!iszero, coefficients(y))-1)
            @test wedge(x,y) == (r+s <= 4 ? grade(x*y,r+s) : zero(x*y))
            @test leftcontract(x,y) == (r <= s ? grade(x*y,s-r) : zero(x*y))
        end
        @test leftcontract(e[1],e[1]∧e[2]) == -e[2]
        @test leftcontract(e[2],e[2]∧e[3]) == e[3]
        @test leftcontract(e[1]∧e[2],e[1]) == zero(a)
        @test wedge(a,a+e[1]) == wedge(a,a)+wedge(a,e[1])
        @test leftcontract(a+e[1],a) == leftcontract(a,a)+leftcontract(e[1],a)
    end
    @testset "Involutions and duals" begin
        a, b = Multivector(Tuple(1:16)), Multivector(Tuple(-7:8))
        @test reverse(reverse(a)) == a
        @test reverse(a*b) == reverse(b)*reverse(a)
        @test gradeinvolution(gradeinvolution(a)) == a
        @test gradeinvolution(a*b) == gradeinvolution(a)*gradeinvolution(b)
        @test cliffordconjugate(cliffordconjugate(a)) == a
        @test cliffordconjugate(a*b) == cliffordconjugate(b)*cliffordconjugate(a)
        @test basisblade(15)^2 == scalar(-1)
        @test dual(dual(a)) == -a
        @test undual(dual(a)) == a == dual(undual(a))
        @test reverse(e[1]∧e[2]) == -(e[1]∧e[2])
        @test cliffordconjugate(scalar(1+2im)) == scalar(1+2im)
    end
    @testset "Restricted inverse and reflection" begin
        v = clvector(5//1,1//1,2//1,3//1)
        @test inv(v)*v == one(v) == v*inv(v)
        @test inv(v) == checked_inverse(v)
        @test v^-2 == inv(v)*inv(v)
        @test e[1]^-1 == -e[1]
        @test copy(v) === v
        @test v^(-3) == inv(v)*inv(v)*inv(v)
        @test e[1]^(typemin(Int8)) == one(e[1])
        r = scalar(5//3) + (4//3)*(e[1]*e[4])
        @test r*inv(r) == one(r)
        @test inv(r) == reverse(r)
        @test_throws DomainError inv(e[1]+e[4])
        @test_throws DomainError inv(zero(v))
        @test_throws ArgumentError inv(scalar(1)+e[2])
        @test iszero((scalar(1)+e[2])*(scalar(1)-e[2]))
        n = clvector(0//1,1//1,2//1,0//1)
        reflect(x) = -n*x*inv(n)
        expected = v - 2*scalarpart(leftcontract(v,n))/scalarpart(n*n)*n
        @test reflect(v) == expected
        @test reflect(reflect(v)) == v
        @test reflect(v)^2 == v^2
        @test reflect(e[4]) == e[4]
        @test reflect(n) == -n
        vf = clvector(5.0,1.0,2.0,3.0)
        @test vf*inv(vf) ≈ one(vf)
    end
end
