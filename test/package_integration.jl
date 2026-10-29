using Test, Kleinian, Poincare, Riley, Random
import Farey
@testset "Separated packages and Farey adapter" begin
    @test !isdefined(Kleinian,:limitset)
    @test !isdefined(Kleinian,:HyperbolicPolyhedron)
    @test !isdefined(Kleinian,:LimitSetHull)
    @test Poincare.IdealPoint === Kleinian.IdealPoint === Riley.IdealPoint
    @test parentmodule(typeof(complexidealpoint(1im))) === Kleinian
    sample = collect(randomlimitset(rileygroup(2.0im);samples=128,rng=Xoshiro(42)))
    hull = limitsethull(sample)
    @test hull isa Poincare.LimitSetHull
    @test hull.dimension == 3
    @test all(p -> p isa Kleinian.IdealPoint,vertices(hull))
    @test Base.get_extension(Riley,:RileyFareyExt) !== nothing
    f,g = Farey.Monoids.FreeGen("f g")
    X,Y = rileygroup(2.0im)
    images = Dict(f=>X,g=>Y)
    @test wordisometry(f*inv(g),images) ≈ X*inv(Y)
    @test isidentity(wordisometry(one(f),images))
    @test isidentity(wordisometry(one(f),Dict();identity=one(X)))
    word = Farey.farey_word(2//5,f,g)
    value = wordisometry(word,images)
    @test wordisometry(inv(word),images) ≈ inv(value)
    @test wordisometry(word*inv(word),images) ≈ one(X)
    p = complexidealpoint(.3+.4im)
    direct = foldr((letter,q) -> (haskey(images,letter) ? images[letter] : inv(images[inv(letter)]))(q),
                   collect(Farey.Monoids.eachgen(word));init=p)
    @test value(p) ≈ direct
    @test_throws KeyError wordisometry(g,Dict(f=>X))
    @test_throws ArgumentError wordisometry(f,Dict(f=>1))
    @test_throws ArgumentError wordisometry(one(f),Dict())
end
