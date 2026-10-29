@testset "Complex ideal boundary and reduced orbits" begin
    for z in (0,1,im,2//3+im*(4//5),0.3+0.7im,1e150+2e150im)
        p = complexidealpoint(z)
        @test complexcoordinate(p) ≈ z
        @test complexcoordinate(IdealPoint(3representative(p))) ≈ z
    end
    @test complexcoordinate(complexidealpoint(nothing)) === nothing
    @test_throws ArgumentError complexidealpoint(Inf)
    @test_throws ArgumentError complexidealpoint(NaN)
    z = 2//3+im*(4//5)
    p = complexidealpoint(z)
    g = MobiusIsometry([1+im 2; 1 3-im])
    A = conj.(matrix(g))
    @test complexcoordinate(g(p)) ≈ (A[1,1]*z+A[1,2])/(A[2,1]*z+A[2,2])
    @test crossratio(complexidealpoint(nothing),complexidealpoint(0),complexidealpoint(1),p) ≈ z
    h = MobiusIsometry([2.0 0; 0 1])
    seed = complexidealpoint(1.0)
    @test limitset([h];seed,depth=0) == [seed]
    @test sort(real.(complexcoordinate.(limitset([h];seed,depth=3)))) ≈ [1/8,8]
    @test length(limitset([g,h];depth=3)) == 4*3^2
    @test_throws ArgumentError limitset([g,h];depth=20,maxpoints=100)
    @test_throws ArgumentError limitset([h];depth=-1)
    @test_throws ArgumentError limitset([])
end

using Random

# Observe chosen letters without involving numerical group relations.
struct RecordedLetter <: AbstractIsometry
    label::Int
    log::Vector{Int}
end
Base.inv(g::RecordedLetter) = RecordedLetter(-g.label,g.log)
function (g::RecordedLetter)(p::IdealPoint)
    push!(g.log,g.label)
    p
end

@testset "Streaming and random reduced words" begin
    gs = rileygroup(2.0im)
    seed = complexidealpoint(0.3+0.2im)
    # Independent breadth-first oracle also checks preservation of public order.
    level = [(seed,0)]
    actions = [gs...,inv.(gs)...]
    for depth in 0:4
        stream = limitsetstream(gs;seed,depth)
        @test all(isapprox.(collect(stream),first.(level)))
        @test collect(stream) == collect(stream)
        @test length(stream) == length(level)
        level = [(actions[j](p),j) for (p,last) in level for j in 1:4
                 if last == 0 || abs(last-j) != 2]
    end
    @test length(limitsetstream(gs;depth=12)) == 4*3^11
    @test Base.summarysize(limitsetstream(gs;depth=12)) < 10_000
    @test_throws ArgumentError limitsetstream(gs;depth=13)
    @test_throws ArgumentError limitsetstream(gs;depth=-1)
    @test_throws ArgumentError limitsetstream(gs;maxpoints=0)
    @test isempty(collect(randomlimitset(gs;samples=0)))
    rng = Xoshiro(123)
    original = copy(rng)
    s = randomlimitset(gs;rng,seed,samples=100,burnin=20,thin=3)
    @test rand(rng,UInt64) == rand(original,UInt64)
    @test collect(s) == collect(s)
    @test collect(s) == collect(randomlimitset(gs;rng=Xoshiro(123),seed,samples=100,burnin=20,thin=3))
    @test collect(s) != collect(randomlimitset(gs;rng=Xoshiro(124),seed,samples=100,burnin=20,thin=3))
    @test length(s) == 100
    @test all(p -> p isa IdealPoint,s)
    log = Int[]
    recorded = [RecordedLetter(1,log),RecordedLetter(2,log)]
    collect(randomlimitset(recorded;rng=Xoshiro(7),samples=50,burnin=17,thin=3))
    @test length(log) == 17+50*3
    @test all(i -> log[i] != -log[i-1],2:length(log))
    @test Set(log) == Set([-2,-1,1,2])
    @test_throws ArgumentError randomlimitset(gs;samples=-1)
    @test_throws ArgumentError randomlimitset(gs;burnin=-1)
    @test_throws ArgumentError randomlimitset(gs;thin=0)
    @test_throws ArgumentError randomlimitset([])
end

@testset "Riley family and cusp words" begin
    for ρ in (2im,4,2.2im,complex(big"0.5",big"2.0"))
        X,Y = rileygroup(ρ)
        z = real(ρ) isa BigFloat ? complex(big"0.4",big"0.3") : 0.4+0.3im
        p = complexidealpoint(z)
        @test complexcoordinate(X(p)) ≈ z+1
        @test complexcoordinate(Y(p)) ≈ z/(ρ*z+1)
        @test X(complexidealpoint(nothing)) ≈ complexidealpoint(nothing)
        @test Y(complexidealpoint(0)) ≈ complexidealpoint(0)
    end
    X,Y = rileygroup(2im)
    C = matrix(X*Y*inv(X)*inv(Y))
    @test C[1,1]+C[2,2] == -2
    @test !isidentity(X*Y*inv(X)*inv(Y))
    X,Y = rileygroup(4)
    C = matrix(X*inv(Y))
    @test C[1,1]+C[2,2] == -2
    @test !isidentity(X*inv(Y))
    @test_throws ArgumentError rileygroup(0)
    @test_throws ArgumentError rileygroup(Inf)
    @test_throws ArgumentError rileygroup(NaN)
end
