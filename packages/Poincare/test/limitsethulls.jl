using LinearAlgebra: norm, cross, dot
using Random

@testset "Sampled limit-set convex hulls" begin
    tetra = [IdealPoint(3,x,y,z) for (x,y,z) in ((1,2,2),(1,-2,-2),(-1,2,-2),(-1,-2,2))]
    h = limitsethull(tetra)
    @test h.dimension == 3
    @test length(vertices(h)) == 4
    @test length(faceindices(h)) == 4
    @test length(edges(h)) == 6
    @test all(p -> any(q -> p === q,tetra),vertices(h))
    @test limitsethull([tetra;tetra]).samplecount == 8
    @test length(vertices(limitsethull([tetra;tetra]))) == 4
    # Independent support, orientation, and Euler characteristic checks.
    function checkhull(h,input)
        ks = h.klein
        @test length(ks)-length(edges(h))+length(faceindices(h)) == 2
        center = sum(ks)/length(ks)
        for (a,b,c) in faceindices(h)
            n = cross(ks[b]-ks[a],ks[c]-ks[a]); n /= norm(n)
            @test dot(n,center-ks[a]) < 0
            @test all(p -> dot(n,modelcoordinates(p;model=:klein)-ks[a]) <= 8h.atol,input)
        end
    end
    checkhull(h,tetra)
    rng = Xoshiro(16)
    points = [modelidealpoint((v=randn(rng,3);v/norm(v))) for _ in 1:100]
    h = limitsethull(points)
    @test length(vertices(h)) == 100
    checkhull(h,points)
    h2 = limitsethull(reverse(points))
    @test length(faceindices(h2)) == length(faceindices(h))
    checkhull(h2,points)
    # Fuchsian circles, including a plane which does not pass through the origin.
    for height in (0.0,0.4)
        ring = [modelidealpoint((sqrt(1-height^2)*cos(t),sqrt(1-height^2)*sin(t),height))
                for t in range(0,2pi;length=25)[1:end-1]]
        h = limitsethull(shuffle(rng,ring))
        @test h.dimension == 2
        @test length(vertices(h)) == 24
        @test length(faceindices(h)) == 22
        @test length(edges(h)) == 24 # no triangulation diagonals
        areas = [norm(cross(h.klein[b]-h.klein[a],h.klein[c]-h.klein[a]))/2 for (a,b,c) in faceindices(h)]
        @test sum(areas) ≈ 12*(1-height^2)*sin(2pi/24)
    end
    p,q = complexidealpoint(0),complexidealpoint(nothing)
    @test limitsethull([p,p]).dimension == 0
    @test isempty(edges(limitsethull([p])))
    @test limitsethull([p,q]).dimension == 1
    @test edges(limitsethull([p,q])) == [(1,2)]
    @test limitsethull((p for p in points[1:20])).samplecount == 20
    @test_throws ArgumentError limitsethull([])
    @test_throws ArgumentError limitsethull(tetra;maxpoints=3)
    @test_throws ArgumentError limitsethull(tetra;maxpoints=0)
    @test_throws ArgumentError limitsethull(tetra;atol=0)
    @test_throws ArgumentError limitsethull(tetra;atol=NaN)
    @test_throws ArgumentError limitsethull([HyperbolicPoint(1,0,0,0)])
end

@testset "Hyperbolic hull bending" begin
    # The ideal cube has six coplanar square faces triangulated by the hull.
    cube = [IdealPoint(3,x,y,z) for x in (-1,1), y in (-2,2), z in (-2,2)] |> vec
    h = limitsethull(cube)
    @test length(edges(h)) == 18
    angles = hullbendingangles(h)
    @test count(r -> r.angle == 0,angles) == 6
    @test length(pleatedges(h)) == 12
    @test isempty(pleatedges(h;minangle=pi))
    # A regular ideal tetrahedron has interior π/3, exterior 2π/3.
    tetra = [modelidealpoint((x,y,z)./sqrt(3)) for (x,y,z) in ((1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1))]
    t = limitsethull(tetra)
    @test all(r -> isapprox(r.angle,2pi/3;atol=1e-12),hullbendingangles(t))
    @test length(pleatedges(t;minangle=pi/2)) == 6
    @test isempty(pleatedges(t;minangle=3pi/4))
    change = MobiusIsometry([1.0 .2+.3im; 0 1])
    transformed = limitsethull(change.(tetra))
    @test all(r -> isapprox(r.angle,2pi/3;atol=1e-10),hullbendingangles(transformed))
    planar = limitsethull([complexidealpoint(x) for x in (-2,-1,1,2)])
    @test isempty(pleatedges(planar))
    @test length(edges(planar)) == 4
    @test_throws ArgumentError pleatedges(t;minangle=-1)
    @test_throws ArgumentError pleatedges(t;minangle=NaN)
    @test_throws ArgumentError hullbendingangles(t;planetol=-1)
end
