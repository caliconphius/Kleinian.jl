@testset "Geodesic polygons and polyhedra" begin
    ps = (IdealPoint(1,0,0,1),IdealPoint(1,0,0,-1),IdealPoint(1,1,0,0),IdealPoint(1,0,1,0))
    t = HyperbolicTetrahedron(ps)
    @test vertices(t) === ps
    @test length(edges(t))==6
    @test all(i -> !(i in faceindices(t)[i]),1:4)
    @test interiorpoint(t) in t
    @test all(p -> p in t,ps)
    @test !(IdealPoint(1,-1,0,0) in t)
    @test all(i -> polygonarea(face(t,i))≈pi,1:4)
    @test all(i -> all(p -> isincident(p,faceplane(t,i)),vertices(face(t,i))),1:4)
    @test all(e -> isinf(edgelength(t,e)),edges(t))
    @test idealshape(t)==im
    G = grammatrix(t)
    @test G≈transpose(G)
    @test all(i -> G[i,i]==1,1:4)
    for e in edges(t)
        @test angle(edgeparameter(t,e))≈dihedralangle(t,e)
        opposite = Tuple(i for i in 1:4 if !(i in e))
        @test edgeparameter(t,e)≈edgeparameter(t,opposite)
    end
    @test prod(edgeparameter(t,(1,j)) for j in 2:4)≈-1
    @test dihedralangle(t,(1,2))≈pi/2
    @test dihedralangle(t,(1,3))≈pi/4
    @test_throws ArgumentError dihedralangle(t,(1,1))
    @test_throws ArgumentError HyperbolicTetrahedron(ps[1:3])
    @test_throws ArgumentError HyperbolicTetrahedron((ps[1],ps[2],ps[3],ps[3]))
    @test_throws ArgumentError HyperbolicPolygon((ps[1],ps[2],ps[1]))
    @test_throws ArgumentError HyperbolicPolygon(ps)

    square = HyperbolicPolygon([modelpoint((x,y,0);model=:klein) for (x,y) in ((-.2,-.2),(.2,-.2),(.2,.2),(-.2,.2))])
    @test 0<polygonarea(square)<2pi
    @test interiorpoint(square) in square
    @test !(modelpoint((0,0,.1);model=:klein) in square)
    @test !(modelpoint((.5,0,0);model=:klein) in square)
    @test_throws ArgumentError HyperbolicPolygon(vertices(square)[[1,3,2,4]])
    @test_throws ArgumentError HyperbolicPolygon(vertices(square);atol=-1)
    @test edgelength(square,(1,2))≈distance(vertices(square)[1],vertices(square)[2])
    @test_throws ArgumentError edgelength(square,(1,3))
    segment = GeodesicSegment(vertices(square)[1],vertices(square)[2])
    @test edgelength(segment)≈edgelength(square,(1,2))
    @test_throws ArgumentError GeodesicSegment(ps[1],ps[1])

    cubevertices = [modelpoint((x,y,z);model=:klein) for z in (-.2,.2) for y in (-.2,.2) for x in (-.2,.2)]
    cubefaces = ((1,2,4,3),(5,6,8,7),(1,2,6,5),(3,4,8,7),(1,3,7,5),(2,4,8,6))
    cube = HyperbolicPolyhedron(cubevertices,cubefaces)
    @test length(edges(cube))==12
    @test length(faceindices(cube))==6
    @test interiorpoint(cube) in cube
    @test_throws ArgumentError HyperbolicPolyhedron(cubevertices,cubefaces[1:5])
    @test_throws ArgumentError HyperbolicTetrahedron(cube)
    @test_throws ArgumentError HyperbolicPolyhedron(cubevertices,((1,2,4,9),cubefaces[2:end]...))
    @test_throws ArgumentError idealshape(HyperbolicTetrahedron(cubevertices[[1,2,3,5]]))
    shift = MobiusIsometry([1 1;0 1])
    @test polygonarea(shift(square))≈polygonarea(square)
    @test grammatrix(shift(t))≈grammatrix(t)
    @test idealshape(shift(t))==idealshape(t)
    @test grammatrix(mapcoefficients(x -> complex(Float64(x)),t))≈grammatrix(t)
    @test shift(interiorpoint(t)) in shift(t)
    @test edgelength(shift(segment))≈edgelength(segment)
end
