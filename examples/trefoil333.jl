# Load GLMakie (interactive/high-quality transparency) or CairoMakie first.
include("trefoil333_data.jl")

const TREFOIL_EDGE_COLORS = ["#7359b8","#168f89","#d59a28","#c95368"]
const TREFOIL_EDGE_NAMES = ["W cycle · order 3","b cycle · order 3","Regular cycle · order 1","a cycle · order 3"]

function _trefoil_edges(d,report)
    globalids = ((1,2,3,4),(1,2,3,5))
    seen = Set{Tuple{Int,Int}}()
    result = []
    for (class,r) in enumerate(report.edge_cycles),(cell,(i,j)) in r.cycle.sectors
        key = minmax(globalids[cell][i],globalids[cell][j])
        key in seen && continue
        push!(seen,key)
        push!(result,(segment=GeodesicSegment(vertices(d.cells[cell])[i],vertices(d.cells[cell])[j]),class=class))
    end
    result
end

"""Comparison of the same verified trefoil domain in the three models."""
function trefoil333_models(;resolution=28)
    original = trefoil333()
    # Recenter at the midpoint of the reflex edge 0--1: (z,h) ↦ (2z-1,2h).
    change = MobiusIsometry([2.0 -1.0;0.0 1.0])
    d = change(original.domain)
    report = verify_edge_gluings(d;atol=1e-9)
    @assert report.valid
    points = map(change,values(original.points))
    edgeobjects = _trefoil_edges(d,report)
    fig = Figure(size=(1800,760),figure_padding=28,backgroundcolor="#fbfcfe")
    Label(fig[1,1:3],"The trefoil (3,3,3) domain",fontsize=30,font=:bold,color="#17243a")
    Label(fig[2,1:3],"Two ideal tetrahedra · paired faces a, b, c · edge colours identify quotient cycles",
          fontsize=16,color="#52627a")
    for (column,(model,title)) in enumerate(zip((:klein,:ball,:halfspace),
                    ("Klein  |  straight edges","Poincaré ball  |  conformal","Upper half-space")))
        ax = Axis3(fig[3,column];title,aspect=:data,azimuth=1.04,elevation=0.32,
                   perspectiveness=0.15,titlecolor="#17243a",titlesize=21)
        hidedecorations!(ax); hidespines!(ax)
        if model != :halfspace
            modelboundary!(ax;model,resolution=16,color=(:slategray,0.08),linewidth=0.6)
        else
            modelboundary!(ax;model,extent=1.25,color=(:slategray,0.045))
        end
        for (cell,fi) in boundaryfaces(d)
            poly = face(d.cells[cell],fi)
            color = cell==1 ? "#4b9dbb" : "#de9d68"
            hyperbolicplot!(ax,poly;model,resolution,cutoff=1e-5,
                            color=(color,0.22),transparency=true)
        end
        for item in edgeobjects
            hyperbolicplot!(ax,item.segment;model,resolution=180,cutoff=1e-5,
                            color=TREFOIL_EDGE_COLORS[item.class],linewidth=3.8)
        end
        for p in d.pairings
            p.internal && continue
            for ((cell,fi),label) in ((p.source,p.label),(p.target,p.label*"⁻¹"))
                pos = modelcoordinates(interiorpoint(face(d.cells[cell],fi));model)
                text!(ax,[Point3f(pos)];text=[label],fontsize=20,font=:bold,
                      color="#23354c",strokecolor=(:white,0.9),strokewidth=1.5,align=(:center,:center))
            end
        end
        for (i,p) in enumerate(points)
            model==:halfspace && i==1 && continue
            pos = modelcoordinates(p;model)
            scatter!(ax,[Point3f(pos)];color="#23354c",markersize=8)
            text!(ax,[Point3f(pos)];text=[("∞","0","1","u","d")[i]],
                  color="#23354c",fontsize=19,offset=(7,5),font=:bold)
        end
        model==:halfspace ? limits!(ax,-1.3,1.3,-1.05,1.05,0,2.7) : limits!(ax,-1.12,1.12,-1.12,1.12,-1.12,1.12)
    end
    Legend(fig[4,1:3],[LineElement(color=c,linewidth=4) for c in TREFOIL_EDGE_COLORS],
           TREFOIL_EDGE_NAMES;orientation=:horizontal,framevisible=false,labelsize=15)
    Label(fig[5,1:3],"The 0–1 edge has interior angle 4π/3. The upper half-space view is clipped at finite height.",
          color="#52627a",fontsize=14)
    fig
end

"""Computed edge equations and a conformal transverse view of the reflex edge."""
function trefoil333_edge_figure()
    data = trefoil333()
    report = verify_edge_gluings(data.domain;atol=1e-10)
    @assert report.valid
    fig = Figure(size=(1350,720),figure_padding=35,backgroundcolor="#fbfcfe")
    Label(fig[1,1:2],"Four edge cycles, verified in the library",fontsize=28,font=:bold,color="#17243a")
    table = GridLayout(fig[2,1])
    headers = ("Cycle word","Angle sum","Holonomy order","Shape product")
    for (j,title) in enumerate(headers)
        Label(table[1,j],title;font=:bold,fontsize=16,color="#52627a",padding=(10,10,15,15))
    end
    for (i,r) in enumerate(report.edge_cycles)
        vals = (r.word,r.order==1 ? "2π" : "2π/3",string(r.order),r.order==1 ? "1" : "exp(2πi/3)")
        for (j,value) in enumerate(vals)
            Label(table[i+1,j],value;fontsize=17,color=TREFOIL_EDGE_COLORS[i],padding=(10,10,18,18))
        end
    end
    Label(table[6,1:4],"Each cycle fixes its endpoints; the complete isometry product has the stated order.",
          fontsize=14,color="#52627a",word_wrap=true,width=650)
    ax = Axis(fig[2,2];aspect=DataAspect(),title="A transverse section at edge 0–1",titlesize=19)
    hidedecorations!(ax); hidespines!(ax)
    a,b = (dihedralangle(t,(2,3)) for t in data.domain.cells)
    for (lo,hi,col) in ((pi/2-a,pi/2,"#4b9dbb"),(pi/2,pi/2+b,"#de9d68"))
        arc = [Point2f(cos(t),sin(t)) for t in range(lo,hi;length=120)]
        poly!(ax,vcat([Point2f(0,0)],arc);color=(col,0.3))
        lines!(ax,arc;color=col,linewidth=2.5)
    end
    for theta in (pi/2-a,pi/2,pi/2+b)
        lines!(ax,[Point2f(0,0),Point2f(cos(theta),sin(theta))];color="#52627a",linewidth=1.6)
    end
    scatter!(ax,[Point2f(0,0)];color=TREFOIL_EDGE_COLORS[3],markersize=12)
    text!(ax,[Point2f(.36,.27),Point2f(-.36,.27)];text=["2π/3","2π/3"],align=(:center,:center),fontsize=20,color="#23354c")
    text!(ax,[Point2f(0,-.43)];text=["missing sector\n2π/3"],align=(:center,:center),fontsize=15,color="#79869a")
    text!(ax,[Point2f(0,-1.18)];text=["Interior angle: 4π/3"],align=(:center,:center),fontsize=20,color=TREFOIL_EDGE_COLORS[3])
    limits!(ax,-1.15,1.15,-1.35,1.15)
    Label(fig[3,1:2],"Regular edge: 4π/3 + 4(π/6) = 2π    ·    Three singular edge classes: cone angle 2π/3",
          fontsize=18,color="#17243a",padding=(0,0,18,15))
    fig
end

"""Open a standard Makie interactive figure with a model selector (use GLMakie)."""
function trefoil333_explorer(;resolution=24)
    data = trefoil333()
    d = MobiusIsometry([2.0 -1.0;0.0 1.0])(data.domain)
    report = verify_edge_gluings(d;atol=1e-9)
    fig = Figure(size=(1000,850),figure_padding=25)
    Label(fig[1,1],"Trefoil (3,3,3) · drag to rotate, scroll to zoom",fontsize=22)
    menu = Menu(fig[2,1];options=[("Poincaré ball",:ball),("Klein",:klein),("Upper half-space",:halfspace)])
    ax = Axis3(fig[3,1];aspect=:data,azimuth=1.04,elevation=.32)
    hidedecorations!(ax); hidespines!(ax)
    function redraw(model)
        empty!(ax)
        hyperbolicplot!(ax,d;model,resolution,cutoff=1e-5,color=(:steelblue,.2),transparency=true)
        for item in _trefoil_edges(d,report)
            hyperbolicplot!(ax,item.segment;model,color=TREFOIL_EDGE_COLORS[item.class],linewidth=4,resolution=160)
        end
        model==:halfspace ? limits!(ax,-1.3,1.3,-1.05,1.05,0,2.7) : limits!(ax,-1.1,1.1,-1.1,1.1,-1.1,1.1)
    end
    on(redraw,menu.selection)
    redraw(menu.selection[])
    (figure=fig,axis=ax,menu=menu)
end
