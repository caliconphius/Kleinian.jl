# Load GLMakie (interactive) or CairoMakie (headless) before including this file.
using Kleinian, Poincare
isdefined(@__MODULE__, :SakumaWeeks) || include("SakumaWeeks.jl")
using .SakumaWeeks

function sakuma_weeks_viewer(q=2//5; parameter=nothing, resolution=24)
    s = PolynomialSystem()
    data = isnothing(parameter) ? developed_domain(s,q) : developed_domain(s,q;parameter)
    fig = Makie.Figure(size=(1250,850),backgroundcolor=:white)
    z = SakumaWeeks.complexfloat(data.parameter)
    Makie.Label(fig[1,1:2],"Sakuma–Weeks • $(data.slope) • $(length(data.tetrahedra)) tetrahedra";
        fontsize=25,tellwidth=false)
    Makie.Label(fig[2,1:2],"x = $(round(z;digits=7))  |  " *
        (data.complete ? "Q(x) = 0" : "deformed positive-shape domain") *
        "  |  developed geometry; face-pairing certificate pending";fontsize=14)
    controls = Makie.GridLayout(fig[3,1:2])
    Makie.Label(controls[1,1],"Model")
    modelmenu = Makie.Menu(controls[1,2];options=["Poincaré ball"=>:ball,"Klein"=>:klein,"Upper half-space"=>:halfspace])
    Makie.Label(controls[1,3],"Transform domain")
    transformmenu = Makie.Menu(controls[1,4];options=["Original"=>:identity,"U₀"=>:U,"W₀"=>:W,"S₀"=>:S])
    Makie.Label(controls[1,5],"Layers shown")
    layers = Makie.Slider(controls[1,6];range=1:length(data.edges),startvalue=length(data.edges))
    axf = Makie.Axis(fig[4,1];title="Funnel and tetrahedron pairs",xlabel="slope",ylabel="1 / denominator")
    for (j,e) in enumerate(funnel(q))
        c = j in (1,length(funnel(q))) ? :gray : :steelblue
        Makie.lines!(axf,Float64.([e.qL,e.qR]),[1/Float64(denominator(e.qL)),1/Float64(denominator(e.qR))];color=c,linewidth=2)
    end
    points = unique(vcat([[e.qL,e.qR] for e in funnel(q)]...))
    Makie.scatter!(axf,Float64.(points),1 ./ Float64.(denominator.(points));color=:steelblue,markersize=8)
    Makie.text!(axf,Float64.(points),1 ./ Float64.(denominator.(points));text=string.(points),offset=(5,5),fontsize=12)
    ax = Makie.Axis3(fig[4,2];aspect=:data,title="Developed tetrahedra",azimuth=1.15,elevation=.45)
    palette = Makie.cgrad(:viridis)
    function redraw(_=nothing)
        Makie.empty!(ax)
        model = modelmenu.selection[]
        which = transformmenu.selection[]
        g = if which == :identity
            nothing
        else
            A = which == :U ? Umat(s,0) : which == :W ? Wmat(s,0) : Smat(s,0)
            boundary_isometry(map(f -> SakumaWeeks.complexfloat(SakumaWeeks.value(f,data.parameter)),A))
        end
        for (i,t0) in enumerate(data.tetrahedra[1:2*layers.value[]])
            t = isnothing(g) ? t0 : g(t0)
            c = palette[(cld(i,2)-.5)/length(data.edges)]
            hyperbolicplot!(ax,t;model,resolution,color=(c,.28))
            for (a,b) in edges(t)
                hyperbolicplot!(ax,GeodesicSegment(vertices(t)[a],vertices(t)[b]);model,
                    resolution=max(64,resolution),color=(:black,.75),linewidth=1.3)
            end
        end
        if model == :halfspace
            Makie.limits!(ax,-1.5,2.5,-2,2,0,3)
        else
            Makie.limits!(ax,-1.05,1.05,-1.05,1.05,-1.05,1.05)
        end
    end
    Makie.on(redraw,modelmenu.selection)
    Makie.on(redraw,transformmenu.selection)
    Makie.on(redraw,layers.value)
    redraw()
    (;figure=fig,axis=ax,funnel_axis=axf,modelmenu,transformmenu,layers,data)
end
