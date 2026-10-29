# Load GLMakie for interaction or CairoMakie for headless exports first.
using Kleinian, Riley, Poincare
using Random

function riley_hull_sample(ρ;samples=30_000,hullsamples=256,rngseed=2026)
    chart = MobiusIsometry(conj.([0.0 1.0;1.0 im]))
    generators = map(g -> MobiusIsometry(matrix(chart*g*inv(chart))),rileygroup(ρ))
    stream = randomlimitset(generators;samples=max(samples,hullsamples),burnin=500,
        seed=chart(complexidealpoint(0.0)),rng=Xoshiro(rngseed))
    boundary = Point3d[]
    complexpoints = Point2d[]
    hullpoints = IdealPoint[]
    for (i,p) in enumerate(stream)
        i <= hullsamples && push!(hullpoints,p)
        if i <= samples
            push!(boundary,Point3d(modelcoordinates(p;model=:klein)))
            z = complexcoordinate(p)
            if !isnothing(z)
                q = Point2d(real(z),imag(z))
                all(isfinite,q) && push!(complexpoints,q)
            end
        end
    end
    (;hull=limitsethull(hullpoints),boundary,complexpoints)
end

"""
    riley_hull_explorer(; samples=30_000, hullsamples=256, rngseed=2026)

Click/drag the parameter plane to change ρ. During a drag use up to 3,000
orbit points and 64 hull samples, throttled to at most 8 updates/second;
release refines to the selected sample budget. Fixed RNG seed across updates.
Returns controls and `set_parameter!(ρ; preview=false)` for scripted updates.
This explores the two-parabolic family, not a certified map of the Riley slice.
"""
function riley_hull_explorer(;samples=30_000,hullsamples=256,rngseed=2026)
    samples >= 1 || throw(ArgumentError("samples must be positive"))
    budgets = [32,64,128,256,512,1024]
    hullsamples in budgets || throw(ArgumentError("hullsamples must be one of $budgets"))
    initial = riley_hull_sample(2.0im;samples,hullsamples,rngseed)
    fig = Figure(size=(1600,780),figure_padding=22)
    Label(fig[1,1:3],"Riley family · drag the parameter to explore",fontsize=25)
    pax = Axis(fig[2,1];title="Parameter ρ · click or drag",xlabel="Re(ρ)",ylabel="Im(ρ)",
               aspect=DataAspect(),limits=(0,4.5,0,3.5),xpanlock=true,ypanlock=true,xzoomlock=true,yzoomlock=true)
    hlines!(pax,[0];color=(:gray,.3)); vlines!(pax,[0];color=(:gray,.3))
    scatter!(pax,[Point2f(4,0),Point2f(0,2)];color=:gray,markersize=8)
    text!(pax,[Point2f(2.8,.3),Point2f(.2,2.3)];
          text=["real cusp","2i cusp"],fontsize=13)
    parameter = Observable(2.0im)
    marker = lift(ρ -> [Point2f(real(ρ),imag(ρ))],parameter)
    scatter!(pax,marker;color=:crimson,markersize=18)
    cax = Axis(fig[2,2];title="Boundary orbit · w = 1/(z+i)",xlabel="Re(w)",ylabel="Im(w)",
               aspect=DataAspect(),limits=(-1.1,1.1,-2.1,.1))
    hax = Axis3(fig[2,3];title="Sampled convex hull · Klein ball",aspect=:equal,
                azimuth=1.0,elevation=.35,limits=(-1.1,1.1,-1.1,1.1,-1.1,1.1))
    hidedecorations!(hax); hidespines!(hax)
    modelmenu = Menu(fig[3,1:3];options=[("Klein",:klein),("Poincaré ball",:ball),("Upper half-space",:halfspace)])
    xy = Observable(initial.complexpoints)
    xyz = Observable(initial.boundary)
    meshdata = Observable(plotgeometry(initial.hull))
    minbend = Observable(0.0)
    edgepoints(h) = hulledgesgeometry(h;model=modelmenu.selection[],resolution=32,
                                    pleats=true,minangle=deg2rad(minbend[]))
    linepoints = Observable(edgepoints(initial.hull))
    valid = Observable(true)
    cloud2 = datashader!(cax,xy;async=false,colormap=:inferno,local_operation=log1p,visible=valid)
    cloud3 = scatter!(hax,xyz;color=(:darkorange,.55),markersize=2,visible=valid)
    controls = GridLayout(fig[4,1:3])
    Label(controls[1,1],"Hull sample budget")
    slider = Slider(controls[1,2];range=budgets,startvalue=hullsamples)
    budgetlabel = Label(controls[1,3],lift(n -> string(n),slider.value))
    faces = Toggle(controls[1,4];active=true)
    Label(controls[1,5],"Faces")
    wire = Toggle(controls[1,6];active=true)
    Label(controls[1,7],"Pleats")
    Label(controls[2,1],"Minimum bend")
    bendslider = Slider(controls[2,2:6];range=0:1:180,startvalue=0)
    Label(controls[2,7],lift(x -> "$(Int(x))°",bendslider.value))
    on(bendslider.value) do angle
        minbend[] = Float64(angle)
        linepoints[] = edgepoints(state[].hull)
    end
    surface = mesh!(hax,meshdata;color=(:steelblue,.25),transparency=true,
                    visible=lift((v,s)->v&&s,valid,faces.active))
    lines = lines!(hax,linepoints;color=:black,linewidth=1.2,
                         visible=lift((v,s)->v&&s,valid,wire.active))
    state = Ref(initial)
    sampledparameter = Ref(2.0im)
    status = Observable("ρ = 2i · $(initial.hull.dimension)D hull · $(length(vertices(initial.hull))) vertices")
    Label(fig[5,1:3],status;fontsize=14)
    Label(fig[6,1:3],"Finite orbit samples; arbitrary parameters need not be discrete. Release to refine. Hulls use a smaller sample.";
          fontsize=13,color=:gray40)
    function updategeometry(data)
        model = modelmenu.selection[]
        geometry = plotgeometry(data.hull;model,resolution=8)
        boundary = Point3d[]
        for k in data.boundary
            p = modelidealpoint(k;model=:klein)
            try
                push!(boundary,Point3d(modelcoordinates(p;model)))
            catch err
                # Half-space infinity has no finite point marker.
                model === :halfspace && err isa DomainError || rethrow()
            end
        end
        meshdata[] = geometry
        xyz[] = boundary
        linepoints[] = edgepoints(data.hull)
        titles = Dict(:klein=>"Klein ball",:ball=>"Poincaré ball",:halfspace=>"Upper half-space")
        hax.title[] = "Sampled convex hull · $(titles[model])"
        model === :halfspace ? limits!(hax,-1.2,1.2,-2.2,.2,0,2.4) :
                              limits!(hax,-1.1,1.1,-1.1,1.1,-1.1,1.1)
    end
    on(modelmenu.selection) do _
        parameter[] == sampledparameter[] || return
        try
            updategeometry(state[])
            valid[] = true
        catch err
            err isa Union{ArgumentError,DomainError} || rethrow()
            valid[] = false
            status[] = sprint(showerror,err)
        end
    end
    function refresh(ρ;preview=false)
        isfinite(real(ρ)) && isfinite(imag(ρ)) && real(ρ) >= 0 && imag(ρ) >= 0 ||
            throw(ArgumentError("viewer parameter must have finite nonnegative real and imaginary parts"))
        parameter[] = ComplexF64(ρ)
        n = preview ? min(samples,3000) : samples
        budget = preview ? min(Int(slider.value[]),64) : Int(slider.value[])
        try
            data = riley_hull_sample(parameter[];samples=n,hullsamples=budget,rngseed)
            data.hull.dimension >= 2 || throw(ArgumentError("sample has no two-dimensional hull surface"))
            updategeometry(data)
            xy[] = data.complexpoints
            state[] = data
            sampledparameter[] = parameter[]
            valid[] = true
            phase = preview ? "preview" : "refined"
            status[] = "ρ = $(round(parameter[];digits=4)) · $phase · $n orbit points · $(data.hull.dimension)D hull · $(length(vertices(data.hull))) vertices"
            return true
        catch err
            err isa Union{ArgumentError,DomainError} || rethrow()
            valid[] = false
            status[] = "ρ = $(round(parameter[];digits=4)) · no current geometry: $(sprint(showerror,err))"
            return false
        end
    end
    on(slider.value) do _
        refresh(parameter[])
    end
    dragging = Ref(false)
    lastupdate = Ref(0.0)
    function mouseupdate(preview)
        q = mouseposition(pax)
        refresh(complex(max(0.0,q[1]),max(0.0,q[2]));preview)
    end
    on(events(fig).mousebutton;priority=100) do event
        if event.button == Mouse.left
            if event.action == Mouse.press && is_mouseinside(pax.scene)
                dragging[] = true
                lastupdate[] = time()
                mouseupdate(true)
                return Consume(true)
            elseif event.action == Mouse.release && dragging[]
                dragging[] = false
                is_mouseinside(pax.scene) ? mouseupdate(false) : refresh(parameter[])
                return Consume(true)
            end
        end
        Consume(false)
    end
    on(events(fig).mouseposition;priority=100) do _
        if dragging[]
            if is_mouseinside(pax.scene) && time()-lastupdate[] >= .125
                mouseupdate(true)
                lastupdate[] = time()
            end
            return Consume(true)
        end
        Consume(false)
    end
    (;figure=fig,parameteraxis=pax,complexaxis=cax,hullaxis=hax,modelmenu,slider,bendslider,faces,wire,
      parameter,status,valid,state,set_parameter! = refresh,plots=(cloud2,cloud3,surface,lines))
end
