# Load GLMakie for interaction, or CairoMakie for static exports, first.
using Kleinian, Riley
using Random

"""
    riley_fractal_viewer(; samples=1_000_000, previewpoints=10_000,
                         rngseed=2026, background=true)

Fractal-only Riley explorer with 100k–10m sample presets. Click/drag the
parameter; release refines. Generation yields in batches so newer requests
can cancel obsolete work. Only compact Point2d coordinates are retained.
`background=false` makes updates synchronous for scripts/tests. For background
exports, `wait(view.task[])` before saving. Returns controls and
`set_parameter!(ρ; preview=false)` for scripted updates.
"""
function riley_fractal_viewer(;samples::Integer=1_000_000,previewpoints::Integer=10_000,
                              rngseed=2026,background=true)
    1 <= samples <= 10_000_000 || throw(ArgumentError("samples must be between 1 and 10 million"))
    previewpoints > 0 || throw(ArgumentError("previewpoints must be positive"))
    budgets = sort!(unique!([100_000,250_000,500_000,1_000_000,2_000_000,5_000_000,10_000_000,Int(samples)]))
    fig = Figure(size=(1600,950),figure_padding=22)
    Label(fig[1,1:2],"Riley fractals · drag the parameter, release for detail",fontsize=25)
    pax = Axis(fig[2,1];title="Parameter ρ",xlabel="Re(ρ)",ylabel="Im(ρ)",
               aspect=DataAspect(),limits=(0,4.5,0,3.5),xpanlock=true,ypanlock=true,xzoomlock=true,yzoomlock=true)
    hlines!(pax,[0];color=(:gray,.3)); vlines!(pax,[0];color=(:gray,.3))
    scatter!(pax,[Point2f(4,0),Point2f(0,2)];color=:gray,markersize=8)
    text!(pax,[Point2f(2.5,.3),Point2f(.2,2.3)];
          text=["real cusp","2i cusp"],fontsize=12)
    parameter = Observable(2.0im)
    scatter!(pax,lift(ρ -> [Point2f(real(ρ),imag(ρ))],parameter);color=:crimson,markersize=18)
    ax = Axis(fig[2,2];title="Boundary orbit · w = 1/(z+i)",aspect=DataAspect(),
              backgroundcolor=:black,limits=(-1.1,1.1,-2.1,.1))
    hidedecorations!(ax); hidespines!(ax)
    colsize!(fig.layout,1,Relative(.27))
    controls = GridLayout(fig[3,1:2])
    Label(controls[1,1],"Refined sample count")
    countmenu = Menu(controls[1,2];options=[(string(n),n) for n in budgets],default=findfirst(==(samples),budgets))
    Label(controls[1,3],"Colours")
    colors = Menu(controls[1,4];options=[("Inferno",:inferno),("Magma",:magma),("Viridis",:viridis),("Gray",:grays)])
    reset = Button(controls[1,5];label="Reset view")
    xy = Observable(Point2d[])
    valid = Observable(false)
    plot = datashader!(ax,xy;async=false,colormap=colors.selection,local_operation=log1p,visible=valid)
    status = Observable("Preparing preview…")
    Label(fig[4,1:2],status;fontsize=14)
    Label(fig[5,1:2],"Scroll to zoom, drag to pan the fractal. Finite orbit samples; arbitrary parameters need not be discrete.";
          color=:gray40,fontsize=13)
    generation = Ref(0)
    task = Ref{Union{Nothing,Task}}(nothing)
    rendered = Ref((parameter=0.0im,samples=0,preview=true))
    function refresh(ρ;preview=false)
        isfinite(real(ρ)) && isfinite(imag(ρ)) && real(ρ) >= 0 && imag(ρ) >= 0 ||
            throw(ArgumentError("viewer parameter must have finite nonnegative real and imaginary parts"))
        generation[] += 1
        ticket = generation[]
        parameter[] = ComplexF64(ρ)
        target = parameter[]
        n = preview ? min(Int(countmenu.selection[]),Int(previewpoints)) : Int(countmenu.selection[])
        # Hide a previous parameter's geometry until this request completes.
        valid[] = rendered[].parameter == target && !isempty(xy[])
        status[] = "ρ = $(round(target;digits=4)) · generating $n samples…"
        function generate()
            ticket == generation[] || return false
            try
                chart = MobiusIsometry(conj.([0.0 1.0;1.0 im]))
                generators = map(g -> MobiusIsometry(matrix(chart*g*inv(chart))),rileygroup(target))
                stream = randomlimitset(generators;samples=n,burnin=1000,
                    seed=chart(complexidealpoint(0.0)),rng=Xoshiro(rngseed))
                points = Vector{Point2d}(undef,n)
                kept = 0
                for (i,p) in enumerate(stream)
                    z = complexcoordinate(p)
                    if !isnothing(z)
                        q = Point2d(real(z),imag(z))
                        if all(isfinite,q)
                            kept += 1
                            points[kept] = q
                        end
                    end
                    if i % 10_000 == 0
                        background && yield()
                        ticket == generation[] || return false
                        if i % 100_000 == 0
                            status[] = "ρ = $(round(target;digits=4)) · generating $i / $n samples…"
                        end
                    end
                end
                ticket == generation[] || return false
                kept > 0 || throw(ArgumentError("no finite points in the sample"))
                resize!(points,kept)
                xy[] = points
                ticket == generation[] || return false
                rendered[] = (parameter=target,samples=n,preview=preview)
                valid[] = true
                phase = preview ? "preview" : "refined"
                status[] = "ρ = $(round(target;digits=4)) · $phase · $n samples · $kept finite points"
                true
            catch err
                err isa Union{ArgumentError,DomainError} || rethrow()
                ticket == generation[] || return false
                valid[] = false
                xy[] = Point2d[]
                status[] = "ρ = $(round(target;digits=4)) · $(sprint(showerror,err))"
                false
            end
        end
        if background
            task[] = @async generate()
            return task[]
        end
        generate()
    end
    on(countmenu.selection) do _
        refresh(parameter[])
    end
    on(reset.clicks) do _
        limits!(ax,-1.1,1.1,-2.1,.1)
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
    # Closing the window cancels any unfinished sample generation.
    on(events(fig).window_open) do open
        open || (generation[] += 1)
    end
    refresh(2.0im)
    (;figure=fig,parameteraxis=pax,axis=ax,countmenu,colors,reset,parameter,status,
      valid,points=xy,rendered,task,plot,set_parameter! = refresh)
end
