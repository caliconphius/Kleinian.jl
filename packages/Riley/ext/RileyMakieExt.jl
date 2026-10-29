module RileyMakieExt
import Kleinian as K
import Riley as R
import Makie
import GeometryBasics as GB
function R.limitsetplot!(axis, points::Union{AbstractVector{<:K.IdealPoint},R.LimitSetStream,R.RandomLimitSetStream};
                         renderer=:datashader, kwargs...)
    renderer in (:datashader,:scatter) ||
        throw(ArgumentError("renderer must be :datashader or :scatter"))
    vertices = GB.Point2d[]
    for p in points
        z = K.complexcoordinate(p)
        isnothing(z) && continue
        v = GB.Point2d(Float64(real(z)),Float64(imag(z)))
        all(isfinite,v) || throw(DomainError(z,"point is outside rendering precision"))
        push!(vertices,v)
    end
    isempty(vertices) && throw(ArgumentError("no finite boundary points to plot"))
    if renderer === :datashader
        options = merge((;async=false,colormap=:inferno), (;kwargs...))
        return Makie.datashader!(axis,vertices;options...)
    end
    Makie.scatter!(axis,vertices;kwargs...)
end

function R.limitsetplot(points::Union{AbstractVector{<:K.IdealPoint},R.LimitSetStream,R.RandomLimitSetStream};
                        figure=(;),axis=(;),kwargs...)
    fig = Makie.Figure(;figure...)
    options = merge((;aspect=Makie.DataAspect(),xlabel="Re(z)",ylabel="Im(z)"),axis)
    ax = Makie.Axis(fig[1,1];options...)
    plot = R.limitsetplot!(ax,points;kwargs...)
    (figure=fig,axis=ax,plot=plot)
end


end
