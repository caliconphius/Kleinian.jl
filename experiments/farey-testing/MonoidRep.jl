module MonoidReps
  using Farey
  using Farey:AbstractElement, AbstractMonoidWord, AbstractMonoidElement, AbstractMonoidGen
  using Farey.Monoids:MonoidWord
  using Farey:eachgen

  struct MonoidRep{T}<:AbstractElement
    map::AbstractDict
    id::T
  end

  function MonoidRep(m::Dict{K,V}; identity = nothing, withinvs=false) where {K,V}
    identity = isnothing(identity) ? one(V) : identity

    if withinvs
      merge!(m, Dict([-k=>inv(v) for (k,v) in pairs(m)]))
    end



    MonoidRep{V}(m, identity)

  end

  struct RepElement{T}<:AbstractElement
    word::AbstractMonoidWord
    parent::MonoidRep{T}
    rep::T
  end


  function monoideval(x::AbstractMonoidElement, di::Dict{K,V}; init=one(first(values(di))))::V where {K,V}
    mapreduce(*, eachgen(x); init=init) do x
      di[x]
    end
  end


  function Base.getindex(x::RepElement{T}, i::S) where {T,S}
    x.parent(x.word[i])
  end

  function (L::MonoidRep{T})(w::S)::RepElement{T} where {T, S<:AbstractMonoidElement}
    rep = monoideval(w, L.map; init=L.id)
    w = Farey.Monoids.MonoidWord(w)

    RepElement{T}(w, L, rep)
  end

  Base.one(s::RepElement{T}) where T = s.parent(one(s.word))

  function Base.:*(a::RepElement{T}, b::RepElement{T})::RepElement{T} where T
    a.parent(a.word*b.word)
  end

  Base.inv(a::RepElement{T}) where {T} = a.parent(-a.word)

  function Base.show(io::IO, x::RepElement{T}) where T
    print(io, "$(x.word) ↦ $(x.rep)")

  end

  function Base.show(io::IO, mime::MIME"text/plain", x::RepElement{T}) where T
    println(io, "Monoid Represenation")
    show(io, x)
  end

end
