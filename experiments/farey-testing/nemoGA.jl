module NemoGroupAlgebras
  using Farey
  using Farey:AbstractMonoidElement, AbstractElement, DEFAULT_INT
  using Farey.Monoids
  using Nemo

  const CoefficientType = Union{Number, RingElement}

  export NemoGroupAlgebra, ∂, ϵ, invt, mapt

  abstract type AbstractGroupRing<:AbstractElement end

  struct NemoGroupAlgebra{T<:AbstractMonoidElement,S<:RingElement}<:AbstractGroupRing
    terms::Dict{T,S}
    coeff_ring::Nemo.Ring
    function NemoGroupAlgebra(kv::Dict{K,V}) where {K,V<:Union{RingElement, Number}}
      kv0 = filter(kv) do x
        !iszero(x[2])
      end
      dtype = typeof(kv0)
      K0, V0 = dtype.parameters
      if V0<:Number
        ctype = parent_type(poly_type(QQ))
      else
        ctype = first(values(kv0))|>parent
      end
      kv0 = Dict([l=>ctype(r) for (l,r) in  pairs(kv0)])
      new{K0,V0}(kv0, ctype)
    end
  end

  (::Type{T})(x...) where T<:AbstractGroupRing = NemoGroupAlgebra(x...)

  NemoGroupAlgebra(x::T) where T<:AbstractArray = NemoGroupAlgebra(Dict(x))
  NemoGroupAlgebra(x::S) where S<:AbstractMonoidElement = NemoGroupAlgebra(Dict([x=>ZZ(1)]))
  Base.zero(::Type{T}) where T<:AbstractGroupRing = NemoGroupAlgebra(Dict{T.parameters...}([]))
  Base.iszero(x::T) where T<:AbstractGroupRing = isempty(x.terms)

  Base.zero(x::T) where T<:AbstractGroupRing = NemoGroupAlgebra(empty(x.terms))

  Base.promote_rule(::Type{NemoGroupAlgebra{S, I1}}, ::Type{NemoGroupAlgebra{T, I2}}) where {T,S,I1, I2} = NemoGroupAlgebra{promote_type(T, S), promote_type(I1, I2)}

  Base.promote_rule(::Type{S}, ::Type{NemoGroupAlgebra{T, I1}}) where {T,S<:Number,I1} = NemoGroupAlgebra{T, promote_type(I1, S)}
  Base.promote_rule(::Type{S}, ::Type{NemoGroupAlgebra{T, I1}}) where {T,S<:RingElement,I1} = NemoGroupAlgebra{T, S}
  Base.promote_rule(::Type{S}, ::Type{NemoGroupAlgebra{T, I1}}) where {T,S<:AbstractElement,I1} = NemoGroupAlgebra{promote_type(T, S), I1}

  Base.convert(::Type{NemoGroupAlgebra{T,S}}, x::R) where {T, S<:RingElement, R<:AbstractGroupRing} = NemoGroupAlgebra([T(l)=>x.coeff_ring(r) for (l,r) in pairs(x.terms)]|>Dict)

  Base.convert(::Type{NemoGroupAlgebra{T,S}}, x::R) where {T, S<:Number, R<:AbstractMonoidElement} = NemoGroupAlgebra([T(x)=>ZZ(1)]|>Dict)
  Base.convert(::Type{NemoGroupAlgebra{T,S}}, x::R) where {T, S<:RingElement, R<:AbstractMonoidElement} = NemoGroupAlgebra([T(x)=>one(S)]|>Dict)

  Base.convert(::Type{NemoGroupAlgebra{T,S}}, x::R) where {T, S, R} = NemoGroupAlgebra([one(T)=>x]|>Dict)

  function _plus(x::T, y::S) where {T<:AbstractGroupRing, S<:AbstractGroupRing}
    T(mergewith(+, x.terms, y.terms))
  end

  function _plus(x0::T, y0::S) where {T<:AbstractMonoidElement, S<:AbstractMonoidElement}
    x, y = NemoGroupAlgebra.([x0,y0])
    NemoGroupAlgebra(mergewith(+, x.terms, y.terms))
  end

  function Base.:-(x::T) where T<:AbstractGroupRing
    T([l=>-r for (l,r) in pairs(x.terms)])
  end

  function _subt(x::T, y::S) where {T, S}
    T(mergewith(+, x.terms, (-y).terms))
  end

  function _mult(x::T, y::S) where {T, S}
    (iszero(x) || iszero(y)) && return zero(T)
    multerms = foldl(
      mergewith!(+),
      [Dict(first(a)*first(b)=>last(a)*last(b)) for (a,b) in Iterators.product(pairs(x.terms), pairs(y.terms))]
    )

    T(multerms)
  end

  binary_ops = (:+, :-, :*)
  named_ops = (_plus, _subt, _mult)

  for (op, nm) in zip(binary_ops, named_ops)
    @eval Base.$op(c::$AbstractElement, d::$AbstractElement) = $nm(promote(c, d)...)
    @eval Base.$op(c::$AbstractElement, d::RingElement) = $nm(promote(c, d)...)
    @eval Base.$op(c::RingElement, d::$AbstractElement) = $nm(promote(c, d)...)
  end


  function Base.show(io::IO, x::T) where T<: AbstractGroupRing
    isempty(x.terms) && return print(io, "0")
    rep = join(["$(y==1 ? "" : y==-1 ? "-" : "($y)*")$x" for (x,y) in pairs(x.terms)], " + ")
    print(io, rep)
  end

  function _foxdx(f::S, g::S) where S<:AbstractMonoidGen
    g.exp !== 1 && error("Fox derivative not implemented for $g")
    id(f)==id(g) ? ispositive(f) ? NemoGroupAlgebra(one(f)) : -NemoGroupAlgebra(f) : zero(NemoGroupAlgebra(f))
  end

  function ∂(w::T, g::S) where {T<:AbstractMonoidElement, S<:AbstractMonoidGen}
    # dx = one(x)
    y = MonoidWord(w)
    foldl(reverse(y), init=zero(NemoGroupAlgebra(y))) do dx, x
      _foxdx(x, g) + x*dx
    end

  end


  function ∂(w::T, g::S) where {T<:AbstractGroupRing, S<:AbstractMonoidGen}
    mapreduce(l->l[2]*∂(l[1],g), +, pairs(w.terms), init=zero(w))
  end


  function ∂(w::T, g::Vector{S}) where {T, S<:AbstractMonoidGen}
    [∂(w, x) for x in g]
  end

  function ϵ(w::T) where T<:AbstractGroupRing
    sum(values(w.terms))
  end

  function mapt(f::F, w::T) where {F, T<:AbstractGroupRing}
    NemoGroupAlgebra([f(l)=>r for (l,r) in pairs(w.terms)])
  end

  invt(w) = mapt(inv, w)

  function Base.show(io::IO, m::MIME"text/plain", x::T) where T<:AbstractGroupRing
    println(io, "$(typeof(x))")
    show(io, x)
  end

end
