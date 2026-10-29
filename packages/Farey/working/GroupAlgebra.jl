  module GroupAlgebras
  using Farey: DEFAULT_INT
  using Farey
  using Farey:AbstractMonoidElement, AbstractElement
  using Farey.Monoids
  import Farey.Monoids:MonoidHom
  export GroupAlgebra, ∂, ϵ, invt, mapt

  abstract type AbstractGroupRing<:AbstractElement end

  struct GroupAlgebra{T,S}<:AbstractGroupRing
    terms::Dict{T,S}
    function GroupAlgebra(kv::Dict{K,V}) where {K,V}
      kv0 = filter(kv) do x
        !iszero(x[2])
      end
      dtype = typeof(kv0)
      K0, V0 = dtype.parameters
      if V0<: Number
        V0 = big(V0)
      end
      new{K0,V0}(kv0)
    end
  end

  function _plus(x::T, y::S) where {T<:AbstractGroupRing, S<:AbstractGroupRing}
    T(mergewith(+, x.terms, y.terms))
  end

  function _plus(x0::T, y0::S) where {T<:AbstractMonoidElement, S<:AbstractMonoidElement}
    x, y = GroupAlgebra.([x0,y0])
    GroupAlgebra(mergewith(+, x.terms, y.terms))
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
    @eval Base.$op(c::$AbstractElement, d::$AbstractGroupRing) = $nm(promote(c, d)...)
    @eval Base.$op(c::$AbstractGroupRing, d::$AbstractElement) = $nm(promote(c, d)...)
    @eval Base.$op(c::$AbstractGroupRing, d::$AbstractGroupRing) = $nm(promote(c, d)...)
    @eval Base.$op(c::$AbstractElement, d::Number) = $nm(promote(c, d)...)
    @eval Base.$op(c::Number, d::$AbstractElement) = $nm(promote(c, d)...)
  end

  Base.:^(c::AbstractGroupRing, n::Integer) = prod([c for i in 1:n])

  (::Type{T})(x...) where T<:AbstractGroupRing = GroupAlgebra(x...)

  GroupAlgebra(x::T) where T<:AbstractArray = GroupAlgebra(Dict(x))
  GroupAlgebra(x::S) where S<:AbstractMonoidElement = GroupAlgebra(Dict{S, DEFAULT_INT}([x=>1]))
  Base.zero(::Type{T}) where T<:AbstractGroupRing = GroupAlgebra(Dict{T.parameters...}([]))
  Base.iszero(x::T) where T<:AbstractGroupRing = isempty(x.terms)

  Base.zero(x::T) where T<:AbstractGroupRing = GroupAlgebra(empty(x.terms))

  Base.promote_rule(::Type{GroupAlgebra{S, I1}}, ::Type{GroupAlgebra{T, I2}}) where {T,S,I1, I2} = GroupAlgebra{promote_type(T, S), promote_type(I1, I2)}

  Base.promote_rule(::Type{S}, ::Type{GroupAlgebra{T, I1}}) where {T,S<:Number,I1} = GroupAlgebra{T, promote_type(I1, S)}
  Base.promote_rule(::Type{S}, ::Type{GroupAlgebra{T, I1}}) where {T,S,I1} = GroupAlgebra{promote_type(T, S), I1}

  Base.convert(::Type{GroupAlgebra{T,S}}, x::R) where {T, S, R<:AbstractGroupRing} = GroupAlgebra([T(l)=>S(r) for (l,r) in pairs(x.terms)]|>Dict{T,S})
  Base.convert(::Type{GroupAlgebra{T,S}}, x::R) where {T, S, R<:AbstractMonoidElement} = GroupAlgebra([T(x)=>S(1)]|>Dict)

  Base.convert(::Type{GroupAlgebra{T,S}}, x::R) where {T, S, R<:Number} = GroupAlgebra([one(T)=>x]|>Dict{T,R})


  function Base.show(io::IO, x::T) where T<: AbstractGroupRing
    isempty(x.terms) && return print(io, "0")
    rep = join(["$(y==1 ? "" : y==-1 ? "-" : "($y)*")$x" for (x,y) in sort(pairs(x.terms), by=x->length(x))], " + ")
    print(io, rep)
  end

  function _foxdx(f::S, g::S) where S<:AbstractMonoidGen
    g.exp != 1 && error("Fox derivative not implemented for $g")
    id(f)==id(g) ? !xor(ispositive(f), ispositive(g)) ? GroupAlgebra(one(f)) : -GroupAlgebra(f) : zero(GroupAlgebra(f))
  end

  function ∂(w::T, g::S) where {T<:AbstractMonoidElement, S<:AbstractMonoidGen}
    # dx = one(x)
    y = MonoidWord(w)
    foldl(reverse(y), init=zero(GroupAlgebra(y))) do dx, x
      _foxdx(x, g) + x*dx
    end

  end


  function ∂(w::T, g::S) where {T<:AbstractGroupRing, S<:AbstractMonoidGen}
    mapreduce(l->l[2]*∂(l[1],g), +, pairs(w.terms), init=zero(w))
  end


  function MagnusCoeff(w::T, g::Vector{S}) where {T, S<:AbstractMonoidGen}
    reduce((x, l)->∂(x,l), g, init=w)
  end

  function MagnusHessian(w::T, g::Vector{S}; n = 2) where {T, S<:AbstractMonoidGen}
    Dict([(l)=>MagnusCoeff(w,l)])

  end

  function ∂(w::T, g::Vector{S}) where {T, S<:AbstractMonoidGen}
    [∂(w, x) for x in g]
  end

  function ∂(w::T) where T
    Dict([l=>∂(w, l) for l in unique(w)∪inv.(unique(w))])
  end

  function ϵ(w::T) where T<:AbstractGroupRing
    sum(values(w.terms))
  end

  function mapt(f::MonoidHom{K,V}, w::GroupAlgebra{T,S}) where {K,V, T, S}
    GroupAlgebra([V(f(l))=>r for (l,r) in pairs(w.terms)]|>Dict{promote_type(V,T), S})
  end

  function mapt(f::F, w::T) where {F, T<:AbstractGroupRing}
    GroupAlgebra([f(l)=>r for (l,r) in pairs(w.terms)])
  end

  invt(w) = mapt(inv, w)

  function Base.show(io::IO, m::MIME"text/plain", x::T) where T<:AbstractGroupRing
    println(io, "$(typeof(x))")
    show(io, x)
  end



end

# import .GroupAlgebras as ZG
# using Farey
# using Farey.Monoids
# f, g = FreeGen("f g")
# cof = Dict([f*g => 1, g/f => 2])
# cof1 = Dict([f => 1, g => 1])
# l = ZG.GroupAlgebra(cof1)
# l - l
# ll = l*l - l
# l
# ZG.GroupAlgebra(f)

# Fw(Q) = farey_word(Q, f, g)
# Θ(Q) = -((f\Fw(Q))|>conj_prefix)[1]
# Θ(3//7)

# il = l * im
# il^2

# Q = @cf0 [2]
# W = Fw(Q)
# # ZG.∂(W, [f, g])
# dWx, dWy = ZG.∂(Fw(Q), [f, g])
# ZG.∂(dWy, [f,g])
# dWy
# ZG.∂(W, g)

# ll.terms
# using StyledStrings
# for n in 1:5
#   qq = CF(0, [n])
#   global W = farey_word(qq, f, g);
#   # global ZG.∂(W);
#   # pairs(ZG.∂(W))|>x->join(x, '\n')
#   dWs = join([styled"{code:$(l)-->$r}\n" for (l, r) in pairs(ZG.∂(W))])
#   #
#   styled"{red:∂(W), W = $W} \n$(dWs)"|>printstyled
# end

# dWy*ZG.invt(dWy)
# dWy * ZG.mapt(reverse, ZG.invt(dWy))
