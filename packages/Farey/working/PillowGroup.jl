module PillowGroup
export π1O, Klein4, (⋊), pi1OGroup
using Farey
F2 = FreeGroup(:f, :g)
E = Farey.Hom(F2, F2, (
    1 => 2,
    2 => 1,
))
O = Farey.Hom(F2, F2, (
    1 => inv(F2(1)),
    2 => F2(2),
))

vhom = E
hhom = (E*O)^2
phom = vhom*hhom

struct pi1OGroup<:Farey.AbstractGroup end

@enum Klein4 begin
    id = 0b00
    h = 0b01
    p = 0b11
    v= 0b10
end

applyK4 = Dict([
    Klein4(0) => E^2,
    Klein4(1)=> hhom,
    Klein4(3)=> phom,
    Klein4(2)=> vhom,
])

struct PillowElement <: Farey.AbstractGroupElement
    elem::Farey.FreeGroupElement
    hom::Klein4
    parent::pi1OGroup
    
    function PillowElement(x::Farey.FreeGroupElement, hom::Klein4)
        new(x, hom, pi1OGroup())        
    end
end

function π1O(x)
    x==0 && return one(pi1OGroup())
    PillowElement(F2(x), Klein4(0))
end

function π1O(x, z::Klein4)
    x==0 && return PillowElement(one(F2), z)
    PillowElement(F2(x), z)
end

(::pi1OGroup)(x...) = π1O(x...)

(z::Klein4)(g::Farey.FreeGroupElement) = applyK4[z](g)

Base.one(::pi1OGroup)::PillowElement = PillowElement(one(F2), Klein4(0))
Base.one(x::PillowElement)::PillowElement = one(x.parent)
Base.similar(x::PillowElement)::PillowElement = one(x.parent)


GPC.isfinite(::pi1OGroup) = false

GPC.gens(::pi1OGroup) = [
    [ π1O(i, Klein4(0)) for i in 1:2]...,
    [ π1O(0, Klein4(i)) for i in 1:3]...
]

GPC.parent(c::PillowElement) = c.parent

# function Base.getindex(g::PillowElement, I::Union{UnitRange, Integer})
#     return g.parent(g.word[I])
# end

function Base.inv(g::PillowElement)
    return (z = g.hom; x=g.elem; π1O(z(inv(x)), z))
end

function Base.inv(g::Klein4)
    return g
end

function Base.:(*)(g::Klein4, h::Klein4)
    Klein4(xor(Int(g), Int(h)))
end

function Base.:(*)(g::PillowElement, h::PillowElement)
    return π1O(g.elem * g.hom(h.elem), g.hom * h.hom)
end

function ⋊(x::Farey.FreeGroupElement, z::Klein4)
    π1O(x, z)
end

function Base.:(^)(m::PillowElement, n::Integer)::PillowElement
    n < 0 && return inv(m)^-n
    return Base.power_by_squaring(m, n)
end


function Base.show(io::IO,  c::PillowElement)
    # if get(io, :compact, false)::Bool
    #     repr_string =  pretty_rep(c)
    #     elems = split(repr_string, ".")
    #     if length(elems)>8
    #         repr_string = join(elems[1:4], ".")*"..."*join(elems[end-4:end], ".")
    #     end
    # else
    #     repr_string =  pretty_rep(c)
    # end
    return print(io, "$(c.elem)⋊$(c.hom)")
end

end