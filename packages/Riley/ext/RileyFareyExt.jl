module RileyFareyExt
import Riley as R
import Kleinian as K
import Farey as F

function R.wordisometry(word::Union{F.AbstractMonoidGen,F.AbstractMonoidWord},images::AbstractDict;
                        identity=nothing)
    all(g -> g isa K.AbstractIsometry,values(images)) ||
        throw(ArgumentError("generator images must be Kleinian isometries"))
    if isnothing(identity)
        isempty(images) && throw(ArgumentError("provide generator images or an explicit identity"))
        identity = one(first(values(images)))
    end
    identity isa K.AbstractIsometry && K.isidentity(identity) ||
        throw(ArgumentError("identity must be an identity isometry"))
    result = K.MobiusIsometry(K.matrix(identity))
    for letter in F.Monoids.eachgen(word)
        g = if haskey(images,letter)
            images[letter]
        elseif haskey(images,inv(letter))
            inv(images[inv(letter)])
        else
            throw(KeyError(letter))
        end
        # Evaluate each step rather than retaining a growing IsometryProduct.
        result = K.MobiusIsometry(K.matrix(result*g))
    end
    result
end
end
