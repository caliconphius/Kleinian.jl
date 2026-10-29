include("sl2algebra.jl")
import AbstractAlgebra as ALG
g = rand(ALG.SymmetricGroup(Int8(5)));
h = rand(ALG.SymmetricGroup(UInt32(5)));
typeof(g*h)

function Rset(a::Perm{T}) where T
  L = Set{Pair{T,T}}([])
  n = parent(a).n
  for i in 1:n
    for j in (i+1):n
      if a[j] < a[i]
        push!(L,i=>j)
      end
    end
  end
  L
end

function geodesic_length(a::Perm{T}) where T
  Rset(a)|>length
end

function perm_element(L::Set{Pair{T,T}}, S::ALG.Generic.SymmetricGroup{T}) where T
  R = deepcopy(L)
  filter!((l->(first(l)-last(l)==1)), R)
  isempty(R) && return one(S)
  maxδ = maximum(l->first(l)-last(l), L)
  for i in 2:maxδ


  end

end
R2 = Rset(sw(1,5))
[last(l)-first(l) for l in R2]
filter((l->(last(l)-first(l)==1)), R2)

S = ALG.SymmetricGroup(UInt32(10))
R1 = Rset(sw(1,4))
R2 = Rset(sw(1,5))
R2∩ R1
R2-R1

longest_element(S::ALG.Generic.SymmetricGroup{T}) where T = S(collect(S.n:-1:1))
longest_element(S)|>geodesic_length
sw = (i,j) -> S("($i,$j)")
sw(1,2)^(S("(5,4,3,2)"))
parent(g).n
l = S("(1,5)")
# parent(l)
l[1]
symgens = S->[S("($i,$(i+1))") for i in 1:S.n-1]
sympairs = S->[i=>i+1 for i in 1:S.n-1]
sympairs(S) ∩ Rset(sw(1,5))
sw(3,4)*sw(2,3)*sw(4,5)*sw(1,2)*sw(1,5)
