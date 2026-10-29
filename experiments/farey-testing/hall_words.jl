using Base: isgreater
import Transducers as TD
struct HallTree
  leaf::Union{Pair{HallTree, HallTree}, Int}
  weight::Int
end
# using Farey
# Threads.@spawn include("working_setup.jl")
function HallTree(n::Int)
  HallTree(n, 1)
end

function ∧(h1::HallTree, h2::HallTree)
  HallTree(h1=>h2, h1.weight + h2.weight)
end

h1 = HallTree(1)
h2 = HallTree(2)
h1∧h2∧h2∧h1∧h2

L = [(h1).leaf...]
function hall_foliage(w::HallTree; f = FreeGen(0))
  word = f
  L = [w]
  while length(L) > 0
    h = popfirst!(L)
    if h.weight == 1
      word*= FreeGen(h.leaf)
    else
      pushfirst!(L, h.leaf[2])
      pushfirst!(L, h.leaf[1])
    end
  end
  return word
end
function nextword(n, w, alphabet)
    x = (w ^ ((n ÷ length(w)) + 1))[begin:n]
    while x != "" && x[end] == alphabet[end]
        x = x[begin:end-1]
    end
    if x != ""
        last_char = x[end]
        next_char_index = something(findfirst(==(last_char), alphabet), 0) + 1
        x = x[begin:end-1] * alphabet[next_char_index]
    end
    return x
end

function generate_lyndon_words(n, alphabet)
    lwords = String[]
    w = string(alphabet[begin])
    while length(w) <= n
        push!(lwords, w)
        w = nextword(n, w, alphabet)
        w == "" && break
    end
    return lwords
end
f,g,h,k = FreeGen("f g h k")
const lyndon_words = generate_lyndon_words(4, "fghk")
wrds = generate_lyndon_words(3, "ab")



ϕ = MonoidHom([
  h=>-f,k=>-g
  ])
M = 6

ϕ.(sort( lyndon_words, by=length).|>element)|>s->filter(x->length(x)>M-2, s).|>ab|>unique.|>println
ϕ.(sort( lyndon_words, by=length).|>element)|>
s->filter(x->abs.(ab(x))∈[[0,0],[1,1]], s)|>
s->filter(x->length(x)>M-2, s)|>
s->filter(x->length(conj_prefix(x)[1])==0, s).|>(x->compute_magnus(x,2).|>(l->(l[1],l[3])))




function duval(s; by = isless)
  n = length(s)
  i = 1
  factorization = []
  while (i <= n)
    j, k = i+1,i
    while ((j <= n) && (by(s[k],s[j])||s[k]==s[j]))
      k = (s[k] == s[j]) ? k+1 : i
      j += 1
    end
    while (i <= k)
      push!(factorization, s[i:i+j-k-1])
      i += j-k
    end
  end

  return factorization
end
Farey.id(-f)

duval("fgffgg"^2)
duval("ffgffgg"; by = (x,y)->x>y)

function standard_factor(w)

  length(w) == 1 && return w
  last = duval(w[2:end]^2)[end]
  first = w[1:end-length(last)]
  [first, last]
end

function min_cycle(w)
  lset = duval(w^2)
  conj = lset|>m->TD.foldl(*, TD.ReduceIf(x->length(x)<length(strn)),m)
  w[length(conj)+1:end], w[length(conj)+1:end] * conj
end
# duval(wordrep(Fw(1//3)^2))
strn = "fggf"

Int(true)
function Base.isless(x::T,y::T) where T<:Farey.AbstractMonoidGen
  (Farey.id(x)<Farey.id(y)) && return true
  ((Farey.id(x)==Farey.id(y)) && (Int(x.inv) < Int(y.inv)))
end
f < -f < g < -g


a,A,b,B = FreeGen("a A b1 B1")
b
a < A < b < B
duval([-1,10,3,7,5,1,9])
w = -[2,3,2,1,1,3,1,2,1]

duval(("fffg")^2)



# q = 3//5
# duval(Fw(2//5)^2)
# Fw(2//3)
# -Fw(3//5)|>duval
# f<-f<g<-g

# duval(Fw(24//103))[1]|>signs_seq|>x->filter(==(5), x)
# duval(Fw(24//103))[2]|>signs_seq|>x->filter(==(4), x)
# duval(Fw(24//103))[3]|>signs_seq|>x->filter(==(4), x)
# duval(-(Fw(24//103)))
# duval(Fw(24//103))[3]|>ab
# @cf 24//103
# 1101/1900 -  11/19
# 4491/7750
# 85691/147750
# 572891/987750
# 20636/35611

# L = [1,1,2,2,-2,-1,-1,2, -1, -1, -2, 2, 2, 1, 1,2,1,1,2,2,-2,-1,-1,2, -1, -1, -2, 2, 2, 1, 1]
# zq = CF0(vcat([L for i in 1:10]...))
# ql = zq + 0 |> CF
# # qll = @cf 70074988645719856485435076741543952473492446955961670165653024482384764616949625525032//120926556450635346995745990570818386212639103542661863619232176319046501420921001943585
# ql.L |> unique
# qnl = @cf0 [1,1,2,2,-2,-1,-1,2, -1, -1, -2, 8, 2, 1, 1,4,1,1,2,2,-2,-1,-1,2, -1, -1, -2, 8, 2, 1, 1]
# qnl + 0.0
# q = @cf0 [2,2,2,4,-2,-2,-2]
# q = [1,1,2,2,-2,-1,-1,2, -1, -1, -2, 8, 2, 1, 1] |> CF0
# q = @cf0 [1,1,2,2,-2,-1,-1,2, -1, -1, -2, 8, 2, 1, 1,4]
# [CF0([2,4,2,2m,-2, -4, -2]) for m in 0:6]
# @cf [2,4,2]
# decomp = farey_word(q,f,g)^2|>eachgen|>l->map(x->x.inv ? "z$(uppercase(repr(x))[1])" : repr(x),l) |>duval
# decomp.|>length
# decomp.|>length|>l->accumulate(+,l)
# denominator(q)*2
# q +0.0

# loong = [1, 1, 2, 1, 1, 1, 4, 1, 1, 7, 2, 1, 1, 4, 1, 1, 2, 1, 1, 1, 4, 1, 1, 7, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 1, 1, 7, 2, 1, 1, 4, 1, 1, 2, 1, 1, 1, 4, 1, 1, 7, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 1, 1, 7, 2, 1, 1, 4, 1, 1, 2, 1, 1, 1, 4, 1, 1, 7, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 1, 1, 7, 2, 1, 1, 4, 1, 1, 2, 1, 1, 1, 4, 1, 1, 7, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 1, 1, 7, 2, 1, 1, 4, 1, 1, 2, 1, 1, 1, 4, 1, 1, 7, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 1, 1, 7, 2, 1, 1, 4, 1, 1, 2, 1, 1, 1, 4, 1, 1, 7, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 1, 1, 7, 2, 1, 1, 4, 1, 1, 2, 1, 1, 1, 4, 1, 1, 7, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 1, 1, 7, 2, 1, 1, 4, 1, 1, 2, 1, 1, 1, 4, 1, 1, 7, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 1, 1, 7, 2, 1, 1, 4, 1, 1, 2, 1, 1, 1, 4, 1, 1, 7, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 1, 1, 7, 2, 1, 1, 4, 1, 1, 2, 1, 1, 1, 4, 1, 1, 7, 2, 2]|>
# duval;
# loong.|> length
# q|>farey_neighbours

# q|>FareyDecomp!
# "f" < "zF"
# "g" < "zF"
