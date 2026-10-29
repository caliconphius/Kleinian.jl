using Memoize
include("hall_words.jl")
@memoize function standard_comm(s::String)
  length(s)==1 && return FreeGen(s)
  length(duval(s))==1|| return standard_comm.(duval(s))
  first, last = standard_factor(s)
  commu(standard_comm(first), standard_comm(last))
end

macro fgroup(expression)

end



~(f,g) = commu(f,g)
# Memoize.memoize_cache(standard_comm)|>empty!
# 1//3
# fḡ[g,f][[f,g],g]
Fw(1//3)
(f/g)*(g~f)*((f~g)~g)
# 2//3
# fg[f̄,g]
Fw(2//3)
f*g*(-f~g)
(-f~g) # =
(g~f)*((f~g)~-f)# =
(g~f)*(f~(f~g))*(((f~g)~f)~-f)# =
f*g*(g~f)*
(f~(f~g))*
(((f~g)~f)~-f)# =
(((f~g)~f)~-f)# =

-standard_comm("fg")*standard_comm("ffg")*(((f~g)~f)~-f)
f/g*-standard_comm("fg")*standard_comm("fgg")
standard_comm("fgg")*g*f
# 1//4
(f~g)^2 *
((g~f)~g) * ((g~f)~f) *
((f~(g~f))~g)

# 3//4
w = (f~g)^2 * -standard_comm("fgg") * -standard_comm("ffg") * standard_comm("ffgg")
# compute_magnus(Fw(3//4)/w,4)
#
# ffhi(1//1)(Fw(2//3))|>x->magnus_lyndon(x, g)
# ffh(1)(Fw(1//3))|>x->magnus_lyndon(x, g)
# import Transducers as TD
# function magnus_lyndon(w, x)
#   head = 0
#   word = one(x)
#   for c in w
#     if Farey.id(c)==Farey.id(x)
#       head += (c.inv ? -1 : 1)
#     else
#       gen = head >= 0 ? FreeGen("$(string(c)[1])$head") : FreeGen("$(uppercase(string(c)[1]))$(-head)")
#       word *= (c.inv ? inv(gen) : gen)
#     end
#   end
#   word * (x^head)
# end
#
# function expand_conjugates(w; x = FreeGen("g"))
#   word = one(w)
#   for c in w
#     idnum = (c.id - 1) % (2 * 26)
#     name =
#       (idnum ) in 0:25 ? "$('a' + idnum)" :
#       "$('A' + idnum - 26)"
#     conj_p = (-1 + Int((c.id - 1) ÷ (26 * 2)))
#     conj_p = conj_p < 0 ? 0 : conj_p
#     if isuppercase(only(name))
#       conj_p = -conj_p
#       name = lowercase(name)
#     end
#     gen = FreeGen("$(name)")^(inv(x^conj_p))
#     word *= c.inv ? inv(gen) : gen
#   end
#   word
# end
# (Fw(3//8))
#
# magnus_lyndon(Fw(3//8),f)
#
# c = magnus_lyndon(Fw(2//5),g)[1]
# -1+Int((c.id-1)÷(26*2))
# magnus_lyndon(Fw(1//4), g)
# magnus_lyndon(Fw(5//12), f)
# ffh(1)
# u,v = FreeGen("u v")
#
# newvar = (ffhi(1))|(MonoidHom([f=>element("u"), g=>element("v")]))
# newvar1 = (ffh(1))|(MonoidHom([f=>element("u"), g=>element("v")]))
#
# newvar1(Fw(1//3))|>x->(magnus_lyndon(x, FreeGen("v")),x)
#
# newvar|(MonoidHom([u=>element("f g"), v=>element("g")]))
# newvar1|(MonoidHom([u=>element("f g -1"), v=>element("g")]))
