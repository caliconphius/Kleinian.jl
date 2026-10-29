function magnus_lyndon(w, x)
  head = 0
  word = one(x)
  for c in w
    if Farey.id(c)==Farey.id(x)
      head += (c.inv ? -1 : 1)
    else
      gen = head >= 0 ? FreeGen("$(string(c)[1])$head") : FreeGen("$(uppercase(string(c)[1]))$(-head)")
      word *= (c.inv ? inv(gen) : gen)
    end
  end
  word * (x^head)
end

function expand_conjugates(w; x = FreeGen("g"))
  word = one(w)
  for c in w
    idnum = (c.id - 1) % (2 * 26)
    name =
      (idnum ) in 0:25 ? "$('a' + idnum)" :
      "$('A' + idnum - 26)"
    conj_p = (-1 + Int((c.id - 1) ÷ (26 * 2)))
    conj_p = conj_p < 0 ? 0 : conj_p
    if isuppercase(only(name))
      conj_p = -conj_p
      name = lowercase(name)
    end
    gen = FreeGen("$(name)")^(inv(x^conj_p))
    word *= c.inv ? inv(gen) : gen
  end
  word
end

newvar = (ee*(ffhi(1))*ee)|(MonoidHom([f=>element("u"), g=>element("v")]))
newvar1 = (ee*(ffh(1))*ee)|
(MonoidHom([f=>element("u"), g=>element("v")]))

u,v = FreeGen("u v")

magnus_lyndon(Fw(2//5)|newvar,u)
magnus_lyndon(Fw(3//8)|newvar,v)
magnus_lyndon(Fw(2//5)|newvar,v)
magnus_lyndon(Fw(1//5)|newvar,u)
magnus_lyndon(Fw(1//2)|newvar,v)
