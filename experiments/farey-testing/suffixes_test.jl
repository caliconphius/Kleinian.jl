    include("working_setup.jl")
    import .SuffixTrees as sftr
    const TREE = sftr.SuffixTree
    FINAL_CHAR = sftr.FINAL_CHAR

    # function isprefix(tr, pw)
    #   tkns = sftr.text(pw)
    #   all()  
    # end
    # wrd = Γ(10, 71)
    twrd = ((f*g)^3*f)
    pal_conj(twrd)
    # pal_conj(Fw(1//2))  

    tr = let strg = twrd
      # Sw = replace(strg, " "=>"_")
      # Sw = strg[2:end-1]
      Sw = strg
      TREE(Sw)
    end


    function prefix_tree(tr::TREE{T}, idx::Int)::TREE{T} where T
      # nodes = filter(tr.tree) do x
      #   x.position <= idx
      # end
      # tr_prfx = TREE{T}(tr.word[1:idx], nodes, idx)
      # sftr.extend_suffix!(tr_prfx, FINAL_CHAR, idx+1)
      # tr_prfx
      TREE(tr.word[1:idx])
    end
  tr1 = sftr.SuffixTree("abacabaabaca")

    # sftr.printpath(tr, 31)|>print
    # sftr.printpath(tr, 31)|>print
    # tr = let strg = Γ(3,5)("a","b")
    #   Sw = replace(strg, " "=>"_")

    #   TREE(Sw)
    # end

    # [sftr.subword(tr, i) for i in 1:3]
    # tr
    # prefix_tree(tr, 5)

    ###
      # Qa, Qb, Qc = farey_neighbours(Qw)
      # ww=Γ(Qw)
      # G0i(ww)
      # D1i(ww)
      # D1
      # D1i(ww)
      # tr = SuffixTrees.SuffixTree(ww)
      # # SuffixTrees.reprtree(tr)
      # (ffi(Qa)(ww))
      # (G0*ffi(Qa)*G0)(ww)
      # (G0*ffi(Qa)*G0)(ww)
      # (ww)
      # G1(ffi(Qb)(ww))
      # G1(ffi(Qc)(ww))

      # Qw = @cf0 [2,1,2]
      # let q = Qw
      # Qs = Q_path(q)
      # Qs1 = [q1⊖q2 for (q1,q2) in zip(Q_path(q)[1:end-1],Q_path(q)[2:end])]
      # ϕ = d0/g1
      # ϕh = hom(ϕ^7)

      # [ab.([D0(Γ(q)), Γ(q)]) for q in Qs]
      # [([ϕh(Γ(q)), Γ(q)]) for q in Qs]
      # [ϕh(Γ(q)) for q in Qs].|>conj_prefix.|>println
      # [ϕh(Γ(q)) for q in Qs1].|>conj_prefix.|>println
      # [signs_seq.(conj_prefix(ϕh(Γ(q)))) for q in Qs]
      # [(ϕh(Γ(q))) for q in Qs]
      # w1 = ϕh(Γ(Q[2]))
      # w1 = ϕh(Γ(Q[3]))
      # SuffixTrees.SuffixTree(w1)
      # end
      # Γ(Qw)
      # # tr for


      # ϕss(Q) = prod([(g1*e)^l for l in Q.L])
      # ϕss(Qw)|>hom|>x->x.([a,b]).|>conj_prefix
      # ϕss(Qw)^-1|>hom|>x->x.([a,b]).|>conj_prefix
      # # ϕss(Qw)|>hom|>x->x.([a,b]).|>conj_prefix
      # Qw = @cf0 [3,1,2]
      # SuffixTrees.SuffixTree(Fw(Qw))


      # # ((G0i*D1))(

      # (ff(qwd))((Γ(qwd)))
      # ff(3//7)(a)
      # Γ(qwd)

      # qw = @cf0 [4, 2, 2, 8]
      # # qw = @cf0 [1, 2, 2, 3]
      # HH = ffi(qw)
      # Qs = Q_path(qw)
      # HH = (G0*ffi(Qs[3]))(Γ(qw))
      # Γ(qw)|>SuffixTrees.SuffixTree
      # qa, qb, qc = farey_neighbours(qw)
      # HH = ffi(qw)
      # (HH)(Γ(qw))
      # (HH)(Γ(qa))
      # (HH)(Γ(qb))
      # (HH)(Γ(qc))
      # #|>SuffixTrees.SuffixTree



      # let wd = (f*g)^3*g^4
      # ab(wd)
      # f*Pal(wd)*g

      # abl = f*Pal(wd)*g|>ab
      # qwd = CF(abl[2]//abl[1])\
      # wchr = f*Pal(wd)*g
      # qa, qb, qc = farey_neighbours(qwd)
      # HH = G0i*(ff(qwd))
      # # HH = (ffi(qwd))
      # end
#       print("""
#        ⨁-⋅-⋅   ⋅
#       ╱ ╲  
           
#      ⋅   ⋅ 


#        ⋅ ⋅⋅ ╲╱ ⋅-⋅   
#       """)


# print(raw"""
#  ╲ ╱  ⎹⎸ ⎹⎸ ⎹⎸  ╲ ╱  ⎹⎸ ⎹⎸ ⎹⎸ ⎹.⎸
#   ╱   ⎹⎸ ⎹⎸ ⎹⎸   ╱   ⎹⎸ ⎹⎸ ⎹⎸ ⎹.⎸
#  ╱ ╲  ⎹⎸ ⎹⎸ ⎹⎸  ╱ ╲  ⎹⎸ ⎹⎸ ⎹⎸ ⎹.⎸
# ⎹⎸  ╲ ╱   ╲ ╱  ⎹⎸  ╲ ╱   ╲ ╱  ⎹.⎸
# ⎹⎸   ╲     ╲   ⎹⎸   ╲     ╲   ⎹.⎸
# ⎹⎸  ╱ ╲   ╱ ╲  ⎹⎸  ╱ ╲   ╱ ╲  ⎹.⎸
# ⎹⎸ ⎹⎸  ╲ ╱  ⎹⎸ ⎹⎸ ⎹⎸  ╲ ╱  ⎹⎸ ⎹.⎸
# ⎹⎸ ⎹⎸   ╱   ⎹⎸ ⎹⎸ ⎹⎸   ╱   ⎹⎸ ⎹.⎸
# ⎹⎸ ⎹⎸  ╱  ╲ ⎹⎸ ⎹⎸ ⎹⎸  ╱ ╲  ⎹⎸ ⎹.⎸""")
# join([
# "⎹⎸  ╲ ╱  ⎹⎸ ⎹⎸ ⎹⎸  ╲ ╱  ⎹⎸ ⎹⎸ ⎹.⎸" ,
# "⎹⎸   ╲   ⎹⎸ ⎹⎸ ⎹⎸   ╲   ⎹⎸ ⎹⎸ ⎹.⎸",
# "⎹⎸  ╱ ╲  ⎹⎸ ⎹⎸ ⎹⎸  ╱ ╲  ⎹⎸ ⎹⎸ ⎹.⎸",
# "⎹⎸  ╲ ╱  ⎹⎸ ⎹⎸ ⎹⎸  ╲ ╱  ⎹⎸ ⎹⎸ ⎹.⎸",
# "⎹⎸   ╲   ⎹⎸ ⎹⎸ ⎹⎸   ╲   ⎹⎸ ⎹⎸ ⎹.⎸",
# "⎹⎸  ╱ ╲  ⎹⎸ ⎹⎸ ⎹⎸  ╱ ╲  ⎹⎸ ⎹⎸ ⎹.⎸",
# "⎹⎸  ╲ ╱  ⎹⎸ ⎹⎸ ⎹⎸  ╲ ╱  ⎹⎸ ⎹⎸ ⎹.⎸",
# "⎹⎸   ╲   ⎹⎸ ⎹⎸ ⎹⎸   ╲   ⎹⎸ ⎹⎸ ⎹.⎸",
# "⎹⎸  ╱ ╲  ⎹⎸ ⎹⎸ ⎹⎸  ╱ ╲  ⎹⎸ ⎹⎸ ⎹.⎸",
# "⎹⎸  ╲ ╱  ⎹⎸ ⎹⎸ ⎹⎸  ╲ ╱  ⎹⎸ ⎹⎸ ⎹.⎸",
# "⎹⎸   ╲   ⎹⎸ ⎹⎸ ⎹⎸   ╲   ⎹⎸ ⎹⎸ ⎹.⎸",
# "⎹⎸  ╱ ╲  ⎹⎸ ⎹⎸ ⎹⎸  ╱ ╲  ⎹⎸ ⎹⎸ ⎹.⎸"
# ]
# , "\n")|>print


