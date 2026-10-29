function decomposeWord(qq)
  f, g = FreeGen("f g")

  qq = CF(qq)
  qa, qb, qc = farey_neighbours(qq)
  if qq[end] % 2 == 0

    qu, l = descendQ(qq, qa; l=1)
    k = l ÷ 2
    qr = qa
    qt = qu ⊖ qa
  else
    qu, l = descendQ(qq, qa; l=0)
    k = l ÷ 2
    qr = qa
    qt = qu ⊕ qa

  end

  fwr = Fw(qr)
  fwt = Fw(qt)

  if qr<qt
    @info (qr,qt,k)
    (Fw(qr)^k, Fw(qt), -(Fw(qr)^f)^k)
  else
    @info (qr,qt,k)
    # (qr,qt,k,-1)
    ((-Fw(qr)^-f)^k, Fw(qt),(Fw(qr))^k)
  end
end

  function FareyDecomp!(qq)
    f, g = FreeGen("f g")

    qq = CF(qq)
    qa, qb, qc = farey_neighbours(qq)
    if qq[end] % 2 == 0

      qu, l = descendQ(qq, qa; l=1)
      k = l ÷ 2
      qr = qa
      qt = qu ⊖ qa
    else
      qu, l = descendQ(qq, qa; l=0)
      k = l ÷ 2
      qr = qa
      qt = qu ⊕ qa

    end
    pow(k) = to_superscript("$k")
    npow(k) = to_superscript("-$k")
    # qr = addK(qv, qu, l1); qt = addK(qv, qu, l1+1);
    printstyled(styled"{yellow:$(addK(qt, qr, k*2))}   {code:qq = $(qq)}\n{tip:k = $k}; {bright_blue:(k-$(l%2))/2 = $(k÷2)}; {magenta:l = $l}\n\n")
    printstyled(styled"{red:qr = $(qr)}   {green:qt = $(qt)}\n{bright_blue:qc = $(qc)} {yellow:qu = $(qu)} {magenta:qa = $(qa)}\n\n")

    fwr = Fw(qr)
    fwt = Fw(qt)
    nt = 2 * numerator(qt)
    dt = 2 * denominator(qt)
    nr = 2 * numerator(qr)
    dr = 2 * denominator(qr)
    F = inv(f)
    if qr < qt
      printstyled("" *
                  (styled"{bright_magenta:$(fwr)}" *
                   styled"{bright_blue:$(fwr)}")^(k ÷ 2) *
                  styled"{red:$(fwr)}"^(k % 2) *
                  styled"{green:$(fwt)}" *
                  styled"{yellow:$(inv(f))}" *
                  styled"{bright_red:$(inv(fwr))}"^(k % 2) *
                  (styled"{magenta:$(inv(fwr))}" *
                   styled"{code:$(inv(fwr))}")^(k ÷ 2) *
                  styled"{yellow:f}\n"
      )
      printstyled(styled"{emphasis:$(Fw(qq))}\n")
      printstyled(styled"\n{red:Wr[$(Rational(qq))]} \t= \t...{emphasis:$(Fw(qq)[k*dr+1:k*dr+dt+1])}...\n")
      printstyled(styled"{red:[$(k*dr+1)...$(k*dr+dt+1)]}\n")
      printstyled(styled"{bright_blue:Wr[$(Rational(qt))]/f} \t= \t   {magenta:$(fwt)}{yellow:$(F)}\n")

      printstyled(styled"\n{red:Wr[$(Rational(qq))]} \t= \t{emphasis:$(Fw(qq)[1:dr])}...\n")
      printstyled(styled"{red:[1...$(dr+1)]}\n")
      printstyled(styled"{bright_blue:Wr[$(Rational(qr))]} \t= \t{magenta:$(fwr)}\n")
      printstyled(styled"\n{code:Wr[$(Rational(qq))]} \t= {bright_red:Wr}{bright_blue:$(pow(k))}{bright_red:[$(Rational(qr))]}.{green:Wr[$(Rational(qt))]}.{yellow:$(inv(f))}.{bright_red:Wr}{bright_blue:$(npow(k))}{bright_red:[$(Rational(qr))]}.{yellow:$(f)}\n")

      printstyled(styled"\n \t\t= ({bright_red:$(fwr)})$(pow(k)).{green:$(fwt)}.{yellow:$(inv(f))}.({bright_red:$(inv(fwr)*f)})$(pow(k))\n")
    else

      printstyled("" *
                  styled"{yellow:$(f)}" *
                  styled"{bright_red:$(inv(fwr))}"^(k % 2) *
                  (styled"{magenta:$(inv(fwr))}" *
                   styled"{code:$(inv(fwr))}")^(k ÷ 2) *
                  styled"{yellow:$(inv(f))}" *
                  styled"{green:$(fwt)}" *
                  (styled"{bright_magenta:$(fwr)}" *
                   styled"{bright_blue:$(fwr)}")^(k ÷ 2) *
                  styled"{red:$(fwr)}\n"^(k % 2)
                  * "\n"
      )
      printstyled(" ↓\n" *
                  styled"{yellow:$(f)}" *
                  styled"{bright_red:$(inv(fwr))}"^(k % 2) *
                  (styled"{magenta:$(inv(fwr))}" *
                   styled"{code:$(inv(fwr))}")^(k ÷ 2) *
                  styled"{green:$(fwt[2:end])}" *
                  (styled"{bright_magenta:$(fwr)}" *
                   styled"{bright_blue:$(fwr)}")^(k ÷ 2) *
                  styled"{red:$(fwr)}"^(k % 2)
                  * "\n"
      )

      printstyled(styled"{emphasis:$(Fw(qq))}\n")


      printstyled(styled"\n{red:Wr[$(Rational(qq))]} \t= \t  ...{emphasis:$(Fw(qq)[k*dr+2:k*dr+dt])}...\n")
      printstyled(styled"{red:[$(k*dr+2)...$(k*dr+dt)]}\n")
      printstyled(styled"{bright_blue:$(F).Wr[$(Rational(qt))]} \t= \t{yellow:$(F)}{magenta:$(fwt)}\n")

      printstyled(styled"\n{red:Wr[$(Rational(qq))]} \t= \t{emphasis:$(Fw(qq)[1:dr+1])}...\n")
      printstyled(styled"{red:[1...$(dr+1)]}\n")
      printstyled(styled"{bright_blue:f/Wr[$(Rational(qr))]} \t= \t{yellow:$(f)}{magenta:$(inv(fwr))}\n")
      printstyled(styled"\n{code:Wr[$(Rational(qq))]} \t= {yellow:$(f)}.{bright_red:Wr}{bright_blue:$(npow(k))}{bright_red:[$(Rational(qr))]}.{yellow:$(inv(f))}.{green:Wr[$(Rational(qt))]}.{bright_red:Wr}{bright_blue:$(pow(k))}{bright_red:[$(Rational(qr))]}\n")

      printstyled(styled"\n \t\t= ({yellow:$((f))}.{bright_red:$(inv(fwr))})$(pow(k)).{green:$(F*fwt)}.({bright_red:$(fwr)})$(pow(k))\n")
    end
  end
