using Farey
using Farey.CHV:Garside as GAR
import Chevie as CHV

using Chevie:Garside as GAR

W = CHV.coxgroup(:A, 2)
# W = CHV.coxgroup(:A,3)
B = CHV.BraidMonoid(W)
BKL = CHV.DualBraidMonoid(W)
B(BKL.δ)
BKL(B.δ)

B(B.δ) |>GAR.word
BKL(BKL.δ) |>GAR.left_divisors


Δ = B(B.δ)
δ = BKL(1,2)
∂(s) = inv(s)*Δ

s0 = BKL(1)
t0 = BKL(2)
S0 = BKL(-1)
T0 = BKL(-2)

s = B(1)
t = B(2)
S = B(-1)
T = B(-2)

Δ* S

els0 = GAR.elements(BKL, 5);
els = GAR.elements(B, 5);

val = W |> typeof
CHV.generators(W)
ww = W(2)
# W(3)

CHV.PermGroups.Perm(fd.elm) |> B

BKL.δword
chvWord = CHV.Groups.word
weird = Δ * s*Δ*T*t*S*S
ppp = GAR.preferred_prefix(weird)
BKL(ppp)

GAR.left_divisors(Δ^3)

BKL.((weird/Δ*Δ^5)|>chvWord)|>prod
BKL.(weird.elm)|>prod
B(1,2,1,-2,-1,-2)
B(1,2,-1,-2)
B(1,-2)




(Δ^5 * s *T * Δ |> chvWord ) .|> B |> prod

s0*t0*s0
B.(BKL(3,1,2,2, 3).elm)|>prod

s*t

GAR.elements(BKL, 5)

BKL(1,2,1)

CHV.Garside.rightgcd(S, T)

s*t*s

s0*δ*s0*t0|>CHV.Garside.Brieskorn_normal_form

GAR.decycle((t*s)^4)
GAR.cycle((t*s)^4)
[Δ^x * S for x in 1:4]
GAR.leftgcd(s/t,t*s*t/s)
GAR.leftlcm(s/t*s,t*s*t/s)

S0 = inv(s0)
T0 = inv(t0)

s0 \ δ
t0 \ δ



www = (s*t*s)^4
www.pd

CHV

p = CHV.Presentation(Symbol.(collect("abcd")), W)
p |> CHV.simplify
p |> CHV.display_balanced






# Brieskorn Saito
