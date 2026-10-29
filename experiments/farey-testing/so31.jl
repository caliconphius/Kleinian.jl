include("erics_functions.jl")
zlx
Qx, xx = QQ[:x]

Qi, iq = number_field(xx^2+1, :i)
L, Lzs = Qi[:t,:x1,:x2,:x3]
K = fraction_field(L)
t,x1,x2,x3 = Nemo.gens(K)
X = K[
  t+x3 x1-iq*x2;
  x1+iq*x2 t-x3;
  ]

Msp = parent(X)
inv2x(X) = Nemo.tr(X)*one(X) - X
inv2x(X)

σ1 = K[
  0 1 ;
  1 0
  ]

σ2 = K[
  0 -iq ;
  iq 0
  ]
σ3 = K[
  1 0 ;
  0 -1
  ]

τ1 = iq*σ1
τ2 = iq*σ2
τ3 = iq*σ3
U = t*one(σ1)+ x1*τ1 + x2*τ2 + x3*τ3
using Nemo


function frac_conj(x)
  map_coefficients(conj,numerator(x))//
  map_coefficients(conj, denominator(x))
end

function adjt(A)
M = parent(A)
M([
  frac_conj(A[1,1]) frac_conj(A[2,1]);
  frac_conj(A[1,2]) frac_conj(A[2,2])
  ])
end

adjt(σ1)*σ2
