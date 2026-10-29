import Transducers as TD
using GLMakie

τ = MonoidHom([
  f => -f,
  g => -g,
])
Q = 2 // 5
Fw(Q)
τ(Fw(Q)) |> t -> println.([t, inv(t)]);
lr.(funnel(5 // 12))
SternBrocot(N) = 1:N |> TD.Iterated(next_minimal, QQ(2 // 5)) |> TD.Filter(s -> s <= 1) |> TD.Filter(s -> numerator(s) % denominator(s) ∉ [denominator(s) - 1, 1])


# function plotroots(N)
#   rts = SternBrocot(N) |> (TD.Map(t -> (t, geometric_root(t)))) |> collect

#   xys = (rts) .|> t -> Point(real(t[2]), imag(t[2]), Float64(t[1]))

#   fig, ax, pl = datashader(xys)
#   fig
# end

function plotroots(N)
  rts = SternBrocot(N) |> (TD.Map(t -> (t, geometric_root(t)))) |> collect

  xys = last.(rts) .|> t -> (real(t), imag(t))
  xs = first.(xys) .|> Float64
  ys = last.(xys) .|> Float64

  fig, ax, pl = scatter(xs, ys; color=Float64.(first.(rts)))
  @info "finished"
  (fig)
end

plotroots(100)
# task.result



function farey_rec(n; init=1 // 2)
  init = ContinuedFraction(init)
  dfs = Queue{ContinuedFraction{BigInt}}()
  out = Vector{ContinuedFraction{BigInt}}()
  push!(dfs, init)
  for i in 1:n
    tnew = popfirst!(dfs)
    ql, qr, _ = farey_neighbours(tnew)
    push!(dfs, ql ⊕ tnew)
    push!(dfs, qr ⊕ tnew)
    # (numerator(tnew) % denominator(tnew) ∉ [denominator(tnew) - 1, 1]) &&
    push!(out, tnew)
  end
  out
end

farey_rec(50; init=1 // 3)
