using Transducers
import Transducers as TC
using Transducers:
    @next, R_, Transducer, combine, complete, inner, next, start, unwrap, wrap, wrapping
struct Chunk
    s::String
end
struct Vacant
    l::String
    r::String
end
vacant_or_chunk(c::Char) = c == ' ' ? Vacant("", "") : Chunk(string(c))
extract(x::Chunk, y::Chunk) = nothing, Chunk(x.s * y.s)
extract(x::Chunk, y::Vacant) = nothing, Vacant(x.s * y.l, y.r)
extract(x::Vacant, y::Chunk) = nothing, Vacant(x.l, x.r * y.s)
extract(x::Vacant, y::Vacant) = maybeword(x.r * y.l), Vacant(x.l, y.r)

maybeword(s) = isempty(s) ? nothing : s
struct WordsXF <: Transducer end
  Transducers.start(rf::R_{WordsXF}, init) = wrap(rf, Chunk(""), start(inner(rf), init))
  function Transducers.next(rf::R_{WordsXF}, acc, x)
      wrapping(rf, acc) do state, iacc
          word, state = extract(state, x)
          iacc = next(inner(rf), iacc, word)
          return state, iacc
      end
  end
  function Transducers.complete(rf::R_{WordsXF}, acc)
      state, iacc = unwrap(rf, acc)
      if state isa Vacant
          pre = @next(inner(rf), start(inner(rf), Init), maybeword(state.l))
          iacc = combine(inner(rf), pre, iacc)  # prepending `state.l`
          iacc = @next(inner(rf), iacc, maybeword(state.r))  # appending `state.r`
      else
          @assert state isa Chunk
          iacc = @next(inner(rf), iacc, maybeword(state.s))
      end
      return complete(inner(rf), iacc)
  end
  function Transducers.combine(rf::R_{WordsXF}, a, b)
      ua, ira = unwrap(rf, a)
      ub, irb = unwrap(rf, b)
      word, uc = extract(ua, ub)
      ira = @next(inner(rf), ira, word)
      irc = combine(inner(rf), ira, irb)
      return wrap(rf, uc, irc)
  end

  wordsxf = opcompose(Map(vacant_or_chunk), WordsXF(), Filter(!isnothing))
  collect(wordsxf, "This is a sample")

  using TypedTables
  astable(xs) = copy(Table, xs)  # using `TypedTables` for a nice display

  nested_objects = [
      (a = 1,  b = [(c = 2,  d = 3),  (c = 4,  d = 5)]),
      (a = 10, b = [(c = 20, d = 30), (c = 40, d = 50)]),
  ];
  more_nested_objects = [
      (a = 1,  b = [(c = 2,  d = [(e = 3,  f = 4),  (e = 4,  f = 5)]),
                    (c = 6,  d = [])]),
      (a = 10, b = [(c = 20, d = [(e = 30, f = 40), (e = 40, f = 50)])]),
  ];
  more_nested_objects |>
  MapCat() do x
      a2 = x.a * 2
      x.b |> MapCat() do b
          a2_plus_c = a2 + b.c
          b.d |> Map() do d
              c_plus_e = b.c + d.e
              c_plus_f = b.c + d.f
              (q1 = a2_plus_c, q2 = c_plus_e, q3 = c_plus_f)
          end
      end
  end |> astable

  using Transducers
  using BangBang  # for `push!!`
  using Farey
using Farey.Monoids
Fw(1)
1:9 |> PartitionBy(x -> (x + 1) ÷ 3) |> TC.Map(copy) |> collect
ψf(2//5)
foldl(
  TeeRF(push!!, TC.Map(ψf(2//5).rep)'(push!!), commu),
  Replace(Dict([a=>Fw(1//2)])),
  -Γ(1//2)
)
foldl(
  TeeRF(*, push!!),
  opcompose(
    TC.MapCat(Replace(Dict("f"=>f, "g"=>g,"h"=>h, "k"=>k)))
  ),
  lyndon_words
)

xf = Filter(x->length(x)==3) |>
  MapCat() do x
    x|>(Replace(Dict("f"=>f, "g"=>g,"h"=>h, "k"=>k)))
  end
lyndon_words |> MapCat() do x
  x |> TC.Map() do
    x|> Replace(Dict(["f"=>f, "g"=>g,"h"=>h, "k"=>k]))
  end
end|>collect
xf = TC.Map() do x
  x|>
  TC.Map(string) |>
  TC.Map(c->Dict(["f"=>f, "g"=>g,"h"=>-f, "k"=>-g])[c]) |>
  prod
end

lyndon_words|> Filter(x->length(x)==4)|>xf |> unique |>collect
