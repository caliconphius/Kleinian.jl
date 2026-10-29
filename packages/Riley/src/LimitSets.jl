"""
    limitset(generators; seed=complexidealpoint(0.0), depth=8, maxpoints=1_000_000)

Sample the orbit of an ideal seed under freely reduced words of exactly
`depth`. Inverses are added automatically; supply each generator once.
Returns a vector of `IdealPoint`s, including infinity if encountered.
Depth zero returns the seed. A point-count bound is checked before sampling.

For suitable discrete groups (e.g. classical Schottky groups), increasing
depth approximates the limit set. This is finite orbit sampling, not a
discreteness test or a certified error bound; arbitrary groups and seeds
need not produce a fractal. Relations and repeated points are not removed.
Actions use the supplied isometries and their existing coefficient arithmetic.
"""
function limitset(generators;seed::IdealPoint=complexidealpoint(0.0),depth::Integer=8,
                  maxpoints::Integer=1_000_000)
    collect(limitsetstream(generators;seed,depth,maxpoints))
end

function _limitsetactions(generators)
    gs = collect(generators)
    !isempty(gs) && all(g -> g isa AbstractIsometry,gs) ||
        throw(ArgumentError("provide at least one isometry generator"))
    vcat(gs,inv.(gs))
end
_inverseletter(j,n) = j <= n ? j+n : j-n

struct LimitSetStream{A,P}
    actions::A
    seed::P
    depth::Int
    count::Int
end
Base.IteratorSize(::Type{<:LimitSetStream}) = Base.HasLength()
Base.eltype(::Type{<:LimitSetStream}) = IdealPoint
Base.length(s::LimitSetStream) = s.count

"""
    limitsetstream(generators; seed=complexidealpoint(0.0), depth=8,
                   maxpoints=1_000_000)

Lazy depth-first enumeration of freely reduced words of exactly `depth`.
Retains O(depth) points rather than whole word levels. Each iteration yields
an `IdealPoint`; `collect` materializes it, and `limitsetplot` accepts it directly.
Traversal order is the same as `limitset`. Independent traversals are repeatable.
`maxpoints` bounds the number of yielded points, not bytes or total work.
"""
function limitsetstream(generators;seed::IdealPoint=complexidealpoint(0.0),
                        depth::Integer=8,maxpoints::Integer=1_000_000)
    0 <= depth <= typemax(Int)-1 || throw(ArgumentError("depth must be a nonnegative Int"))
    maxpoints >= 1 || throw(ArgumentError("maxpoints must be positive"))
    actions = _limitsetactions(generators)
    count = big(1)
    for level in 1:depth
        count *= level == 1 ? length(actions) : length(actions)-1
        count <= min(maxpoints,typemax(Int)) ||
            throw(ArgumentError("word count exceeds maxpoints; reduce depth"))
    end
    LimitSetStream(actions,seed,Int(depth),Int(count))
end

function Base.iterate(s::LimitSetStream)
    s.depth == 0 && return (s.seed,nothing)
    points = IdealPoint[s.seed]
    letters = Int[]
    # State is local to a traversal; the iterator itself remains reusable.
    _nextlimitset(s,(points,letters,ones(Int,s.depth)))
end
Base.iterate(::LimitSetStream,::Nothing) = nothing
Base.iterate(s::LimitSetStream,state) = _nextlimitset(s,state)
function _nextlimitset(s::LimitSetStream,state)
    points,letters,next = state
    n = length(s.actions) ÷ 2
    while !isempty(points)
        level = length(points)
        if next[level] > length(s.actions)
            pop!(points)
            isempty(letters) || pop!(letters)
            continue
        end
        j = next[level]
        next[level] += 1
        !isempty(letters) && j == _inverseletter(last(letters),n) && continue
        p = s.actions[j](last(points))
        level == s.depth && return (p,state)
        push!(points,p)
        push!(letters,j)
        next[level+1] = 1
    end
    nothing
end

struct RandomLimitSetStream{A,P,R}
    actions::A
    seed::P
    rng::R
    samples::Int
    burnin::Int
    thin::Int
end
Base.IteratorSize(::Type{<:RandomLimitSetStream}) = Base.HasLength()
Base.eltype(::Type{<:RandomLimitSetStream}) = IdealPoint
Base.length(s::RandomLimitSetStream) = s.samples

"""
    randomlimitset(generators; samples=100_000, burnin=100, thin=1,
                   seed=complexidealpoint(0.0), rng=Random.default_rng())

Lazy non-backtracking random walk on ideal points, with O(1) retained points.
Choose the first letter uniformly from generators and inverses, then uniformly
among letters other than the previous inverse. Discard `burnin` actions and
perform `thin` actions before each emitted sample. Results are correlated
orbit samples, not uniform area samples or independent fixed-length words.

The RNG is copied on construction and per traversal: an explicit seeded RNG
makes repeated traversals reproducible without advancing the caller's RNG.
The default obtains a fresh seed from Julia's default RNG. With one generator,
non-backtracking cannot switch direction, so a walk samples only one end.
Use floating generators/seeds for long walks; exact coefficients may grow.
"""
function randomlimitset(generators; samples::Integer=100_000,burnin::Integer=100,
                        thin::Integer=1,seed::IdealPoint=complexidealpoint(0.0),
                        rng::Random.AbstractRNG=Random.default_rng())
    0 <= samples <= typemax(Int) || throw(ArgumentError("samples must be a nonnegative Int"))
    0 <= burnin <= typemax(Int) || throw(ArgumentError("burnin must be a nonnegative Int"))
    1 <= thin <= typemax(Int) || throw(ArgumentError("thin must be a positive Int"))
    actions = _limitsetactions(generators)
    ownedrng = rng === Random.default_rng() ? Random.Xoshiro(rand(rng,UInt64)) : copy(rng)
    RandomLimitSetStream(actions,seed,ownedrng,Int(samples),Int(burnin),Int(thin))
end

function _randomsteps(s,p,last,rng,steps)
    m = length(s.actions)
    for _ in 1:steps
        if last == 0
            j = rand(rng,1:m)
        else
            excluded = _inverseletter(last,m÷2)
            j = rand(rng,1:m-1)
            j >= excluded && (j += 1)
        end
        p = s.actions[j](p)
        last = j
    end
    p,last
end
function Base.iterate(s::RandomLimitSetStream)
    s.samples == 0 && return nothing
    rng = copy(s.rng)
    p,last = _randomsteps(s,s.seed,0,rng,s.burnin)
    Base.iterate(s,(p,last,rng,0))
end
function Base.iterate(s::RandomLimitSetStream,state)
    p,last,rng,emitted = state
    emitted == s.samples && return nothing
    p,last = _randomsteps(s,p,last,rng,s.thin)
    p,(p,last,rng,emitted+1)
end

"""
    rileygroup(ρ)

Two parabolic generators acting in the displayed complex coordinate as
`X(z)=z+1`, `Y(z)=z/(ρ*z+1)`. Returns existing `MobiusIsometry` objects,
with matrix conjugation accounting for the library's spinor convention.
Requires finite nonzero ρ; does not certify membership in the Riley slice.
At ρ=2im the group has the Apollonian gasket limit set; at ρ=4,
`X*inv(Y)` is an additional parabolic. See docs/limitsets.md for references.
"""
function rileygroup(ρ)
    coefficient_isfinite(ρ) && !iszero(ρ) ||
        throw(ArgumentError("Riley parameter must be finite and nonzero"))
    u,z = one(ρ),zero(ρ)
    (MobiusIsometry([u u;z u]),MobiusIsometry(conj.([u z;ρ u])))
end


"""
    limitsetplot(points; renderer=:datashader, figure=(;), axis=(;), kwargs...)

Plot ideal points in a 2D complex plane with Makie. Returns `(figure,axis,plot)`.
Infinity is omitted; nonfinite Float64 conversions throw. Styling keywords
pass to Makie. The default `renderer=:datashader` aggregates pixel density,
with `async=false` for reliable saving; use `async=true` for interaction.
Use `renderer=:scatter` for individual markers. Load a Makie backend first.
Pass a vector of ideal points, `limitsetstream`, or `randomlimitset`.
Streams are projected directly, but Makie retains O(samples) 2D coordinates
for reaggregation; the renderer itself is not constant-memory.
"""
function limitsetplot end

"""`limitsetplot!(axis, points; kwargs...)` adds a complex-boundary plot."""
function limitsetplot! end

"""
    wordisometry(word, images; identity=nothing)

Evaluate a Farey.Monoids generator or word using a dictionary of generator
images in Kleinian. Load Farey to enable this optional adapter. Inverse letters
are inferred from a mapped inverse when necessary. Multiplication follows
Kleinian's convention: `g*h` acts by `g(h(p))`. Returns a `MobiusIsometry`.
An empty word gives the identity; an empty image map requires an explicit
identity isometry. This evaluates words, not cusp roots or discreteness.
"""
function wordisometry end
