export Endgame, Result, allowed_keywords

const allowed_keywords = [:sampling_factor, :egtol, :minradius, :maxnorm, :minimal_maxnorm,
    :maxwindingnumber, :max_extrapolation_samples, :cauchy_loop_closed_tolerance,
    :cauchy_samples_per_loop, PathTracking.allowed_keywords...]

struct Options
    # See Endgame docstring for explanations
    sampling_factor::Float64
    tol::Float64
    minradius::Float64
    maxnorm::Float64
    minimal_maxnorm::Float64
    maxwindingnumber::Float64
    max_extrapolation_samples::Int
    cauchy_loop_closed_tolerance::Float64
    cauchy_samples_per_loop::Int
end

mutable struct State{V<:ProjectiveVectors.AbstractProjectiveVector, C, T}
    # Current value
    x::V
    # prediction and previous prediction
    p::V
    pprev::V
    pbest::V
    pbest_delta::Float64
    # Samples of the solution path obtained from a geometric series
    # defined by x(λ^k R₀).
    # The raw vectors are stored.
    samples::Matrix{C}
    nsamples::Int
    # The data structure is that `logabs_samples` Contains a vector for each
    # coordinate.
    logabs_samples::Matrix{T}
    directions::Matrix{T}

    R::Float64
    npredictions::Int
    iters::Int
    # current status
    status::Symbol
    # current estimate of winding number
    windingnumber::Int
    in_endgame_zone::Bool
end

function State(x, R₀, options)
    n = maxsamples(R₀, options)

    R = real(R₀)

    samples = fill(zero(eltype(x)), length(x), n)
    nsamples = 1
    directions = fill(NaN, length(x), n)
    logabs_samples = copy(directions)

    p = copy(x)
    pprev = copy(x)
    pbest = copy(x)
    pbest_delta = Inf

    npredictions = 0
    iters = 0

    status = :ok
    windingnumber = 1
    in_endgame_zone = false

    State(copy(x), p, pprev, pbest, pbest_delta, samples, nsamples, logabs_samples,
        directions,
        R, npredictions, iters, status,
        windingnumber, in_endgame_zone)
end

function maxsamples(R₀, options)
    ceil(Int, log(options.sampling_factor, options.minradius / real(R₀))) + 2
end

function reset!(state::State, x, R)
    state.nsamples = 1
    state.x .= x
    for i = 1:length(x)
        state.samples[i, 1] = x[i]
        state.logabs_samples[i, 1] = logabs(x[i])
    end

    state.pbest_delta = Inf

    state.npredictions = 0
    state.iters = 0
    state.R = R
    state.status = :ok

    state.windingnumber = 0
    state.in_endgame_zone = false
    nothing
end

struct Cache{V}
    log_sampling_factor::Float64
    windingnumbers::Vector{Int}
    unitroots::Vector{ComplexF64}
    direction_buffer::Matrix{Float64}
    xbuffer::V
    pbuffer::V
end

function Cache(state::State, options::Options)
    log_sampling_factor = log(options.sampling_factor)
    windingnumbers = zeros(Int, length(state.x))
    unitroots = Vector{ComplexF64}()
    direction_buffer = zeros(length(state.x), options.max_extrapolation_samples)
    pbuffer = copy(state.x)
    xbuffer = copy(state.x)
    Cache(log_sampling_factor, windingnumbers, unitroots, direction_buffer,
        pbuffer, xbuffer)
end

"""
    Endgame(H, x; options...)

Construct an `Endgame` to run the endgame for paths of the type of `x` and the homotopy
`H`.

## Options
* `sampling_factor=0.5` During the endgame we approach ``0`` by the geometric series ``h^kR₀``
where ``h`` is `sampling_factor` and `R₀` the endgame start provided in `runendgame`.
* `egtol=1e-10` This is the tolerance necessary to declare the endgame converged.
* `minradius=1e-15` A path is declared false if the endgame didn't finished until then.
* `maxnorm=1e5` If our original problem is affine we declare a path at infinity if the infinity norm
with respect to the standard patch is larger than `maxnorm`.
* `minimal_maxnorm=min(1e3, maxnorm)` A path is **not** declared going to infinity if the infinity norm
with respect to the standard patch is not larger than `minimal_maxnorm`.
* `maxwindingnumber=15` The maximal windingnumber we try to find using Cauchys integral formula.
* `max_extrapolation_samples=4` During the endgame a Richardson extrapolation is used to improve the accuracy
of certain approximations. This is the maximal number of samples used for this.
* `cauchy_loop_closed_tolerance=1e-3` The tolerance for which is used to determine whether a loop is closed.
The distance between endpoints is normalized by the maximal difference between any point in the loop and the starting point.
* `cauchy_samples_per_loop=6` The number of samples used to predict an endpoint. A higher number of samples should result
in a better approximation. Note that the error should be roughly ``t^n`` where ``t`` is the current time of the loop
and ``n`` is `cauchy_samples_per_loop`.
* `pathtrackerkwargs...` During the endgame a [`PathTracking.PathTracker`](@ref) is used. These are all arguments possible
to be supplied to it (with the excemption of `patch` this is always [`AffinePatches.FixedPatch()`](@ref)).
"""
struct Endgame{P<:PathTracking.PathTracker, V}
    tracker::P
    state::State{V}
    cache::Cache
    options::Options
end

function Endgame(H::Homotopies.AbstractHomotopy, x::ProjectiveVectors.AbstractProjectiveVector;
    sampling_factor=0.5, egtol=1e-10, minradius=1e-15, maxnorm=1e5, minimal_maxnorm=min(1e3, maxnorm),
    maxwindingnumber=15,
    max_extrapolation_samples=4,
    cauchy_loop_closed_tolerance=1e-3, cauchy_samples_per_loop=6,
    patch=nothing, pathtrackerkwargs...)

    options = Options(sampling_factor, egtol, minradius, maxnorm, minimal_maxnorm,
        maxwindingnumber, max_extrapolation_samples,
        cauchy_loop_closed_tolerance, cauchy_samples_per_loop)
    H = Homotopies.PatchedHomotopy(H, AffinePatches.state(AffinePatches.FixedPatch(), x))
    tracker = PathTracking.PathTracker(H, x, complex(0.1,0.0), 0.0im; pathtrackerkwargs...)
    state = State(x, complex(1.0,0.0), options)

    Endgame(tracker, state, Cache(state, options), options)
end

"""
    Result(endgame)

## Fields
* `returncode::Symbol`
* `x::Vector{T}`: The solution or last prediction.
* `t::Float64`: The point in time corresponding to `x`.
* `res::Float64`: The residual of the homotopy at `x` and `t`.
* `windingnumber::Int`: The windingnumber estimated by the endgame. It can be 0 if no estimation happened.
* `npredictions`: The number of predictions.
* `iters`: The number of iterations.
"""
struct Result{V}
    returncode::Symbol
    x::V
    t::Float64
    res::Float64
    windingnumber::Int
    npredictions::Int
    iters::Int
end

function Result(endgame::Endgame)
    state = endgame.state

    if state.npredictions == 0 || state.status == :at_infinity
        x = copy(state.x)
    else
        x = copy(state.pbest)
    end
    Result(state.status, x, real(state.R), state.pbest_delta,
        state.windingnumber, state.npredictions, state.iters)
end

Base.show(io::IO, R::Result) = Utilities.print_fieldnames(io, R)
Base.show(io::IO, ::MIME"application/prs.juno.inline", R::Result) = R
