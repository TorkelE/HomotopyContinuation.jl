module Correctors

using ..NewHomotopies
using ..Utilities

export AbstractCorrector,
    AbstractCorrectorCache,
    cache,
    correct!,
    Newton, NewtonCache

# Interface
abstract type AbstractCorrector end
abstract type AbstractCorrectorCache end

"""
    cache(::AbstractCorrector, ::HomotopyWithCache{M, N}, x, t)::AbstractCorrectorCache

Construct a cache to avoid allocations.
"""
function cache end


"""
    correct!(xnext, ::AbstractCorrector, ::AbstractCorrectorCache, H::HomotopyWithCache, x, t, tol)

Perform a correction step such that in the end `H(xnext,t) < tol`. Return `true` if the
correction was successfull otherwise return `false`.
"""
function correct! end


# Newton
"""
    Newton(;maxiters=3)

* `maxiters`: Maximal number of iterations until a correction is declared as failed.
"""
struct Newton <: AbstractCorrector
    maxiters::Int
end
Newton(;maxiters=3) = Newton(maxiters)


struct NewtonCache{T} <: AbstractCorrectorCache
    A::Matrix{T}
    b::Vector{T}
end

function cache(::Newton, H::HomotopyWithCache, x, t)
    A = jacobian(H, x, t)
    b = evaluate(H, x, t)
    NewtonCache(A, b)
end

function correct!(xnext, alg::Newton, cache::NewtonCache, H::HomotopyWithCache{N, N}, x, t, tol) where N
    A, b = cache.A, cache.b
    evaluate_and_jacobian!(b, A, H, x, t)
    solve_with_lu_inplace!(A, b)
    @. xnext = xnext - b

    k = 1
    while true
        evaluate!(b, H, xnext, t)

        if norm(b, Inf) < tol
            return true
        elseif k ≥ alg.maxiters
            return false
        end

        k += 1

        # put jacobian in A
        jacobian!(A, H, xnext, t)

        solve_with_lu_inplace!(A, b)
        @. xnext = xnext - b
    end
end

end
