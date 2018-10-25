export StaticParameterHomotopy, ParameterHomotopyCache, gamma

import LinearAlgebra
import MultivariatePolynomials
const MP = MultivariatePolynomials
import StaticPolynomials
const SP = StaticPolynomials
import StaticArrays: MVector, SVector

import ..Utilities


"""
    ParameterHomotopy(F, parameters;
        variables=setdiff(MP.variables(F), parameters),
        p₁=randn(ComplexF64, length(parameters)),
        p₀=randn(ComplexF64, length(parameters)),
        γ₁=nothing, γ₀=nothing)

Construct the homotopy
```math
H(x, t) = F(x, (tγ₁p₁+(1-t)γ₀p₀) / (tγ₁+(1-t)γ₀))
```,
where `p₁` and `p₀` are a vector of parameter values for ``F`` and
`γ₁` and `γ₀` are complex numbers. If `γ₁` or `γ₀` is `nothing`, it is assumed
that `γ₁` and `γ₀` are ``1``.
The input ``parameters`` specifies the parameter variables of ``F``.
Neccessarily, ``length(parameters) == length(p₁) == length(p₀)``.

Note that `p₁` and `p₀` are stored as a tuple `p` of `SVectors` and `γ₁` and `γ₀`
are stored as a tuple `γ` or as `γ=nothing`

    ParameterHomotopy(F, parameters;
        variables=setdiff(MP.variables(F), parameters),
        startparameters=randn(ComplexF64, length(parameters)),
        targetparameters=randn(ComplexF64, length(parameters)),
        startgamma=nothing, targetgamma=nothing)

This is a non-unicode variant where `γ₁=startparameters`, `γ₀=targetparameters`,
`γ₁=startgamma`, `γ₀=targetgamma`.
"""
mutable struct ParameterHomotopy{Sys<:AbstractSystem, NParams, T<:Number} <: AbstractHomotopy
    F::Sys
    p::NTuple{2, SVector{NParams, T}}
    γ::Union{Nothing, NTuple{2, ComplexF64}}
end

function ParameterHomotopy(F::AbstractSystem, p₁::AbstractVector, p₀::AbstractVector;
    startgamma=nothing, γ₀ = startgamma,
    targetgamma=nothing, γ₁ = targetgamma
    )


    if length(p₁) != length(p₀)
        error("Length of parameters provided doesn't match.")
    end

    if γ₁ === nothing || γ₀ === nothing
        γ = nothing
    else
        γ = (γ₁, γ₀)
    end

    p = promote(SVector{length(p₁)}(p₁), SVector{length(p₁)}(p₀))
    ParameterHomotopy(F, p, γ)
end

function ParameterHomotopy(F::Vector{T},
    parameters::AbstractVector{V};
    variables=setdiff(MP.variables(F), parameters),
    startparameters=randn(ComplexF64, length(parameters)), p₁ = startparameters,
    targetparameters=randn(ComplexF64, length(parameters)), p₀ = targetparameters,
    kwargs...) where {T<:MP.AbstractPolynomialLike, V<:MP.AbstractVariable}
    G = Systems.SPSystem(F; variables=variables, parameters=parameters)
    ParameterHomotopy(G, p₁, p₀; kwargs...)
end

struct ParameterHomotopyCache{C<:AbstractSystemCache, T1<:Number, T2<:Number} <: AbstractHomotopyCache
    F_cache::C
    ∂p∂t::Vector{T1}
    J_p::Matrix{T2}
end

(H::ParameterHomotopy)(x, t, c=cache(H, x, t)) = evaluate(H, x, t, c)

function cache(H::ParameterHomotopy, x, t)
    pt = p(H, t)
    F_cache = Systems.cache(H.F, x, pt)
    ∂p∂t = Vector{eltype(pt)}(undef, length(pt))
    J_p = Matrix(Systems.differentiate_parameters(H.F, x, pt, F_cache))

    ParameterHomotopyCache(F_cache, ∂p∂t, J_p)
end

Base.size(H::ParameterHomotopy) = size(H.F)

@inline function p(H::ParameterHomotopy, t)
    p₁, p₀ = H.p
    if H.γ === nothing
        return t * p₁ + (1 - t) * p₀
    else
        γ₁, γ₀ = H.γ
        tγ₁, γ₀_₁₋t = t * γ₁, (1 - t) * γ₀
        γ = (tγ₁ + γ₀_₁₋t)
        return (@fastmath tγ₁ / γ) * p₁ + (@fastmath γ₀_₁₋t / γ) * p₀
    end
end

@inline function ∂p∂t!(u, H::ParameterHomotopy, t, c::ParameterHomotopyCache)
    p₁, p₀ = H.p
    if H.γ === nothing
        for i in eachindex(p₁)
            u[i] = p₁[i] - p₀[i]
        end
    else
        γ₁, γ₀ = H.γ
        tγ₁, γ₀_₁₋t = t * γ₁, (1 - t) * γ₀
        γ = (tγ₁ + γ₀_₁₋t)
        pt = (@fastmath tγ₁ / γ) * p₁ + (@fastmath γ₀_₁₋t / γ) * p₀
        # quotient rule to get the derivative of p(t)
        λ = @fastmath γ₁ * γ₀ / (γ * γ)
        for i in eachindex(p₁)
            u[i] = λ * (p₁[i] - p₀[i])
        end
    end
end


function evaluate!(u, H::ParameterHomotopy, x, t, c::ParameterHomotopyCache)
    Systems.evaluate!(u, H.F, x, p(H, t), c.F_cache)
end
function evaluate(H::ParameterHomotopy, x, t, c::ParameterHomotopyCache)
    Systems.evaluate(H.F, x, p(H, t), c.F_cache)
end

function jacobian!(u, H::ParameterHomotopy, x, t, c::ParameterHomotopyCache)
    Systems.jacobian!(u, H.F, x, p(H, t), c.F_cache)
end
function jacobian(H::ParameterHomotopy, x, t, c::ParameterHomotopyCache)
    Systems.jacobian(H.F, x, p(H, t), c.F_cache)
end

function evaluate_and_jacobian!(u, H::ParameterHomotopy, x, t, c::ParameterHomotopyCache)
    Systems.evaluate_and_jacobian!(u, H.F, x, p(H, t), c.F_cache)
end

function dt!(u, H::ParameterHomotopy, x, t, c::ParameterHomotopyCache)
    # apply chain rule to H(x, p(t))
    pt = p(H, t)
    ∂p∂t!(c.∂p∂t, H, t, c)
    Systems.differentiate_parameters!(c.J_p, H.F, x, pt, c.F_cache)
    LinearAlgebra.mul!(u, c.J_p, c.∂p∂t)
end
