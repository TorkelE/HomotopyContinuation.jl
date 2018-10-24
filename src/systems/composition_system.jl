export CompositionSystem

"""
    CompositionSystem(composition::Composition, systems_constructor) <: AbstractSystem

A system representing the composition of polynomial maps.
"""
struct CompositionSystem{S1<:AbstractSystem, S2<:AbstractSystem} <: AbstractSystem
    f::S1 # Can be a CompositionSystem again
    g::S2 # Never a composition system
end

function CompositionSystem(C::Composition{<:Vector{<:MP.AbstractPolynomial}}, vars, system_constructor)
    CompositionSystem(system_constructor(C.f, vars), system_constructor(C.g))
end
function CompositionSystem(C::Composition{<:Composition}, vars, system_constructor)
    CompositionSystem(CompositionSystem(C.f, vars, system_constructor), system_constructor(C.g))
end
function CompositionSystem(C::Composition{<:Vector{<:MP.AbstractPolynomial}}, system_constructor)
    CompositionSystem(system_constructor(C.f), system_constructor(C.g))
end
function CompositionSystem(C::Composition{<:Composition}, system_constructor)
    CompositionSystem(CompositionSystem(C.f, system_constructor), system_constructor(C.g))
end

struct CompositionSystemCache{C1<:AbstractSystemCache, C2<:AbstractSystemCache, T} <: AbstractSystemCache
    cache_f::C1
    cache_g::C2

    eval_f::Vector{T}
    J_f::Matrix{T}
    J_g::Matrix{T}
end

function cache(C::CompositionSystem, x)
    cache_f = cache(C.f, x)
    eval_f = evaluate(C.f, x, cache_f)
    cache_g = cache(C.g, eval_f)
    J_f = jacobian(C.f, x, cache_f)
    J_g = jacobian(C.g, eval_f, cache_g)
    # need to bring eval_f, J_f, J_g to the same element type
    T = promote_type(eltype(eval_f), eltype(J_f), eltype(J_g))
    eval_f_T = convert(Vector{T}, eval_f)
    J_f_T = convert(Matrix{T}, J_f)
    J_g_T = convert(Matrix{T}, J_g)

    CompositionSystemCache(cache_f, cache_g, eval_f_T, J_f_T, J_g_T)
end

Base.size(C::CompositionSystem) = (length(C.g), size(C.f)[2])

function evaluate!(u, C::CompositionSystem, x, c::CompositionSystemCache)
    evaluate!(c.eval_f, C.f, x, c.cache_f)
    evaluate!(u, C.g, c.eval_f, c.cache_g)
end
function evaluate(C::CompositionSystem, x, c::CompositionSystemCache)
    evaluate!(c.eval_f, C.f, x, c.cache_f)
    evaluate(C.g, c.eval_f, c.cache_g)
end
function jacobian!(U, C::CompositionSystem, x, c::CompositionSystemCache)
    # chain rule
    evaluate_and_jacobian!(c.eval_f, c.J_f, C.f, x, c.cache_f)
    jacobian!(c.J_g, C.g, c.eval_f, c.cache_g)
    LinearAlgebra.mul!(U, c.J_g, c.J_f)
end
function jacobian(C::CompositionSystem, x, c::CompositionSystemCache)
    # chain rule
    evaluate_and_jacobian!(c.eval_f, c.J_f, C.f, x, c.cache_f)
    jacobian!(c.J_g, C.g, c.eval_f, c.cache_g)
    c.J_g * c.J_f
end
function evaluate_and_jacobian!(u, U, C::CompositionSystem, x, c::CompositionSystemCache)
    evaluate_and_jacobian!(c.eval_f, c.J_f, C.f, x, c.cache_f)
    evaluate_and_jacobian!(u, c.J_g, C.g, c.eval_f, c.cache_g)
    LinearAlgebra.mul!(U, c.J_g, c.J_f)
    nothing
end
