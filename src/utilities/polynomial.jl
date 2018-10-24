export Composition, ishomogenous,
    uniquevar,
    homogenize, check_zero_dimensional, remove_zeros!, ncompositions, expansion, nvariables,
    maxdegrees, variables

const MPPoly = MP.AbstractPolynomialLike
const MPPolys = Vector{<:MP.AbstractPolynomialLike}

const WeightedVariable = Tuple{<:MP.AbstractVariable, Int}


abstract type AbstractComposition end
struct Composition{T<:Union{<:MPPolys, <:AbstractComposition}} <: AbstractComposition
    g::MPPolys
    f::T
end

expansion(F::MPPolys) = F
function expansion(C::Composition)
    vars = MP.variables(C.g)
    map(g -> MP.subs(g, vars => expansion(C.f)), C.g)
end

function Base.show(io::IO, C::Composition)
    println(io, "Composition of $(ncompositions(C)) polynomial systems with expansion:")

    show(io, expansion(C))
end
Base.length(C::Composition) = length(C.g)

ncompositions(C::Composition{<:MPPolys}) = 2
ncompositions(C::Composition{<:Composition}) = 1 + ncompositions(C.f)

compose(g::MPPolys, f::Composition) = Composition(g, f)
compose(g::MPPolys, fs::MPPolys...) = Composition(g, compose(fs...))
compose(f::MPPolys) = f
import Base: ∘
∘(g::MPPolys, f::Union{<:MPPolys, <:Composition}) = compose(g, f)
function ∘(C::Composition, f::Union{<:MPPolys, <:Composition})
    # ∘ is left associative and not right associative ....
    Composition(C.g, C.f ∘ f)
end

nvariables(C::Composition) = nvariables(C.f)
nvariables(F::MPPolys) = MP.nvariables(F)

variables(C::Composition) = variables(C.f)
variables(F::MPPolys) = MP.variables(F)


"""
    degree(term::MP.AbstractTermLike, vars)

Compute the (weighted) degree of `f` in the variables `vars`.
"""
function degree(term::MP.AbstractTermLike, variables::Vector{<:MP.AbstractVariable})
    sum(MP.degree(term, v) for v in variables)
end
function degree(term::MP.AbstractTermLike, weighted_variables)
    deg = 0
    for (v, w) in weighted_variables
        deg += w * MP.degree(term, v)
    end
    deg
end


"""
    minmaxdegree(f::MP.AbstractPolynomialLike, variables)

Compute the minimum and maximum (weighted total) degree of `f` with respect to the given variables.
"""
function minmaxdegree(f::MP.AbstractPolynomialLike, variables)
    d_min, d_max = typemax(Int), 0
    for term in f
        d = degree(term, variables)
        d_min, d_max = min(d, d_min), max(d, d_max)
    end
    d_min, d_max
end
function minmaxdegree(F::MPPolys, variables)
    map(f -> minmaxdegree(f, variables), F)
end


maxdegrees(F::MPPolys) = map(f -> MP.maxdegree(f), F)
maxdegrees(F::MPPolys, variables) = map(f -> minmaxdegree(f, variables)[2], F)

function maxdegrees(C::Composition)
    degrees = maxdegrees(C.f)
    maxdegrees(C.g, zip(MP.variables(C.g), degrees))
end



###############
# ishomogenous
##############

"""
    ishomogenous(f::MP.AbstractPolynomialLike, vars)

Checks whether `f` is homogenous in the variables `vars` with possible weights.

    ishomogenous(F::Vector{MP.AbstractPolynomialLike}, variables; parameters=nothing)

Checks whether each polynomial in `F` is homogenous in the variables `variables`.
"""
function ishomogenous(f::MP.AbstractPolynomialLike, variables)
    d_min, d_max = minmaxdegree(f, variables)
    d_min == d_max
end

function ishomogenous(F::MPPolys, vars=MP.variables(F); parameters=nothing)
    if parameters === nothing
        all(f -> ishomogenous(f, vars), F)
    else
        variables = setdiff(vars, parameters)
        all(f -> ishomogenous(f, variables), F)
    end
end

function ishomogenous(C::Composition; kwargs...)
    homogenous_degrees(C; kwargs...) !== nothing
end

function homogenous_degrees(C::Composition{<:MPPolys}; parameters=nothing)
    variables = MP.variables(C.f)
    if parameters !== nothing
        setdiff!(variables, parameters)
    end
    degrees = homogenous_degrees(C.f, variables)
    degrees === nothing && return nothing

    homogenous_degrees(C.g, zip(MP.variables(C.g), degrees))
end

function homogenous_degrees(C::Composition{<:Composition}; kwargs...)
    degrees = homogenous_degrees(C.f; kwargs...)
    degrees === nothing && return nothing
    homogenous_degrees(C.g, zip(MP.variables(C.g), degrees))
end

function homogenous_degrees(F::MPPolys, variables)
    weights = Int[]
    for f in F
        mindeg, maxdeg = minmaxdegree(f, variables)
        if mindeg != maxdeg
            return nothing
        else
            push!(weights, mindeg)
        end
    end
    weights
end



#############
# homogenize
#############
"""
    homogenize(f::MP.AbstractPolynomial, v::Vector{<:MP.AbstractVariable}, variable=uniquevar(f))

Homogenize the variables `v` in the polynomial `f` by using the given variable `variable`.

    homogenize(F::Vector{<:MP.AbstractPolynomial}, v::Vector{<:MP.AbstractVariable}, variable=uniquevar(F))

Homogenize the variables `v` in each polynomial in `F` by using the given variable `variable`.

    homogenize(f::MP.AbstractPolynomial, variable=uniquevar(f))

Homogenize the polynomial `f` by using the given variable `variable`.

    homogenize(F::Vector{<:MP.AbstractPolynomial}, variable=uniquevar(F))

Homogenize each polynomial in `F` by using the given variable `variable`.
"""
function homogenize(F::MPPolys, var::MP.AbstractVariable=uniquevar(F); parameters=nothing)
    variables = MP.variables(F)
    if parameters !== nothing
        homogenize(F, setdiff(variables, parameters), var)
    else
        homogenize(F, variables, var)
    end
end
function homogenize(F::MPPolys, variables, var::MP.AbstractVariable=uniquevar(F))
    map(f -> homogenize(f, variables, var), F)
end

function homogenize(C::Composition, var::MP.AbstractVariable=uniquevar(C.g); parameters=nothing)
    C̄, _ = homogenize_degrees(C, var; parameters=parameters)
    C̄
end


function homogenize(f::MPPoly, variables, var::MP.AbstractVariable)
    _, d_max = minmaxdegree(f, variables)
    MP.polynomial(map(f) do t
        d = d_max - degree(t, variables)
        t * var^d
    end)
end

function homogenize_degrees(F::MPPolys, variables, var::MP.AbstractVariable)
    F̄ = similar(F)
    degrees = Int[]
    for (i,f) in enumerate(F)
        f̄, d = homogenize_degree(f, variables, var)
        F̄[i] = f̄
        push!(degrees, d)
    end

    F̄, degrees
end

function homogenize_degree(f::MPPoly, variables, var::MP.AbstractVariable)
    _, d = minmaxdegree(f, variables)
    MP.polynomial(map(t -> t * var^(d - degree(t, variables)), f)), d
end

function homogenize_degrees(C::Composition{<:MPPolys}, var::MP.AbstractVariable; parameters=nothing)
    variables = MP.variables(C.f)
    if parameters !== nothing
        variables = setdiff(variables, parameters)
    end
    f̄, degrees = homogenize_degrees(C.f, variables, var)
    push!(f̄, var)
    ḡ, degrees = homogenize_degrees(C.g, zip(MP.variables(C.g), degrees), var)
    Composition(ḡ, f̄), degrees
end

function homogenize_degrees(C::Composition{<:Composition}, var::MP.AbstractVariable; kwargs...)
    f̄, degrees = homogenize_degrees(C.f, var; kwargs...)
    push!(f̄.g, var)
    ḡ, degrees = homogenize_degrees(C.g, zip(MP.variables(C.g), degrees), var)
    Composition(ḡ, f̄), degrees
end


"""
    remove_zeros!(F::Vector{<:MP.AbstractPolynomialLike})

Remove zero polynomials from the given system `F`.
"""
function remove_zeros!(F::MPPolys)
    filter!(!iszero, F)
    F
end
function remove_zeros!(C::Composition)
    filter!(!iszero, C.g)
    C
end

"""
    check_zero_dimensional(F::Vector{<:MP.AbstractPolynomial})

Check that the given polynomial system can have zero dimensional components.
"""
function check_zero_dimensional(F::Union{MPPolys, Composition})
    N = nvariables(F)
    n = length(F)

    if n ≥ N || (n == N - 1 && ishomogenous(F))
        return nothing
    end
    error("The input system will not result in a finite number of solutions.")
end


function homogenous_weights(F::MPPolys, variables=MP.variables(F))
    weights = Int[]
    for f in F
        mindeg, maxdeg = minmaxdegree(f, variables)
        if mindeg != maxdeg
            return nothing
        else
            push!(weights, mindeg)
        end
    end
    weights
end
function homogenous_weights(C::Composition)
    f_weights = homogenous_weights(C.f)
    f_weights === nothing && return nothing
    homogenous_weights(C.g, zip(variables(C.g), weights))
end






"""
    uniquevar(f::MP.AbstractPolynomialLike, tag=:x0)
    uniquevar(F::Vector{<:MP.AbstractPolynomialLike}, tag=:x0)

Creates a unique variable.
"""
uniquevar(f::MP.AbstractPolynomialLike, tag=:x0) = MP.similarvariable(f, gensym(tag))
uniquevar(F::Vector{<:MP.AbstractPolynomialLike}, tag=:x0) = uniquevar(F[1], tag)
uniquevar(C::Composition, tag=:x0) = uniquevar(C.g, tag)
