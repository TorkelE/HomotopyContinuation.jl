export ishomogenous,
    uniquevar,
    homogenize

const WeightedVariable = Tuple{<:MP.AbstractVariable, Int}

"""
    ishomogenous(f::MP.AbstractPolynomialLike)

Checks whether `f` is homogenous.

    ishomogenous(f::MP.AbstractPolynomialLike, vars)

Checks whether `f` is homogenous in the variables `vars` with possible weights.
"""
ishomogenous(f::MP.AbstractPolynomialLike) = MP.mindegree(f) == MP.maxdegree(f)
function ishomogenous(f::MP.AbstractPolynomialLike, variables)
    d_min, d_max = minmaxdegree(f, variables)
    d_min == d_max
end

"""
    ishomogenous(F::Vector{MP.AbstractPolynomialLike}, variables)

Checks whether each polynomial in `F` is homogenous in the variables `variables`.
"""
function ishomogenous(F::Vector{<:MP.AbstractPolynomialLike}, variables)
    all(f -> ishomogenous(f, variables), F)
end
function ishomogenous(F::Vector{<:MP.AbstractPolynomialLike}; parameters=nothing)
    if parameters !== nothing
        ishomogenous(F, setdiff(MP.variables(F), parameters))
    else
        all(ishomogenous, F)
    end
end


"""
    degree(term::MP.AbstractTermLike, vars)

Compute the (weighted) degree of `f` in the variables `vars`.
"""
function degree(term::MP.AbstractTermLike, variables::Vector{<:MP.AbstractVariable})
    sum(MP.degree(term, v) for v in variables)
end
function degree(term::MP.AbstractTermLike, weighted_variables::Vector{<:WeightedVariable})
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

"""
    uniquevar(f::MP.AbstractPolynomialLike, tag=:x0)
    uniquevar(F::Vector{<:MP.AbstractPolynomialLike}, tag=:x0)

Creates a unique variable.
"""
uniquevar(f::MP.AbstractPolynomialLike, tag=:x0) = MP.similarvariable(f, gensym(tag))
uniquevar(F::Vector{<:MP.AbstractPolynomialLike}, tag=:x0) = uniquevar(F[1], tag)

"""
    homogenize(f::MP.AbstractPolynomial, variable=uniquevar(f))

Homogenize the polynomial `f` by using the given variable `variable`.

    homogenize(F::Vector{<:MP.AbstractPolynomial}, variable=uniquevar(F))

Homogenize each polynomial in `F` by using the given variable `variable`.
"""
function homogenize(f::MP.AbstractPolynomialLike, var::MP.AbstractVariable=uniquevar(f))
    d = MP.maxdegree(f)
    MP.polynomial(map(t -> var^(d - MP.degree(t)) * t, MP.terms(f)))
end
function homogenize(F::Vector{<:MP.AbstractPolynomialLike}, var::MP.AbstractVariable=uniquevar(F); parameters=nothing)
    if parameters !== nothing
        homogenize(F, setdiff(MP.variables(F), parameters), var)
    else
        homogenize.(F, Ref(var))
    end
end

"""
    homogenize(f::MP.AbstractPolynomial, v::Vector{<:MP.AbstractVariable}, variable=uniquevar(f))

Homogenize the variables `v` in the polynomial `f` by using the given variable `variable`.

    homogenize(F::Vector{<:MP.AbstractPolynomial}, v::Vector{<:MP.AbstractVariable}, variable=uniquevar(F))

Homogenize the variables `v` in each polynomial in `F` by using the given variable `variable`.
"""
function homogenize(f::MP.AbstractPolynomialLike, variables::Vector, var::MP.AbstractVariable=uniquevar(f))
    _, d_max = minmaxdegree(f, variables)
    MP.polynomial(map(f) do t
        d = degree(t, variables)
        var^(d_max - d)*t
    end)
end
function homogenize(F::Vector{<:MP.AbstractPolynomialLike}, variables::Vector, var::MP.AbstractVariable=uniquevar(F))
    map(f -> homogenize(f, variables, var), F)
end
