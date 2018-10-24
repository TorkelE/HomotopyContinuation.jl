module Utilities

import LinearAlgebra
import MultivariatePolynomials
const MP = MultivariatePolynomials

export ldiv_lu!, blas_ldiv_lu!,
    infinity_norm,
    unsafe_infinity_norm,
    fubini_study,
    logabs,
    batches,
    randomish_gamma,
    filterkwargs,
    splitkwargs,
    solve!,
    set_num_BLAS_threads,
    get_num_BLAS_threads,
    randseed,
    check_kwargs_empty,
    start_solution_sample

include("utilities/polynomial.jl")

"""
    check_kwargs_empty(kwargs, [allowed_kwargs])

Chack that the list of `kwargs` is empty. If not, print all unsupported keywords
with their arguments.
"""
function check_kwargs_empty(kwargs, allowed_kwargs=[])
    if !isempty(kwargs)
        msg = "Unexpected keyword argument(s): "
        first_el = true
        for kwarg in kwargs
            if !first_el
                msg *= ", "
            end
            msg *= "$(first(kwarg))=$(last(kwarg))"
            first_el = false
        end
        if !isempty(allowed_kwargs)
            msg *= "\nAllowed keywords are\n"
            msg *= join(allowed_kwargs, ", ")
        end
        throw(ErrorException(msg))
    end
end



"""
    randseed(range=1_000:1_000_000)

Return a random seed in the range `range`.
"""
randseed(range=1_000:1_000_000) = rand(range)

"""
    solve!(A, b)

Solve ``Ax=b`` inplace. This overwrites `A` and `b`
and stores the result in `b`.
"""
function solve!(A::StridedMatrix, b::StridedVecOrMat)
    m, n = size(A)
    if m == n
        lusolve!(A, b)
    else
        LinearAlgebra.ldiv!(LinearAlgebra.qr!(A), b)
    end
    b
end

# This is an adoption of LinearAlgebra.generic_lufact!
# with 3 changes:
# 1) For choosing the pivot we use abs2 instead of abs
# 2) Instead of using the robust complex division we
#    just use the naive division. This shouldn't
#    lead to problems since numerical diffuculties only
#    arise for very large or small exponents
# 3) We fold lu! and ldiv! into one routine
#    this has the effect that we do not need to allocate
#    the pivot vector anymore and also avoid the allocations
#    coming from the LU wrapper
function lusolve!(A::AbstractMatrix{T}, b::Vector{T}, ::Val{Pivot} = Val(true)) where {T,Pivot}
    m, n = size(A)
    minmn = min(m,n)
    # LU Factorization
    @inbounds begin
        for k = 1:minmn
            # find index max
            kp = k
            if Pivot
                amax = zero(real(T))
                for i = k:m
                    absi = abs2(A[i,k])
                    if absi > amax
                        kp = i
                        amax = absi
                    end
                end
            end
            if !iszero(A[kp,k])
                if k != kp
                    # Interchange
                    for i = 1:n
                        tmp = A[k,i]
                        A[k,i] = A[kp,i]
                        A[kp,i] = tmp
                    end
                    b[k], b[kp] = b[kp], b[k]
                end
                # Scale first column
                Akkinv = @fastmath inv(A[k,k])
                for i = k+1:m
                    A[i,k] *= Akkinv
                end
            end
            # Update the rest
            for j = k+1:n
                for i = k+1:m
                    A[i,j] -= A[i,k]*A[k,j]
                end
            end
        end
    end
    # now forward and backward substitution
    ldiv_unit_lower!(A, b)
    ldiv_upper!(A, b)
    b
end
@inline function ldiv_upper!(A::AbstractMatrix, b::AbstractVector, x::AbstractVector = b)
    n = size(A, 2)
    for j in n:-1:1
        @inbounds iszero(A[j,j]) && throw(LinearAlgebra.SingularException(j))
        @inbounds xj = x[j] = (@fastmath A[j,j] \ b[j])
        for i in 1:(j-1)
            @inbounds b[i] -= A[i,j] * xj
        end
    end
    b
end
@inline function ldiv_unit_lower!(A::AbstractMatrix, b::AbstractVector, x::AbstractVector = b)
    n = size(A, 2)
    @inbounds for j in 1:n
        xj = x[j] = b[j]
        for i in j+1:n
            b[i] -= A[i,j] * xj
        end
    end
    x
end


"""
    infinity_norm(z)

Compute the ∞-norm of `z`. If `z` is a complex vector this is more efficient
than `norm(z, Inf)`.

    infinity_norm(z₁, z₂)

Compute the ∞-norm of `z₁-z₂`.
"""
infinity_norm(z::AbstractVector{<:Complex}) = sqrt(maximum(abs2, z))
function infinity_norm(z₁::AbstractVector{<:Complex}, z₂::AbstractVector{<:Complex})
    m = abs2(z₁[1] - z₂[1])
    n₁, n₂ = length(z₁), length(z₂)
    if n₁ ≠ n₂
        return convert(typeof(m), Inf)
    end
    @inbounds for k=2:n₁
        m = max(m, abs2(z₁[k] - z₂[k]))
    end
    sqrt(m)
end
unsafe_infinity_norm(v, w) = infinity_norm(v, w)


"""
    fubini_study(x, y)

Computes the Fubini-Study norm of `x` and `y`.
"""
fubini_study(x,y) = acos(min(1.0, abs(LinearAlgebra.dot(x,y))))

"""
    logabs(z)

The log absolute map `log(abs(z))`.
"""
logabs(z::Complex) = 0.5 * log(abs2(z))
logabs(x) = log(abs(x))


function randomish_gamma()
    # Usually values near 1, i, -i, -1 are not good randomization
    # Therefore we artificially constrain the choices
    theta = rand() * 0.30 + 0.075 + (rand(Bool) ? 0.0 : 0.5)
    cis(2π * theta)
end

"""
    filterkwargs(kwargs, allowed_kwargs)

Remove all keyword arguments out of `kwargs` where the keyword is not contained
in `allowed_kwargs`.
"""
function filterkwargs(kwargs, allowed_kwargs)
    [kwarg for kwarg in kwargs if any(kw -> kw == first(kwarg), allowed_kwargs)]
end

"""
    splitkwargs(kwargs, supported_keywords)

Split the vector of `kwargs` in two vectors, the first contains all `kwargs`
whose keywords appear in `supported_keywords` and the rest the other one.
"""
function splitkwargs(kwargs, supported_keywords)
    supported = []
    rest = []
    for kwarg in kwargs
        if any(kw -> kw == first(kwarg), supported_keywords)
            push!(supported, kwarg)
        else
            push!(rest, kwarg)
        end
    end
    supported, rest
end


start_solution_sample(xs) = first(xs) |> promote_start_solution
start_solution_sample(x::AbstractVector{<:Number}) = promote_start_solution(x )

promote_start_solution(x::AbstractVector{ComplexF64}) = x
function promote_start_solution(x)
    x_new =similar(x, promote_type(eltype(x), ComplexF64), length(x))
    copyto!(x_new, x)
    x_new
end


# Parallelization

set_num_BLAS_threads(n) = LinearAlgebra.BLAS.set_num_threads(n)
get_num_BLAS_threads() = convert(Int, _get_num_BLAS_threads())
# This is into 0.7 but we need it for 0.6 as well
const _get_num_BLAS_threads = function() # anonymous so it will be serialized when called
    blas = LinearAlgebra.BLAS.vendor()
    # Wrap in a try to catch unsupported blas versions
    try
        if blas == :openblas
            return ccall((:openblas_get_num_threads, Base.libblas_name), Cint, ())
        elseif blas == :openblas64
            return ccall((:openblas_get_num_threads64_, Base.libblas_name), Cint, ())
        elseif blas == :mkl
            return ccall((:MKL_Get_Max_Num_Threads, Base.libblas_name), Cint, ())
        end

        # OSX BLAS looks at an environment variable
        if Sys.isapple()
            return ENV["VECLIB_MAXIMUM_THREADS"]
        end
    catch
    end

    return nothing
end

end
