import StaticPolynomials
const SP = StaticPolynomials

export SPSystem

"""
    SPSystem(polynomials, vars) <: AbstractSystem

Create a system using the [`StaticPolynomials`](https://github.com/JuliaAlgebra/StaticPolynomials.jl) package.
Note that `StaticPolynomials` leverages Julias metaprogramming capabilities to automatically generate
functions to evaluate the system and its Jacobian. These generated functions are *very fast* but
at the cost of possibly large compile times. The compile time depends on the size of the support of the polynomial system.
If you intend to solve a large system or you need to solve a system
with the *same support* but different coefficients even large compile times can be worthwile.
As a general rule of thumb this usually is twice as fast as solving the same system using [`FPSystem`](@ref).

## Example
You can use `SPSystem` as follows with solve
```julia
@polyvar x y
F = [x^2+3y^4-2, 2y^2+3x*y+4]
solve(F, system=SPSystem)
```
"""
struct SPSystem{S<:SP.PolynomialSystem} <: AbstractSystem
    system::S
end

SPSystem(polys::Vector{<:MP.AbstractPolynomial}, vars) = SPSystem(SP.PolynomialSystem(polys, variables=vars))
SPSystem(polys::Vector{<:MP.AbstractPolynomial}; kwargs...) = SPSystem(SP.PolynomialSystem(polys; kwargs...))

Base.size(F::SPSystem) = (SP.npolynomials(F.system), SP.nvariables(F.system))

cache(F::SPSystem, x, p=nothing) = NullCache()
Base.@propagate_inbounds evaluate!(u, F::SPSystem, x, ::NullCache) = SP.evaluate!(u, F.system, x)
Base.@propagate_inbounds evaluate!(u, F::SPSystem, x, p, ::NullCache) = SP.evaluate!(u, F.system, x, p)
evaluate(F::SPSystem, x, ::NullCache) = SP.evaluate(F.system, x)
evaluate(F::SPSystem, x, p, ::NullCache) = SP.evaluate(F.system, x, p)
Base.@propagate_inbounds jacobian!(U, F::SPSystem, x, ::NullCache) = SP.jacobian!(U, F.system, x)
Base.@propagate_inbounds jacobian!(U, F::SPSystem, x, p, ::NullCache) = SP.jacobian!(U, F.system, x, p)
jacobian(F::SPSystem, x, ::NullCache) = SP.jacobian(F.system, x)
jacobian(F::SPSystem, x, p, ::NullCache) = SP.jacobian(F.system, x, p)
Base.@propagate_inbounds evaluate_and_jacobian!(u, U, F::SPSystem, x, ::NullCache) = SP.evaluate_and_jacobian!(u, U, F.system, x)
Base.@propagate_inbounds evaluate_and_jacobian!(u, U, F::SPSystem, x, p, ::NullCache) = SP.evaluate_and_jacobian!(u, U, F.system, x, p)
evaluate_and_jacobian(F::SPSystem, x, ::NullCache) = SP.evaluate_and_jacobian(F.system, x)
evaluate_and_jacobian(F::SPSystem, x, p, ::NullCache) = SP.evaluate_and_jacobian(F.system, x, p)
Base.@propagate_inbounds differentiate_parameters!(U, F::SPSystem, x, p, ::NullCache) = SP.differentiate_parameters!(U, F.system, x, p)
differentiate_parameters(F::SPSystem, x, p, ::NullCache) = SP.differentiate_parameters(F.system, x, p)
