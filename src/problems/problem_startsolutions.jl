export problem_startsolutions

const supported_keywords = [[:seed, :homvar, :homotopy, :system]; Input.supported_keywords]
const DEFAULT_SYSTEM = FPSystem
const DEFAULT_HOMOTOPY = StraightLineHomotopy

"""
    problem_startsolutions(F; options...)
    problem_startsolutions(G, F, startsolutions; options...)
    problem_startsolutions(H::AbstractHomotopy, startsolutions; options...)
    problem_startsolutions(prob::TotalDegreeProblem; options...)
    problem_startsolutions(prob::StartTargetProblem; options...)

    Construct the problem and if necessary the startsolutions. This steps
    constructs a homotopy and homogenizes the systems if necessary.

    The `options` are
    * `seed::Int`: Random seed used in the construction.
    * `system=FPSystem`: A constructor to assemble a [`Systems.AbstractSystem`](@ref). The constructor
    is called with `system(polys, variables)` where `variables` determines the variable ordering.
    * `homotopy=StraightLineHomotopy`: A constructor to construct a [`Homotopies.AbstractHomotopy`](@ref) an `Systems.AbstractSystem`. The constructor
    is called with `homotopy(start, target)` where `start` and `target` are systems constructed
    with `system`.
"""
function problem_startsolutions end

function problem_startsolutions(args...; seed=randseed(), kwargs...)
    Random.seed!(seed)
    supported, rest = Utilities.splitkwargs(kwargs, Input.supported_keywords)
    input = Input.input(args...; supported...)
    problem_startsolutions(input, seed; rest...)
end

function problem_startsolutions(input::AbstractInput; seed=randseed(), kwargs...)
    problem_startsolutions(input, seed; kwargs...)
end
function problem_startsolutions(input::AbstractInput, seed::Int;
	homvar::Union{Nothing, Int, MP.AbstractVariable}=nothing, kwargs...)
    problem_startsolutions(input, homvar, seed; kwargs...)
end

function problem_startsolutions(input::Input.Homotopy, homvar, seed; kwargs...)
    Projective(input.H, homogenization(homvar), seed), input.startsolutions
end


const overdetermined_error_msg = """
The input system is overdetermined. Therefore it is necessary to provide an explicit start system.
See
    https://www.JuliaHomotopyContinuation.org/guides/latest/overdetermined_tracking/
for details.
"""

##############
# TOTALDEGREE
##############

function problem_startsolutions(prob::TotalDegree{<:MPPolyInputs}, _homvar::Nothing, seed::Int; system=DEFAULT_SYSTEM, kwargs...)
    F, vars, homogenization = homogenize_if_necessary(prob.system)
    if ishomogenized(homogenization)
		# Check overdetermined case
		length(F) ≥ length(vars) && error(overdetermined_error_msg)

        proj = Projective(
            Systems.TotalDegreeSystem(prob.degrees, vars, vars[homogenization.homvaridx]),
            construct_system(system, F, vars), homogenization, seed; kwargs...)
        proj, totaldegree_solutions(prob.degrees, homogenization)
    else
		# Check overdetermined case
		length(F) > length(vars) && error(overdetermined_error_msg)

        G = Systems.TotalDegreeSystem(prob.degrees, vars, vars[1])
        start = totaldegree_solutions(prob.degrees, NullHomogenization())
        Projective(G, construct_system(system, F), NullHomogenization(), seed; kwargs...), start
    end
end

function problem_startsolutions(prob::TotalDegree{<:MPPolyInputs},
    homvar::MP.AbstractVariable, seed; system=DEFAULT_SYSTEM, kwargs...)

    if !ishomogenous(prob.system)
        error("Input system is not homogenous although `homvar=$(homvar)` was passed.")
    end
    F, vars, homogenization = homogenize_if_necessary(prob.system; homvar=homvar)
	# Check overdetermined case
	length(F) > length(vars) && error(overdetermined_error_msg)

    start = totaldegree_solutions(prob.degrees, homogenization)
    proj = Projective(
        Systems.TotalDegreeSystem(prob.degrees, vars, homvar),
        construct_system(system, F, vars), homogenization, seed; kwargs...)
    proj, start
end

function problem_startsolutions(prob::TotalDegree{<:AbstractSystem}, homvaridx::Nothing, seed; system=DEFAULT_SYSTEM, kwargs...)
    n, N = size(prob.system)
    G = Systems.TotalDegreeSystem(prob.degrees, collect(2:N), 1)
	# Check overdetermined case
	n > N && error(overdetermined_error_msg)

    (Projective(G, prob.system, NullHomogenization(), seed; kwargs...),
     totaldegree_solutions(prob.degrees, NullHomogenization()))
end

function problem_startsolutions(prob::TotalDegree{<:AbstractSystem}, homvaridx::Int, seed; system=DEFAULT_SYSTEM, kwargs...)
    n, N = size(prob.system)

    homogenization = Homogenization(homvaridx)
    G = Systems.TotalDegreeSystem(prob.degrees, [1:homvaridx-1;homvaridx+1:N], homvaridx)

    (Projective(G, prob.system, homogenization, seed; kwargs...),
     totaldegree_solutions(prob.degrees, homogenization))
end


###############
# START TARGET
###############

function problem_startsolutions(prob::StartTarget{<:MPPolyInputs, <:MPPolyInputs}, homvar, seed; system=DEFAULT_SYSTEM, kwargs...)
    F, G = prob.target, prob.start
    F_ishom, G_ishom = ishomogenous.((F, G))
    if F_ishom && G_ishom && homvar !== nothing
        Projective(construct_system(system, G),
			construct_system(system, F),
			Homogenization(homvar, variables(F)), seed; kwargs...), prob.startsolutions
    elseif F_ishom && G_ishom && homvar === nothing
        Projective(construct_system(system, G),
			construct_system(system, F),
			NullHomogenization(), seed; kwargs...), prob.startsolutions
    elseif F_ishom || G_ishom
        error("One of the input polynomials is homogenous and the other not!")
    else
        if homvar !== nothing
            error("Input system is not homogenous although `homvar` was passed.")
        end
        homvar = uniquevar(F)
        homogenization = Homogenization(1)
        var_ordering = [homvar; variables(F)]
        Gₕ = construct_system(system, homogenize(G, homvar), var_ordering)
        Fₕ = construct_system(system, homogenize(F, homvar), var_ordering)
        Projective(Gₕ, Fₕ, homogenization, seed; kwargs...), prob.startsolutions
    end
end

#####################
# Parameter homotopy
#####################

function problem_startsolutions(prob::ParameterSystem, homvar, seed; system=FPSystem, kwargs...)
    F, vars, homogenization = homogenize_if_necessary(prob.system, homvar=homvar, parameters=prob.parameters)
    H = ParameterHomotopy(construct_system(system, SPSystem, F, variables=vars, parameters=prob.parameters),
						  p₁=prob.p₁, p₀=prob.p₀, γ₁=prob.γ₁, γ₀=prob.γ₀)

    Projective(H, homogenization, seed), prob.startsolutions
end

##########
# HELPERS
##########

construct_system(constructor, F::MPPolys; kwargs...) = constructor(F; kwargs...)
construct_system(constructor, F::MPPolys, vars) = constructor(F, vars)
construct_system(constructor, F::Composition, vars) = Systems.CompositionSystem(F, vars, constructor)
construct_system(constructor, F::Composition; kwargs...) = Systems.CompositionSystem(F, constructor; kwargs...)

construct_system(constructor, final_constructor, F::MPPolys; kwargs...) = final_constructor(F; kwargs...)
construct_system(constructor, final_constructor, F::Composition; kwargs...) = Systems.CompositionSystem(F, vars, constructor, final_constructor; kwargs...)


"""
    homogenize_if_necessary(F::Vector{<:MP.AbstractPolynomialLike})

Homogenizes the system `F` if necessary and returns the (new) system `F` its variables
and a subtype of [`AbstractHomogenization`] indicating whether it was homegenized.
If it was homogenized and no then the new variable is the **first one**.
"""
function homogenize_if_necessary(F::MPPolyInputs; homvar=nothing, parameters=nothing)
    vars = variables(F)
    if parameters !== nothing
        setdiff!(vars, parameters)
    end

    n, N = length(F), length(vars)
    if ishomogenous(F; parameters=parameters)
        # N = n+1 is the only valid size configuration
        if n + 1 > N
            error(overdetermined_error_msg)
        end
        if homvar === nothing
            F, vars, NullHomogenization()
        else
            F, vars, Homogenization(homvar, vars)
        end
    else
        if homvar !== nothing
            error("Input system is not homogenous although `homvar` was passed.")
        end
        # We create a new variable to homogenize the system
        homvar = uniquevar(F)
        push!(vars, homvar)
        sort!(vars, rev=true)

        homogenize(F, homvar; parameters=parameters), vars, Homogenization(homvar, vars)
    end
end
