using DynamicPolynomials, HomotopyContinuation
import MultivariatePolynomials
const MP = MultivariatePolynomials;
using Bertini

function sclgen(F)
    A = assemble_lhs(F)
    b = assemble_rhs(F)
    w = exp10.(A \ b)
    n = div(length(w), 2)
    w[1:n], w[n+1:end]
end

function assemble_lhs(F)
    vars = MP.variables(F)
    n = length(vars)
    A = zeros(2n, 2n)
    T = map(MP.terms, F)
    for (r, f) in enumerate(F)
        A[r, r] = length(T[r])
        for s=1:n
            A[r, n+s] = sum(t -> MP.degree(t, vars[s]), T[r])
        end
        for s=1:n
            A[n+r, s] = sum(t -> MP.degree(t, vars[r]), T[s])
        end
        for s=1:n
            A[n+r, n+s] = sum(T) do T_i
                sum(t -> MP.degree(t, vars[s]) * MP.degree(t, vars[r]), T_i)
            end
        end
    end
    A
end

function assemble_rhs(F)
    vars = MP.variables(F)
    n = length(vars)
    b = zeros(2n)
    for r in 1:length(F)
        b[r] = -sum(log10.(abs.(MP.coefficients(F[r]))))
        b[n+r] = -sum(F) do f
            sum(t -> log10(abs(MP.coefficient(t))) * MP.degree(t, vars[r]), MP.terms(f))
        end
    end
    b
end

vars = @polyvar w w2 w2v v w2v2 vP sB w2sB vPp phos
F = [(-1 * 0.7 * w + -2 * 3600.0 * (w ^ 2 / 2) + 2 * 18.0 * w2)*(0.2 + sB) + 4.0 * 0.4 * (1 + 30.0sB),
 -1 * 0.7 * w2 + 3600.0 * (w ^ 2 / 2) + -1 * 18.0 * w2 + -1 * 3600.0 * w2 * v + 18.0w2v + 36.0w2v + -1 * 3600.0 * w2 * sB + 18.0w2sB,
 -1 * 0.7 * w2v + 3600.0 * w2 * v + -1 * 18.0 * w2v + -1 * 3600.0 * w2v * v + 18.0w2v2 + -1 * 36.0 * w2v + 36.0w2v2 + 1800.0 * w2sB * v + -1 * 1800.0 * w2v * sB,
 (-1 * 0.7 * v + -1 * 3600.0 * w2 * v + 18.0w2v + -1 * 3600.0 * w2v * v + 18.0w2v2 + -1 * 1800.0 * w2sB * v + 1800.0 * w2v * sB + 180.0vPp)*(0.2 + sB) + 4.5 * 0.4 * (1 + 30.0sB),
 -1 * 0.7 * w2v2 + 3600.0 * w2v * v + -1 * 18.0 * w2v2 + -1 * 36.0 * w2v2,
 -1 * 0.7 * vP + 36.0w2v + 36.0w2v2 + -1 * 3600.0 * vP * phos + 18.0vPp,
 (-1 * 0.7 * sB + -1 * 3600.0 * w2 * sB + 18.0w2sB + 1800.0 * w2sB * v + -1 * 1800.0 * w2v * sB)*(0.2 + sB) + 0.4 * (1 + 30.0sB),
 -1 * 0.7 * w2sB + 3600.0 * w2 * sB + -1 * 18.0 * w2sB + -1 * 1800.0 * w2sB * v + 1800.0 * w2v * sB,
 -1 * 0.7 * vPp + 3600.0 * vP * phos + -1 * 18.0 * vPp + -1 * 180.0 * vPp,
 (phos + vPp) - 2.0]

scales, substitutions = sclgen(F)

G = DynamicPolynomials.subs(F, vars => substitutions .* vars)

prob, starts = Problems.problem_startsolutions(F, seed=12345)
S = collect(starts)

solver = Solving.Solver(prob, S[1], 1.0, 0.0, maxiters=100, tol=1e-6)

solver
@time solve(solver, S, threading=false)

solver












tracker, starts = pathtracker_startsolutions(F; seed=12345)

S = collect(starts)
r = PathTracking.track(tracker, S[1224], 1.0, 0.0)

y = r.x.data

 y[1:10] ./ (substitutions .* y[11])
abs2.(y)

using LinearAlgebra

H = tracker.homotopy
Homotopies.evaluate(H, y ./ y[11], 0.0)
y

solve(F, threading=false)
