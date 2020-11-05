using Optim, LinearAlgebra

include("utils.jl")

include("../process.jl")
include("../markov-types.jl")
include("../simulation.jl")


"""
Computes the parametric (p, q) -> P matrix in a recursive manner.
FIXME: This is very inefficient both in memory and allocation time.
"""
function P_N(N)
    if N == 2
        return (p, q) -> [p (1 - p); q (1 - q)]
    end

    P_p = P_N(N - 1)  

    function nested_P(p, q)
        Ψ = P_p(p, q)
        n, n = size(Ψ)

        Ω = zeros((n + 1, n + 1))

        # TODO: Find a nicer way to do this
        Ω[1:n, 1:n] += p * Ψ
        Ω[1:n, end - n + 1:end] += (1 - p) * Ψ
        Ω[end - n + 1:end, 1:n] +=  (1 - q) * Ψ
        Ω[end - n + 1:end, end - n + 1:end] +=  q * Ψ

        return Ω
    end


    return nested_P
end

function rouwenhorst(proc::Process, N::Int; x0=[3., 0.5, 0.5])

    if N > 100 @warn "Current implementation is too slow for big matrices" end

    P = P_N(N)

    function sim_moments(ψ::Float64, p::Float64, q::Float64)
        markov = MarkovDiscrete(
            colnormalized(P(p, q)), 
            makepartition(ψ, N)
        )

        sim = discrete_sim(markov)

        return [mean(sim), var(sim)]
    end

    moments = [mean(proc), var(proc)]

    function minimize(params::Array{Float64}) 
        ψ, p, q = params
        return norm(sim_moments(ψ, p, q) - moments)
    end

    res = optimize(minimize, x0, LBFGS())

    if Optim.converged(res)
        minimum = Optim.minimum(res)
        if minimum > 0.05 @warn "Matching the moments led to a minimum of $minimum > .05" end

        ψ, p, q = Optim.minimizer(res)

        return colnormalized(P(p, q)), makepartition(ψ, N)
    
    else 
        error("Method of moments matching did not converge")
    end
end