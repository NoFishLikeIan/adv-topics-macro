using Distributions, Random

include("../process.jl")
include("../markov-types.jl")
include("utils.jl")

"""
Discretize a Process based on equidistant points (mode = "equi") or equal probabilities (mode = "imp").
"""
function tauchen(proc::Process, N::Int, m::Int; mode="equi")
    F = cdf(proc)

    ψ = m * sqrt(var(proc))
    
    partition = makepartition(ψ, N)
    d = distance(partition)

    function scale_f(z::Real) 
        forward = z_p::Real -> F((z + d - proc.evol(z_p)) / var(proc.error))
        backward = z_p::Real -> F((z - d - proc.evol(z_p)) / var(proc.error))

        return forward, backward
    end
    

    P = zeros((N, N))

    for i in 1:N
        f, b = scale_f(partition[i])
        P[:, i] = f.(partition) - b.(partition)
    end

    return P, partition
end


