using LinearAlgebra, Distributions, Random, StatsBase

include("markov-types.jl")

# TODO: Use stationary distribution to compute the moments

mutable struct MarkovDiscrete
    P::Matrix{Float64}
    S::Partition
    transformation::Function
end

"""
Construct the psuedo cdf, based on P matrix
"""
function pseudoF(P::Matrix{Float64})
    N, N = size(P)
    F = copy(P)

    for n in 2:N
        F[:, n] = F[:, n] + F[:, n - 1]
    end

    return F
end

"""
Constructs the inverse of the pseudo F matrix.
"""
function make_inv_F(pseudoF::Matrix{Float64})

    """
    Takes the partition index of z and a probability value. 
    Returns the originating z' partition index.
    """
    function inv(u::Float64, z_idx::Int)::Int
        cond_F = pseudoF[z_idx, :]

        for (jdx, prob) in enumerate(cond_F)
            if u < prob
                return jdx
            end
        end

        return length(cond_F)
    end
    
    return inv
end

function discrete_sim(markov::MarkovDiscrete; 
    T::Int=2000,
    drop::Int=500)

    if drop > T error("Variables to drop ($drop) > total ($total)") end

    inverse_F = make_inv_F(pseudoF(markov.P))

    evolution_jdx = zeros(Int, T)
    evolution_jdx[1] = rand(1:length(markov.S))

    for (t, u) in enumerate(rand(Uniform(), T - 1))
        evolution_jdx[t + 1] = inverse_F(u, evolution_jdx[t])
    end

    sim = currygetindex(markov.S).(evolution_jdx[drop:end])
    
    return markov.transformation.(sim)
end

function summary_stats(markov::MarkovDiscrete; T=10_000)
    path = discrete_sim(markov; T=T)

    est = Dict(
        "μ" => mean(path),
        "ν" => var(path),
        "ρ" => autocor(path)
    )

    return est
end