include("../../week-two/markov/simulation.jl")

using Distributions, Parameters
using Base.Threads

State = Union{Any,Tuple{Vararg{Any}}} # make variable type

struct StateMarkov
    S::Vector{State}
    P::Matrix{Float64}
end

function simulation(
    markov::StateMarkov, T::Int; 
    drop::Int=0, s0::State=nothing)::Vector{State}

    if drop > T throw("Variables to drop ($drop) > total ($total)") end
    if !isnothing(s0) && s0 ∉ markov.S throw("Initial state not in state space") end
    F_inv = make_inv_F(pseudoF(markov.P))
    
    evolution_jdx = zeros(Int, T)
    evolution_jdx[1] = isnothing(s0) ? 
        rand(1:length(markov.S)) : 
        findfirst(==(s0), markov.S)

    for (t, u) in enumerate(rand(Uniform(), T - 1))
        evolution_jdx[t + 1] = F_inv(u, evolution_jdx[t])
    end

    start = drop + 1
    sim = [markov.S[j] for j in evolution_jdx[start:end]]

    return sim
end

"""
Simulate conditionally on an aggregate vector of shock
"""
function conditional_simulation(
    jointmarkov::StateMarkov, aggshock::Vector{State}, N::Int;
    drop=0)::Matrix{Float64}

    @unpack S, P = jointmarkov

    S_z = unique([st[1] for st in S])

    T = length(aggshock)

    # Construct an array of π(z′ | ϵ′, z, ϵ) conditional transition matrices
    rzeros = findall(st -> st[2] == 0, S) 
    rones = findall(st -> st[2] == 1, S)

    Ns = length(rzeros)
    shock_Ps = Array{Float64}(undef, T - 1, Ns, Ns)

    # TODO: Check the order of the simulation
    @threads for t in 2:T
        from = aggshock[t - 1] == 1 ? rones : rzeros
        to = aggshock[t] == 1 ? rones : rzeros

        shock_Ps[t - 1, :, :] = rownormal(P[from, to])
    end

    # Simulate a process for z, long N
    sim_jdx = Array{Int}(undef, T, N)
    sim_jdx[1, :] = rand(1:Ns, N)

    for (t, sampled_u) in enumerate(eachrow(rand(Uniform(), T - 1, N)))
        thisP = shock_Ps[t, :, :]
        F_inv = make_inv_F(pseudoF(thisP))

        sim_jdx[t + 1, :] = @. F_inv(sampled_u, sim_jdx[t, :])
    end

    # FIXME: Drop should happen before
    sim = (j -> S_z[j]).(sim_jdx[drop + 1:end, :]) 

    return sim
    
end