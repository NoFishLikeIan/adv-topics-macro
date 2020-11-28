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
    model::DenHaanModel, aggshock::Vector{State}, N::Int;
    drop=0)::Matrix{Float64}

    @unpack S_ϵ, S_z, ζ = model
    @unpack S, P = ζ

    T = length(aggshock)

    # Construct an array of π(ϵ′ | z, ϵ, z′) conditional transition matrices
    recession = findall(st -> st[1] == S_z[1], S) 
    boom = findall(st -> st[1] == S_z[2], S)

    Ns = length(recession)
    shock_Ps = Array{Float64}(undef, T - 1, Ns, Ns)

    @threads for t in 2:T
        from = aggshock[t - 1] == S_z[1] ? recession : boom
        to = aggshock[t] == S_z[1] ? recession : boom

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

    sim = sim_jdx[drop + 1:end, :] .- 1

    return sim
    
end