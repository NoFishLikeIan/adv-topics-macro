include("../../week-two/markov/simulation.jl")

using Distributions

State = Union{Any,Tuple{Vararg{Any}}} # make variable type

struct StateMarkov
    S::Vector{State}
    P::Matrix{Float64}
end

function markovstate_simulation(
    markov::StateMarkov, T::Int; 
    drop::Int=0, s0::State=nothing
)::Vector{State}

    if drop > T throw("Variables to drop ($drop) > total ($total)") end
    if s0 ∉ markov.S throw("Initial state not in state space") end
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
