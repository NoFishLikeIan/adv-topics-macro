using Distributions

include("../week-two/markov/simulation.jl")
include("../week-two/markov/process.jl")


struct Aiyagari 
    y::MarkovDiscrete
    a_::Float64
    β::Float64
    α::Float64
    δ::Float64
    σ_u::Float64
    
    # Inner constructor for the Markov discretization of an AR(1)
    function Aiyagari(ρ::Float64, σ_inn::Float64, N::Int, params...)
        σ_y = σ_inn / √(1 - ρ^2)
        ar = Process(Normal(0, σ_y), ρ, Normal(0, σ_inn))

        P, S = tauchen(ar, N)
        y = MarkovDiscrete(P, S, z -> z)

        return new(y, params...)
    end

end

"""
Returns utility and first derivative based on Aiyagari parameters
"""
function make_u(ai::Aiyagari)::Tuple{Function,Function,Function}
    esp = 1 - ai.σ_u

    u(c::Float64)::Float64 = (c^esp - 1) / esp
    u′(c::Float64)::Float64 = c^(-ai.σ_u)
    invu′(x::Float64)::Float64 = x^(-1 / ai.σ_u)

    return u, u′, invu′
end
