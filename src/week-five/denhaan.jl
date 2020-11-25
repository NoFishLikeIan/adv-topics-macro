include("markov/state.jl")

using Parameters


const Π = [
    0.525 0.35 0.03125 0.09375;
    0.038889 0.836111 0.002083 0.122917;
    0.09375 0.03125 0.291667 0.583333; 
    0.009115 0.115885 0.024306 0.850694
] # exogenous transition matrix as specified in the paper


cartesianstate(s1, s2) = vec(collect(Iterators.product(s1, s2)))


"""
Specification for the Aggregate Uncertainty model from
    Den Haan, Judd, and Juillard (2010)
"""
struct DenHaanModel

    ζ::StateMarkov
    Δ::Float64

    β::Float64 
    γ::Float64  
    α::Float64  
    δ::Float64  
    l::Float64 
    μ::Float64

    function DenHaanModel(
        Δ::Float64=0.01,
        β::Float64=0.99,
        γ::Float64=1.,
        α::Float64=0.36,
        δ::Float64=0.025,
        l::Float64=10 / 9,
        μ::Float64=0.15)
        
        S = cartesianstate((1 - Δ, 1 + Δ), (0, 1)) # Cartesian product of ϵ and s

        process = StateMarkov(S, Π)

        return new(process, Δ, β, γ, α, δ, l, μ)
    end
end

"""
Construct the production functions
"""
function makeproduction(model::DenHaanModel)
    @unpack α, l = model

    w(K, L, z) = (1 - α) * z * (K / (L * l))^α
    r(K, L, z) = α * z * (K / (L * l))^(α - 1) # FIXME: Does this require R = 1 + f′ - δ?

    return w, r
end


"""
Construct the utility functions 
"""
function makeutility(model::DenHaanModel)
    @unpack γ = model
    
    u(c) = γ == 1 ? log(c) : (c^(1 - γ) - 1) / (1 - γ)
    u′(c) = c^(-γ)
    inv_u′(x) = x^(- 1 / γ)

    return u, u′, inv_u′
end