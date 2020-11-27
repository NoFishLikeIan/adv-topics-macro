using Parameters, Kronecker

# exogenous transition matrices as specified in the paper
const P_z = [
    0.875 (1 - 0.875); 
    0.15625 (1 - 0.15625)]

const P_ϵ = [
    0.6 0.4; 
    2 / 45 43 / 45
] 

const Π = P_z ⊗ P_ϵ

cartesianstate(s1, s2) = vec(collect(Iterators.product(s1, s2)))
rownormal(M::Matrix{Float64}) = M ./ sum(M, dims=2)


"""
Specification for the Aggregate Uncertainty model from
    Den Haan, Judd, and Juillard (2010)
"""
struct DenHaanModel
    Ε::StateMarkov
    ζ::StateMarkov
    S_z::Tuple{Vararg{Float64}}
    S_ϵ::Tuple{Vararg{Int}}

    β::Float64 
    γ::Float64  
    α::Float64  
    δ::Float64  
    l::Float64 
    μ::Float64

    function DenHaanModel(
        Δ=0.01,
        β=0.99,
        γ=1.,
        α=0.36,
        δ=0.025,
        l=10 / 9,
        μ=0.15,
        S_ϵ=(0, 1))

        S_z = (1 - Δ, 1 + Δ)
        
        S = cartesianstate(S_z, S_ϵ) # Cartesian product of ϵ and s

        joint_process = StateMarkov(S, Π)
        ϵ_process = StateMarkov(S_ϵ, P_ϵ)

        return new(
            ϵ_process, joint_process, 
            S_z, S_ϵ,
            β, γ, α, δ, l, μ)
    end
end


"""
Construct the production functions
"""
function makeproduction(model::DenHaanModel)
    @unpack δ, α, l, ζ, S_ϵ = model

    ϵ = collect(S_ϵ)
    
    H(z) = (ϵ' * π_ϵ(z, model))[1]

    R(z, K) = 1 + α * z * (K / (H(z) * l))^(α - 1) - δ
    w(z, K) = (1 - α) * z * (K / (H(z) * l))^α

    return R, w
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

"""
Get the Markov matrix for ϵ conditional on z
"""
function π_ϵ(z::State, model::DenHaanModel)
    @unpack S, P = model.ζ
    cond = findall(st -> st[1] == z, S)

    if isempty(cond) throw("State $z not in space $S") end

    P_ϵ = sum(P[cond, cond], dims=1)

    return rownormal(P_ϵ)'
end