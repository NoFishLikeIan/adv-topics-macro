using Parameters, Kronecker

# exogenous transition matrices as specified in the paper
P_ϵ = [
    0.6 0.4; 
    2 / 45 43 / 45
] 

Π = [
    0.525 0.35 0.03125 0.09375;
    0.038889 0.836111 0.002083 0.122917;
    0.09375 0.03125 0.291667 0.583333; 
    0.009115 0.115885 0.024306 0.850694
]


cartesianstate(s1, s2) = vec(collect(Iterators.product(s1, s2)))
rownormal(M::Matrix{Float64}) = M ./ sum(M, dims=2)

P_z = rownormal([
    sum(Π[1:2, 1:2]) sum(Π[1:2, 3:4]);
    sum(Π[3:4, 1:2]) sum(Π[3:4, 3:4])
])



"""
Specification for the Aggregate Uncertainty model from
    Den Haan, Judd, and Juillard (2010)
"""
struct DenHaanModel
    Ε::StateMarkov
    Z::StateMarkov
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
        ϵ_process = StateMarkov(collect(S_ϵ), P_ϵ)
        z_process = StateMarkov(collect(S_z), P_z)

        return new(
            ϵ_process, z_process, joint_process, 
            S_z, S_ϵ,
            β, γ, α, δ, l, μ)
    end
end


"""
Construct the production functions
"""
function makeproduction(model::DenHaanModel)
    @unpack δ, α, l, μ = model
    @unpack ζ, S_ϵ = model

    ϵ = collect(S_ϵ)
    
    H(z) = (ϵ' * π_ϵ(z, model))[1]

    R(z, K) = 1 + α * z * (K / (H(z) * l))^(α - 1) - δ
    w(z, K) = (1 - α) * z * (K / (H(z) * l))^α

    τ(z) = (μ / l) * (-1 + 1 / H(z))

    return R, w, τ
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

    P_ϵ = sum(P[:, cond], dims=1)

    return rownormal(P_ϵ)'
end


function makeconsumption(model::DenHaanModel)

    @unpack μ, l = model
    R, w, τ = makeproduction(model)

    tax(z, ϵ) = (1 - τ(z)) * l * ϵ + μ * (1 - ϵ)

    """
    Default consumption function based on parameters
    """
    function consumption(a::Float64, K::Float64, z::Float64, ϵ::Float64, g::Function)
        return R(z, K) * a + tax(z, ϵ) * w(z, K) - g(a, K, ϵ, z)
    end
    function consumption(a′::Float64, a::Float64, m::Float64, z::Float64, ϵ::Float64)
        return consumption(a, m, z, ϵ, (_...) -> a′)
    end
    function consumption(x::Vector{Float64})
        return consumption(x...)
    end

    function invconsumption(x, K, z, ϵ, a′)
        return (x - tax(z, ϵ) * w(z, K) + a′) / R(z, K)
    end

    return consumption, invconsumption
end  