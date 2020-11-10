include("../markov/simulation.jl")
include("../markov/process.jl")
include("../markov/discretization/rouwenhorst.jl")

struct StochGrowthModel
    β::Float64
    α::Float64
    z::MarkovDiscrete
    f::Function
    f_prime::Function
    u::Function
    u_c::Function
    inv_u_c::Function
end

function make_z_proc(μ::Float64, σ::Float64, N::Int, ρ::Float64)::MarkovDiscrete
    μ_y, σ_y = 0., σ / √(1 - ρ^2)

    z = Process(
        Normal(μ_y, σ_y),
        ρ,
        Normal(μ, σ)
    )

    P, S = rouwenhorst(z, N; numerical=false)
    markov = MarkovDiscrete(P, S, exp)

    return markov
end