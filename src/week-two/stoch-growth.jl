using Distributions, Statistics, Plots, Roots

include("markov/process.jl")
include("markov/simulation.jl")
include("markov/discretization/rouwenhorst.jl")
include("../week-one/model/detgrowth.jl")

function make_z_proc(μ::Float64, σ::Float64, N::Int, ρ::Float64)::Tuple{MarkovDiscrete,Process}
    μ_y, σ_y = 0., σ / √(1 - ρ^2)

    z = Process(
        Normal(μ_y, σ_y),
        ρ,
        Normal(μ, σ)
    )

    P, S = rouwenhorst(z, N; numerical=false)
    markov = MarkovDiscrete(P, S, z -> exp(z))

    return markov, z
end

markov, z = make_z_proc(0., 1., 3, 0.7)
model = StochGrowthModel(.9, .5, markov, .01, 1_000)

u_c(x::Real) = 1 / x

function solve(m::StochGrowthModel, u_c::Function; max_iter=10_000)

    k = Partition(collect(range(1e-4, 5, length=m.k_size)))

    support_z = Partition(m.z.transformation.(m.z.S))

    zs = discrete_sim(m.z; T=max_iter + 500, drop=500)

    f(k, z) = z * k^m.α
    f_prime(k, z) = z * m.α * k^(m.α - 1)

    cond_exp(z) = P[get_row(support_z, z), :] .* support_z

    function euler_diff(c::Float64, k::Float64, z::Float64)
        exp_z = cond_exp(z)
        vals = (1. .+ f_prime.(k, exp_z)) * u_c(c)
        
        return u_c(c) - sum(vals)
    end

    policy_k = zeros(m.k_size) 

    for i in 1:max_iter
        z = zs[i]

        current_policy = zeros(m.k_size)

        for (j, k_row) in enumerate(k)
            opt_c = find_zero(c -> euler_diff(c, k_row, z), .5)
            opt_k = opt_c - f(k_row, z) 

            current_policy[j] = opt_k
        end

        if maximum(abs.(policy_k - current_policy)) > m.ε return current_policy end 

        policy_k = current_policy
    end
end

policy = solve(
    model, u_c
)