using Distributions, Statistics, Plots, Roots

include("ogm/stochgrowth.jl")
include("ogm/policy_fi.jl")


plot_path = "src/week-two/solutions/plots/"

function run_policy(β::Float64, α::Float64, μ::Float64, σ::Float64, ρ::Float64, N::Int; do_plot=false, verbose=false)

    z = make_z_proc(μ, σ, N, ρ)

    f(k, prod) = prod * k^α
    f_prime(k, prod) = prod * α * k^(α - 1)
    u_c(c) = 1 / c

    grid_N = 1_000

    model = StochGrowthModel(
        β, α, z,
        f, f_prime, u_c
    )

    policy = policy_solve(model; grid_N=grid_N, tol=1e-3, verbose=verbose)

    if do_plot
        plot(
            collect(range(1e-1, 5, length=grid_N)),
            policy, xaxis="k", yaxis="k'(k)", dpi=1000, legend=false
        )

        savefig("$plot_path/policy-$β-$α.png")
    end
end

β, α = .9, .5

run_policy(
    β, α, 0., 1., 0.7, 3; do_plot=true
)