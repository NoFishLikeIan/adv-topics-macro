using Distributions, Statistics, Plots, Roots

include("ogm/stochgrowth.jl")
include("ogm/policy_fi.jl")
include("ogm/value_fi.jl")


plot_path = "src/week-two/solutions/plots/"

function run_solver(β::Float64, α::Float64, μ::Float64, σ::Float64, ρ::Float64, N::Int; do_plot=false, verbose=false)

    z = make_z_proc(μ, σ, N, ρ)

    f(k, prod) = prod * k^α
    f_prime(k, prod) = prod * α * k^(α - 1)
    u(c) = log(c)
    u_c(c) = 1 / c

    grid_N = 1_000

    model = StochGrowthModel(
        β, α, z,
        f, f_prime, u, u_c, u_c
    )

    policy = policy_solve(model; grid_N=grid_N, tol=1e-3, verbose=verbose)
    V = value_solve(model; grid_N=grid_N, tol=1e-3, verbose=verbose)

    if do_plot
        print("Plotting... \n")
        k_x = collect(range(1e-1, 5, length=grid_N))
        plot(k_x, support(z), policy', xaxis="k", yaxis="z", zaxis="k'", dpi=1000)

        savefig("$plot_path/policy-$β-$α-$ρ.png")

        plot(k_x, support(z), V', xaxis="k", yaxis="z", zaxis="V", dpi=1000)
        savefig("$plot_path/value-$β-$α-$ρ.png")

        plot(xaxis="k", yaxis="k'(k)", dpi=1000, legend=true, title="Policy function")

        for (h, z_h) in enumerate(support(z))
            rounded = round(z_h, digits=2)
            plot!(k_x, policy[:, h], label="z = $rounded")
        end
        savefig("$plot_path/cross-policy-$β-$α-$ρ.png")

        plot(xaxis="k", yaxis="V(k)", dpi=1000, legend=true, title="Value function")
        for (h, z_h) in enumerate(support(z))
            rounded = round(z_h, digits=2)
            plot!(k_x, V[:, h], label="z = $rounded")
        end
        savefig("$plot_path/cross-value-$β-$α-$ρ.png")
        print("...done!\n")
    end

    return V, policy
end

β, α = .9, .5
μ, σ, ρ = 0., 1., .7
N = 3

V, policy = run_solver(
    β, α, μ, σ, ρ, N; do_plot=true, verbose=true
)
