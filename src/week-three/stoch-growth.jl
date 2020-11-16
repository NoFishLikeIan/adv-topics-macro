include("solve-sgm/analytical.jl")
include("solve-sgm/log.jl")
include("solve-sgm/perturbation.jl")

using Dolo, Distances, Distributions, Random, Plots

Random.seed!(123)

model = yaml_import("src/week-three/models/sgm.yaml")
plot_path = "src/week-three/solutions/plots"

function construct_shock(model; T=200)
    
    ρ, σ = model.calibration.flat[:rho], model.calibration.flat[:sigma]
    
    ϵ = rand(Normal(0, σ), T)

    z = zeros(T)
    z[1] = 0.1

    for t in 2:T
        z[t] = z[t - 1] * ρ + ϵ[t]
    end

    return z
end


"""
Compares policy function and Euler Equation errors for the three methods,
    log-linearization, numerical quadratic approximation, analytical
"""
function compare_methods(model; bounds=[0.01, 1.], n_steps=500)
    quad_eee, tab = pert_eee(model; bounds=bounds, n_steps=n_steps)
    k_space = tab[:k]
    quad_policy = tab[:c]
    
    y_ss, k_ss, c_ss = analytical_ss(model)
    analy_policy = analytical_policy(model)[1].(k_space, y_ss)
    analy_eee = analytical_eee(model).(k_space, y_ss)

    plot(
        title="Euler equation errors", dpi=800,
        xaxis="k", yaxis="EEE(k)"
    )

    plot!(k_space, analy_eee, label="Analytical")
    plot!(k_space, quad_eee, label="Quadratic perturbation")

    savefig("$plot_path/eee_comparison.png")

    plot(
        title="Policy", dpi=800,
        xaxis="k", yaxis="c(k)"
    )
    plot!(k_space, analy_policy, label="Analytical")
    plot!(k_space, quad_policy, label="Quadratic perturbation")

    savefig("$plot_path/policy_comparison.png")

end

function perfect_foresigh_simulation(model; T=200)
    z = construct_shock(model; T=T)
    shocks = Dict(:z => z)

    sol_ref = perfect_foresight(model, shocks, T=T, complementarities=false)

end

if false

    plot(collect(1:T), z, title="AR(1) prod. simulation", label="z")
    savefig("$plot_path/ar1.png")

    shocks = Dict(:z => z)

    sol_ref = perfect_foresight(
        model, shocks, 
        T=T, complementarities=false)

    k = sol_ref[2, :][1:end - 1]

    plot(collect(0:T), k, title="Simulation", label="k(t)")

    savefig("$plot_path/k.png")

    x = loglin_poicy(model)

    x_policy = (k -> x(c_ss, k, y_ss)[1]).(k_space)


    eee = analytical_eee(model)
    plot(k_space, ((k) -> eee(k, y_ss)).(k_space))
    savefig("eee.png")
end


compare_methods(model)