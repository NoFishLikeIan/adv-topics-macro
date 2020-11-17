include("solve-sgm/analytical.jl")
include("solve-sgm/perturbation.jl")
include("solve-sgm/log.jl")

include("solve-sgm/eee.jl")
include("solve-sgm/simulation.jl")

using Dolo, Distances, Distributions, Random, Plots

Random.seed!(11148705)

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
function compare_methods(model; bounds=[0.01, 0.7], n_steps=500)
    quad_eee, tab = pert_eee(model; bounds=bounds, n_steps=n_steps)
    k_space = tab[:k]
    quad_policy = tab[:c]
    
    y_ss, k_ss, c_ss = analytical_ss(model)
    analy_policy = analytical_policy(model)[1].(k_space, y_ss)
    analy_eee = analytical_eee(model).(k_space, y_ss)

    c_imp_pol, k_imp_pol = implicit_policy(model)
    c_imp_eee = eee(c_imp_pol, k_imp_pol, model).(k_space, y_ss)

    plot(title="Euler equation errors", dpi=800, xaxis="k", yaxis="EEE(k)")

    plot!(k_space, analy_eee, label="Analytical")
    plot!(k_space, quad_eee, label="Quad. pert.")
    plot!(k_space, c_imp_eee, label="Imp. Fn. Th.")


    savefig("$plot_path/sgm_comp/eee_comparison.png")

    plot(title="Policy", dpi=800, xaxis="k", yaxis="c(k)")
    plot!(k_space, analy_policy, label="Analytical")
    plot!(k_space, quad_policy, label="Quad. pert.")
    plot!(k_space, c_imp_pol.(k_space, y_ss), label="Imp. Fn. Th.")

    savefig("$plot_path/sgm_comp/policy_comparison.png")

end

function perfect_foresight_simulation(model; T=200)
    ts = collect(1:T)
    shocks = Dict(:z => construct_shock(model; T=T))

    # Dolo.jl simulations
    sol_ref = perfect_foresight(model, shocks, T=T, complementarities=false, verbose=false)
    z = sol_ref[1, 1:end - 1] # take the row for k excluding the last entry (undetermined)
    k = sol_ref[2, 1:end - 1]
    c = sol_ref[3, 1:end - 1]

    k1 = k[1]

    # Analytical simulation
    c_p_analytical, k_p_analytical = analytical_policy(model)
    tab_analytical = simulate_shock(c_p_analytical, k_p_analytical, k1, z)
    k_analytical = tab_analytical[:, 2]

    # Implicit policy simulation
    c_p_imp, k_p_imp = loglin_policy(model)
    tab_imp = simulate_shock(c_p_imp, k_p_imp, k1, z)
    k_imp = tab_imp[:, 2]

    plot(title="Simulation", xaxis="T", dpi=800)
    plot!(ts, k[2:end], yaxis="k", label="Quad. pert.", color=:blue, linewidth=1)
    plot!(ts, k_analytical[2:end], label="Analytical", color=:red, linewidth=1)


    plot!(ts, NaN .* ts, label="z", color=:black, alpha=0.5) # :(
    plt = twinx()
    plot!(plt, ts, z[2:end], label="z", color=:black, alpha=0.5, legend=false)

    savefig("$plot_path/simulation/simulation.png")

end


compare_methods(model)
perfect_foresight_simulation(model)