include("solve-sgm/analytical.jl")
include("solve-sgm/log.jl")
include("solve-sgm/perturbation.jl")

using Dolo, Distances, Distributions, Random, Plots

Random.seed!(11148705)

# TODO: Euler equation should not depend on model but only on policy

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

    #  c_imp_pol = implicit_policy(model)[1].(k_space, y_ss) 
    #  c_imp_eee = implicit_eee(model).(k_space, y_ss)
#  
    #  log_pol = loglin_poicy(model)[1].(k_space, y_ss)
    #  log_eee = loglin_eee(model).(k_space, y_ss)


    plot(
        title="Euler equation errors", dpi=800,
        xaxis="k", yaxis="EEE(k)"
    )

    plot!(k_space, analy_eee, label="Analytical")
    plot!(k_space, quad_eee, label="Quadratic perturbation")
    # plot!(k_space, c_imp_eee, label="Implicit function theorem")
    # plot!(k_space, log_eee, label="Log-linearization")


    savefig("$plot_path/sgm_comp/eee_comparison.png")

    plot(
        title="Policy", dpi=800,
        xaxis="k", yaxis="c(k)"
    )
    plot!(k_space, analy_policy, label="Analytical")
    plot!(k_space, quad_policy, label="Quadratic perturbation")
    # plot!(k_space, c_imp_pol, label="Implicit function theorem")
    # plot!(k_space, log_pol, label="Log-linearization")

    savefig("$plot_path/sgm_comp/policy_comparison.png")

end

function perfect_foresigh_simulation(model; T=200)
    ts = collect(0:T)
    shocks = Dict(:z => construct_shock(model; T=T))

    sol_ref = perfect_foresight(model, shocks, T=T, complementarities=false, verbose=false)
    k = sol_ref[2, 1:end - 1] # take the row for k excluding the last entry (undetermined)
    z = sol_ref[1, 1:end - 1]


    plot(title="Simulation", xaxis="T", dpi=800)
    plot!(ts, k, yaxis="k", label="k", color=:blue)
    plot!(ts, NaN .* ts, label="z", color=:green, alpha=0.5) # :(

    plt = twinx()
    plot!(plt, ts, z, label="z", color=:green, alpha=0.5,  linewidth=2,
        legend=false)

    savefig("$plot_path/simulation/simulation.png")

end


compare_methods(model)
perfect_foresigh_simulation(model)