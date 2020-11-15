include("solve-sgm/analytical.jl")

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
end

c_policy = analytical_policy(model)

y_ss, k_ss, c_ss = analytical_ss(model)

k_space = range(0.01, 1., length=500)
policies = c_policy.(k_space, y_ss)

plot(k_space, policies)
savefig("test.png")