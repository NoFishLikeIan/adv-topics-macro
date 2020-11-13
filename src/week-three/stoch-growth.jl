using Dolo, Distances, Distributions, Random, Plots

Random.seed!(123)

model = yaml_import("src/week-three/models/sgm.yaml")
plot_path = "src/week-three/solutions/plots"

T = 200

parameters = model.calibration.flat

ρ, σ = parameters[:rho], parameters[:sigma]
ϵ = rand(Normal(0, σ), T)

z = zeros(T)
z[1] = 0.1

for t in 2:T
    z[t] = z[t - 1] * ρ + ϵ[t]
end

plot(collect(1:T), z, title="AR(1) prod. simulation", label="z")
savefig("$plot_path/ar1.png")

shocks = Dict(:z => z)

sol_ref = perfect_foresight(
    model, shocks, 
    T=T, complementarities=false)

k = sol_ref[2, :][1:end - 1]

plot(collect(0:T), k, title="Simulation", label="k(t)")

savefig("$plot_path/k.png")
