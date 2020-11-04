using Plots, ColorSchemes

include("markov/discretization.jl")
include("markov/simulation.jl")
include("markov/markov-types.jl")

plot_path = "src/week-two/solutions/plots/"

# -- Process variables
σ = 1
ρ = 0.7

ar = Process(
    Normal(0, σ^2 / (1 - ρ^2)),
    z -> ρ * z,
    Normal(0, (1 - ρ^2)^2)
)

# -- Partition variables
N = 3_000
m = 3


P, partition = tauchen(ar, N, m)

# -- Plot a number of run simulations
runs = 4
T = 400

plot()

for i in 1:runs
    @time local z = discrete_sim(P, partition, T; drop=0)

    c = get(ColorSchemes.rainbow, i ./ runs)
    plot!(z, linewidth=2, label="It-$i")
end

savefig("$plot_path/ar-simulation.png")


