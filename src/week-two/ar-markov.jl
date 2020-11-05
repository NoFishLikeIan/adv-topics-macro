using Plots, ColorSchemes

include("markov/discretization/tauchen.jl")
include("markov/simulation.jl")
include("markov/markov-types.jl")

plot_path = "src/week-two/solutions/plots/"
do_plot = false

# -- Process variables
σ = 1
ρ = 0.7

σ_ϵ = (1 - ρ^2)^2

ar = Process(
    Normal(0, σ^2 / (1 - ρ^2)),
    z -> ρ * z,
    Normal(0, σ_ϵ)
)

# -- Partition variables
N = 3_000
m = 3


P, S = tauchen(ar, N, m)
markov = MarkovDiscrete(P, S)

stats = summary_stats(markov)

print("μ: ", stats["μ"], "\n")
print("ν: ", stats["ν"], "\n")

plot(stats["ρ"], linewidth=2, label="Autocorrelation")
savefig("$plot_path/autocorr.png")

if do_plot
    # -- Plot a number of run simulations
    runs = 4
    T = 400

    plot()

    for i in 1:runs
        @time local z = discrete_sim(markov; T=T, drop=0)

        c = get(ColorSchemes.rainbow, i ./ runs)
        plot!(z, linewidth=2, label="It-$i")
    end

    savefig("$plot_path/ar-simulation_$N.png")
end
