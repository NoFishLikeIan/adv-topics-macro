using Plots, ColorSchemes, StatsBase

include("markov/discretization/tauchen.jl")
include("markov/discretization/rouwenhorst.jl")
include("markov/simulation.jl")
include("markov/markov-types.jl")
include("../comm_utility.jl")

do_plot = true
plot_path = "src/week-two/solutions/plots/"
MoM_attempts = 3


m = 3
σ = 1

N_sp = [5, 100]
ρ_sp = [0.7, 0.99]
do_tauchen_sp =  [false, true]

param_space = Base.product(N_sp, ρ_sp, do_tauchen_sp)

for (N, ρ, do_tauchen) in param_space
    method_class = do_tauchen ? "tauchen" : "rouwenhorst"

    σ_ϵ = (1 - ρ^2)^2

    ar = Process(
        Normal(0, σ^2 / (1 - ρ^2)),
        z -> ρ * z,
        Normal(0, σ_ϵ)
    )

    if do_tauchen
        P, S = tauchen(ar, N, m)
    else
        P, S = try_n(() -> rouwenhorst(ar, N),
            MoM_attempts, StatsBase.ConvergenceException;
            verbose=true)
    end

    markov = MarkovDiscrete(P, S)

    stats = summary_stats(markov)
    print("Summary for ρ = $ρ and N = $N, with method = $method_class:\n")
    print("μ: ", stats["μ"], "\n")
    print("ν: ", stats["ν"], "\n")

    plot(stats["ρ"], linewidth=2, label="ACF", title="Autocorrelation with ($N, $ρ) for $method_class")
    savefig("$plot_path/$method_class/autocorr_N$N-rho$ρ.png")

    if do_plot
        # -- Plot a number of run simulations
        runs = 4
        T = 400

        plot(title="Markov simulation with ($N, $ρ) for $method_class")

        for i in 1:runs

            local z = @time discrete_sim(markov; T=T, drop=0)

            c = get(ColorSchemes.rainbow, i ./ runs)
            plot!(z, linewidth=2, label="It-$i")

            method_class = do_tauchen ? "tauchen" : "rouwenhorst"

            savefig("$plot_path/$method_class/ar-simulation_N$N-rho$ρ.png")

        end
    end
end