using Plots, ColorSchemes, StatsBase

include("markov/discretization/tauchen.jl")
include("markov/discretization/rouwenhorst.jl")
include("markov/simulation.jl")
include("markov/markov-types.jl")
include("../commons/routine.jl")

plot_path = "src/week-two/solutions/plots/"
MoM_attempts = 3


function run_markov(N::Int, ρ::Float64, do_tauchen::Bool; do_plot=false)
    method_class = do_tauchen ? "tauchen" : "rouwenhorst"

    std_err = √(1 - ρ^2)
    
    ar = Process(
        Normal(0, std_err),
        ρ,
        Normal(0, 1 - ρ^2)
    )

    if do_tauchen
        P, S = tauchen(ar, N; m=m)
    else
        P, S = try_n(() -> rouwenhorst(ar, N; numerical=false),
            MoM_attempts, StatsBase.ConvergenceException;
            verbose=true)
    end

    markov = MarkovDiscrete(P, S, z -> z)

    stats = summary_stats(markov)
    print("\nSummary for ρ = $ρ and N = $N, with method = $method_class:\n")
    print("μ: ", stats["μ"], "\n")
    print("ν: ", stats["ν"], "\n")
    print("ρ: ",  stats["ρ"][2], "\n")

    plot(stats["ρ"], linewidth=2, label="ACF", title="Autocorrelation with ($N, $ρ) for $method_class")
    savefig("$plot_path/$method_class/autocorr_N$N-rho$ρ.png")

    if do_plot
        # -- Plot a number of run simulations
        runs = 1
        T = 1500

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

m = 3
N_sp = [5, 100]
ρ_sp = [0.7, 0.99]
do_tauchen_sp =  [false, true]

for (N, ρ, do_tauchen) in Base.product(N_sp, ρ_sp, do_tauchen_sp)
    run_markov(N, ρ, do_tauchen; do_plot=false)
end
