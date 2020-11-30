include("src/commons/interpolation.jl")
include("src/commons/matrix.jl")

include("src/week-five/markov/state.jl")
include("src/week-five/denhaan.jl")

include("src/week-five/policy/endgrid.jl")
include("src/week-five/policy/vfi.jl")
include("src/week-five/policy/stationary.jl")

include("src/week-five/krusell-smith/simulate.jl")
include("src/week-five/krusell-smith/det_simulation.jl")
include("src/week-five/krusell-smith/main.jl")

Random.seed!(11148706)

using Plots, KernelDensity

# default(dpi=600)

plot_path = "src/week-five/solutions/plots/"
model = DenHaanModel()

# Test algorithm

function run_krusselsmith(model_kwargs...)
    model = DenHaanModel(model_kwargs...)
    g = krusellsmith(
        model;
        verbose=true,
        ρ=0.5
    )

    a_grid = range(0.01, 5., length=500)

    plot(title="Policy function", legend=:left, xlabel="a", ylabel="a′(a)")

    for (z_f, ϵ_f) in model.ζ.S
        ys_pol = g.(as, 1., z_f, collect(Float64, ϵ_f))
        plot!(as, ys_pol, label="a′(a | $z_f, $ϵ_f)")
    end

    savefig("$plot_path/policy_krusell.png")

    as, zs = economysim(g, model; drop=0, T=5_000)

    boom_as = vec(mean(as[boom, :], dims=1))
    recession_as = vec(mean(as[.!boom, :], dims=1))

    a_grid = range(0.01, 5., length=500)
    plot(title="Steady state distribution", xaxis="a", yaxis="probability")

    plot!(a_grid, pdf(kde(boom_as), a_grid), label="boom")

    plot!(a_grid, pdf(kde(recession_as), a_grid), label="recession")

    savefig("$plot_path/stationary_distribution.png")

    return model, g, as, zs
end

