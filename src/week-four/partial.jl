include("aiyagari.jl")
include("algos/endgrid.jl")
include("algos/piecewise.jl")
include("algos/vfi.jl")

include("algos/eigenmethod.jl")

include("../commons/matrix.jl")

include("../week-two/markov/discretization/tauchen.jl")

using Plots

function solvepartial(
    model::Aiyagari, r::Float64, w::Float64;
    grid_N=1_000, gth=true, verbose=false, cache=false, end_grid=true)

    a′, a_grid = (end_grid ? endgrid : value_solve)(
    model, r, w; n_steps=grid_N, upperbound=10., verbose=verbose)

    Φ = distribution_eigenvector(a′, a_grid, model, verbose=verbose, gth=gth)

    return Φ, a′, a_grid 
end


# Run only if directly called, not included
if abspath(PROGRAM_FILE) == @__FILE__

    model = Aiyagari(
    0.9, 0.1, 7, # AR process parameters ρ, σ, N
    0, 0.95, 0.33, 0.1, 2 # Model parameters a, β, α, δ, σ_u
    )

    Φ, a′, a_grid = solvepartial(model, 0.05, 1.)

    plot(
        a_grid,
        (a -> a′(a, 0.)).(a_grid),
        label="a′(a)",
        title="Policy at y = 0",
        xaxis="a"
    )

    savefig("$plotpath/policy.png")

    surface(
        model.y.S, a_grid, 
        Φ, 
        xaxis="a", yaxis="y",
        title="Stationary distribution", 
        linealpha=0.3)

    savefig("$plotpath/dist$(gth ? "-gth" : "").png")

    heatmap(
        Φ,
        xaxis="a", yaxis="y",
        title="Stationary distribution")

    savefig("$plotpath/heat$(gth ? "-gth" : "").png")
end