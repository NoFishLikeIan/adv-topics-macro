include("aiyagari.jl")

include("algos/policy/endgrid.jl")
include("algos/policy/vfi.jl")

include("algos/stationary/eigenmethod.jl")
include("algos/stationary/montecarlo.jl")

include("../week-two/markov/discretization/tauchen.jl")

plot_dens = false

using Plots

function solvepartial(
    model::Aiyagari, r::Float64, w::Float64;
    mc=true, grid_N=1_000, end_grid=false,
    verbose=false, cache=false)

    a′, a_grid = value_solve(model, r, w; n_steps=grid_N, upperbound=10., verbose=verbose)

    if mc
        Φ = distribution_mc(a′, a_grid, model; 
            verbose=verbose, inits=10_000, tol=1e-1, max_iter=5000)

    else 

        stable_Q = distribution_eigenvector(a′, a_grid, model, verbose=verbose, gth=true)
    
        center_y = get_row(model.y.S, 0.0)
        Φ = collect(stable_Q[:, center_y])
        Φ = Φ ./ sum(Φ)

    end

    
    return Φ, a′, a_grid 
end

if plot_dens

    model = Aiyagari(
        0.9, 0.1, 7, # AR process parameters ρ, σ, N
        0, 0.95, 0.33, 0.1, 2 # Model parameters a, β, α, δ, σ_u
    )


    # FIXME: Doesn't work, Φ is a dist
    Φ_mc, Φ_ei, a′, a_grid = solvepartial(model, .05, 1., verbose=true)


    plot(
        a_grid,
        (a -> a′(a, 0.)).(a_grid),
        label="a′(a)",
        title="Policy at y = 0",
        xaxis="a"
    )

    plot(a_grid, Φ_mc, title="Stable distribution", label="MC", xaxis="a", color=:red)
    plot!(a_grid, Φ_ei, label="Ei", color=:black, alpha=0.5)


    savefig("src/week-four/solutions/plots/policy.png")
end