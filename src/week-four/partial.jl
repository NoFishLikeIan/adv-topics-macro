include("aiyagari.jl")

include("algos/policy/solve.jl")
include("algos/policy/pfi.jl")

include("algos/stationary/eigenmethod.jl")
include("algos/stationary/montecarlo.jl")

include("../week-two/markov/discretization/tauchen.jl")

plot_dens = true

using Plots

function solvepartial(
    model::Aiyagari, r::Float64, w::Float64;
    mc=true, grid_N=1_000, end_grid=false,
    verbose=false, cache=false)

    a′, a_grid = policysolve(model, r, w; n_steps=grid_N, upperbound=10., verbose=verbose)

    if mc
        λ = distribution_mc(a′, a_grid, model; 
            verbose=verbose, inits=10_000, tol=1e-1, max_iter=5000)

    else 

        stable_Q = distribution_eigenvector(a′, a_grid, model, verbose=verbose, gth=true)
    
        center_y = get_row(model.y.S, 0.0)
        λ = collect(stable_Q[:, center_y])
        λ = λ ./ sum(λ)

    end

    
    return λ, a′, a_grid 
end

if plot_dens
    model = Aiyagari(
        0.9, 0.1, 7, # AR process parameters ρ, σ, N
        0, 0.95, 0.33, 0.1, 2 # Model parameters a, β, α, δ, σ_u
    )
    
    λ, a′, a_grid = solvepartial(
        model, .05, 1., grid_N=50;
        verbose=true, mc=false, end_grid=true)

    plot(title="Policy", xaxis="a")

    for y in model.y.transformation.(model.y.S)
        plot!(
            a_grid,
            (a -> a′(a, y)).(a_grid),
            label="a′(a | y = $(@sprintf("%.2f", y)))",
            xaxis="a"
        )
    end

    savefig("src/week-four/solutions/plots/policy.png")


    # plot(a_grid, λ, title="Stable distribution", xaxis="a", color=:red)

    # savefig("src/week-four/solutions/plots/stat_dist.png")
end