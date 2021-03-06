include("aiyagari.jl")

include("algos/policy/solve.jl")
include("algos/policy/pfi.jl")
include("algos/policy/endgrid.jl")

include("algos/stationary/eigenmethod.jl")
include("algos/stationary/montecarlo.jl")

include("../commons/matrix.jl")

include("../week-two/markov/discretization/tauchen.jl")

using Plots

"""
Construct the wealth distribution pdf based on the MonteCarlo 
or Eigenvector method
"""
function distribution_pdf(a′, a_grid, model; mc=false, kwargs...)

    if mc
        kde = distribution_mc(a′, a_grid, model; 
            verbose=verbose, inits=10_000)

        λ(x) = pdf(kde, x)

    else 

        stable_Q, finer_grid = distribution_eigenvector(a′, a_grid, model, verbose=verbose, gth=true)
    
        γ = QuantEcon.gth_solve(model.y.P)
        cent_dist = stable_Q * γ # get the equilibrium weighted density
        λ_grid = cent_dist / sum(cent_dist) # normalize

        λ = LinearInterpolation(finer_grid, λ_grid, extrapolation_bc=Line())

    end

    return λ

end

function solvepartial(
    model::Aiyagari, r::Float64, w::Float64;
    mc=true, verbose=false, kwargs...)

    a′, a_grid = policysolve(model, r, w; verbose=verbose, kwargs...)

    λ = distribution_pdf(a′, a_grid, model; mc=mc, verbose=verbose, kwargs...)

    
    return λ, a′, a_grid 
end