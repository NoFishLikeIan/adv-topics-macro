include("aiyagari.jl")
include("algos/endgrid.jl")
include("algos/piecewise.jl")
include("algos/vfi.jl")

include("algos/eigenmethod.jl")

include("../week-two/markov/discretization/tauchen.jl")

using Plots

r = .0
w = .5

model = Aiyagari(
    0.9, 0.1, 150, # AR process parameters ρ, σ, N
    0, 0.95, 0.33, 0.1, 2 # Model parameters a, β, α, δ, σ_u
)

a′, a_grid = value_solve(model, r, w; grid_N=150, upperbound=10)

# a′, a_grid = endgrid(model, r, w; n_steps=150, upperbound=50)

Φ = distribution_eigenvector(a′, a_grid, model)