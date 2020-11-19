include("aiyagari.jl")
include("algos/endgrid.jl")

include("../week-two/markov/discretization/tauchen.jl")

r = .04
w = 1.

model = Aiyagari(
    0.9, 0.1, 7, # AR process parameters ρ, σ, N
    0, 0.95, 0.33, 0.1, 2 # Model parameters a, β, α, δ, σ_u
)

a′, a_grid = endgrid(model, r, w; n_steps=10)
# Q_a = computeQ_a(a′, grid)
