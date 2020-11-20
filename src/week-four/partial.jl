include("aiyagari.jl")
include("algos/endgrid.jl")
include("algos/piecewise.jl")
include("algos/vfi.jl")

include("algos/eigenmethod.jl")

include("../week-two/markov/discretization/tauchen.jl")

using Plots

default(size=(600, 600), dpi=800)

grid_N = 1_000
markov_N = 15

r = 0.05
w = 1.

gth = true
end_grid = true

model = Aiyagari(
    0.9, 0.1, markov_N, # AR process parameters ρ, σ, N
    0, 0.95, 0.33, 0.1, 2 # Model parameters a, β, α, δ, σ_u
)

a′, a_grid = (end_grid ? endgrid : value_solve)(model, r, w; n_steps=grid_N, upperbound=10.)


Φ = distribution_eigenvector(a′, a_grid, model, verbose=true, gth=gth)

plot(
    a_grid,
    (a -> a′(a, 0.)).(a_grid),
    label="a′(a)",
    title="Policy at y = 0",
    xaxis="a"
)

savefig("src/week-four/solutions/plots/policy.png")

surface(
    model.y.S, a_grid, 
    Φ, 
    xaxis="a", yaxis="y",
    title="Stationary distribution", 
    linealpha=0.3)

savefig("src/week-four/solutions/plots/dist$(gth ? "-gth" : "").png")

heatmap(
    Φ,
    xaxis="a", yaxis="y",
    title="Stationary distribution")

savefig("src/week-four/solutions/plots/heat$(gth ? "-gth" : "").png")
