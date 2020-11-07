using Distributions, Statistics

include("markov/process.jl")
include("markov/simulation.jl")
include("markov/discretization/rouwenhorst.jl")

N = 50
ρ = 0.7
μ = 0
σ = 1

autocov = exp(0.5 * ((1 + ρ)^2 * σ_y^2 + σ^2)) - exp(σ_y^2)

z = Process(
    LogNormal(μ, σ / sqrt(1 - ρ^2)),
    z -> z^ρ,
    Normal(μ, σ),
    autocov)

P, S = rouwenhorst(z, N)
markov = MarkovDiscrete(P, S)

stats = summary_stats(markov)
print("Summary for ρ = $ρ and N = $N\n")
print("μ: ", stats["μ"], "\n")
print("ν: ", stats["ν"], "\n")
