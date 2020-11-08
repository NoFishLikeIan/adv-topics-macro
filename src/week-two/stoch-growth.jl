using Distributions, Statistics

include("markov/process.jl")
include("markov/simulation.jl")
include("markov/discretization/rouwenhorst.jl")

N = 3
ρ = 0.7
μ = 0
σ = 1

z = Process(
    LogNormal(μ, σ / √(1 - ρ^2)),
    z -> z^ρ,
    Normal(μ, σ))

P, S = rouwenhorst(z, N)
markov = MarkovDiscrete(P, S)

stats = summary_stats(markov)
print("Summary for ρ = $ρ and N = $N\n")
print("μ: ", stats["μ"], "\n")
print("ν: ", stats["ν"], "\n")
