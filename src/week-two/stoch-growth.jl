using Distributions, Statistics

include("markov/process.jl")
include("markov/simulation.jl")
include("markov/discretization/rouwenhorst.jl")

include("../week-one/model/detgrowth.jl")

function make_z_proc(μ::Float64, σ::Float64, N::Int, ρ::Float64)::MarkovDiscrete
    μ_y, σ_y = 0., σ / √(1 - ρ^2)

    z = Process(
        LogNormal(μ_y, σ_y),
        z -> z^ρ,
        Normal(μ, σ)
    )

    P, S = rouwenhorst(z, N; numerical=false)
    markov = MarkovDiscrete(P, S)

    return markov
end

markov = make_z_proc(0, 1, 3, 0.7)
model = StochGrowthModel(.9, .5, markov, .01, 1_000)

