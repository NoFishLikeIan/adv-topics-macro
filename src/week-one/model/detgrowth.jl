include("../../week-two/markov/simulation.jl")

struct DetGrowthModel
    β::Float64
    α::Float64
    z::Float64
    ε::Float64
    k_size::Int
end

struct StochGrowthModel
    β::Float64
    α::Float64
    z::MarkovDiscrete
    ε::Float64
    k_size::Int
end