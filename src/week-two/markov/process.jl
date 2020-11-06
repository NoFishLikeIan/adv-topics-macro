using Distributions, Statistics

mutable struct Process
    dist::UnivariateDistribution
    evol::Function
    error::UnivariateDistribution
end


Statistics.mean(p::Process) = Statistics.mean(p.dist)
Statistics.var(p::Process) = Statistics.var(p.dist)

σ_ϵ(p::Process) = sqrt(var(p.error))

cdf(p::Process) = x::Real -> Distributions.cdf(p.dist, x)
quantile(p::Process) = q::Real -> Distributions.quantile(p.dist, q)

step(p::Process) = z -> p.evol(z) + rand(p.error)