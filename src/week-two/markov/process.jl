using Distributions, Statistics

mutable struct Process
    dist::UnivariateDistribution
    evol::Function
    error::UnivariateDistribution
end


Statistics.mean(p::Process) = Statistics.mean(p.dist)
Statistics.var(p::Process) = Statistics.var(p.dist)
cdf(p::Process) = x::Real -> Distributions.cdf(p.dist, x)
quantile(p::Process) = q::Real -> Distribution.quantile(p.dist, q)
step(p::Process) = z -> p.evol(z) + rand(p.error)