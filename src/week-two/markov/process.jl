using Distributions, Statistics, StatsBase

struct Process{T <: UnivariateDistribution}
    dist::T
    evol::Function
    error::UnivariateDistribution
end


Statistics.mean(p::Process) = Statistics.mean(p.dist)
Statistics.var(p::Process) = Statistics.var(p.dist)

σ_ϵ(p::Process) = sqrt(var(p.error))

cdf(p::Process) = x::Real -> Distributions.cdf(p.dist, x)
quantile(p::Process) = q::Real -> Distributions.quantile(p.dist, q)

steplinear(p::Process) = z -> p.evol(z) + rand(p.error)
stepmult(p::Process) = z -> p.evol(z) * rand(p.error)

function autocov(p::Process{LogNormal{Float64}})
    ρ = log(p.evol(3)) / log(3)
    μ, σ_proc = params(p.dist)

    return exp(((1 + ρ)^2 * σ_proc^2 + σ^2) / 2) - exp(σ_proc^2)
end

function autocov(p::Process{Normal{Float64}})
    ρ = p.evol(4) / 4
    μ, σ_proc = params(p.dist)

    return ρ * σ_proc^2
end

StatsBase.autocor(p::Process) = autocov(p) / var(p)