using Distributions, Statistics, StatsBase

struct Process{T <: UnivariateDistribution}
    dist::T
    ρ::Float64
    error::UnivariateDistribution
end


Statistics.mean(p::Process) = Statistics.mean(p.dist)
Statistics.var(p::Process) = Statistics.var(p.dist)

σ_ϵ(p::Process) = √(var(p.error))

cdf(p::Process) = x::Real -> Distributions.cdf(p.dist, x)
quantile(p::Process) = q::Real -> Distributions.quantile(p.dist, q)

function evol(p::Process{LogNormal{Float64}}, z::Float64)
    return z^p.ρ
end

function evol(p::Process{Normal{Float64}}, z::Float64)
    return z*p.ρ
end

steplinear(p::Process) = z -> evol(p, z) + rand(p.error)
stepmult(p::Process) = z -> evol(p, z) * rand(p.error)

function autocov(p::Process{LogNormal{Float64}})
    μ, σ_proc = params(p.dist)
    σ_ϵ(p) 

    return exp(((1 + p.ρ)^2 * σ_proc^2 + + σ_ϵ(p)^2) / 2) - exp(σ_proc^2)
end

function autocov(p::Process{Normal{Float64}})
    μ, σ_proc = params(p.dist)

    return p.ρ * σ_proc^2
end

StatsBase.autocor(p::Process) = autocov(p) / var(p)

ispositive(p::Process) = !insupport(p.dist, -1)
