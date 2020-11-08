using Distributions, Random, StatsBase

ρ = 0.7
σ = 1
μ = 0

σ_y = σ / sqrt(1 - ρ^2)

N = 10_000_000

y = zeros((N))

y[1] = 20

for t in 2:N
    y[t] = ρ * y[t - 1] + rand(Normal(μ, σ))
end

z = exp.(y[5000:N])

sample_mean = mean(z)
sample_var = var(z)
sample_cov = StatsBase.autocov(z, [1])[1]

theoretical = LogNormal(μ, σ_y)
mean_th = mean(theoretical)
var_th = var(theoretical)
cov_th = exp(0.5 * ((1 + ρ)^2 * σ_y^2 + σ^2)) - exp(σ_y^2)

print("

μ: $mean_th - $sample_mean
σ^2: $var_th - $sample_var
cov: $cov_th - $sample_cov
")
