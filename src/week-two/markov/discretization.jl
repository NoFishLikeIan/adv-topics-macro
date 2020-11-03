using Distributions, Random

include("process.jl")

"""
Discretize a Process based on equidistant points (mode = "equi") or equal probabilities (mode = "imp").
"""
function tauchen(proc::Process, N::Int, m::Int; mode="equi")
    F = cdf(proc)

    upper = m * sqrt(var(proc))
    lower = -upper
    
    partitions = range(lower, upper, length=N)
    d = (partitions[end] - partitions[end - 1]) / 2

    function scale_f(z::Real) 
        forward = z_p::Real -> F((z + d - proc.evol(z_p)) / var(proc.error))
        backward = z_p::Real -> F((z - d - proc.evol(z_p)) / var(proc.error))

        return forward, backward
    end
    

    P = zeros((N, N))

    z_grid = collect(Iterators.product(partitions, partitions))

    P[:, 1] = scale_f(partitions[1])[1].(partitions)

    for i in 2:(N - 1)
        f, b = scale_f(partitions[i])
        P[:, i] = f.(partitions) - b.(partitions)
    end

    P[:, N] = ones(N) - scale_f(partitions[N])[2].(partitions)

    return P
end

σ = 1
ρ = 0.7
N = 5
m = 3

# TODO: Note that this is not the distribution of the excercise, only to practice
ar = Process(
    Normal(0, σ^2 / (1 - ρ^2)),
    z -> ρ * z,
    Normal(0, 1)
)

P = tauchen(ar, N, m)