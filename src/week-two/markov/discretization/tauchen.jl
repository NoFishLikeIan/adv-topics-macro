using Distributions, Random

include("../process.jl")
include("../markov-types.jl")
include("utils.jl")

function imp_tauch(proc::Process, N::Int)
    Π = quantile(proc)
    F = cdf(proc)

    b = Π.(collect(1:N - 1) / N)

    # FIXME: Need to compute the z
    partition = Partition(b)

    transition(z) = (bj, bj_1) -> F((bj - evol(proc, z)) / σ_ϵ(proc)) - F((bj_1 - evol(proc, z)) / σ_ϵ(proc))

    P = zeros((N, N))

    for i in 1:N
        f = transition(partition[i])
        P[:, i] = [f(partition[i], partition[i - 1]) for i in 1:N]
    end

    return P, partition

end

"""
Discretize a Process based on equidistant points (mode = "equi") or equal probabilities (mode = "imp").
"""
function tauchen(proc::Process, N::Int; mode="equi", m::Int=3)
    F = cdf(proc)

    ψ = m * √(var(proc))
    
    partition = makepartition(ψ, N)
    d = distance(partition)


    f(z_row) = z_col -> F((z_row + d - evol(proc, z_col)) / σ_ϵ(proc))
    b(z_row) = z_col -> F((z_row - d - evol(proc, z_col)) / σ_ϵ(proc))

    P = zeros((N, N))

    P[:, 1] = f(partition[1]).(partition)

    for i in 2:(N - 1)
        P[:, i] = f(partition[i]).(partition) - b(partition[i]).(partition)
    end

    P[:, N] = 1 .- b(partition[N]).(partition)

    return P, partition
end