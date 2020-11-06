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

    transition(z) = (bj, bj_1) -> F((bj - proc.evol(z)) / σ_ϵ(proc)) - F((bj_1 - proc.evol(z)) / σ_ϵ(proc))

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

    ψ = m * sqrt(var(proc))
    
    partition = makepartition(ψ, N)
    d = distance(partition)

    transition(z) = (z_p ->
    F((z + d - proc.evol(z_p)) / σ_ϵ(proc)) 
    - F((z - d - proc.evol(z_p)) / σ_ϵ(proc)))
    
    P = zeros((N, N))

    for i in 1:N
        f = transition(partition[i])
        P[:, i] = f.(partition)
    end

    return P, partition
end


N = 5
σ = 1
ρ = 0.7

proc = Process(
    Normal(0, σ^2 / (1 - ρ^2)),
    z -> ρ * z,
    Normal(0, (1 - ρ^2)^2)
)

P, S = imp_tauch(proc, 5)