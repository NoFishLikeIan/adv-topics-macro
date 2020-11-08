include("../markov-types.jl")

function makepartition(ψ::Float64, N::Int)::Partition 
    Partition(range(-ψ, ψ, length=N))
end

function makepositivepartition(ψ::Float64, N::Int)::Partition 
    Partition(range(0, ψ, length=N))
end

colnormalized(A) = A ./ sum(A, dims=2)