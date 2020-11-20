using SparseArrays, JLD

# A unit vector
struct e{T} <: AbstractVector{T}
    i::Int
    n::Int
end
e(i::Int, n::Int, ::Type{T}=Int) where T = e{T}(i, n)
function Base.getindex(e::e{T}, i) where T
    @boundscheck @assert i <= e.n
    ifelse(i == e.i, one(T), zero(T))
end
Base.size(e::e{T}) = (e.n,)

function matrix_distance(A::Array{Float64,2}, B::Array{Float64,2})::Float64
    N, M = size(A)

    d = 0.

    for col in 1:M
        m = maximum(abs.(A[:, col] - B[:, col]))
        d = m > d ? m : d
    end

    return d
end