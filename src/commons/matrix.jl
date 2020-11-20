using SparseArrays, JLD2, FileIO


struct e{T} <: AbstractVector{T}
    i::Int
    n::Int
end
"""
Constructs a unit vector
"""
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


"""
Caches a matrix and returns a loading JLD path
"""
function cachematrix(
    A::Union{SparseMatrixCSC,Matrix}, 
    path::String, filename::String)::Tuple{String,String}

    destpath = joinpath(pwd(), path)
    isdir(destpath) || mkdir(destpath) # creates directory if it does not exist

    filepath = joinpath(destpath, filename, "$filename.jld2")

    save(filepath, filename, A)

    return filepath, filename
end
