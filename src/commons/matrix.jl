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
Base.size(e::e) = (e.n,)

function matrix_distance(A::Array{Float64,2}, B::Array{Float64,2})::Float64
    
    d = abs.(vec(A) - vec(B))

    return maximum(d)
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
