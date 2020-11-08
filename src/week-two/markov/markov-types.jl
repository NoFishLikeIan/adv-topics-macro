using Base

struct Partition{T,N,AA <: AbstractArray{T,N}} <: AbstractArray{T,N} 
    steps::AA
end

Partition(steps::AbstractArray{T,N}) where {T,N} = Partition{T,N,typeof(steps)}(steps)

Base.length(P::Partition) = Base.length(P.steps)
Base.size(P::Partition) = Base.size(P.steps)
Base.size(P::Partition, d) = Base.size(P.steps, d)


function Base.getindex(P::Partition, i::Int)
    if i == 0
        return -Inf
    elseif i == length(P) + 1
        return Inf
    else 
        return P.steps[i]
    end

end

currygetindex(P::Partition) = (z -> Base.getindex(P, z))

distance(P::Partition)::Real = (P.steps[end] - P.steps[end - 1]) / 2

"""
Computes the row of the exact partition. If partition not found, returns 0
"""
function get_row(P::Partition, z::Real)
    for (i, row_z) in enumerate(P.steps)
        if isapprox(row_z, z) return i end
    end

    return 0
end

function get_closest(P::Partition, z::Real)
    
    for (i, boundary) in enumerate(P.steps)
        if z < boundary
            return i
        end
    end
end