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

distance(P::Partition)::Real = (P.steps[end] - P.steps[end - 1]) / 2
