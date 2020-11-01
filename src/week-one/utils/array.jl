include("types.jl")

function unknown_fn(array::RealArray)
    max::Float64 = -Inf
    idx::Int = -1

    for (index, element) in enumerate(array)
        if element >= max
            max = element
            idx = index
        end
    end

    return idx, max
end

function concave_fn(array::RealArray)
    low::Int = 1
    high::Int = length(array) 

    while low != high

        # TODO: Better way to do this
        if high == low + 1
            high = low
        end

        mid = floor(Int, low + (high - low) / 2)
        
        if array[mid + 1] < array[mid]
            high = mid
        end

        if array[mid - 1] < array[mid]
            low = mid
        end
    end

    return low, array[low]
end

"""
Find the maximum value and index of an array. 
"""
function find_maximum(array::RealArray; mode="unknown")::Tuple{Int,Real}  
    fn = mode == "concave" ? concave_fn : unknown_fn
    
    return fn(array)
end

"""
Computes the sup norm between two arrays
"""
function distance(vec1::RealArray, vec2::RealArray)::Real
    maximum(abs.(vec1 - vec2))
end