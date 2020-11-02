include("types.jl")

function naive_max(array::RealArray)
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

function concave_max(array::RealArray)
    low::Int = 1
    high::Int = length(array) 

    while high - low > 2

        mid = floor(Int, low + (high - low) / 2)

        while array[mid] == -Inf
            high = mid
            mid = floor(Int, low + (high - low) / 2)
        end
        
        if array[mid + 1] < array[mid]
            high = mid
        end

        if array[mid - 1] < array[mid]
            low = mid
        end
    end

    arg_max = array[low] > array[high] ? low : high

    return arg_max, array[arg_max]
end

"""
Computes the sup norm between two arrays
"""
function distance(vec1::RealArray, vec2::RealArray)::Real
    maximum(abs.(vec1 - vec2))
end
