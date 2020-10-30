RealArray{T <: Real} = AbstractVector{T}

function intervalbounds(partition::RealArray, x::Real)::Tuple{Int,Int}

    for (index, bound) in enumerate(partition)

        if x <= bound
            return index - 1, index
        end
    end

    return length(partition), 0
end

# TODO: Force length of array at type level
function constructlinear(arguments::RealArray, evaluations::RealArray)
    
    interp = function (x::Real)  
        lower, upper = intervalbounds(evaluations, x)
        
        # TODO: This behaviour might not be the best
        if lower == 0 return -Inf end
        if upper == 0 return Inf end

        base = evaluations[lower]
        pull = evaluations[upper]

        gap = arguments[upper] - arguments[lower]
        step = x - arguments[lower]

        return base + step * (pull - base) / gap
        
    end

    return interp
end
