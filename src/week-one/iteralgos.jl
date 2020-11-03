include("model/detgrowth.jl")
include("utils/types.jl")
include("utils/array.jl")


asint(x::Real)::Int = trunc(Int, x)

function naive_iteration!(model::DetGrowthModel, u_matrix, policy_vec::RealArray, concave=false)

    maximizer = concave ? concave_max : naive_max
    
    function iterate(V_i::RealArray)::RealArray
        V_iter = copy(V_i)

        H = u_matrix .+ model.β * V_i'
        
        for (k_idx, h_row) in enumerate(eachrow(H))

            k_prime, v_max = maximizer(h_row)
            

            V_iter[k_idx] = v_max
            policy_vec[k_idx] = k_prime
        end

        return V_iter
    end

    return iterate
end

function monotone_iteration!(model::DetGrowthModel, u_matrix, policy_vec::RealArray, concave=false)

    maximizer = concave ? concave_max : naive_max
        
    function iterate(V_i::RealArray)::RealArray
        V_iter = copy(V_i)

        prev = 1

        H = u_matrix .+ model.β * V_i'
        
        for (k_idx, h_row) in enumerate(eachrow(H))

            k_prime, v_max = maximizer(h_row[prev:end])

            V_iter[k_idx] = v_max

            prev = asint(k_prime + prev - 1)

            policy_vec[k_idx] = prev
        end

        return V_iter
    end

    return iterate
end
