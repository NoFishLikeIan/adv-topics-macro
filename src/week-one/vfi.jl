include("utils/types.jl")
include("utils/array.jl")
include("utils/interp.jl")

include("detgrowth.jl")

using Plots

model = DetGrowthModel(0.9, 0.5, 0.01, 2_000)

function construct_utility(model::DetGrowthModel)
    function utility(k_cur::Real, k_prime::Real)::Real
        cap = k_cur^model.α
        return cap > k_prime ? log(cap - k_prime) : -Inf
    end

    function utility(k_tuple::Tuple{Real,Real})::Real
        utility(k_tuple[1], k_tuple[2])
    end
    
    return utility
end

function solve_value_function(model::DetGrowthModel)
    utility = construct_utility(model)

    k_star = equil_k(model)
    k = collect(range(0.01, k_star + 2, length=model.k_size))
    k_grid = collect(Iterators.product(k, k))

    V_0 = utility(k_star, k_star) / (1 - model.β) * ones(model.k_size)
    u_matrix = utility.(k_grid)

    V_i = copy(V_0)
    policy_vec = -1 * ones(model.k_size)

    while true
        V_iter = copy(V_i)
        for (k_idx, utility_row) in enumerate(eachrow(u_matrix))
            H = utility_row + model.β * V_i
            k_prime, v_max = find_maximum(H)

            V_iter[k_idx] = v_max
            policy_vec[k_idx] = k_prime
        end

        if distance(V_iter, V_i) < model.ε break end

        V_i = V_iter
    end

    evaluations = [k[trunc(Int, p)] for p in policy_vec]

    
    policy = constructlinear(k, evaluations)
    return V_i, policy

end


V, policy = solve_value_function(model)
