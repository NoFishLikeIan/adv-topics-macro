include("utils/types.jl")
include("utils/array.jl")
include("utils/interp.jl")

include("iteralgos.jl")
include("model/detgrowth.jl")
include("model/dynamics.jl")

using LinearAlgebra, Plots, Printf

model = DetGrowthModel(0.9, 0.5, 1., 0.01, 1_000)

function construct_utility(model::DetGrowthModel)
    function utility(k_cur::Real, k_prime::Real)::Real
        cap = model.z * k_cur^model.α
        return cap > k_prime ? log(cap - k_prime) : -Inf
    end

    function utility(k_tuple::Tuple{Real,Real})::Real
        utility(k_tuple[1], k_tuple[2])
    end
    
    return utility
end


function solve_value_function(
        model::DetGrowthModel;
        monotone=false, concave=false, howard=false
    )
    utility = construct_utility(model)

    k_star = equil_k(model)
    k = collect(range(0.01, k_star + 1, length=model.k_size))
    k_grid = collect(Iterators.product(k, k))

    V_0 = utility(k_star, k_star) / (1 - model.β) * ones(model.k_size)
    u_matrix = utility.(k_grid)

    V_i = copy(V_0)
    policy_vec = -1 * ones(model.k_size)

    iter_fn! = monotone ?  monotone_iteration!(model, u_matrix, policy_vec, concave) : naive_iteration!(model, u_matrix, policy_vec, concave)

    print("Starting iterations...\n")

    while true
        V_iter = iter_fn!(V_i)

        if distance(V_iter, V_i) < model.ε  break end

        if howard
            # As implemented by Heer and Maubner (2008)
            # TODO: Not as fast
            
            u = [u_matrix[i, asint(policy_vec[i])] for i in 1:model.k_size]
            
            Q = zeros((model.k_size, model.k_size))
            
            for (i, j) in enumerate(policy_vec)
                Q[i, asint(j)] = 1.
            end

            V_iter = inv(I - model.β * Q) * u

        end

        V_i = V_iter
    end

    print("...done!\n")

    evaluations = [k[asint(p)] for p in policy_vec]

    
    policy = constructlinear(k, evaluations)

    return V_i, policy

end


print("Simple\n")
@time V, policy = solve_value_function(model)

print("Monotone\n")
@time V, policy = solve_value_function(model, monotone=true)

print("Concave\n")
@time V, policy = solve_value_function(model, concave=true)

print("Howard\n")
@time V, policy = solve_value_function(model, howard=true)

