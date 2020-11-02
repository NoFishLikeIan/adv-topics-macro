include("utils/types.jl")
include("utils/array.jl")
include("utils/interp.jl")

include("iteralgos.jl")
include("model/detgrowth.jl")

using Plots

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

    iter_fn = monotone ?  monotone_iteration(model, u_matrix, policy_vec, concave, howard) : naive_iteration(model, u_matrix, policy_vec, concave, howard)

    print("Starting iterations...\n")

    while true
        V_iter = iter_fn(V_i)

        if distance(V_iter, V_i) < model.ε  break end

        V_i = V_iter
    end

    evaluations = [k[trunc(Int, p)] for p in policy_vec]

    
    policy = constructlinear(k, evaluations)

    return V_i, policy

end

printpolicy(V, policy) = print(V[2:5], "-", policy(0.5), "-", policy(1), "\n\n")


print("Simple\n")
@time V, policy = solve_value_function(model)
printpolicy(V, policy)


print("Monotone\n")
@time V, policy = solve_value_function(model, monotone=true)
printpolicy(V, policy)


print("Concave\n")
@time V, policy = solve_value_function(model, concave=true)
printpolicy(V, policy)


# FIXME: Doesn't work
print("Howard\n")
@time V, policy = solve_value_function(model, howard=true)
printpolicy(V, policy)
