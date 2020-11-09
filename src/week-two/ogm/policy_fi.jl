include("stochgrowth.jl")
include("../markov/simulation.jl")

include("../../comm_utility.jl")

function policy_solve(
    model::StochGrowthModel; 
    max_iter=10_000, tol=1e-3, grid_N=1_000, verbose=false)

    f, f_prime, u_c = model.f, model.f_prime, model.u_c 
    
    support_z = support(model.z)
    zs = discrete_sim(model.z; T=max_iter + 500, drop=500)
    
    k_space = Partition(collect(range(0.1, 5, length=grid_N)))

    function euler_diff(c::Float64, k::Float64, prod::Float64)::Float64

        vals = 1. .+ f_prime.(k, support_z)
        p_z = p_cond(model.z, prod)
        
        return u_c(c) -  model.β * E(vals, p_z) * u_c(f(k, prod) - k)
    end

    policy_k_matrix = ones(grid_N, length(support_z))
    
    for (h, z) in enumerate(support_z)

        policy_k = ones(grid_N)

        for i in 1:max_iter
            #z = zs[i]

            current_policy = 0.5 * ones(grid_N) # fixes the k = 1, z = 1 bug

            for (j, k_row) in enumerate(k_space)
                k_prime = current_policy[j]
                        
                opt_c = find_zero(c -> euler_diff(c, k_prime, z), .01, verbose=verbose)

                opt_k = f(k_row, z) - opt_c 
                current_policy[j] = opt_k

                
            end

            distance = maximum(abs.(policy_k - current_policy))

            if distance < tol 
                if verbose print("Found policy in $i iterations (|x - x'| = $distance)") end
                return current_policy 
            end 

            policy_k = current_policy
        end

        policy_k_matrix[:, h] = policy_k

    end

    throw(ConvergenceException(max_iter))
    
end