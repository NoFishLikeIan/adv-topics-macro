include("stochgrowth.jl")
include("../markov/simulation.jl")

include("../../comm_utility.jl")

function value_solve(
    model::StochGrowthModel; 
    max_iter=10_000, tol=1e-3, grid_N=1_000, verbose=false)

    f, f_prime, u_c, inv_u_c, u = model.f, model.f_prime, model.u_c, model.inv_u_c, model.u


    function u_k(k, k_prime, z) 
        k_next = f(k, z)
        if k_next < k_prime return -Inf end

        return u(k_next - k_prime) 
    end
    function u_k(param::Tuple{Float64,Float64,Float64})
        k, k_p, z = param
        return u_k(k, k_p, z)
    end

    support_z = support(model.z)    
    k_space = Partition(collect(range(1e-5, 5, length=grid_N)))

    kz_grid = collect(Iterators.product(k_space, k_space, support_z))

    V_i = ones(grid_N, length(support_z))

    utility = u_k.(kz_grid)

    get_k = (k::Real) -> get_closest(k_space, k)

    function euler_diff(k::Float64, prod::Float64, k_prime::Float64)::Float64

        vals = 1. .+ f_prime.(k, support_z)
        p_z = p_cond(model.z, prod)

        return inv_u_c(model.β * E(vals, p_z) * u_c(f(k, prod) - k_prime))
    end
    
    H = zeros(size(utility))

    for i in 1:max_iter

        if verbose print("Iteration $i / $max_iter \r") end

        EV = V_i * model.z.P'

        for h in 1:length(support_z)
            H[:, :, h] = utility[:, :, h] .+ model.β * EV[:, h]
        end

        values, argmax = findmax(H, dims=2)

        next_V = reshape(values, size(V_i))

        distance = matrix_distance(V_i, next_V)

        if distance < tol 
            if verbose print("Found policy in $i iterations (|x - x'| = $distance)") end

            current_policy = reshape(kz_grid[argmax], size(V_i))

            euler_error = zeros(size(V_i))

            for (h, ζ) in enumerate(support_z)
                k_p = [tup[1] for tup in current_policy[:, h]]
                k_double_prime = k_p[get_k.(k_p)]
                
                euler_eq = euler_diff.(k_p, ζ, k_double_prime)
                current_c = f.(k_space, ζ) .- k_p 

                euler_error[:, h] =  log10.(abs.(1. .- (euler_eq ./ current_c)))
            end
            
            
            return V_i, euler_error
        end 
        
        V_i = next_V

    end
    
end