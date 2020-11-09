include("stochgrowth.jl")
include("../markov/simulation.jl")

include("../../comm_utility.jl")

function matrix_distance(A, B)
    N, M = size(A)

    d = 0.

    for col in 1:M
        m = maximum(abs.(A[:, col] - B[:, col]))
        d = m > d ? m : d
    end

    return d
end

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

    for i in 1:max_iter

        if verbose print("Iteration $i / $max_iter \r") end

        EV = V_i * model.z.P'

        H = copy(utility)

        for h in 1:length(support_z)
            H[:, :, h] .+= model.β * EV[:, h]
        end

        next_V = reshape(maximum(H, dims=2), size(V_i))

        distance = matrix_distance(V_i, next_V)

        if distance < tol 
            if verbose print("Found policy in $i iterations (|x - x'| = $distance)") end
            break
        end 
        
        V_i = next_V

    end

    return V_i
    
end