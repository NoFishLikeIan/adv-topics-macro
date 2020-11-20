using Parameters

using Interpolations

include("../../comm_utility.jl")


function value_solve(
    model::Aiyagari, r::Float64, w::Float64; 
    upperbound=2.,
    max_iter=10_000, tol=1e-3, grid_N=1_000, verbose=false)

    @unpack β, y, a_ = model

    Γ = model.y.P
    R = 1 + r

    u, u′, invu′ = make_u(model)

    binding_u = (v) -> u(R * v[1] + w * v[3] - v[2])

    support_y = support(model.y)    
    a_space = Partition(collect(range(1e-5, upperbound, length=grid_N)))

    ay_grid = collect(Iterators.product(a_space, a_space, support_y))

    V_i = ones(grid_N, length(support_y))
    
    get_a = (a::Real) -> get_closest(a_space, a)
    
    utility = binding_u.(ay_grid)

    H = zeros(size(utility))
    

    for i in 1:max_iter

        verbose && print("Iteration $i / $max_iter \r")

        EV = V_i * Γ'

        for h in 1:length(support_y)
            H[:, :, h] = utility[:, :, h] .+ β * EV[:, h]
        end

        values, argmax = findmax(H, dims=2)

        next_V = reshape(values, size(V_i))

        distance = matrix_distance(V_i, next_V)


        if distance < tol 
            verbose && print("Found policy in $i iterations (|x - x'| = $distance)\n")

            current_policy = reshape(ay_grid[argmax], size(V_i))

            ϕ = (v -> v[2]).(current_policy)

            intp = LinearInterpolation((a_space, support_y,), ϕ, extrapolation_bc=(Linear(), Linear()))
            
            function policy(a::Float64, y::Float64)::Float64
                max(intp(a, y), a_)
            end
            policy(v::Array{Float64}) = policy(v...)

            return policy, a_space.steps
        end 
        
        V_i = next_V

    end
    
end