using Printf

using Interpolations, Parameters

include("../../../commons/matrix.jl")


function value_solve(
    model::Aiyagari, r::Float64, w::Float64; 
    upperbound=2.,
    max_iter=10_000, tol=1e-3, n_steps=1_000, verbose=false)

    @unpack β, y, a_, y = model

    ϕ = y.transformation
    Γ = y.P
    R = 1 + r

    u, u′, invu′ = make_u(model)

    function binding_u(v)
        a, a_p, y = v
        return u(R * a + w * ϕ(y) - a_p)
    end

    support_y = support(model.y)    
    a_space = Partition(collect(range(1e-5, upperbound, length=n_steps)))

    ay_grid = collect(Iterators.product(a_space, a_space, support_y))

    V_i = ones(n_steps, length(support_y))
    
    get_a = (a::Real) -> get_closest(a_space, a)
    
    utility = binding_u.(ay_grid)

    H = zeros(size(utility))
    

    for i in 1:max_iter

        verbose && print("Iteration $i / $max_iter \r")

        EV = V_i * Γ'

        for h in 1:length(support_y)
            H[:, :, h] = utility[:, :, h] .+ β * EV[:, h]
        end

        values, argmx = findmax(H, dims=2)

        next_V = reshape(values, size(V_i))

        distance_error = matrix_distance(V_i, next_V)


        if distance_error < tol 
            verbose && print("Found policy in $i iterations (|x - x'| = $(@sprintf("%.4f", distance_error))\n")

            current_policy = reshape(ay_grid[argmx], size(V_i))

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