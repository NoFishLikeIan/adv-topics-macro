using Parameters, Interpolations
using Printf

using Logging

"""
Compute a policy using the endogenous grid method
"""
function endgrid(
    as::Vector{Float64}, ai::Aiyagari, R::Float64, w::Float64;
    n_steps=1_000, verbose=false, tol=1e-3, max_iter=1_000)

    u, ∂u∂c, inv_∂u∂c = make_u(ai)
    @unpack β, a_, y = ai

    support_y = y.S
    ys = collect(support_y)
    ϕ = y.transformation
    Γ = y.P

    N, T = length(as), length(ys)

    ap_grid = copy(as)
    domain = cartesian(as, ys)

    function construct_policy(a_matrix::Matrix{Float64})

        function a(a′::Float64, y::Float64) 
            return fromMtoFn(as, ys, a_matrix)(a′, y)
        end
    
        function a(v::Vector{Float64})
            return a(v...)
        end
    
        function a(a′::Matrix{Float64}, y::Vector{Float64})
            next_a = similar(a′)
            for j in 1:length(y)
                col_a = a′[:, j]
                next_a[:, j] = a.(col_a, y[j])
            end
    
            return next_a
        end
    
        return a
    end    

    a_dprime = repeat(ap_grid, 1, T)

    for iter in 1:max_iter

        a′′ = construct_policy(a_dprime)

        function a(a′::Float64, y::Float64)
            u_values = @. ∂u∂c(R * a′ + w * ϕ(ys) - a′′(a′, ϕ(ys)))

            rhs = β * R *  Γ[get_closest(support_y, y), :]' * u_values
            return (inv_∂u∂c(rhs) - w * ϕ(y) + a′) / R
        end

        next_adp = zeros(size(a_dprime))

        for j in 1:T
            
            reverted_a = a.(ap_grid, ys[j])
            i = sortperm(reverted_a)

            a_pp_intp = fromVtoFn(reverted_a[i], ap_grid[i])

            next_adp[:, j] = max.(a_pp_intp.(ap_grid), a_)

        end

        err_distance = matrix_distance(a_dprime, next_adp)

        verbose && print("Iteration $iter / $max_iter -> $err_distance\r")

        if err_distance < tol 
            verbose && print("Found policy in $iter iterations (|x - x'| = $(@sprintf("%.4f", err_distance))\n")
            return a′′
        end 

        a_dprime = next_adp
    end

    @warn "Could not find policy in $max_iter iterations with tolerance $tol"

    return a′′

end 
