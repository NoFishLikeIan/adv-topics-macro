using Parameters, Interpolations
using Base.Threads


using Logging, Printf
"""
Compute a policy using the endogenous grid method
"""
function endgrid(
    a_grid::Vector{Float64}, ai::Aiyagari, R::Float64, w::Float64;
    n_steps=1_000, verbose=false, tol=1e-3, max_iter=1_000)

    u, u′, inv_u′ = make_u(ai)
    @unpack β, a_, y = ai

    Γ = y.P
    support_y = y.S
    cond_dens(x) = Γ[get_closest(support_y, x), :]
    
    ys = collect(support_y)
    shocks = y.transformation.(ys)
    
    T, N = length(ys), length(a_grid)

    policy = (a, y) -> a + y
    
    as_origin = zeros(N, T)
    prev_policy = zeros(N, T)

    for iter in 1:max_iter

        @threads for i in 1:N
            for j in 1:T
                a_p = a_grid[i]

                vals = @. u′(R * a_p + w * shocks - policy(a_p, ys))

                c = β * R * E(vals, cond_dens(ys[j]))

                as_origin[i, j] = (inv_u′(c) - w * ys[j] + a_p) / R
            end
        end

        eval_policy = zeros(N, T)

        for j in 1:T
            itps = LinearInterpolation(as_origin[:, j], a_grid, extrapolation_bc=Line())
            eval_policy[:, j] = @. max(itps(a_grid), a_)
        end

        err_distance = matrix_distance(eval_policy, prev_policy)

        if err_distance < tol 
            verbose && print("Found policy in $iter iterations (|x - x'| = $(@sprintf("%.4f", err_distance))\n")
            return (a, y) -> max(policy(a, y), a_)
        end

        policy = fromMtoFn(a_grid, ys, eval_policy)
        prev_policy = eval_policy

    end

    @warn "Could not find policy in $max_iter iterations with tolerance $tol"

    return policy
end 
