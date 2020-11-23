using Interpolations, Roots
using Base.Threads

using Logging, Printf

function iterate_pfi(
    as::Vector{Float64}, ai::Aiyagari, R::Float64, w::Float64;
    max_iter=1_000, tol=1e-3, verbose=false
)

    u, u′, inv_u′ = make_u(ai)
    @unpack β, a_, y = ai

    support_y = y.S
    ys = collect(support_y)
    shocks = y.transformation.(ys)
    cond_dens(x) = Γ[get_closest(support_y, x), :]
    Γ = y.P

    N, T = length(as), length(ys)

    a_grid = copy(as)
    domain = cartesian(a_grid, ys)

    a_prime = ones(N, T)

    a′ = fromMtoFn(a_grid, ys, a_prime)
    
    for iter in 1:max_iter

        next_a_prime = copy(a_prime)

        @threads for i in 1:N
            t_a = a_grid[i]
            for j in 1:T
                t_y = ys[j]

                vals(a_p) = @. u′(R * a_p + w * shocks - a′(a_p, shocks))
                rhs(a_p) = β * R * E(vals(a_p), cond_dens(t_y))
                lhs(a_p) = u′(R * t_a + w * t_y - a_p)

                a0 = a′(t_a, t_y)
                try
                    solution = find_zero(z -> rhs(z) - lhs(z), a0)
                    next_a_prime[i, j] = max(solution, a_)
                catch e end
        
            end
        end
        
        err_distance = matrix_distance(a_prime, next_a_prime)

        verbose && print("Iteration $iter / $max_iter: $(@sprintf("%.4f", err_distance)) \r")

        if err_distance < tol 
            verbose && print("Found policy in $iter iterations (|x - x'| = $(@sprintf("%.4f", err_distance))\n")
            return a′
        end

        a_prime = next_a_prime
        a′ = fromMtoFn(a_grid, ys, a_prime)
    end

    @warn "Could not find policy in $max_iter iterations with tolerance $tol"

    return policy
end
