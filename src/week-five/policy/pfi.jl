using Parameters, Interpolations
using Base.Threads
using Roots

using Logging, Printf

cartesianfromsize(I, J, L, M) = collect(Iterators.product(1:I, 1:J, 1:L, 1:M))

"""
Compute a policy using policy function iteration, 
assuming exogenous forecasting rule Ψ
"""
function pfi(
    prices::NTuple{2,Function}, utils::NTuple{3,Function}, 
    Ψ::Function, model::DenHaanModel, grids_sizes::NTuple{2,Int};
    bounds=[.01, 10.],
    max_iter=1_000, tol=1e-3, ρ=0.7, 
    verbose=false)

    @unpack S_ϵ, S_z, ζ, β = model
    u, u′, invu′ = utils
    R, w = prices

    ϵ_grid = collect(Float64, S_ϵ)
    z_grid = collect(S_z)
    D_z, D_ϵ = length(S_ϵ), length(S_z)
    P_cond(state) = ζ.P[findfirst(==(state), ζ.S), :]
    
    N_a, N_m = grids_sizes

    m_grid = range(bounds..., length=N_m)
    a_grid = range(bounds..., length=N_a)

    policy = ones(N_a, N_m, D_z, D_ϵ)
    
    for iter in 1:max_iter
        g = fromMtoFn(policy, a_grid, m_grid, ϵ_grid, z_grid)
        new_policy = similar(policy)

        for (i, j, k, l) in cartesianfromsize(N_a, N_m, D_z, D_ϵ)
            # Solve Euler for ā

            a, m = a_grid[i], m_grid[j]
            z, ϵ = z_grid[k], ϵ_grid[l]

            Ψ′ = Ψ(z, m)

            function righteuler(a′)
                function inner(state::Tuple{State,State})
                    z′, ϵ′ = state
                    ϵ′ = convert(Float64, ϵ′)

                    rate = R(z′, Ψ′)
                    wage = w(z′, Ψ′)

                    return rate * u′(rate * a′ + wage * ϵ′) - g(a′, Ψ′, ϵ′, z′)
                end

                return inner
            end

            lefteuler(a′) = u′(R(z, m) * a + w(z, m) * ϵ - a′)
            euler(a′) = lefteuler(a′) - β * P_cond((z, ϵ))' * righteuler(a′).(ζ.S)

            f1, f2 = euler.(bounds) # Evaluate the function at the boundaries
            
            try
                a′ = f1 * f2 < 0 ? 
                    find_zero(euler, bounds) : 
                    find_zero(euler, bounds[1])

                new_policy[i, j, k, l] = a′

            catch e
                print("$a, $m, $z, $ϵ\n")
                throw(e)
            end
        end

        d = new_policy - policy
        err_distance = maximum(abs.(d))
        verbose && print("Iteration $iter / $max_iter: $(@sprintf("%.4f", err_distance)) \r")

        if err_distance < tol 
            verbose && print("Found policy in $iter iterations (|x - x'| = $(@sprintf("%.4f", err_distance))\n")
            return g
        end

        policy += ρ * d # Update with dumping parameter
    end

    throw(ConvergenceException(max_iter))
end
