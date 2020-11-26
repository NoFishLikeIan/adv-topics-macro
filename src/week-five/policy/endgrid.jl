using Parameters, Interpolations
using Base.Threads
using Roots

using Logging, Printf

cartesianfromsize(Ns...) = collect(
    Iterators.product((1:N for N in Ns)...)
)


"""
Compute a policy using the endogenous grid method, 
assuming exogenous forecasting rule Ψ
"""
function endgrid_method(
    prices::NTuple{2,Function}, utils::NTuple{3,Function}, 
    Ψ::Function, model::DenHaanModel, grids_sizes::NTuple{2,Int};
    grid_bounds=[.01, 10.],
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

    m_grid = range(grid_bounds..., length=N_m)
    a_grid = range(grid_bounds..., length=N_a)

    policy = ones(N_a, N_m, D_z, D_ϵ)
    space = collect.(Iterators.product(a_grid, m_grid, ϵ_grid, z_grid)) 

    for iter in 1:max_iter
        g = (x -> max(x, 0.)) ∘ fromMtoFn(policy, a_grid, m_grid, ϵ_grid, z_grid)
        
        """
        Given the current policy function compute
            next period consumption.
        """
        function next_consumption(
            state::Tuple{Float64,Int}, 
            Ψ′::Float64, a′::Float64)
            z′, ϵ′ = state
            ϵ′ = convert(Float64, ϵ′)
        
            rate = R(z′, Ψ′)
            wage = w(z′, Ψ′)
        
            return rate * u′(rate * a′ + wage * ϵ′) - g(a′, Ψ′, ϵ′, z′)
        end

        # TODO: Nested loop or two loops?

        inv_policy = similar(policy)

        @threads for (i, j, k, l) in cartesianfromsize(N_a, N_m, D_z, D_ϵ)
            # Endogenous grid method for each a′
            a′, m = a_grid[i], m_grid[j]
            z, ϵ = z_grid[k], ϵ_grid[l]

            Ψ′ = Ψ(z, m)

            rhs = β * P_cond((z, ϵ))' * next_consumption.(ζ.S, Ψ′, a′)
            a = (invu′(rhs) + a′ - w(z, m) * ϵ) / R(z, m)

            inv_policy[i, j, k, l] = a
        end

        new_policy = similar(policy)

        @threads for (j, k, l) in cartesianfromsize(N_m, D_z, D_ϵ)
            origin_a = inv_policy[:, j, k, l]
            ixs = sortperm(origin_a)

            forward_policy = LinearInterpolation(origin_a[ixs], a_grid[ixs], extrapolation_bc=Line())
            new_policy[:, j, k, l] = @. max(forward_policy(a_grid), 0.0)
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
