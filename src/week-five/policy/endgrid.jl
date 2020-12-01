using Parameters, Interpolations
using Base.Threads

using Logging, Printf

"""
Compute a policy using the endogenous grid method, 
assuming exogenous forecasting rule Ψ
"""
function endgrid_method( 
    Ψ::Function, model::DenHaanModel, grids_sizes::NTuple{2,Int};
    ρ0=0.8,
    grid_bounds=[.01, 10.],
    max_iter=1_000, tol=1e-3, 
    verbose=false)

    @unpack S_ϵ, S_z, ζ = model
    @unpack β, l, δ, μ = model

    u, u′, invu′ = makeutility(model)
    R, w, τ = makeproduction(model)
    cons, invcons = makeconsumption(model)

    c = positive ∘ cons
    invc = positive ∘ invcons

    ϵ_grid = collect(Float64, S_ϵ)
    z_grid = collect(S_z)
    D_z, D_ϵ = length(S_ϵ), length(S_z)
    P_cond(state) = ζ.P[findfirst(==(state), ζ.S), :]

    N_a, N_m = grids_sizes

    m_grid = range(grid_bounds..., length=N_m)
    a_grid = range(grid_bounds..., length=N_a)

    policy = repeat(a_grid, 1, N_m, D_z, D_ϵ)
    space = collect.(Iterators.product(a_grid, m_grid, z_grid, ϵ_grid)) 

    for iter in 1:max_iter
        g = positive ∘ fromMtoFn(policy, a_grid, m_grid, z_grid, ϵ_grid)
        
        function next_value(state::Tuple{Float64,Int}, Ψ′, a′)
            z′, ϵ′ = state
            ϵ′ = convert(Float64, ϵ′)
            return R(z′, Ψ′) * u′(c(a′, Ψ′, z′, ϵ′, g))
        end

        inv_policy = similar(policy)

        @threads for (i, j, k, l) in cartesianfromsize(N_a, N_m, D_z, D_ϵ)
            # Endogenous grid method for each a′
            a′, m = a_grid[i], m_grid[j]
            z, ϵ = z_grid[k], ϵ_grid[l]

            Ψ′ = Ψ(z, m)

            rhs = β * P_cond((z, ϵ))' * next_value.(ζ.S, Ψ′, a′)
            c′ = invu′(rhs)
            a = invc(c′, m, z, ϵ, a′)

            inv_policy[i, j, k, l] = a
        end

        new_policy = similar(policy)

        @threads for (j, k, l) in cartesianfromsize(N_m, D_z, D_ϵ)
            origin_a = inv_policy[:, j, k, l]
            ixs = sortperm(origin_a)

            forward_policy = LinearInterpolation(origin_a[ixs], a_grid[ixs], extrapolation_bc=Line())
            new_policy[:, j, k, l] = forward_policy.(a_grid)
        end

        d = new_policy - policy
        err_distance = maximum(abs.(d))
        verbose && print("Iteration $iter / $max_iter: $(@sprintf("%.4f", err_distance)) \r")

        if err_distance < tol 
            verbose && print("Found policy in $iter iterations (|x - x'| = $(@sprintf("%.4f", err_distance)))\n")
            return g
        end

        ρ = ρ0 - iter / max_iter # Dynamic dumping parameter
        policy += ρ * d # Update with dumping parameter
    end

    throw(ConvergenceException(max_iter))
end
