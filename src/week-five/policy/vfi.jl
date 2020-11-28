using Parameters, Interpolations
using Base.Threads

using Logging, Printf

first(x) = x[1]

"""
Compute a policy using value function iteration, 
assuming exogenous forecasting rule Ψ
"""
function vfi( 
    Ψ::Function, model::DenHaanModel, grids_sizes::NTuple{2,Int};
    grid_bounds=[.01, 10.],
    max_iter=1_000, tol=1e-3, ρ=0.7, 
    verbose=false)

    @unpack S_ϵ, S_z, ζ = model
    @unpack β, l, δ, μ = model

    u, u′, invu′ = makeutility(model)
    R, w, τ = makeproduction(model)
    c, invc = makeconsumption(model)

    ϵ_grid = collect(Float64, S_ϵ)
    z_grid = collect(S_z)
    D_z, D_ϵ = length(S_ϵ), length(S_z)
    P_cond(state) = ζ.P[findfirst(==(state), ζ.S), :]

    N_a, N_m = grids_sizes

    m_grid = range(grid_bounds..., length=N_m)
    a_grid = range(grid_bounds..., length=N_a)

    # Space for a′, a, m, z, ϵ
    space = collect.(Iterators.product(a_grid, a_grid, m_grid, z_grid, ϵ_grid)) 

    utility = (u ∘ positive ∘ c).(space)
    V = ones(N_a, N_m, D_z, D_ϵ)

    E(V::Array{Float64,4}) = reshape(reshape(V, N_a * N_m, D_z * D_ϵ) * ζ.P, size(V))
    euler_step(V::Array{Float64,4}) = reshape(reshape(utility, (N_a, N_a * N_m * D_z * D_ϵ)) .+ reshape(β * E(V), 1, N_a * N_m * D_z * D_ϵ), size(utility))

    for iter in 1:max_iter
        H = euler_step(V)
        
        H_max, policy_idx = findmax(H, dims=1)

        V′ = reshape(H_max, size(V))

        d = V′ - V
        err_distance = maximum(abs.(d))
        verbose && print("Iteration $iter / $max_iter: $(@sprintf("%.4f", err_distance)) \r")

        if err_distance < tol 
            verbose && print("Found policy in $iter iterations (|x - x'| = $(@sprintf("%.4f", err_distance))\n")

            # get the a′ of the maximized matrix
            policy = (first).(reshape(space[policy_idx], size(V)))
            g = positive ∘ fromMtoFn(policy, a_grid, m_grid, z_grid, ϵ_grid)
            return (a, m, z, ϵ) -> g(a, Ψ(z, ϵ), z, ϵ)
        end

        V += ρ * d
    end

    throw(ConvergenceException(max_iter))
end
