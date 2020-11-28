using StatsBase

"""
Implements the Krusell-Smith Algorithm over a cross-sectional distribution λ
for a Den Haan et al. model.
"""
function krusellsmith(
    model::DenHaanModel;
    B0=[0., 0., 0., 1.],
    ρ=0.7, N_a=100, N_m=10,
    ϵ_m=1e-3, max_iter=1_000,
    kwargs...
)
    u, u′, inv_u′ = makeutility(model)
    R, w, τ = makeproduction(model)

    function ho_Ψ(B)
        function Ψ(z, K)
            boom = z == 1.01 ? 1. : 0.

            return exp(B' * [1., boom, boom * log(K), log(K)])
        end

        return Ψ
    end

    B = B0 # Initial guess for forecasting rule Ψ

    X = Array{Float64}(undef, T - 1, 4) # regressors

    for iter in 1:max_iter
        Ψ = ho_Ψ(B)

        verbose && print("Computing policy...\n")
        policy = endgrid_method(
            Ψ, model, (N_a, N_m);
            tol=3e-1,
            grid_bounds=[0.01, 10.],
            ρ=ρ, kwargs...)

        verbose && print("simulating economy...\n")
        as, zs = economysim(policy, model; kwargs...)
        ms = log.(mean(as, dims=2))

        X[:, 1] .= 1.
        X[:, 2] = [z == 1.01 ? 1. : 0. for z in zs[1:end - 1]]
        X[:, 3] = ms[1:end - 1] .* X[:, 2]
        X[:, 4] = ms[1:end - 1]

        y = ms[2:end]

        B′ = inv(X' * X) * X'y
        verbose && print("OLSing, new coefficients $(B′)...\n")

        d = B′ - B

        err_distance = maximum(abs.(d))

        verbose && print("Ψ iteration: $iter / $max_iter: $(@sprintf("%.4f", err_distance))...\n")

        if err_distance < ϵ_m 
            verbose && print("Found policy in $iter iterations (|x - x'| = $(@sprintf("%.4f", err_distance))\n")
            return policy
        end

        B = B′

    end
end

