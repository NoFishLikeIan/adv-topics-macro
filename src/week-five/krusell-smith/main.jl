using StatsBase

"""
Implements the Krusell-Smith Algorithm over a cross-sectional distribution λ
for a Den Haan et al. model.
"""
function krusellsmith(
    model::DenHaanModel;
    b0=[1., 1.],
    ρ=0.7, N_a=100, N_m=10,
    ϵ_m=1e-3, ϵ_g=1e-3,
    kwargs...
)
    u, u′, inv_u′ = makeutility(model)
    R, w, τ = makeproduction(model)

    ho_Ψ(b0, b1) = (z, K) -> exp(b0 + b1 * log(K)) # FIXME: How is this depending on z?
    Ψ = ho_Ψ(b0...) # Initial guess for forecasting rule Ψ

    policy = endgrid_method(
        makeproduction(model),
        makeutility(model),
        Ψ, model, (N_a, N_m);
        grid_bounds=[0.01, 20.],
        ρ=ρ, kwargs...)


    simulated_a = economysim(policy, model; kwargs...)

end

