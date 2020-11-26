using StatsBase

"""
Implements the Krusell-Smith Algorithm over a cross-sectional distribution λ
for a Den Haan et al. model.
"""
function krusellsmith(
    model::DenHaanModel;
    b0=[0., 1.],
    ρ=0.7, ϵ_m=1e-3, ϵ_g=1e-3,
    kwargs...
)
    u, u′, inv_u′ = makeutility(model)
    R, w = makeproduction(model)

    ho_Ψ(b0, b1) = (z, K) -> exp(b0 + b1 * log(K)) # FIXME: How is this depending on z?
    Ψ = ho_Ψ(b0...) # Initial guess for forecasting rule Ψ

    policy = pfi(
        makeproduction(model),
        makeutility(model),
        Ψ, model, (5, 5);
        ρ=ρ, kwargs...)

end