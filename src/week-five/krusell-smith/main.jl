using StatsBase

"""
Implements the Krusell-Smith Algorithm over a cross-sectional distribution λ
for a Den Haan et al. model.
"""
function krusellsmith(
    model::DenHaanModel;
    ρ=0.7, ϵ_m=1e-3, ϵ_g=1e-3
)

    u, u′, inv_u′ = makeutility(model)
    r, w = makeproduction(model)

    ho_Ψ(b0, b1) = (x) -> exp(b0 + b1 * log(x))
    Ψ = ho_Ψ(0., 1.)

    function R(z::Float64, m̄::Vector{Float64})
        1 + r(m̅[1]) - δ
    end

end