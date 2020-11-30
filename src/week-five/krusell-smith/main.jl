using StatsBase

function regress_m(logm::Vector{Float64}, boom::BitArray{1})
    
    afterboom = circshift(boom, 1) # All periods that follow a boom

    gX = hcat(ones(size(logm[boom])), logm[boom])
    gY = logm[afterboom]

    β_g = inv(gX'gX) * gX'gY

    bX = hcat(ones(size(logm[.!boom])), logm[.!boom])
    bY = logm[.!afterboom]

    β_b = inv(bX'bX) * bX'bY
    
    return β_b, β_g
end

"""
Implements the Krusell-Smith Algorithm over a cross-sectional distribution λ
for a Den Haan et al. model.
"""
function krusellsmith(
    model::DenHaanModel;
    ρ=0.7, N_a=100, N_m=10,
    ϵ_m=1e-1, ϵ_a=1e-2,
    max_iter=1_000,
    kwargs...
)   

    verbose = kwargs[:verbose]

    @unpack S_z = model
    u, u′, inv_u′ = makeutility(model)
    R, w, τ = makeproduction(model)

    function ho_Ψ(B_g, B_b)
        function Ψ(z, K)
            boom = (z == S_z[end])
            b0, b1 = boom ? B_g : B_b
            Ψ′ = exp(b0 + b1 * log(K))
            return Ψ′
        end

        return Ψ
    end

    B_g = B_b = [1., 2.] # Initial guess for forecasting rule Ψ

    for iter in 1:max_iter
        Ψ = ho_Ψ(B_g, B_b)

        policy = endgrid_method(Ψ, model, (N_a, N_m); tol=ϵ_a, ρ0=1., kwargs...)

        as, zs = economysim(policy, model; kwargs...)
        ms = log.(mean(as, dims=2))

        booms = (zs .== S_z[end])

        B_b′, B_g′ = regress_m(vec(ms), booms)
        
        d_g = B_g′ - B_g
        d_b = B_b′ - B_b

        err_distance = maximum(abs.(vcat(d_g, d_b)))

        verbose && print("Ψ iteration: $iter / $max_iter: $(@sprintf("%.4f", err_distance))...\n\n")

        if err_distance < ϵ_m 
            verbose && print("Found policy in $iter iterations (|x - x'| = $(@sprintf("%.4f", err_distance))\n\n")
            return policy
        end

        B_g += ρ * d_g
        B_b += ρ * d_b

    end
end

