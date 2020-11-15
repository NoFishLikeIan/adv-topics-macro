function analytical_policy(model)
    parameters = model.calibration.flat

    α, β = parameters[:alpha], parameters[:beta]
    ρ, σ = parameters[:rho], parameters[:sigma]

    k_policy(k, y) = α * β * y * k^α


    function c_policy(k, y)
        k_next = k_policy(k, y)

        return y * k^α - k_policy(k, y)
    end

    return c_policy
end

function analytical_ss(model)
    parameters = model.calibration.flat

    α, β = parameters[:alpha], parameters[:beta]
    ρ, σ = parameters[:rho], parameters[:sigma]

    y = exp(0.5 * σ^2 / (1 - ρ^2))
    k = (α * β * y)^(1 / (1 - α))
    c = y * k^α - k

    return y, k, c

end