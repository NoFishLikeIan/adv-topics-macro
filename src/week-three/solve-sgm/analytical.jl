function analytical_policy(model)
    parameters = model.calibration.flat

    α, β = parameters[:alpha], parameters[:beta]
    ρ, σ = parameters[:rho], parameters[:sigma]

    k_policy(k, y) = α * β * y * k^α
    c_policy(k, y) = y * k^α - k_policy(k, y)

    return c_policy, k_policy
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

function analytical_eee(model)
    parameters = model.calibration.flat

    α, β = parameters[:alpha], parameters[:beta]
    ρ, σ = parameters[:rho], parameters[:sigma]

    c, k_p = analytical_policy(model)


    f(k, y) = β * (α * y * k_p(k, y)^(α - 1) ) / c(k_p(k, y), y)

    return (k, y) -> log10(abs(1 - f(k, y)))
end