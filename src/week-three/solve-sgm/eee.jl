function eee(c_p, k_p, model)

    parameters = model.calibration.flat

    α, β = parameters[:alpha], parameters[:beta]
    ρ, σ = parameters[:rho], parameters[:sigma]

    f(k, y) = β * (α * y * k_p(k, y)^(α - 1) ) / c_p(k_p(k, y), y)

    return (k, y) -> log10(abs(1 - f(k, y)))
end