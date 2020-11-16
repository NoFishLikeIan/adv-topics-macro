function loglin_poicy(model)

    parameters = model.calibration.flat

    α, β = parameters[:alpha], parameters[:beta]
    ρ, σ = parameters[:rho], parameters[:sigma]

    y_ss, k_ss, c_ss = analytical_ss(model)

    c_h(c) = (c - c_ss) / c_ss
    k_h(k) = (k - k_ss) / k_ss
    y_h(y) = (y - y_ss) / y_ss
 
    yk_ss = y_ss * k_ss^α

    A = [1 (1 - α); 0 k_ss]
    B = [1 0; -c_ss yk_ss * α]
    C = [1; yk_ss]

    ℵ = A \ B
    ℶ = A \ C

    function x(c, k, y) 
        x_h = -[c_h(c); k_h(k)]
        next_x = ℵ * x_h  + y_h(y) .* ℶ

        return @. [c_ss; k_ss] * exp(next_x)
    end

    c(k, y) = x(c_ss, k, y)[1]
    k(k, y) = x(c_ss, k, y)[2]

    return c, k
end

function implicit_policy(model)
    parameters = model.calibration.flat

    α, β = parameters[:alpha], parameters[:beta]
    ρ, σ = parameters[:rho], parameters[:sigma]

    y_ss, k_ss, c_ss = analytical_ss(model)

    akss = α * k_ss^(α - 1)

    k_k = 0.5 * (√(4 * c_ss - akss^2) + akss)
    k_y = k_ss^(α) 

    c_k = 0.5 * (√((k_ss * akss^2 - 4 * c_ss * k_ss^2) / k_ss^2) - akss)
    
    c_policy(k, _) = c_ss + c_k * (k - k_ss) / k_ss
    k_policy(k, y) = k_ss + k_k * (k - k_ss) / k_ss + k_y * (y - y_ss) / y_ss


    return c_policy, k_policy
end

function loglin_eee(model)

    parameters = model.calibration.flat

    α, β = parameters[:alpha], parameters[:beta]
    ρ, σ = parameters[:rho], parameters[:sigma]


    c_p, k_p = loglin_poicy(model)


    f(k, y) = β * (α * y * k_p(k, y)^(α - 1) ) / c_p(k_p(k, y), y)

    return (k, y) -> log10(abs(1 - f(k, y)))
end

function implicit_eee(model)

    parameters = model.calibration.flat

    α, β = parameters[:alpha], parameters[:beta]
    ρ, σ = parameters[:rho], parameters[:sigma]


    c_p, k_p = implicit_policy(model)

    f(k, y) = β * (α * y * k_p(k, y)^(α - 1) ) / c_p(k_p(k, y), y)

    return (k, y) -> log10(abs(1 - f(k, y)))
end