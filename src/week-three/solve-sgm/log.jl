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

    """
    System:
    xy + (ρ - 1)*z - c_ss = 0
    z + w - k_ss^α = 0
    x + y - α * k_ss^(α - 1) = 0
    (y - 1) * x * k_ss^(α - 1) - y^(α - 2)*(α - 1)*c_ss = 0
    """
    function f!(F, x)
        x, y, z, w = x
        F[1] = x * y + (ρ - 1) * z - c_ss
        F[2] = z + w - k_ss^α
        F[3] = x + y - α * k_ss^(α - 1)
        F[4] = (y - 1) * x * k_ss^(α - 1) - y^(α - 2) * (α - 1) * c_ss
    end

    init_x = [0.5, 0.5, 0.5, 0.5]
    res = mcpsolve(f!, [0., 0., 0., 0.], [Inf, Inf, Inf, Inf], init_x)
    c_k, k_k, c_y, k_y = res.zero

    c_policy(k, y) = max(c_ss + c_k * (k - k_ss) / k_ss + c_y * (y - y_ss) / y_ss, 0.01)
    k_policy(k, y) = max(k_ss + k_k * (k - k_ss) / k_ss + k_y * (y - y_ss) / y_ss, 0.01) # Avoid negative values for k too far from k_ss


    return c_policy, k_policy
end