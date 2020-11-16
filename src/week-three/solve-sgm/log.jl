function loglin_poicy(model)

    parameters = model.calibration.flat

    α, β = parameters[:alpha], parameters[:beta]
    ρ, σ = parameters[:rho], parameters[:sigma]

    y_ss, k_ss, c_ss = analytical_ss(model)

    c_h(c) = log(c / c_ss)
    k_h(k) = log(k / k_ss)
    y_h(y) = log(y / y_ss)
 
    yk_ss = y_ss * k_ss^α

    A = [1 (1 - α); 0 k_ss]
    B = [1 0; -c_ss yk_ss * α]
    C = [1; yk_ss]

    ℵ = A \ B
    ℶ = A \ C

    function x(c, k, y) 
        x_h = [c_h(c); k_h(k)]
        next_x = ℵ * x_h  + y_h(y) .* ℶ

        return @. [c_ss; k_ss] * exp(next_x)
    end

    return x
end