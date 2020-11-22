using Interpolations, NLsolve
using Logging

function fromMtoFn(
    x::Vector{Float64}, y::Vector{Float64}, m::Matrix{Float64})


    intp = LinearInterpolation((x, y,), m, extrapolation_bc=(Linear(), Linear()))

    function fn(x::Float64, y::Float64) intp(x, y) end

    function fn(v::Vector{Float64}) return fn(v...) end

    return fn
end

function cartesian(x::Vector{Real}, y::Vector{Real})
    return collect.(Iterators.product(x, y))
end

function iterate_pfi(
    as::Vector{Float64}, ai::Aiyagari, R::Float64, w::Float64;
    max_iter=1_000, tol=1e-3, verbose=false
)

    u, ∂u∂c, inv_∂u∂c = make_u(ai)
    @unpack β, a_, y = ai

    support_y = y.S
    ys = collect(support_y)
    ϕ = y.transformation
    Γ = y.P

    N, T = length(as), length(ys)

    a_grid = copy(as)
    domain = cartesian(a_grid, ys)

    function construct_ap(a_prime::Matrix{Float64})

        function a_n(a::Float64, y::Float64) 
            return fromMtoFn(as, ys, a_prime)(a, y)
        end

        function a_n(v::Vector{Float64})
            return a_n(v...)
        end

        function a_n(a::Matrix{Float64}, y::Vector{Float64})
            next_a = similar(a)
            for j in 1:length(y)
                col_a = a[:, j]
                next_a[:, j] = a_n.(col_a, y[j])
            end

            return next_a
        end

        return a_n
    end

    a_prime = 0.5 * ones(N, T)

    function factory_euler(a::Float64, y::Float64, a′′::Function)
        function euler(a_p::Float64)
            rhs = β * R * (Γ[get_row(support_y, y), :]' * (@. ∂u∂c(R * a_p + w * ϕ(ys) - a′′(a_p, ϕ(ys)))))
            lhs =  ∂u∂c(R * a + w * ϕ(y) - a_p)
            return rhs - lhs
        end
    end

    for iter in 1:max_iter
        a′′ = construct_ap(a_prime)
        next_ap = copy(a_prime)
        
        for (i, j) in Iterators.product(1:N, 1:T)
            a, y = domain[i, j]

            function tosolve!(F, x) F[1] = factory_euler(a, y, a′′)(x[1]) end
            root = mcpsolve(tosolve!, [a_], [Inf], [a_prime[i, j]]) 

            next_ap[i, j] = max(root.zero[1], a_)
        end

        err_distance = matrix_distance(a_prime, next_ap)

        verbose && print("Iteration $iter / $max_iter -> $err_distance\r")

        if err_distance < tol 
            verbose && print("Found policy in $iter iterations (|x - x'| = $(@sprintf("%.4f", err_distance))\n")
            return a′′
        end 

        a_prime = next_ap
    end

    @warn "Could not find policy in $max_iter iterations with tolerance $tol"

    return a′′

end    