# Get the square root of machine epsilon
sq_ϵ = √eps(Float64)

function constructtwosided(f::Function)::Function

    function f_prime(x::Float64)::Float64
        h = max(abs(x), 1.) * sq_ϵ
        return (f(x + h) - f(x - h)) / 2h
    end

    function f_prime(x::Array{Float64,1})::Array{Float64,1}
        h = max(abs.(x)..., 1)
        n = length(x)
        return  [ f(x + h * e(i, n)) - f(x - h * e(i, n)) for i in 1:n ] ./ 2h # Is this the best way?
    end

    return f_prime
end

function constructonesided(f::Function)::Function

    function f_prime(x::Float64)::Float64
        h = max(x, 1.) * sq_ϵ
        return (f(x + h) - f(x)) / h
    end

    function f_prime(x::Array{Float64,1})::Array{Float64,1}
        h = max(abs.(x)..., 1)
        n = length(x)
        return  [ f(x + h * e(i, n)) - f(x) for i in 1:n ] ./ h
    end

    return f_prime
end 