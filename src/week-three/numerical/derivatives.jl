# Get the square root of machine epsilon
sq_ϵ = √eps(Float64)

function constructtwosided(f::Function)::Function

    function f_prime(x::Real)::Real

        h = max(x, 1.) * sq_ϵ

        return (f(x + h) - f(x - h)) / 2h
    end

    return f_prime
end

function constructonesided(f::Function)::Function

    function f_prime(x::Real)::Real

        h = max(x, 1.) * sq_ϵ

        return (f(x + h) - f(x)) / h
    end

    return f_prime
end 