using Printf, ForwardDiff

mid(a, b) = 0.5 * (a + b)

notdone(max_iter) = throw(ErrorException("Could not find a solution in $max_iter iterations"))

function bisection(f::Function, a::Real, b::Real; 
    ε::Real=.01, max_iter::Int=1000)::Real

    inter = b - a

    while sign(f(a)) == sign(f(b)) 
        a += ε * inter
        b -= ε * inter
    end

    c = mid(a, b)
    i = 0

    while abs(f(c)) > ε
        i += 1

        if f(a) * f(c) < 0
            b = c
        else
            a = c
        end

        c = mid(a, b)

        if i > max_iter notdone(max_iter) end
    end

    return c
end

function newton(f::Function, a::Real, b::Real;
    ε::Real=.01, max_iter::Int=1000)::Real

    x = rand() * (b - a) + a
    i = 0

    g(x) = ForwardDiff.derivative(f, x)

    while abs(f(x)) > ε
        x = x - f(x) / g(x)
        i += 1
        if i > max_iter notdone(max_iter) end
    end

    return x
end

function secant(f::Function, a::Real, b::Real;
    ε::Real=.01, max_iter::Int=1000)
    prev_x = rand() * (b - a) + a
    prev_fx = f(prev_x)

    x = rand() * (b - a) + a
    i = 0

    while abs(f(x)) > ε
        f_x = f(x)
        Δf = f_x - prev_fx
        Δx = x - prev_x 

        new_x = x - f_x * (Δx / Δf)

        prev_x = x
        prev_fx = f_x

        x = new_x

        i += 1
        if i > max_iter notdone(max_iter) end
    end

    return x
end
