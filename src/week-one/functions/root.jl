mid(a, b) = 0.5 * (a + b)

function bisection(f, a::Real, b::Real, ε::Real=.01, max_iter::Int=1000)
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

        if i > max_iter break end
    end

    return c
end

function newton(f, a::Real, b::Real, ε::Real=.01, max_iter::Int=1000)

end


fn(x) =  sin(2π * x) - 2 * x

c = bisection(fn, 5, 6)
