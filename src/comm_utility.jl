function try_n(routine, n, exception; verbose=false)
    for _ in 1:n
        try
            return routine()
        catch exc
            if exc isa exception 
                if verbose
                    print("Trying again...")
                end
            else
                throw(exception)
            end 
        end
    end

    error("Routine failed $n times with $exception")
end

function matrix_distance(A::Array{Float64,2}, B::Array{Float64,2})::Float64
    N, M = size(A)

    d = 0.

    for col in 1:M
        m = maximum(abs.(A[:, col] - B[:, col]))
        d = m > d ? m : d
    end

    return d
end

toupper(x, i) = x[1:i - 1] * uppercase(x[i:i]) * x[i + 1:end]
