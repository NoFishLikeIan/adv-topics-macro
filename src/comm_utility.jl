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