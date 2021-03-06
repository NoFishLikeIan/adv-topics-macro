
using Distributions, KernelDensity
using Logging

function moments(a::Vector{Float64})
    return [mean(a), var(a)]
end

function distribution_mc(
    a′, a_grid::Vector{Float64}, model::Aiyagari;
    inits=1_000, verbose=false, max_iter=500, tol=1e-2)

    markov = model.y    
    sim = (y0) -> discrete_sim(markov, T=2, drop=1, y0=y0)[2]

    y0s = sample(markov.S, inits)
    a0s = sample(a_grid, inits)

    min_err = Inf

    for i in 1:max_iter
 
        y_next = sim.(y0s)
        a_next = a′.(a0s, y0s)
        
        err = norm(moments(a_next) - moments(a0s))

        if isnan(err)
            print(a_next, " - ", a0s, " - ", y0s, "\n\n")
        end
        if verbose
            min_err = min(err, min_err)
            print("$i / $max_iter : $err - $min_err \r")
        end

        if err < tol 
            print("Found stationary!\n\n")
            return kde(a_next)
        end

    end

    @warn "Algorithm did non converge with tol=$tol in $max_iter iterations"

    return kde(a0s)
    

end


