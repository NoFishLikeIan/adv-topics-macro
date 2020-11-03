using Plots, Printf, NLsolve

include("model/detgrowth.jl")
include("model/dynamics.jl")
include("utils/types.jl")

function shooting(
    k_start::Real, k_end::Real, fn_k2::Function;
    size_grid::Int=5000, steps::Int=200, ε::Real=0.01, pert::Real=1.
    )::Union{RealArray,Nothing}

    k_path = zeros(steps + 1)
    k_path[1] = k_start

    k_grid = range(max(0, k_start - pert), k_start + pert, length=size_grid)
    print("Searching from $k_start to $k_end, around $k_grid\n")

    for k_t1 in k_grid
        
        k_path[2] = k_t1

        for i in 3:length(k_path)

            k_t2 = fn_k2(k_path[i - 2], k_path[i - 1])

            if k_t2 < 0 break end

            k_path[i] = k_t2 
        end

        if abs(k_path[end] - k_end) < ε
            return k_path
        end

    end

    return k_path

end 

function system(k_start::Real, k_end::Real, fn_k2::Function; steps::Int=200)
    iota = ones(steps)
    mid = (k_end - k_start) / 2

    lower = 0 * iota
    upper = 100 * iota
    init_vec = mid * iota

    init_F = similar(init_vec)
    
    function f!(F, K)
        F[1] = k_start - K[1]
        F[2] = fn_k2(k_start, K[1]) - K[2]

        for i in 3:length(steps) - 1
            F[i] = fn_k2(K[i - 2], K[i - 1]) - K[i]
        end

        F[end] = fn_k2(K[end - 2], K[end - 1]) - K[end]

    end

    df = OnceDifferentiable(f!, init_vec, init_F)

    return nlsolve(df, init_vec)

end


"""
Computes the transition between to steady states of capital using one of three techniques,

- shooting
- inverse shooting
- a system of 3 equations
"""
function compute_transition(first::DetGrowthModel, last::DetGrowthModel, mode::String)
    k_start = equil_k(first)

    if mode == "shooting"

        construct_path = shooting
        fn_k2 = construct_2dk(last)

    elseif mode == "inverse shooting"

        construct_path = (f(start, last, fn) = shooting(last, start, fn))
        fn_k2 = construct_inv2dk(last)

    elseif mode == "system"
        fn_k2 = construct_2dk(last)

        function construct_path(k_start::Real, k_end::Real, fn_k2::Function; size_grid::Int=5000, steps::Int=200, ε::Real=0.01, pert::Real=1.)
            sol = system(k_start, k_end, fn_k2; steps=steps)
            
            return sol.zero
        end
        
    else 
        throw(ErrorException("Mode $mode not implemented"))
    
    end

    equils = equil_k(last, unique=false)

    print("Found ", length(equils), " equilibria! $equils\n
        Solving with $mode\n")

    return [construct_path(k_start, k_end, fn_k2) for k_end in equils]
end

alpha = 0.5
beta = 0.9

first = DetGrowthModel(beta, alpha, 1., 0.01, 1_000)
last = DetGrowthModel(beta, alpha, 2., 0.01, 1_000)


for method in ["shooting", "inverse shooting", "system"]
    path = compute_transition(first, last, method)

    plot(path[2], title="Transition path", label="k", lw=2)

    savefig("src/week-one/plots/transition_$method.png")
end