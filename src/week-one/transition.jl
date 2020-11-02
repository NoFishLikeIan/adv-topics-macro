using Random, Distributions

Random.seed!(11148705)

include("model/detgrowth.jl")
include("model/dynamics.jl")
include("utils/types.jl")

function shooting(k_start::Real, k_end::Real, fn_k2::Function, ε::Real, steps::Int)::RealArray
    k_path = zeros(steps + 1)
    k_path[1] = k_start

    while abs(k_path[end] - k_end) > ε

        k_t1 = rand(TruncatedNormal(k_start, 1, 0, Inf))
        
        k_path[2] = k_t1

        for i in 3:steps
            k_t2 = fn_k2(k_path[i - 2], k_path[i - 1])

            if k_t2 < 0 break end

            k_path[i] = k_t2 
        end

    end

    return k_path
end 

function compute_transition(first::DetGrowthModel, last::DetGrowthModel, steps::Int, mode::String)
    k_start = equil_k(first)
    k_end = equil_k(last)

    fn_k2 = construct_2dk(last)

    if mode == "shooting"
        return shooting(k_start, k_end, fn_k2, 0.01, steps)

    else
        throw(ErrorException("Not implemented mode, $mode"))
    end

end

first = DetGrowthModel(0.9, 0.5, 1., 0.01, 1_000)
last = DetGrowthModel(0.9, 0.5, 2., 0.01, 1_000)

path = compute_transition(first, last, 200, "shooting")