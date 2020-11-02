using Random, Distributions, Plots

# Random.seed!(11148705)

include("model/detgrowth.jl")
include("model/dynamics.jl")
include("utils/types.jl")

function shooting(
    k_start::Real, k_end::Real, fn_k2::Function;
    size_grid::Int=5000, steps::Int=200, ε::Real=0.01)::Union{RealArray,Nothing}

    k_path = zeros(steps + 1)
    k_path[1] = k_start

    k_grid = range(0, k_start * 2, length=size_grid)

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

    return nothing

end 

function compute_transition(first::DetGrowthModel, last::DetGrowthModel, mode::String)
    k_start = equil_k(first)

    fn_k2 = construct_2dk(last)

    construct_path = mode == "shooting" ? shooting : shooting

    equils = equil_k(last, unique=false)

    print("Found ", length(equils), " equilibria! $equils\n\n")

    return [construct_path(k_start, k_end, fn_k2) for k_end in equils]


end

alpha = 0.5
beta = 0.9

first = DetGrowthModel(beta, alpha, 1., 0.01, 1_000)
last = DetGrowthModel(beta, alpha, 2., 0.01, 1_000)

path = compute_transition(first, last, "shooting")

# plot(path[1])
plot(path[2])

savefig("plots/transition.png")