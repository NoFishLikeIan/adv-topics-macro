using SparseArrays
using Base.Threads

function columnormal(M::Matrix{Float64})
    M ./ sum(M, dims=1)
end

function normal(M::Matrix{Float64})
    M ./ sum(M)
end

"""
Simulate the economy for N individuals and T periods
    given a policy function
"""
function economysim_det(
g::Function, model::DenHaanModel;
T::Int=1_500, J::Int=1000, verbose=false)

    S_ϵ = collect(Float64, model.S_ϵ)
    S = length(S_ϵ)

    zs = simulation(model.Z, T)

    a_grid = collect(range(.01, 100., length=J)) 

    λ0 = rand(Uniform(a_grid[1], a_grid[end]), (J, length(model.S_ϵ)))

    λ = Array{Float64}(undef, T, J, S)
    λ[1, :, :] = λ0

    for (t, z) in enumerate(zs[1:end - 1])
        verbose && print("Simulating economy $t / $T\r")

        λt = λ[t, :, :]
        m = mean(λt)
        a′ = g.(λt, m, z, S_ϵ')
        
        λt′ = zeros(J, S)

        @threads for (j, s) in cartesianfromsize(J, S)
            ϵ = S_ϵ[s]
            x = a′[j, s]

            k = findfirst(a_k -> x < a_k, a_grid)
            dens = π_ϵ′(ϵ, zs[t + 1], z, model)
    
            if isnothing(k)
                λt′[end, :] += dens * x
            elseif k == 1
                λt′[1, :] += dens * x
            else
                ω = 1. - (x - a_grid[k - 1]) / (a_grid[k] - a_grid[k - 1])
                λt′[k - 1, :] += ω * dens * x
                λt′[k, :] += (1 - ω) * dens * x
            end
        end

        λ[t + 1,:, :] = λt′

    end

    verbose && print("\n")

    as = reshape(mean(λ, dims=3), (T, J))

    return as, zs
end

