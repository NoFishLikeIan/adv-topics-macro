using SparseArrays
using Base.Threads

function density_policy(a_grid::Vector{Float64})

    function populate(a::Float64) 

        dens_vec = zeros(length(a_grid))
        for (k, ak) in enumerate(a_grid) # :(
            if a < ak
                if k == 1
                    dens_vec[k] = 1.
                    return dens_vec
                else
                    ak_1 = a_grid[k - 1]
                    dens_vec[k - 1] = (ak - a) / (ak - ak_1)
                    dens_vec[k] = (a - ak_1) / (ak - ak_1)
                        
                    return dens_vec
                end
            end
        end

        dens_vec[end] = 1.
        return dens_vec
    end

end

"""
Simulate the economy for N individuals and T periods
    given a policy function
"""
function economysim_det(
    g::Function, model::DenHaanModel;
    N::Int=10_000, T::Int=1_500, J::Int=100, verbose=false)

    S_ϵ = collect(Float64, model.S_ϵ)
    S = length(S_ϵ)
    
    zs = simulation(model.Z, T)
    
    a_grid = collect(range(.01, 10., length=J)) 
    function findgrid(x) 
        for (k, a_k) in enumerate(a_grid)
            if x < a_k return k end
        end

        return length(a_grid)
    end

    λ0 = rand(Uniform(), (J, length(model.S_ϵ)))

    λ = Array{Float64}(undef, T, J, S)
    λ[1, :, :] = λ0 ./ sum(λ0)

    for (t, z) in enumerate(zs[1:end - 1])
        verbose && print("Simulating economy $t / $T\r")
        z′ = zs[t + 1]
        π_dens = π_ϵ′(z′, z, model)
        λ_t = λ[t, :, :]

        λ′ = zeros(size(λ_t))
        m = sum(λ_t .* a_grid)

        @threads for i in 1:J
            a = a_grid[i]
            
            a′ = g.(a, S_ϵ, z, m)
            k = findgrid.(a′) 
            prev_k = k .- 1

            for i in 1:2 
                upper_prop = k[i] == J ?  1. : (a_k′[i] - a′[i]) / (a_k′[i] - a_k[i])
                lower_prop = k[i] == J ?  0. : (a′[i] - a_k[i]) / (a_k′[i] - a_k[i])

                λ′[k[i], i] += ((π_dens'upper_prop) .* λ_t[i, :])[i]
                λ′[prev_k[i], i] += ((π_dens'lower_prop) .* λ_t[i, :])[i]
            end
        end

        λ[t + 1, :, :] = λ′ / sum(λ′)

    end 

    verbose && print("\n")

    return as, zs
end

