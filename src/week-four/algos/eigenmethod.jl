using LinearAlgebra, Arpack, Kronecker, SparseArrays


function a_todensity(
    Q_aprime::Array{Float64,2}, 
    a_grid::Array{Float64,1})

    function populate(a::Float64) 
        dens_vec = spzeros(length(a_grid))
        for (k, ak) in enumerate(a_grid)
            # :(
            if a < ak
                if k == 1
                    dens_vec[k] = 1
                    return dens_vec
                else
                    ak_1 = a_grid[k - 1]
                    dens_vec[k - 1] = (ak - a) / (ak - ak_1)
                    dens_vec[k] = (a - ak_1) / (ak - ak_1)
                    
                    return dens_vec
                end
            end
        end

        dens_vec[end] = 1
        return dens_vec
    end

    return hcat(populate.(Q_aprime)...)
end

function computeQ(a′::Function, a_grid::Array{Float64}, ai::Aiyagari)
    ys = ai.y.S
    Γ = ai.y.P

    double_grid = collect.(Iterators.product(a_grid, ys))
    N, T = size(double_grid)
    
    Q_aprime = a′.(double_grid)

    Q_a = reshape(a_todensity(Q_aprime, a_grid), (N, N, T))

    Q = spzeros(N * T, N * T)

    for j in 1:T
        st = (j - 1) * N + 1
        en = j * N 

        Q[st:en, :] = Γ[j, :]' ⊗ Q_a[:, :, j]'
    end

    droptol!(Q, 1e-10)

    return Q, (N, T)
end 

function distribution_eigenvector(
    a′::Function, a_grid::Array{Float64}, ai::Aiyagari; verbose=true)

    if verbose print("Computing Q matrix...\n\n") end
    Q, dims = computeQ(a′::Function, a_grid::Array{Float64}, ai::Aiyagari)

    if verbose print("Finding eigenvalues...\n\n") end
    λ, Φ = eigs(Q', nev=1, which=:LR)

    if verbose print("Associated eigenvalue $(real(λ))\n\n") end

    return reshape(real(Φ) ./ sum(real(Φ)), dims)

end