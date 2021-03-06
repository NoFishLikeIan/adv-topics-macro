using LinearAlgebra, Arpack, Kronecker, SparseArrays
using QuantEcon
using Interpolations

function eigenstationary(Q::SparseMatrixCSC)
    return [1; (I - Q[2:end,2:end]) \ Vector(Q[2:end,1])]
end

"""
Given the policy function matrix and the grid, 
    maps the a′ values as a density
"""
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

    return hcat(populate.(Q_aprime)...)'
end


"""
Creates the transition matrix Q ∈ R(NT, NT) as a Markov
matrix of (a, y) -> (a, y)′
"""
function computeQ(
        a′::Function, a_grid::Array{Float64}, ai::Aiyagari;
        fact_finer::Int=2
    )
    ys = ai.y.S
    Γ = ai.y.P

    finer_grid = collect(
        range(a_grid[1], a_grid[end], length=length(a_grid) * fact_finer)
    )

    double_grid = collect.(Iterators.product(finer_grid, ys))
    N, T = size(double_grid)
    
    Q_aprime = (v -> a′(v...)).(double_grid)

    Q_a = a_todensity(Q_aprime, finer_grid)

    Q = spzeros(N * T, N * T)

    for j in 1:T
        st = (j - 1) * N + 1
        en = j * N 

        Q[st:en, :] = Q_a[st:en, :] ⊗ Γ[j, :]'
    end

    droptol!(Q, 1e-10)

    return Q, (N, T), finer_grid
end 


"""
Computes the stationary distribution of the Aiyagari
model with the eigenvector method, given a policy function
"""
function distribution_eigenvector(
    a′::Function, a_grid::Array{Float64}, model::Aiyagari; 
    gth=true, verbose=false)

    verbose && print("Computing Q matrix...\n")
    Q, dims, finer_grid = computeQ(a′, a_grid, model)

    verbose && print("Finding stationary distribution...\n")
    
    Φ_bar = gth ? QuantEcon.gth_solve(collect(Q)) : eigenstationary(Q)

    return reshape(sparse(Φ_bar), dims), finer_grid

end