using Parameters, Interpolations

"""
Compute a policy using the endogenous grid method
"""
function endgrid(
    ai::Aiyagari, r::Float64, w::Float64;
    n_steps=1_000, upperbound=10.)

    u, u′, inv_u′ = make_u(ai)
    @unpack β, a_, y = ai 
    ys = y.S
    Γ = y.P

    R = 1 + r
    a_p = collect(range(a_, upperbound, length=n_steps))
    double_grid = collect.(Iterators.product(a_p, ys))

    a_pp = ones((n_steps, length(ys)))

    next_mu = @. u′(R * a_p + w * ys' - a_pp)
    rhs = β * R * next_mu * Γ'
    a = @. (inv_u′(rhs) - w * ys' + a_p) / R
        
    intp = LinearInterpolation((a_p, ys,), a, extrapolation_bc=(Linear(), Linear()))

    function policy(a::Float64, y::Float64)::Float64
        max(intp(a, y), a_)
    end
    
    return policy, a_p
end 


function a_todensity(
    Q::Array{Float64,3}, 
    Q_aprime::Array{Float64,2}, 
    a_grid::Array{Float64,1})

    function populate(a::Float64) 
        dens_vec = zeros(length(a_grid))
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

    return reshape(hcat(populate.(Q_aprime)...), size(Q))
end

function computeQ(a′::Function, a_grid::Array{Float64}, ai::Aiyagari)
    ys = ai.y.S
    Γ = ai.y.P

    double_grid = collect.(Iterators.product(a_grid, ys))
    N, T = size(double_grid)

    Q_a = zeros((N, N, T))
    
    Q_aprime = (v -> a′(v...)).(double_grid)

    Q_a = a_todensity(Q_a, Q_aprime, a_grid)

end 