include("detgrowth.jl")

using Roots

function equil_k(model::DetGrowthModel; unique::Bool=true)

    k_analytic = model.z * (model.β / (1 + model.β))^(1 / (1 - model.α))

    if unique
        return k_analytic
    else        
        ab = model.β * model.α
        constant = model.z * (1 + ab)
        
        fn(k) = k^(1 - model.α) + ab * model.z^2 * k^(model.α - 1) - constant
        
        return find_zeros(fn, 0, k_analytic * 10)
    end
end

"""
Construct function to compute the net capital transition
    and its derivative
"""
function construct_netK(model::DetGrowthModel)
    F(k::Real)::Real = model.z * (k^model.α)
    F_k(k::Real)::Real = model.α * model.z * (k^(model.α - 1))
    F_inv(x::Real)::Real = (x / model.z)^(1 / model.α)

    return F, F_k, F_inv
end


"""
Constructs the difference equation of second order for K:
    k_{t+2} = g(k_t, k_{t+1})
"""
function construct_2dk(model::DetGrowthModel)
    F, F_k, _ = construct_netK(model)

    function k_2(k_t::Real, k_t1::Real)::Real
        F(k_t1) - model.β * F_k(k_t1) * (F(k_t) - k_t1)
    end

    return k_2
end

"""
Constructs the inverse difference equation of second order for K:
    k_t = g(k_{t+2}, k_{t+1})
"""
function construct_inv2dk(model::DetGrowthModel)
    F, F_k, F_inv = construct_netK(model)

    function k_t(k_t2::Real, k_t1::Real)::Real
        x = (F(k_t1) - k_t2) / (model.β * F_k(k_t1)) + k_t1
        return F_inv(x)
    end

    return k_t
end
