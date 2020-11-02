include("detgrowth.jl")


function equil_k(model::DetGrowthModel)
    model.z * (model.β / (1 + model.β))^(1 / (1 - model.α))
end

"""
Construct function to compute the net capital transition
    and its derivative
"""
function construct_netK(model::DetGrowthModel)
    F(k::Real)::Real = model.z * (k^model.α)
    F_k(k::Real)::Real = model.α * model.z * (k^(model.α - 1))

    return F, F_k
end


"""
Constructs the difference equation of second order for K:
    k_{t+2} = g(k_t, k_{t+1})
"""
function construct_2dk(model::DetGrowthModel)
    F, F_k = construct_netK(model)

    function k_2(k_t::Real, k_t1::Real)::Real
        return F(k_t1) - model.β * F_k(k_t1) * (F(k_t) - k_t1)
    end

    return k_2
end

