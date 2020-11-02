include("detgrowth.jl")

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

    function k_2(k_1::Real, k::Real)::Real
        return F(k_1) - model.β * F_k(k_1) * (F(k) - k_1)
    end

    return k_2
end