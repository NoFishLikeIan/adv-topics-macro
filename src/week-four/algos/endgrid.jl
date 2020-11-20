using Parameters, Interpolations

"""
Compute a policy using the endogenous grid method
"""
function endgrid(
    ai::Aiyagari, r::Float64, w::Float64;
    n_steps=1_000, upperbound=10., verbose=false)

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
    policy(v::Array{Float64}) = policy(v...)
    
    return policy, a_p
end 
