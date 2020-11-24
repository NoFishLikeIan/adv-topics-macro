include("partial.jl")
include("aiyagari.jl")

using Parameters
using QuadGK, Roots

ϵ = 1e-3

function solvegeneral(
    model::Aiyagari; verbose=false, n_steps=75)

    @unpack δ, β, a_ = model

    bounds = [ϵ - δ, 1 / β - (1 + ϵ)]

    F, F_k, F_l, invF_k = make_F(model)

    K = (r) -> invF_k(r + δ)
    w = (r) -> F_l(K(r))

    function A(r::Float64)
        
        λ, a′, a_grid = solvepartial(
                model, r, w(r),
                verbose=false, 
                n_steps=n_steps, upperbound=50.,
                mc=false, end_grid=true)
                
        return quadgk(a -> a * λ(a), model.a_, a_grid[end])[1]

    end

    return A, K, w
end