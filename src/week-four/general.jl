include("partial.jl")
include("aiyagari.jl")

using Parameters
using QuadGK, Roots

ϵ = 1e-4

function solve_general(
    model::Aiyagari; verbose=false, n_steps=200, mc=true)

    @unpack δ, β, a_ = model

    bounds = [ϵ - δ, 1 / β - (1 + ϵ)]

    F, F_k, F_l, invF_k = make_F(model)

    K = (r) -> invF_k(r + δ)
    w = (r) -> F_l(K(r))

    function A(r::Float64)
        λ, a′, a_grid = solvepartial(model, r, w(r), verbose=true, n_steps=n_steps, mc=mc)

        return quadgk(a -> a * λ(a), model.a_, a_grid[end])[1]
    end

    clearprices(r) = A(r) - K(r)

    @time interest = find_zero(clearprices, bounds, verbose=verbose)

    return interest
end


model = Aiyagari(
    0.9, 0.1, 7, # AR process parameters ρ, σ, N
    0, 0.95, 0.33, 0.1, 2 # Model parameters a, β, α, δ, σ_u
)

r = solve_general(model, verbose=true, n_steps=100, mc=false)

