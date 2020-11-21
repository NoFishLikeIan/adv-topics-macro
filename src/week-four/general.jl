include("partial.jl")
include("aiyagari.jl")

using Parameters
using QuadGK, KernelDensity
using Roots


function solve_general(
    model::Aiyagari; verbose=false, grid_N=1_000)
    @unpack δ, β, a_ = model

    bounds = [-δ, -1 + 1 / β]

    F, F_k, F_l, invF_k = make_F(model)

    K = (r) -> invF_k(r + δ)
    w = (r) -> F_l(K(r))

    function A(r::Float64)
        Φ, a′, a_grid = solvepartial(model, r, w(r), verbose=verbose, grid_N=grid_N)
        return quadgk(a -> a * pdf(Φ, a), model.a_, Inf)[1]
    end

    clearprices(r) = A(r) - K(r)

    @time interest = find_zero(clearprices, bounds, verbose=verbose)

    return interest
end


model = Aiyagari(
    0.9, 0.1, 7, # AR process parameters ρ, σ, N
    0, 0.95, 0.33, 0.1, 2 # Model parameters a, β, α, δ, σ_u
)

r = solve_general(model, verbose=true, grid_N=500)

