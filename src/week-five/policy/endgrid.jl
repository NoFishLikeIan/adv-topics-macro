include("../../commons/interpolation.jl")
include("../denhaan.jl")

using Parameters, Interpolations
using Base.Threads

using Logging, Printf


"""
Compute a policy using the endogenous grid method, 
assuming exogenous forecasting rule Ψ
"""
function endgrid(
    prices::NTuple{2,Function}, utils::NTuple{3,Function}, 
    Ψ::Function, model::DenHaanModel,
    grids_sizes::NTuple{2,Int})

    @unpack S_ϵ, S_z, ζ = model
    u, u′, invu′ = utils
    R, w = prices

    ϵ_grid = collect(S_ϵ)
    z_grid = collect(S_z)
    D_z, D_ϵ = length(S_ϵ), length(S_z)
    
    N_a, N_m = grids_sizes
    m_grid = range(0.01, 10., length=N_m)
    a_grid = range(0.01, 10., length=N_m)

    policy = ones(N_a, N_m, D_z, D_ϵ)
    g = fromMtoFn(policy, a_grid, m_grid, ϵ_grid, z_grid)

    space = Iterators.product(a_grid, m_grid, ϵ_grid, z_grid)

    @threads for (a, m, z, ϵ) in collect(space)
        # Solve Euler for ā

        Ψ′ = Ψ(z, m)

        function righteuler(a′)
            function inner(state::Tuple{State,State})
                z′, ϵ′ = state
                ϵ′ = convert(Float64, ϵ′)

                rate = R(z′, Ψ′)
                wage = w(z′, Ψ′)

                return rate * u′(rate * a′ + wage * ϵ′) - g(a′, Ψ′, ϵ′, z′)
            end
        end
        
    end

end
