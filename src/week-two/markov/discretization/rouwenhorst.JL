using Optim, LinearAlgebra, StatsBase, NLsolve

include("utils.jl")

include("../process.jl")
include("../markov-types.jl")
include("../simulation.jl")

L_BN = [3., 0., 0.]
U_BN = [Inf, 1., 1.]

"""
Computes the parametric (p, q) -> P matrix in a recursive manner.
FIXME: This is very inefficient both in memory and allocation time.
"""
function P_N(N)
    if N == 2
        return (p, q) -> [p (1 - p); (1 - q) q]
    end

    P_p = P_N(N - 1)  

    function nested_P(p, q)
        Ψ = P_p(p, q)
        n, n = size(Ψ)

        Ω = zeros((n + 1, n + 1))

        # TODO: Find a nicer way to do this
        Ω[1:n, 1:n] += p * Ψ
        Ω[1:n, end - n + 1:end] += (1 - p) * Ψ
        Ω[end - n + 1:end, 1:n] +=  (1 - q) * Ψ
        Ω[end - n + 1:end, end - n + 1:end] +=  q * Ψ

        return Ω
    end


    return nested_P
end

function rouwenhorst_numerical(proc::Process, N::Int, x0)

    if N > 100 @warn "Current implementation is too slow for big matrices" end

    P = P_N(N)

    function sim_moments(ψ::Float64, p::Float64, q::Float64)
        S = makepartition(ψ, N)

        markov = MarkovDiscrete(colnormalized(P(p, q)), S)

        sim = discrete_sim(markov)

        return [mean(sim), var(sim)]
    end

    moments = [mean(proc), var(proc)]

    function minimize(params::Array{Float64}) 
        ψ, p, q = params
        return norm(sim_moments(ψ, p, q) - moments)
    end

    inner_optimizer = GradientDescent()

    res = optimize(minimize, L_BN, U_BN, x0, Fminbox(inner_optimizer))

    if Optim.converged(res)
        minimum = Optim.minimum(res)
        if minimum > 0.05 @warn "Matching the moments led to a minimum of $minimum > .05" end

        ψ, p, q = Optim.minimizer(res)

        return colnormalized(P(p, q)), makepartition(ψ, N)
    
    else 
        throw(
            StatsBase.ConvergenceException(
                Optim.iterations(res)
            )
        )
    end
end



function rouwenhorst_analy(proc::Process, N::Int, x0)
    Ez = mean(proc)
    Vz = var(proc)
    corr = autocor(proc)

    function f!(F, parameters)
        ψ, q, p = parameters

        s = (1 - p) / (2 - (p + q))

        four_s = 4 * s * (1 - s)

        F[1] = (q - p) * ψ / (2 - (p + q)) - Ez
        F[2] = ψ^2 * (1 - four_s - four_s / (N - 1)) - Vz
        F[3] = p + q - 1 - corr

    end

    res = mcpsolve(f!, L_BN, U_BN, x0)

    ψ, q, p = res.zeros

    return colnormalized(P(p, q)), makepartition(ψ, N)

end

function rouwenhorst(proc::Process, N::Int; method="anal", x0=[5., 0.5, 0.5])
    if method == "num"
        return rouwenhorst_numerical(proc, N, x0)

    else
        return rouwenhorst_analy(proc, N, x0)
    end

end


if false
    N = 3
    ρ = 0.7
    μ = 0
    σ = 1

    σ_proc = σ / sqrt(1 - ρ^2)

    autocov = exp(0.5 * ((1 + ρ)^2 * σ_proc^2 + σ^2)) - exp(σ_proc^2)

    z = Process(
        LogNormal(μ, σ_proc),
        z -> z^ρ,
        Normal(μ, σ),
        autocov)

    res = rouwenhorst(z, N; method="anal")

end