using Optim, LinearAlgebra, StatsBase, Roots

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

function rouwenhorst_numerical(proc::Process, N::Int, x0, m)

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

"""
Solves the analytical rouwenhorst process for the special case E[z] = 0
"""
function zero_process(Vz::Float64, corr::Float64, N::Int)
    q = (corr + 1) / 2
    
    # For p = q, s = 0.5, and 4s(1-s) = 1.
    # Write ψ and p as a functional for consistency
    ψ = (_) -> √(Vz * (N - 1))
    p = (_) -> q

    return ψ, p, q

end

function rouwenhorst_analy(proc::Process, N::Int)
    Ez = mean(proc)
    Vz = var(proc)
    corr = autocor(proc)

    q0 = 0.8 * rand(Uniform(.9, 1.1))

    if Ez == 0
        ψ, p, q = zero_process(Vz, corr, N)
    else
        N_scale = -1 + 1 / (N - 1)

        p(q) = corr + 1 - q
        ψ(q) = Ez * (2 - p(q) - q) / (q - p(q))
        s(q) = (1 - p(q)) / (2 - p(q) - q)
    
        four_s(q) = 4 * s(q) * (1 - s(q)) 
    
        function f(q)
            return ψ(q)^2 * (1 + four_s(q) * N_scale) - Vz
        end
    
        q = find_zero(f, q0)
    end
    
    return colnormalized(P_N(N)(p(q), q)), makepartition(ψ(q), N)

end

function rouwenhorst(proc::Process, N::Int; numerical=false, x0=[5., 0.5, 0.5], m=3)
    if numerical
        return rouwenhorst_numerical(proc, N, x0, m)

    else
        return rouwenhorst_analy(proc, N)
    end

end