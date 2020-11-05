include("../process.jl")
include("../markov-types.jl")
include("utils.jl")


"""
Computes the parametric (p, q) -> P matrix in a recursive manner.
FIXME: This is very inefficient both in memory and allocation time.
"""
function P_N(N)
    if N == 2
        return (p, q) -> [p (1 - p); q (1 - q)]
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

function rouwenhorst(proc::Process, N::Int, m::Int)

    if N > 100 @warn "Current implementation is too slow for big matrices" end

end