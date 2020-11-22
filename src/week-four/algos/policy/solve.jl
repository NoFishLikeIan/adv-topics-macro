using Parameters

function policysolve(
    ai::Aiyagari, r::Float64, w::Float64;
    n_steps=1_000, upperbound=10., kwargs...)

    R = 1 + r

    as = collect(range(a_, upperbound, length=n_steps))

    a′ = iterate_pfi(as, ai, R, w; kwargs...)

    return a′, as
end