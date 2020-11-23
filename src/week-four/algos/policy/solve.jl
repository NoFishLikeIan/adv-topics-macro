using Parameters

"""
Computes the discrete expectation
"""
function E(vals::Array{Float64}, density::Array{Float64})::Float64
    return sum(vals .* density)
end


function fromVtoFn(x::Vector{Float64}, f::Vector{Float64})

    intp = LinearInterpolation((x,), f, extrapolation_bc=(Linear(),))

    function fn(x::Float64) intp(x) end

    return fn
    
end

function fromMtoFn(
    x::Vector{Float64}, y::Vector{Float64}, m::Matrix{Float64})


    intp = LinearInterpolation((x, y,), m, extrapolation_bc=(Linear(), Linear()))

    function fn(x::Float64, y::Float64) intp(x, y) end

    function fn(v::Vector{Float64}) return fn(v...) end

    return fn
end

function cartesian(x, y)
    return collect.(Iterators.product(x, y))
end


function policysolve(
    ai::Aiyagari, r::Float64, w::Float64;
    n_steps=100, upperbound=10., end_grid=false, kwargs...)

    R = 1 + r

    as = collect(range(ai.a_, upperbound, length=n_steps))

    solver = end_grid ? endgrid : iterate_pfi

    a′ = solver(as, ai, R, w; kwargs...)

    return a′, as
end