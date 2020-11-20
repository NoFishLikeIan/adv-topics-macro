include("numerical/derivatives.jl")
include("../commons/matrix.jl")

using Plots, StatsPlots, LinearAlgebra

plot_path = "src/week-three/solutions/plots/derivative"

Grid = Union{Array{Float64,1},Array{Array{Float64,1}}}

function test_diff(
    f::Function, f_prime::Function, x_grid::Grid; 
    do_plot=false, filename="errors")

    multivariate = length(size(x_grid)) > 1

    analy_der = f_prime.(x_grid)

    errors = Dict{String,Array{Float64}}()
    evals = Dict()

    for (name, derivative) in techniques

        df = derivative(f)

        eval_deriv = df.(x_grid)

        error = norm.(analy_der .- eval_deriv)

        errors[name] = multivariate ? sum(error, dims=2) : error
        evals[name] = eval_deriv

    end

    if do_plot

        x = multivariate ? [el[1] for el in x_grid[:, 1]] : x_grid

        fn_name = toupper(filename, 1)
        plot(title="$fn_name derivative error", dpi=600)

        for (name, err_vec) in errors
            plot!(x, err_vec, label="Error $name")
        end

        savefig("$plot_path/$filename.png")

    end

    return errors, evals

end

c_grid = collect(range(0, 2, length=100))

techniques = Dict(
    "One sided" => constructonesided,
    "Two sided" => constructtwosided
)

print("Testing u(c)...\n")

u(c) = - 1 / c
u_p(c) = 1 / c^2

test_diff(u, u_p,
    c_grid, do_plot=true, filename="utility")

print("Testing production...\n")

const µ = 0.05 
const σ = 0.6
const λ = 0.25
const ρ = 0.2

G(K, N, H) = μ * N^σ + (1 - μ) * (λ * K^ρ + (1 - λ) * H^ρ)^(σ / ρ)

F(K::Float64, N::Float64, H::Float64) = G(K, N, H)^(1 / σ)
function F(v::Array{Float64})
    K, N, H = v
    return F(K, N, H)
end

"""
The analytical derivative of the production function
"""
function F_p(K::Float64, N::Float64, H::Float64) 
    
    outer = (1 / σ) * G(K, N, H)^((1 - σ) / σ)
    outer_kh = (σ * (1 - μ) / ρ) * (λ * K^ρ + (1 - λ) * H^ρ)^(-1 + σ / ρ)

    return outer * [
        outer_kh * λ * ρ * K^(ρ - 1),
        σ * μ * N^(σ - 1),
        outer_kh * (1 - λ) * H^(ρ - 1)
    ]

end


function F_p(v::Array{Float64})
    K, N, H = v
    return F_p(K, N, H)
end


grid = range(0.01, .99, length=10)
mul_grid = collect.(Iterators.product(grid, grid, grid)) 

errors, evals = test_diff(F, F_p,
    mul_grid, do_plot=false, filename="production")