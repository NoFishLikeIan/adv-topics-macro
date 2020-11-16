include("numerical/derivatives.jl")
include("../comm_utility.jl")

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

c_grid = collect(range(-2π, 2π, length=100))

techniques = Dict(
    "One sided" => constructonesided,
    "Two sided" => constructtwosided
)

print("Testing x^2...\n")

test_diff(
    x -> x^2, 
    x -> 2x, 
    c_grid, do_plot=true, filename="squared")

print("Testing u(c)...\n")

test_diff(
    c -> -1 / c,
    c -> 1 / c^2, 
    c_grid, do_plot=true, filename="utility")

print("Testing multivariate...\n")

mul_grid = collect.(Iterators.product(c_grid, c_grid))

function f(v::Array{Float64})::Float64
    x, y = v
    return sin(x^2) * cos(y)
end

function f_p(v::Array{Float64})::Array{Float64}
    x, y = v
    return [
        2 * x * cos(x^2) * cos(y),
        -1 * sin(x^2) * sin(y)
    ]
end

test_diff(
    f, 
    f_p, 
    mul_grid, do_plot=true, filename="trigonometric")

test_diff(
    x -> x[1]^3 + x[2]^3, 
    x -> [3 * x[1]^2, 3 * x[2]^2 ], 
    mul_grid, do_plot=true, filename="multivariate")


if false

    function f(v::Array{Float64})::Array{Float64}
        x, y, z = v 

        return sin(x^2) * exp(-y) * z
    end

    function f_p(v::Array{Float64})::Array{Float64}
        x, y, z = v

        return [
            2 * x * cos(x^2) * exp(-y) * z,
            - f(v),
            sin(x^2) * exp(-y)
        ]
    end

    test_diff(
        f, f_p, 
        collect.(Iterators.product(c_grid, c_grid, c_grid)), do_plot=true, filename="trivariate"
    )

end