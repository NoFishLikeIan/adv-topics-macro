include("numerical/derivatives.jl")
include("../comm_utility.jl")

using Plots, StatsPlots

plot_path = "src/week-three/solutions/plots"

Grid = Union{Array{Float64,1},Array{Array{Float64,1}}}

function test_diff(
    f::Function, f_prime::Function, x_grid::Grid; 
    do_plot=false, filename="errors")

    analy_der = f_prime.(x_grid)

    errors = Dict{String,Array{Float64}}()

    for (name, derivative) in techniques

        @time eval_deriv = derivative(f).(x_grid)

        error = abs.(analy_der .- eval_deriv)

        errors[name] = error

    end

    if do_plot

        fn_name = toupper(filename, 1)
        plot(title="$fn_name derivative error", dpi=600)

        for (name, error) in errors
            plot!(x_grid, error, label="Error $name")
        end

        savefig("$plot_path/$filename.png")

    end

    return errors

end

c_grid = collect(range(-1., 1., length=100))

techniques = Dict(
    "One sided" => constructonesided,
    "Two sided" => constructtwosided
)

print("Testing x^2...\n")

test_diff(
    x -> x^2, 
    x -> 2x, 
    c_grid, do_plot=false, filename="squared")

print("Testing u(c)...\n")

test_diff(
    c -> -1 / c,
    c -> 1 / c^2, 
    c_grid, do_plot=false, filename="utility")

mul_grid = collect.(Iterators.product(c_grid, c_grid))

test_diff(
    (x) -> x[1] + x[2],
    (x) -> [1, 1], 
    mul_grid, do_plot=false, filename="multivariate")


