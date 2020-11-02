include("utils/interp.jl")
include("utils/array.jl")

using Plots

if false
    # --- Excercise 3

    low = 0
    up = 2π

    space = collect(range(low, up, step=1))
    evaluations = exp.(space)

    thickerspace = collect(range(low + 0.1, up - 0.1, step=.1))

    linear_int = constructlinear(space, evaluations)
    interpolated = linear_int.(thickerspace)


    # --- Exercise 4
    size = 100_000_000

    test = vcat(range(0, 100, length=size), range(99, 0, length=50))

    @time naive_max(test)
    @time concave_max(test)
end

xs = collect(range(-2, 2, step=0.1))
fn(x) = sin(2π * x) - 2 * x