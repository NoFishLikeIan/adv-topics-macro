include("utils/interp.jl")
include("utils/array.jl")

using Plots

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

@time find_maximum(test, mode="unknown")
@time find_maximum(test, mode="concave")

