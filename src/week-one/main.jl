include("utils/interp.jl")
include("utils/array.jl")

using Plots

# --- Excercise 3

space = collect(range(0, 1, step=.1))
evaluations = exp.(space)

thickerspace = collect(range(0, 1, step=.01))

linear_int = constructlinear(space, evaluations)
interpolated = linear_int.(thickerspace)


# --- Exercise 

size = 100000000

test = vcat(range(0, 1, length=size), range(1, 0, length=50))

@time find_maximum(test, mode="unknown")
@time find_maximum(test, mode="concave")