include("utils/interp.jl")
include("utils/array.jl")
include("functions/root.jl")

using Plots

# --- Excercise 3

low = 0
up = 2π

space = collect(range(low, up, step=1))
evaluations = exp.(space)

thickerspace = collect(range(low + 0.1, up - 0.1, step=.1))

linear_int = constructlinear(space, evaluations)
interpolated = linear_int.(thickerspace)

plot(thickerspace, interpolated, title="Interpolated exponential", label="exp(x)", lw=2)
savefig("src/week-one/plots/linear_interpolation.png")

# --- Exercise 4
size = 100_000_000

test = vcat(range(0, 100, length=size), range(99, 0, length=50))

@time naive_max(test)
@time concave_max(test)

# --- Excercise 5

functions = [
    [x -> sin(2π * x) - 2 * x, "sin(2π x) - 2x"],
    [x -> 2π * x - 2 * x, "2π x - 2x"],
    [x -> 2π * x - x / 2, "2π x - x/2"]
]


for (fn, body) in functions

    print("-- Evaluating: ", body, "\n")

    a, b = -2, 2

    @time bis = bisection(fn, a, b)
    @time new = newton(fn, a, b)
    @time sec = secant(fn, a, b)

    @printf(
        "
        Bisection: %.3f, \n
        Newton: %.3f, \n 
        Secant: %.3f \n
        ------ \n
        ", 
        bis, new, sec
    )
end