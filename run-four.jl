# Environment variables
const scriptpath = rsplit(@__FILE__, "/", limit=2)[1]
const do_plot = false
verbose = get!(ENV, "VERBOSE", "false") == "true"
plotpath = joinpath(scriptpath, "src/week-four/solutions/plots/")

include("src/week-four/partial.jl")