# Environment variables
const scriptpath = rsplit(tst, "/", limit=2)[1]
const plot = false
verbose = get!(ENV, "VERBOSE", "false") == "true"
plotpath = joinpath(scriptpath, "solutions/plots/")

include("src/week-four/partial.jl")