using Pkg

Pkg.activate(".")
Pkg.instantiate() # this will install the packages listed in Project.toml

# Runs excercise 1
include("src/week-one/vfi.jl")

# Runs excercise 2
include("src/week-one/transition.jl")

# Runs excercise 3 and 5
include("src/week-one/extra.jl")