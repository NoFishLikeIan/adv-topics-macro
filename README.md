# Advanced Topics in Macro I

The repository for the code of the assignemnts of the course Advanced Topics in Macro I.

## Requirements

Requires `Julia 1.5`. Make sure you initialize the environment based on the `Project.toml` file. 

If one wishes to use the Julia repl, the `init.jl` file runs the boilerplate for you (environment and package installing). So something like,

```julia
include("init.jl") # deals with startup
include("src/$week/$my_file.jl")
```

## Assignments

The assignment of a given week can be found in the `src/` folder, for example, for week two in `src/week-two`. One can run all the code for a given assignment by simply doing,


```
julia run-one.jl

julia run-two.jl
```

and so on.

## LaTex submission

The `.pdf` file with the LaTex submissions can be found under `src/week-$i/solutions/sol.pdf`