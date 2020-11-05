using Pkg

Pkg.activate(".")
Pkg.instantiate() # this will install the packages listed in Project.toml

using DotEnv

DotEnv.config()