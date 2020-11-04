using Pkg, DotEnv

DotEnv.config()

Pkg.activate(".")
Pkg.instantiate() # this will install the packages listed in Project.toml