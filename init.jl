# Run this file to install all packages listed in Project.toml/Manifest.toml
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
Pkg.precompile()