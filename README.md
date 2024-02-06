# Workflow

To add a new package dependency, start a new Julia REPL, switch to package mode (]) and run :
] add PackageName.jl

To try out a piece of code in file.jl, start a new Julia REPL (ALT J + ALT R) and run :
> includet("file.jl")

To profile using nsys :
$ nsys launch --trace=cuda,nvtx julia
And then in julia :
> Cuda.@profile external=true ...
