# Workflow

To add a new package dependency, start a new Julia REPL, switch to package mode (]) and run :
] add PackageName.jl

To try out a piece of code in file.jl, start a new Julia REPL (ALT J + ALT R) and run :
> includet("file.jl")

To profile using nsys :
$ nsys launch --trace=cuda,nvtx julia
And then in julia :
> Cuda.@profile external=true ...

To profile using ncu :
$ ncu --mode=launch julia
This will hang on the first call to the CUDA API. Launch the GUI :
$ ncu-ui
Select an interactive profile, attach to the julia process, and resume the application (run until kernel launch), then profile the kernel.