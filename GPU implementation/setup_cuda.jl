# Run once on the NVIDIA machine: julia --project=. setup_cuda.jl
# Only this GPU implementation project is changed; the original files and global environment are untouched.
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
Pkg.add("CUDA")
Pkg.precompile()
println("CUDA installed. Next: julia --project=. check_cuda.jl")
