include("SchwingerGPU.jl")
include("CUDABackend.jl")
using .SchwingerGPU, .CUDABackend
main(ARGS, opts -> cuda_backend(device=parse(Int,get(opts,"device","0"))))
