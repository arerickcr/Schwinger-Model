include("SchwingerGPU.jl")
using .SchwingerGPU
main(ARGS, _ -> cpu_backend())
