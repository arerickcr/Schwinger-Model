module CUDABackend
using CUDA, Adapt, ITensors, ITensorMPS
using ..SchwingerGPU: Backend
export cuda_backend

function cuda_backend(; device=0)
    CUDA.functional() || error("CUDA is unavailable. An NVIDIA GPU and working driver are required; use run_cpu.jl here.")
    CUDA.device!(device)
    CUDA.allowscalar(false)
    metadata = Dict{String,Any}("device"=>string(CUDA.name(CUDA.device())),
        "device_id"=>device,"CUDA_version"=>string(pkgversion(CUDA)),"precision"=>"Float64")
    # cu(...) would normally demote Float64 to Float32. Adapt preserves precision.
    upload = x -> adapt(CUDA.CuArray,x)
    download = x -> adapt(Array,x)
    return Backend("cuda",upload,download,CUDA.synchronize; metadata)
end
end
