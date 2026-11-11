# Explicit opt-in study. Never invoked by setup or tests.
include("SchwingerGPU.jl")
include("CUDABackend.jl")
using .SchwingerGPU, .CUDABackend
using LinearAlgebra, TOML, Statistics

function benchmark(args)
    opts = parse_options(args)
    BLAS.set_num_threads(parse(Int,get(opts,"blas_threads","4")))
    backend = cuda_backend(device=parse(Int,get(opts,"device","0")))
    parts = parse.(Int,split(get(opts,"sizes","10:20"),':'))
    sizes = length(parts)==1 ? parts : length(parts)==2 ? collect(parts[1]:parts[2]) : error("sizes must be N or first:last")
    isempty(sizes) && error("Empty size range")
    repeats = parse(Int,get(opts,"repeats","3")); repeats>0 || error("repeats must be positive")
    warmup = parse(Bool,get(opts,"warmup","true"))
    directory = get(opts,"directory",joinpath(@__DIR__,"results","benchmark_"*string(time_ns())))
    mkpath(directory)
    for N in sizes
        p,c = options_config(opts;N)
        if warmup
            # Same-size warmup initializes both library specializations and GPU workspaces.
            compare_case(p,c,backend;directory,label="warmup_N$N")
        end
        reports = [compare_case(p,c,backend;directory,label="N$(N)_repeat$(r)") for r in 1:repeats]
        cpu_median = median([r["cpu_dmrg_seconds"] for r in reports])
        gpu_median = median([r["gpu_dmrg_seconds"] for r in reports])
        summary = Dict("N"=>N,"repeats"=>repeats,"warmup"=>warmup,
            "median_cpu_seconds"=>cpu_median,"median_gpu_seconds"=>gpu_median,
            "median_speedup"=>cpu_median/gpu_median,
            "all_energies_agree"=>all(all(s["energy_agrees"] for s in r["states"]) for r in reports),
            "all_energy_converged"=>all(all(s["both_energy_converged"] for s in r["states"]) for r in reports))
        open(joinpath(directory,"N$(N)_summary.toml"),"w") do io
            TOML.print(io,summary;sorted=true)
        end
        println(summary)
    end
end
benchmark(ARGS)
