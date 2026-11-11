# Base-only parse check, including NVIDIA files on machines without CUDA.
root = dirname(@__DIR__)
for file in sort(readdir(root;join=true))
    endswith(file,".jl") || continue
    Meta.parseall(read(file,String);filename=file)
    println("Parsed ",basename(file))
end
