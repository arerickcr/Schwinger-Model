# Unchanged CPU reference sources

These are the two local Julia files against which this GPU port was developed and tested. They are preserved byte-for-byte, with hashes in `../reference`. They are not replacements for the separate examples at the repository root.

`../tests/runtests.jl` includes only `OBC_Schwinger_dmrg.jl`, inside a module. It does not execute the N=20 driver. The driver contains historical user-specific commented paths and is retained as source provenance, not as the recommended entry point. Use `../run_cpu.jl` or, after NVIDIA validation, `../run_cuda.jl`.
