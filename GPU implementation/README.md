# GPU implementation: Schwinger DMRG

> **Experimental implementation - GPU validation pending.** Published on November 11, 2026 for review and future testing. No NVIDIA execution, CPU/GPU agreement, GPU speedup, or N=10 through 20 benchmark is claimed. Complete the checklist in [VALIDATION.md](VALIDATION.md) when hardware results are available.

This folder adds a portable CPU solver and an NVIDIA CUDA backend to the repository. The existing repository Julia scripts are unchanged. The exact local CPU sources used to develop and check this port are preserved unchanged in `reference/`; they are distinct from the repository-root examples.

Read [the detailed PDF guide](docs/Schwinger_GPU_guide.pdf) for the Hamiltonian derivation, explanation of both original programs, every new component, validation, and the planned CPU/GPU study.

## Status and hardware

- Portable CPU solver and backend interface: implemented.
- NVIDIA CUDA adapter: implemented, **not yet executed on an NVIDIA GPU**.
- This computer reports **Intel UHD Graphics 620**, not an NVIDIA device. Do not use its integrated GPU with `run_cuda.jl`.
- A general backend interface is provided, not a claim that all GPU vendors can run ITensor DMRG. CPU is the tested portable fallback. AMD/Metal/Intel adapters are not supplied as validated implementations.
- N=10 through 20 benchmarking is a separate, explicit command, to run after reviewing the guide. Setup and tests do not launch that study.

Historical CPU validation (October 22, 2026): **115 checks passed**. N=4,5,6 ground-state energies agree with independent exact diagonalization within 1.5e-14; original-helper, excited-state, entropy, and report checks also pass. See `results/validation_log.txt` and `VALIDATION.md`. This does not validate CUDA execution.

## Setup and CPU check

From the repository root, enter the folder first. All subsequent commands run inside it:

```text
cd "GPU implementation"
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. tests/runtests.jl
julia --project=. run_cpu.jl --N=4 --maxdim=16 --nsweeps=12
```

The CPU project has a manifest initially seeded from this machine's Julia 1.12.4 environment (ITensors 0.9.31, ITensorMPS 0.4.1). CUDA setup modifies only the new project's dependencies; keep its resulting manifest with future results. Julia may precompile packages on first use.

## NVIDIA machine

Clone the repository or copy this entire folder, including `reference/`, then enter `GPU implementation`. Install Julia and a compatible NVIDIA driver, then:

```text
julia --project=. setup_cuda.jl
julia --project=. check_cuda.jl
julia --project=. run_cuda.jl --N=4 --maxdim=16 --nsweeps=12
```

`check_cuda.jl` tests actual device storage, Float64 preservation, and N=4 CPU/GPU agreement. It must pass before interpreting GPU results. CUDA failures are errors, not silent CPU fallbacks. The CUDA package supplies its runtime dependencies; a standalone system CUDA toolkit is not normally needed.

## Later: the requested N=10,...,20 study

```text
julia --project=. benchmark.jl --sizes=10:20 --repeats=3 --warmup=true
```

Defaults match the original physics and main solver controls: x=10, y=1, mass=0.33, theta=pi, g=0, 10 sweeps, maxdim=100, cutoff=0, noise=1e-12, Krylov dimension=10, eigensolver maxiter=1. Energy convergence uses a fresh, zero-safe observer for each state. All floating-point tensors use Float64.

Both even and odd sizes are included. Initial spin defaults to 0 for even N and +0.5 for odd N. `--spin=0` with odd N is rejected. Dense tensors do **not** enforce a total-spin sector: this option sets only the product-state initialization.

Warmup runs each size on each backend and is excluded from the three measured repeats. Reports include energy differences, fidelity, residuals, truncation errors, convergence flags, and median DMRG time. Matching energy does not prove convergence. Small-system GPU runs can be slower than CPU runs.

Use the same changes on both backends when improving convergence, for example:

```text
julia --project=. benchmark.jl --sizes=10:20 --nsweeps=30 --maxdim=200 --eig_maxiter=5 --noise=0 --repeats=3
```

This is a second convergence experiment, not the original numerical settings. No dense 2^N by 2^N Hamiltonian is constructed in production; the independent exact-matrix tests stop at N=6.

## Options and outputs

Options use `--key=value`. Model: `N,x,y,mass,theta,g,spin`. Solver: `nsweeps,maxdim,cutoff,noise,energy_tol,min_sweeps,krylovdim,eig_maxiter,eig_tol,num_states,seed,outputlevel`. Runtime: `blas_threads` (4), `device` (0), `output`. Benchmark: `sizes`, `repeats`, `warmup`, `directory`. Theta is numeric in radians (default pi); use `--theta=3.141592653589793`, not `--theta=pi`.

Each run writes TOML under `results/` unless an output path is specified. Default single-run filenames are overwritten on repeat; specify `--output=...` to preserve distinct experiments. Benchmark directories are unique by default. The reports save scalars and histories, **not MPS checkpoints**. Julia API return values contain states for further analysis.

`entropy_profile` and `middle_tensors` operate on host states returned by diagnostics. Entropy includes every bond and skips zero probabilities. Excited-state `energy` is the normalized expectation of the physical Hamiltonian; `dmrg_objective` also includes the orthogonality penalty.

See the PDF for function-level explanations, limitations of variance subtraction and energy stopping, backend adaptation, and scientific next steps. No N=10:20 study or GPU speedup is claimed yet.

To rebuild the PDF after editing its LaTeX sources, run `docs/build_pdf.ps1` with MiKTeX/LaTeX installed. The build requires geometry, amsmath, amsfonts, booktabs, tools, xcolor, listings, and hyperref packages.
