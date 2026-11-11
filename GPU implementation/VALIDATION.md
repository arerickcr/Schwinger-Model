# Validation status

**GPU validation: PENDING.** Publication date: November 11, 2026. This is an experimental implementation; the results below validate only the CPU path. No CUDA execution, GPU comparison, GPU timing, or N=10 through 20 study has been completed.

## Historical CPU validation - October 22, 2026

All **115 checks passed**, with no failures or errors:

| Check group | Passed |
|---|---:|
| Initial states and configuration | 8 |
| Full Hamiltonian vs independent spin matrix and unchanged original | 60 |
| CPU DMRG, exact energies, observables, and original solver | 40 |
| Excited states and zero-energy convergence | 5 |
| Generic backend interface and reports | 2 |

| N | New CPU ground energy | Independent exact energy | Absolute difference |
|---|---:|---:|---:|
| 4 | -21.85462741028993 | -21.854627410289932 | 3.55e-15 |
| 5 | -26.7175399777961 | -26.717539977796104 | 3.55e-15 |
| 6 | -33.965949729130365 | -33.96594972913035 | 1.42e-14 |

Physics defaults: x=10, y=1, mass=0.33, theta=pi, g=0. Initial spin is 0 for even N and +0.5 for odd N. The exact-matrix comparisons additionally covered N=2 through 6, theta=0,pi,0.73, g=0.17, y=0.8, mass=0.21.

The original helper file was loaded unchanged inside a module. Its N=4 solver returned -21.854627410289932. Its N=20 driver was not run. SHA-256 values for both original files match `original_files.sha256`.

The main tiny DMRG checks used nsweeps=16, maxdim=32, noise=0, eig_maxiter=8, and Float64. These are correctness checks with generous small-system accuracy, not the original ten-sweep benchmark. The original solver retains its own internal settings when invoked. The N=4 two-state test uses nsweeps=20, maxdim=16, eig_maxiter=10, noise=0. Per-run TOML records contain the exact settings.

Environment: Julia 1.12.4, ITensors 0.9.31, ITensorMPS 0.4.1, Windows, four BLAS threads, one Julia thread. The manifest is synchronized with the new CPU project's declared dependencies. Seven top-level scripts pass `tests/check_syntax.jl`.

**Not run:** NVIDIA execution, CPU/GPU agreement, N=10–20 benchmarking. Windows reported Intel UHD Graphics 620; no NVIDIA GPU was available. First-run CPU timings contain compilation costs and are not performance results.

Full transcript: `results/validation_log.txt`. Run `julia --project=. tests/runtests.jl` from inside `GPU implementation` to reproduce. On the NVIDIA machine, run `setup_cuda.jl` and then `check_cuda.jl` before benchmarking.


## Pending NVIDIA validation - fill in after execution

Leave the following items unchecked until actual NVIDIA runs and their saved reports support them. CPU tests and the `host_adapter` reports do not validate GPU execution.

- [ ] Record GPU model, VRAM, driver, OS, Julia, CUDA.jl, ITensors, and ITensorMPS versions, plus the tested commit and post-setup manifest.
- [ ] Run `julia --project=. setup_cuda.jl` on a compatible NVIDIA machine and retain the resolved environment.
- [ ] Run `julia --project=. check_cuda.jl`; verify Float64 device storage and N=4 CPU/GPU agreement, with no silent CPU fallback.
- [ ] Check nonzero quartic coupling, background angles, odd/even sizes, and the excited-state path on GPU before interpreting those modes.
- [ ] Run `julia --project=. benchmark.jl --sizes=10:20 --repeats=3 --warmup=true` and retain the reports.
- [ ] Check energy differences, fidelity, residuals, orthogonality, measured charge, truncation errors, and convergence on both backends.
- [ ] Establish stability with sweeps, bond dimension, and local eigensolver accuracy before accepting physical results.
- [ ] Report synchronized timings with warmup excluded; state the timed region and make speedup claims only from measured, converged comparisons.

| Future evidence | Status / value |
|---|---|
| Test date and Git commit | Pending |
| GPU model and VRAM | Pending |
| Driver, OS, and software environment | Pending |
| N=4 CUDA smoke-test transcript | Pending |
| N=10 through 20 comparison reports | Pending |
| Convergence study | Pending |
| Timing results and limitations | Pending |

The committed `results/` files are historical CPU evidence, not GPU results. The tests reuse their filenames: archive existing reports or run in a separate checkout when preserving multiple validation sessions. Add future GPU logs, resolved dependencies, and measured results in a follow-up commit; update this status only for checks actually completed.
