# Precision studies of the lattice Schwinger model

**Numerical quantum field theory · Julia / ITensors · DMRG and matrix product states**

## 1. Problem

How precisely can a lattice calculation locate the critical point of the massive Schwinger model? Near a quantum phase transition, finite system size and numerical truncation can obscure the physical answer. This repository contains spectrum calculations with open and periodic boundary conditions, providing computational examples related to a broader precision study.

## 2. Result

The related coauthored paper reports a critical-mass determination to five digits and agreement among four independent criticality criteria.

| Published study | Reported result |
| --- | --- |
| System size | Up to 3,000 qubits in a classical tensor-network simulation |
| Critical mass | Determined to five digits |
| Cross-checks | Four criticality criteria in agreement |
| Extension | Critical mass studied under a four-fermion deformation |

Source: E. Arguello Cruz, G. Tarnopolsky, and Y. Xin, [*Precision study of the massive Schwinger model near quantum criticality*](https://arxiv.org/abs/2412.01902).

The table describes the published study. The scripts here are compact computational examples, not an archived reproduction of every large-system run or figure in the paper.

## Additional GPU implementation (validation pending)

The separate [GPU implementation](GPU%20implementation/README.md) folder contains a shared CPU solver, an NVIDIA CUDA backend, comparison scripts, unchanged local reference sources, and a [detailed PDF guide](GPU%20implementation/docs/Schwinger_GPU_guide.pdf). **The CUDA path has not yet been run or validated on an NVIDIA GPU.** Historical tiny-system CPU checks do not establish GPU correctness or performance; the N=10 through 20 CPU/GPU study remains pending. See the [validation record and checklist](GPU%20implementation/VALIDATION.md). The original Julia scripts and their environment remain unchanged.

## 3. Reproduce

Install Julia 1.12.4 and Git, then instantiate the recorded project-local environment:

```bash
git clone https://github.com/arerickcr/Schwinger-Model.git
cd Schwinger-Model
julia --project=. -e 'using Pkg; Pkg.instantiate(; allow_autoprecomp=false)'
julia --project=. --startup-file=no --compiled-modules=no Schwinger_model_with_OBC.jl
```

The supplied open-boundary example uses **8 sites**, `x = 10`, `theta = 0`, five sweeps, and maximum bond dimension 500. It prints DMRG progress, overlaps between computed states, and four energies labeled `E0`, `E1`, `E2`, and `E3`. The recorded run gave `E0 ≈ -45.5374129880` and `E1 ≈ -37.4066023569`. Later excited states need further convergence checks; see the [recorded output and limitations](docs/reproduction-status.md).

The command disables compiled modules because package precompilation failed in the review environment. Package loading and compilation can be slow; the flag can be omitted in an environment where precompilation succeeds.

| File | Purpose |
| --- | --- |
| [Open-boundary calculation](Schwinger_model_with_OBC.jl) | Small spectrum example and overlap diagnostics |
| [Periodic-boundary calculation](Schwinger_model_with_PBC.jl) | Spectrum with an additional electric-flux degree of freedom |
| [Custom site operators](S5%28modified%29.jl) | Truncated flux-site definitions used in the periodic formulation |

The periodic script is a research reference and uses legacy ITensors imports. Its environment and custom-site integration require separate checking; the open-boundary example is the intended starting point. Record Julia and package versions with `julia --version` and `julia --project=. -e 'using Pkg; Pkg.status()'`.

## 4. Validation

The original open-boundary script prints pairwise overlaps of its computed states and symmetry-related diagnostics. Numerical acceptance should additionally include stability with sweep count and bond dimension, low-energy agreement with exact diagonalization at small sizes, and checks of finite-size extrapolations.

The paper's agreement among four criticality criteria is evidence for the published result; it does not substitute for a convergence study of a new parameter choice. Small solver tolerances alone do not guarantee physical accuracy. A complete automated regression suite and reproduction of the paper's largest runs are outside the current repository scope.

See [reproduction status](docs/reproduction-status.md) for the environment check performed during documentation preparation.

## 5. Methods

The Schwinger model is quantum electrodynamics in one spatial dimension. Staggered fermions and a Jordan–Wigner transformation map the lattice problem to spins. DMRG approximates low-energy states as matrix product states; the electric-flux basis for periodic boundaries also requires truncation.

The dimensionless parameters used in the notes include

$$x=\frac{1}{a^2e^2},\qquad \mu=\frac{2m_{\mathrm{latt}}}{ae^2},\qquad m_{\mathrm{latt}}=m-\frac{e^2a}{8}.$$

The detailed Hamiltonians, boundary conventions, and bibliography are retained in [theory and references](docs/theory.md).

[Paper](https://arxiv.org/abs/2412.01902) · [Open-boundary notes](https://earguelloc.com/assets/pdfs/Schwinger_Model_Notes.pdf) · [Periodic-boundary notes](https://earguelloc.com/assets/pdfs/Schwinger_Model_in_a_Circle.pdf) · [ITensorMPS documentation](https://docs.itensor.org/ITensorMPS/stable/)
