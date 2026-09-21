# Reproduction status

The existing `Schwinger_model_with_OBC.jl` was executed unchanged during portfolio preparation on Windows with Julia 1.12.4, ITensors 0.9.31, ITensorMPS 0.4.1, and HDF5 0.17.4. The resolved `Project.toml` and `Manifest.toml` are included at the repository root. That environment also includes Plots 1.41.7 for the companion tensor-network examples.

## Command

After instantiating the supplied environment:

```bash
julia --project=. --startup-file=no --compiled-modules=no Schwinger_model_with_OBC.jl
```

Package precompilation failed in the review environment. Disabling compiled modules allowed the original example to complete. This is a workaround for loading packages, not a change to the numerical algorithm. Startup and compilation time can be substantial; no runtime performance claim is made.

## Recorded output

The default eight-site example printed:

```text
E0 = -45.537412988049724
E1 = -37.406602356915606
E2 = -31.45619803247753
E3 = -30.23737531817939
```

These are outputs from the supplied settings, not independently verified exact energies. The last ground-state sweep changed the energy by approximately 10^-14, and the first excited state was also stable at that scale in its final sweeps. In contrast, the second excited-state estimate changed by about 0.167 in its final sweep. Several overlaps involving the later states were of order 10^-5.

**Interpretation:** the script executes and produces the documented observables, but the supplied sweep count does not establish convergence of every excited state. Do not use the later energies as precision reference values. No independent exact-diagonalization comparison, residual-based validation, or large-system reproduction was performed in this documentation check.

The periodic-boundary script and custom flux-site integration were not executed. They retain their original library imports and require a separate compatibility review.
