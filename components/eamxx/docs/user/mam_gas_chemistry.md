# MAM4xx gas-chemistry solver selection

The build option `MAM4XX_GAS_CHEMISTRY_SOLVER` selects `legacy` (the default,
iterative solver) or `direct` (a direct backward-Euler solve).

## Select a solver

To use `direct` in an existing CIME case:

```sh
./xmlchange --append SCREAM_CMAKE_OPTIONS="MAM4XX_GAS_CHEMISTRY_SOLVER direct"
./case.build
```

For a standalone EAMxx build, add
`-DMAM4XX_GAS_CHEMISTRY_SOLVER=direct` to the CMake configure command.

Rebuild after changing the selection. To switch back, replace `direct` with
`legacy` in the command or flag above. Specify `legacy` explicitly: an existing
build cache retains its previous selection when the option is omitted.

## Solver failures

The `legacy` solver continues after failure, preserving any accepted internal
substeps. The `direct` solver aborts in debug builds; in release builds it
reports failure, leaves the input state unchanged, and continues the run.

The existing `mam4_imp_sol_*` output fields record solver status and integration
progress when extra MAM4 microphysics diagnostics are enabled. Outcome codes
0-5 retain their legacy meanings; direct failures additionally use
6 (nonfinite result), 7 (unsafe denominator), and 8 (negative result).
