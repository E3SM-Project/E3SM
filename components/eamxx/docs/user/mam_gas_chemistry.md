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

## Failures and output compatibility

If either solver cannot complete a chemistry step, the model stops with an
error message in both debug and release builds.

These fields are no longer available. Remove them from any output YAML that
requests them:

- `mam4_imp_sol_outcome`
- `mam4_imp_sol_failed_attempts`
- `mam4_imp_sol_cut_count`
- `mam4_imp_sol_accepted_steps`
- `mam4_imp_sol_requested_interval`
- `mam4_imp_sol_accepted_interval`
