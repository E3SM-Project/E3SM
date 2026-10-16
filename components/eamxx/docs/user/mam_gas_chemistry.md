# MAM4xx gas-chemistry solver selection

EAMxx uses MAM4xx's compile-time `MAM4XX_GAS_CHEMISTRY_SOLVER` option.
The default, `legacy`, selects the iterative solver with the retry fix from
[MAM4xx PR 586](https://github.com/E3SM-Project/mam4xx/pull/586).
The opt-in `direct` choice solves the existing backward-Euler system directly,
as introduced in [MAM4xx PR 599](https://github.com/E3SM-Project/mam4xx/pull/599).
Both choices retain the chemical mechanism, rate expressions, chemistry
timestep, and process ordering. The direct choice is not an exact
continuous-time integration.

## Select a solver

For an existing CIME case, select the direct solver before building:

```sh
./xmlchange --append SCREAM_CMAKE_OPTIONS="MAM4XX_GAS_CHEMISTRY_SOLVER direct"
./case.build
```

For a matched legacy comparison, use a separate case and append
`MAM4XX_GAS_CHEMISTRY_SOLVER legacy` instead. Keep the other case settings and
inputs identical. The CIME build script accepts option/value pairs and uses
the last occurrence of an option.

For a standalone EAMxx configuration, add
`-DMAM4XX_GAS_CHEMISTRY_SOLVER=direct` or
`-DMAM4XX_GAS_CHEMISTRY_SOLVER=legacy` to the CMake configure command. The
configure log reports the chosen solver; the CMake cache records it.
Only `legacy` and `direct` are accepted.

This choice requires rebuilding. It is not a runtime namelist setting.
An existing build cache retains its explicit choice when the flag is omitted;
use `legacy` explicitly to return to the default.

## Failure handling and output compatibility

For either choice, the MAM4xx public wrapper aborts with a reason-specific
message if it cannot complete the requested chemistry interval. This applies
to debug and release builds. Earlier EAMxx versions issued a debug-only
warning and continued after an incomplete solve.

The wrapper no longer returns its internal solver-statistics record. Remove
these fields from any output YAML that requests them:

- `mam4_imp_sol_outcome`
- `mam4_imp_sol_failed_attempts`
- `mam4_imp_sol_cut_count`
- `mam4_imp_sol_accepted_steps`
- `mam4_imp_sol_requested_interval`
- `mam4_imp_sol_accepted_interval`

The gas-chemistry tendency field
`mam4_microphysics_tendency_gas_phase_chemistry` remains available when extra
MAM4 microphysics diagnostics are enabled. Other microphysics diagnostics are
unchanged.
