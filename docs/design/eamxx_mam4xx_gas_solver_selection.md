# EAMxx integration of MAM4xx gas-chemistry solver selection

Status: approved by the user on 2026-10-15; implemented with focused native CPU
validation. Full EAMxx integration remains pending CI.

## Scope and evidence

Provide the EAMxx companion to [MAM4xx PR 599](https://github.com/E3SM-Project/mam4xx/pull/599),
reviewed at `dc171e8b538a8ba3a5b5cc5380e685ec7bec1f01`.
Use E3SM `master` after [PR 8749](https://github.com/E3SM-Project/E3SM/pull/8749),
which introduced the returned solver status consumed by the current EAMxx caller.
Kai's final approval asks whether a corresponding EAMxx PR is needed.

Pin `externals/mam4xx` to the reviewed PR revision. Retain the upstream
`MAM4XX_GAS_CHEMISTRY_SOLVER=legacy` default and allow `direct` through the
existing EAMxx CMake option mechanism. Do not add an EAMxx solver, duplicate
rate construction, or change the MAM4xx PR.

## Equations, units, and state

The existing reduced mechanism, reaction coefficients, prescribed oxidants,
external forcing, species mapping, and work-array conversions remain in MAM4xx.
For fixed coefficients, both solver choices implement backward Euler:

\[
(I-\Delta t A)q^{n+1}=q^n+\Delta t b.
\]

Here `dt` is in seconds, diagonal and transfer coefficients in `A` are in
inverse seconds, and `b` has each corresponding work-entry unit per second.
The direct scalar update is

\[
q_i^{n+1}=\frac{q_i^n+\Delta t s_i(q^{n+1})}{1-\Delta t a_i}.
\]

The coupled DMS to SO2 to H2SO4 entries are solved in dependency order by
MAM4xx. This remains a backward-Euler discretization, with its existing
truncation error; it is not exact continuous-time integration.

The 31-entry work array contains six gas entries followed by 25 aerosol-mass
and modal-number entries. Gas entries use molar mixing ratios; the remaining
entries retain the host's existing converted work-array basis. O3, entry 0,
is outside the 30-entry implicit class and is unchanged by this solve.
No EAMxx conversion, limiter, boundary treatment, coupling cadence, process
ordering, or tolerance is added or modified.

## Public call and failure handling

Call the reviewed MAM4xx public function directly:

```cpp
template <typename VectorType>
KOKKOS_INLINE_FUNCTION void mam4::microphysics::gas_phase_chemistry(
    Real temperature, Real dt,
    const Real photo_rates[mam4::mo_photo::phtcnt],
    const Real external_forcing[mam4::gas_chemistry::extcnt],
    const Real invariants[mam4::gas_chemistry::nfs],
    const Real heterogeneous_rates[mam4::gas_chemistry::gas_pcnst],
    VectorType& state);
```

| Argument | Role, units, and shape |
| --- | --- |
| `temperature` | Input temperature, K; existing host value. |
| `dt` | Input requested integration interval, finite and positive, s. |
| `photo_rates` | Input photolysis frequencies, inverse seconds; size `phtcnt`. |
| `external_forcing` | Input unnormalized sources; gases in molecules per cubic centimetre per second; size 9. MAM4xx performs the existing air-density normalization. |
| `invariants` | Input prescribed number densities, molecules per cubic centimetre; size 8, ordered M, N2, O2, H2O, OH, NO3, HO2, O3. |
| `heterogeneous_rates` | Input first-order removal coefficients, inverse seconds; size 31. |
| `state` | Input/output chemistry work array, size 31; existing host conversions. |

The signature returns `void` in PR 599. Successful return means that the
requested step completed. On an unsuccessful result, its wrapper calls
`Kokkos::abort` with a reason-specific message in both debug and release builds.
The direct solver retains its transactional state handling; the legacy solver
retains only accepted internal substeps. EAMxx must not replace the abort with
a warning or assume an incomplete step is complete.

Compared with current E3SM master, this adopts PR 599's stricter failure policy:
current EAMxx emits a debug-only warning and continues on an incomplete solve.
This behavior difference must be explicit in the companion PR.

## Diagnostic compatibility

The reviewed public wrapper no longer exposes the status/accounting record.
Remove its consumers and these six optional output fields:

- `mam4_imp_sol_outcome`
- `mam4_imp_sol_failed_attempts`
- `mam4_imp_sol_cut_count`
- `mam4_imp_sol_accepted_steps`
- `mam4_imp_sol_requested_interval`
- `mam4_imp_sol_accepted_interval`

Remove the matching registrations and descriptions. Retain the chemistry
tendency diagnostic and all other microphysics diagnostics. Do not fabricate
legacy retry counts from successful `void` returns. Document that output YAML
files requesting the removed fields must remove those entries. Retaining the
six detailed diagnostics would require a separate MAM4xx API change or an
EAMxx wrapper duplicating MAM4xx's chemistry setup; neither is in this proposal.

## Files and configuration

Implementation files:

- `externals/mam4xx`: pin the reviewed revision.
- `components/eamxx/src/physics/mam/eamxx_mam_microphysics_process_functions.cpp`:
  use the void call; remove obsolete status handling and statistics.
- `components/eamxx/src/physics/mam/eamxx_mam_microphysics_process_interface.cpp`:
  remove unsupported output-field registrations and descriptions.
- `components/eamxx/src/physics/mam/tests/CMakeLists.txt` and
  `mam_gas_chemistry_test.cpp`:
  exercise the configured public wrapper through an EAMxx/Kokkos view.
- `components/eamxx/docs/user/`: document solver selection and compatibility.
- This design document: record approved and executed state before completion.

Use `-DMAM4XX_GAS_CHEMISTRY_SOLVER=direct` for standalone configuration, or
append `MAM4XX_GAS_CHEMISTRY_SOLVER direct` to `SCREAM_CMAKE_OPTIONS` in a CIME
case. Select `legacy` explicitly for a matched comparison. Rebuild after
changing the choice; this is not a runtime namelist setting. Let upstream
MAM4xx validate the option and generate exactly one solver macro.

## STRICT validation and acceptance

1. Confirm the E3SM submodule pin exactly matches the reviewed MAM4xx commit
   and the patch contains no unrelated changes.
2. Configure default, explicit `legacy`, and explicit `direct` selections;
   confirm generated solver macros and reject an invalid choice.
3. Compile and run focused wrapper checks for both selections and precisions
   where local dependencies permit. Use a device-compatible view, finite
   nonnegative fixtures, at least three consecutive steps, and unchanged O3.
4. Check a zero-rate/zero-source reference where applicable and DMS loss sign;
   retain upstream solver tests as evidence for conservation, positivity,
   transactional failure handling, and the backward-Euler algebra.
5. Check that the removed diagnostics have no remaining registrations,
   consumers, or documentation promising their availability.
6. Run `git diff --check`; review the patch, address findings, and perform a
   second brief review. Report exact checks and counts, without claiming a
   complete EAMxx simulation or GPU validation from local wrapper tests.

Full EAMxx/CIME build, restart, GPU, and atmospheric comparison remain CI or
separately authorized supported-machine checks. No remote HPC work or job
submission is authorized by this integration task.

## Executed evidence

Local component checks used Apple Clang 21, C++17, Release (`-DNDEBUG`), and
the real Serial Kokkos/EKAT dependency binaries from an existing build.
EKAT source identity was `007bdcff6cb26948517d2d26f3d7ca005363d119`, the exact
revision requested by MAM4xx. The dependency build was reused read-only;
MAM4xx and the tests were rebuilt in short unsynced scratch directories.
This harness exercised the actual MAM4xx CMake embedding path, generated
configuration, and library rather than substituting solver headers or runtime.

| Solver | Precision | Test cases | Assertions | Result |
| --- | --- | ---: | ---: | --- |
| legacy | double | 7 | 5,166 | Pass |
| legacy | single | 7 | 5,166 | Pass |
| direct | double | 8 | 5,170 | Pass |
| direct | single | 8 | 5,170 | Pass |

Each configuration included the new EAMxx view regression, the upstream direct
solver unit suite, and the upstream multi-step wrapper smoke test. The new
regression checks all 31 entries over three consecutive steps in two columns
and three levels, through a const strided Kokkos view. It checks nonnegative,
finite state, unchanged O3, zero-rate preservation within roundoff, and the
independent scalar backward-Euler removal reference. The upstream suite covers
the coupled system, discrete balance, signed forcing, and failure rollback.

Additional executed checks:

- Four invalid-interval subprocesses aborted with the expected wrapper error
  in Release builds, one per solver/precision combination.
- A fresh configure without a solver flag selected only the legacy macro.
- An invalid solver option failed configuration with the expected error.
- Source inspection found no remaining removed-field registrations or
  consumers. The migration documentation intentionally names the six fields.
- `git diff --check` passed.

Local commands and logs:

```sh
python3 /tmp/eamxx-mam599-validation/run_checks.py
# The driver configures each solver/precision pair with:
cmake -S /tmp/eamxx-mam599-validation -B /tmp/eamxx-mam599-checks/legacy-double \
  -DCMAKE_BUILD_TYPE=Release -DMAM4XX_PRECISION=double \
  -DMAM4XX_GAS_CHEMISTRY_SOLVER=legacy
cmake --build /tmp/eamxx-mam599-checks/legacy-double --parallel 2
/tmp/eamxx-mam599-checks/legacy-double/integration_tests
```

The driver repeats these commands for the four combinations and records
`configure.log`, `build.log`, `tests.log`, `failure.log`, and a final
`results.json` under `/tmp/eamxx-mam599-checks`. This component harness did
not build the complete EAMxx process or its MPI/Fortran dependencies. In a
configured EAMxx build, the committed regression is registered as
`mam_gas_chemistry_test` and can be selected with
`ctest -R mam_gas_chemistry_test --output-on-failure`.

An initial zero-rate fixture demanded bitwise identity for all entries and
exposed roundoff in the legacy solver. The test now uses a 64-epsilon relative
tolerance for evolved entries, while O3 still requires exact equality. No
production-solver arithmetic was changed to satisfy that test.

## Limitations and open questions

- PR 599 is still open; refresh the submodule pin if its reviewed head changes
  or its merge commit is required before the E3SM companion can merge.
- The six output-field removals are an explicit compatibility cost.
- Full atmosphere behavior and performance require paired integration runs;
  standalone results do not establish atmospheric accuracy or speedup.
- The legacy work-array unit conventions are inherited, not reinterpreted.

## Revision History

| Date | Author/agent | Summary |
| --- | --- | --- |
| 2026-10-15 | Codex | Propose the matched EAMxx integration, retaining the legacy default and adopting the reviewed void API, fail-fast behavior, and diagnostic compatibility changes; awaiting user approval. |
| 2026-10-15 | Codex | Record independent plan review with no critical changes and explicit user approval to implement the proposed integration. |
| 2026-10-15 | Codex | Implement the reviewed dependency pin and public void call, remove unsupported output fields, add configuration/migration documentation and the strided-view regression, and record passing local four-configuration checks with full EAMxx/MPI/GPU/restart coverage still pending. |
