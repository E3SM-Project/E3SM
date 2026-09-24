<!-- markdownlint-disable -->
# AGENTS.md

E3SM (Energy Exascale Earth System Model) is a fully coupled Earth System Model
using the CIME case control system. This is a personal AGENTS.md copied into
fresh checkouts along with my own ignore files (`.ignore`, `.opencodeignore`,
`.rgignore`, `.gitignore`) -- not upstream boilerplate. Feel free to propose
edits to this file as we go; it should evolve based on actual usage rather
than staying static.

- **Docs:** https://docs.e3sm.org/E3SM/

## User-specific shell and text preferences

- Use normal ASCII characters for all screen output and for any plain-text
  files created or edited. Avoid Unicode punctuation, symbols, and other
  non-ASCII characters unless the user explicitly requests them.
- The user's interactive shell is `tcsh`. Use `tcsh` syntax for commands
  intended for interactive execution rather than rewriting commands in
  `bash` syntax and requiring adjustments.
- When creating shell scripts, `bash` is acceptable and generally preferred
  unless the script specifically needs to run under `tcsh`.

## Cost/context discipline
- Prefer Grep/Glob over Read for exploring large Fortran/C++ files; only Read
  the specific function/region needed.
- Redirect `case.build`/`case.submit` output to a log file and grep for
  `ERROR`/`FAIL` rather than pasting full build/test logs into context.
- Build is out-of-source, not under the source directory. Only create a case
  or test on a supported machine (`./query_config --machines` in
  `cime/scripts`).

## Daily workflow

Ask the user to specify a CASENAME. Create the case in the top level project
directory unless the case name includes an absolute or relative path.

```bash
cd cime/scripts
./create_newcase --case ../../CASENAME --compset COMPSET --res RESOLUTION --mach MACHINE

cd /path/to/case
./case.setup
./case.build
./case.submit

# Useful lookups
./xmlquery EXEROOT   # build output dir
./xmlquery RUNDIR    # run output dir
```

### Rebuild
```bash
./case.build          # just rebuild, don't clean
./case.build -m       # faster: only components/driver changed

# if pelayout changed:
./case.setup --reset
./case.build
```

### Rerun
```bash
./case.submit
```

### Compsets/grids/machines lookup (on demand, not memorized here)
```bash
cd cime/scripts
./query_config --machines
./query_config --compsets all
./query_config --grids
```
For a created case, the full compset longname is in `README.case` in the case
directory.

### Common compsets/resolutions for testing single-component changes
Use these when a change is isolated to one component:
- elm: COMPSET I1850CNPRDCTCBCTOP, RESOLUTION ne4pg2_ne4pg2
- eam: COMPSET F1850, RESOLUTION ne4pg2_oQU480
- mosart: COMPSET RMOSGPCC, RESOLUTION r05_r05
- eamxx: COMPSET F2010-SCREAMv1, RESOLUTION ne4pg2_ne4pg2
- mpas-ocean: COMPSET CMPASO-NYF, RESOLUTION T62_oQU120
- mpas-seaice: COMPSET DTESTM, RESOLUTION T62_oQU240

If code changed in elm and eam, use F1850. If 2+ other components changed,
use the coupled case: COMPSET WCYCL1850NS, RESOLUTION ne4pg2_r05_oQU480.

## Testing

Most of E3SM does not have unit tests. For component work, create a single
case as above and iterate with `./case.build` / `./case.submit`.

### EAMxx (SCREAM) standalone testing
Code under `components/eamxx` can be tested standalone (no CIME/driver).
See `components/eamxx/AGENTS.md` for details. Quick reference:
```bash
cd components/eamxx
./scripts/test-all-eamxx -m MACHINE                # all tests
./scripts/test-all-eamxx -m MACHINE -t dbg          # one type: sp/dbg/fpe/opt
./scripts/test-all-eamxx --preserve-env -m MACHINE  # keep current env
./scripts/test-all-eamxx -m MACHINE --baseline-dir=LOCAL
```

### CIME system tests
```bash
cd cime/scripts
./create_test TEST_NAME     # single test
./create_test SUITE_NAME    # suite

# Monitor: cd to the test directory, tail TestStatus.log

# Naming: TEST_TYPE.RESOLUTION.COMPSET[.MACHINE][.TESTMOD]
# e.g. SMS_D_Ln5_P4.ne4pg2_oQU480.F2010.pm-cpu_gnu
```
Common test types: SMS (short smoke), SMS_D (debug smoke), ERS (exact
restart), PEM (bit-for-bit across MPI task count), PET (bit-for-bit across
thread count), ERS_Ld (exact restart over d days).

```bash
cd cime/CIME/Tools
./list_e3sm_tests                              # list suites
./list_e3sm_tests -t compsets SUITE_NAME       # compsets in a suite
./list_e3sm_tests -t compsets -l SUITE_NAME    # compset longnames
```
Test suites: `e3sm_land_developer`, `e3sm_mosart_developer`,
`e3sm_land_exeshare` (more via `list_e3sm_tests`).

## Key paths cheat sheet
- `cime_config/machines/config_machines.xml` — supported machines
- `cime_config/machines/cmake_macros/` — machine-specific CMake settings
- `cime_config/tests.py` — test suite definitions
- `cime_config/config_grids.xml` — grid definitions
- `cime_config/config_compilers.xml` — compiler configs
- `cime_config/allactive/config_compsets.xml` + `components/*/cime_config/config_compsets.xml` -- compset definitions

## Submodules
```bash
git submodule update --init --recursive --depth=1
```

## Development workflow
- Feature branch names: `<github username>/<source code area or component>/<feature-description>`
- PR-based workflow with automated testing; coordinate new features via
  E3SM management/science plan.
- Dev guide: https://e3sm.org/model/running-e3sm/developing-e3sm/

## Rules for pull request (PR) description
1. Two parts. First: a brief plain-text description, no markdown. Then a
   blank line. Second: full description with markdown formatting.
2. Use imperative tense, start with a verb like "Fix" or "Add".
