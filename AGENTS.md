<!-- markdownlint-disable -->
# AGENTS.md

E3SM is a coupled Earth system model, configured, built and run through the
CIME case control system (`cime/`, a git submodule). Docs: https://docs.e3sm.org/E3SM/

## Constraints

- Cases and CIME tests only work on the machines in
  `cime_config/machines/config_machines.xml`. Do not create a case or test
  anywhere else.
- Exception: EAMxx has a standalone build that runs anywhere. Read
  `components/eamxx/AGENTS.md` before testing changes under `components/eamxx`.
- Most of E3SM has no unit tests. Test a change by building and running one
  case that uses the modified component.
- `cime/` and `externals/` are empty until
  `git submodule update --init --recursive --depth=1`.

## Testing a change with a case

Ask the user for a CASENAME. Create the case in the repository top level unless
the name includes a path. Builds are out of source.

```bash
cd cime/scripts
./create_newcase --case ../../CASENAME --compset COMPSET --res RESOLUTION --mach MACHINE
cd ../../CASENAME
./case.setup && ./case.build && ./case.submit
```

- Rebuild with `./case.build` (do not clean first). `./case.build -m` is faster
  when only component or driver code changed.
- After changing the PE layout: `./case.setup --reset`, then `./case.build`.
- `./xmlquery EXEROOT` and `./xmlquery RUNDIR` give the build and run directories.

Smallest configuration for each modified component:

| Component | Compset | Resolution |
|-----------|---------|------------|
| elm | I1850CNPRDCTCBCTOP | ne4pg2_ne4pg2 |
| eam, or elm and eam together | F1850 | ne4pg2_oQU480 |
| mosart | RMOSGPCC | r05_r05 |
| eamxx | F2010-SCREAMv1 | ne4pg2_ne4pg2 |
| mpas-ocean | CMPASO-NYF | T62_oQU120 |
| mpas-seaice | DTESTM | T62_oQU240 |
| any other mix of two or more | WCYCL1850NS | ne4pg2_r05_oQU480 |

## CIME system tests

`cd cime/scripts && ./create_test NAME`, where NAME is a single test
(`TYPE.RESOLUTION.COMPSET[.MACHINE_COMPILER][.TESTMOD]`, for example
`SMS_D_Ln5_P4.ne4pg2_oQU480.F2010.pm-cpu_gnu`) or a suite from
`cime_config/tests.py`. Follow progress in `TestStatus.log` in the test
directory. Nightly results: https://my.cdash.org/index.php?project=E3SM

## Branches and pull requests

- Branch names: `<github username>/<component or code area>/<feature-description>`.
- Topic branches merge to `next` for testing, then to `master` (gitworkflows
  without `seen`).
- A PR description has two parts separated by a line: first a brief
  description in plain text with no markdown, then the full description in
  markdown. Write it in the imperative, starting with a verb such as "Fix" or
  "Add".
