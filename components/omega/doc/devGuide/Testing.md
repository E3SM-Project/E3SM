(omega-dev-testing)=

# Testing Code

For instructions on how to build Omega and run CTest unit tests, see
{ref}`omega-dev-quick-start-build-test`.


## Running the Polaris `omega_pr` Test Suite

The CTest unit tests on their own are not sufficient to ensure that changes
to Omega behave as expected.  Developers also need to run the regression tests
in the form of the
[omega_pr](https://docs.e3sm.org/polaris/main/users_guide/ocean/suites.html#omega-pr-suite)
test suite from Polaris.  The tests in the suite are defined in
[omega_pr.txt](https://github.com/E3SM-Project/polaris/blob/main/polaris/suites/ocean/omega_pr.txt)
in the Polaris repository.

You will need to run `omega_pr` once on the Omega
[develop branch](https://github.com/E3SM-Project/Omega)
and once on your development branch.  This means you will need to
{ref}`build Omega <omega-dev-quick-start-build>` for each branch before running
the suite.  You may wish to use the {ref}`omega-dev-quick-start-ctest-util` to
automate the build and run the CTests for each branch.

```{note}
Run both baseline and PR suite on the same machine, CPU/GPU partition,
compiler, and Omega build type (Debug/Release) to avoid false positive/negative
diffs.
```

### Setting up the suite

If you have not yet set up Polaris, you will need to follow its
[quick start](https://docs.e3sm.org/polaris/main/developers_guide/quick_start.html)
including [cloning the repo](https://docs.e3sm.org/polaris/main/developers_guide/quick_start.html#set-up-a-polaris-repository-for-beginners),
[configuring the software environment](https://docs.e3sm.org/polaris/main/developers_guide/quick_start.html#supported-machines),
and [activating the environment](https://docs.e3sm.org/polaris/main/developers_guide/quick_start.html#activating-the-environment).

Then, with the Polaris environment active, you can set up the baseline
`omega_pr` suite as follows:
``` bash
DEVELOP_BUILD=/path/to/omega/develop/build
POLARIS_SCRATCH=/path/to/scratch/polaris_testing/
polaris suite -c ocean -t omega_pr --model omega \
    -w $POLARIS_SCRATCH/baseline_omega_pr \
    -p $DEVELOP_BUILD
```
where `$POLARIS_SCRATCH` is a directory within your scratch space that is
useful for setting up and running the Polaris test suites and
`$DEVELOP_BUILD` is the directory where you built Omega from the `develop`
branch, or where the Polaris CTest utility built if for you (e.g.
`build_omega/build_chrysalis_intel`).

Next, you can set up the `omega_pr` for your development branch, using the
baseline from the `develop` branch:

``` bash
TEST_BRANCH_BUILD=/path/to/omega/test/branch/build
POLARIS_SCRATCH=/path/to/scratch/polaris_testing/
polaris suite -c ocean -t omega_pr --model omega \
    -b $POLARIS_SCRATCH/baseline_omega_pr \
    -w $POLARIS_SCRATCH/my_branch_omega_pr \
    -p $TEST_BRANCH_BUILD
```
Here, `$POLARIS_SCRATCH` is the same scratch directory as above,
`$TEST_BRANCH_BUILD` is the directory where you build Omega from your
development branch (or where the CTest utility built for you), and
`my_branch_omega_pr` can be any name that helps you keep track of what you
were testing.

### Running the suite

First, run the baseline, e.g. for Slurm systems:
``` bash
cd $POLARIS_SCRATCH/baseline_omega_pr
sbatch job_script.omega_pr.sh
```
Note the job ID for the next step.  Next, run the suite for the test branch,
waiting until the baseline is done:
``` bash
cd $POLARIS_SCRATCH/my_branch_omega_pr
sbatch --kill-on-invalid-dep=yes --dependency=afterok:$JOBID job_script.omega_pr.sh
```
where `$JOBID` is the JOB ID of the baseline slurm job.

Monitor the jobs (e.g. `squeue -j $JOBID`) and take a look at combined
stdout/stderr log in `polaris_omega_pr.o$JOBID` in each work directory to see
whether the tests passed.  Report the results as part of any (non-trivial)
Omega PR.

```{note}
For other schedulers (e.g. PBS on Aurora), take a look at the machine's
documentation for equivalent commands.
```

### Handling test failures

What should you do if you see test failures?

First, if you see test execution failures in the baseline run from `develop`,
these need to be reported as [issues](https://github.com/E3SM-Project/Omega/issues)
on the Omega repo, preferably mentioning a relevant Omega developer who may
be able to investigate the cause of the failure.  Please provide the contents
of the relevant log file for the failing test(s) and optionally the location of
the tests (and log files) if other Omega developers will have access.

If you see execution failures on the tests for your development branch, please
take a look at the log files and attempt to debug the issues yourself first,
then reach out to the Omega team if you need help.

If you see errors related to comparisons with the baseline, you will need to
determine if these are expected based on the code changes or not.  If they
are expected, you will need to document how large the differences are in the
tests that show differences.  In addition, you will need to make a pull
request to the Polaris repository to update the Omega submodule in
`e3sm_submodules/Omega` once your branch has been merged so that the changes
become the new baseline in Polaris.

```{note}
We do not update `e3sm_submodules/Omega` after every branch is merged into
`develop` on Omega to reduce maintenance burden.
```

If you see unexpected differences between the results for your development
branch and those from `develop`, you should attempt to debug the cause yourself
and then reach out to the Omega team if you need help.


### What to include in your PR Testing comment

Please add a PR comment titled "Testing" that summarizes what you ran and the
outcomes. Include the following so reviewers can validate apples-to-apples
runs.

- CTest unit tests
    - Machine(s) and partition/queue (if applicable)
    - Compiler(s), and build type(s) (Debug/Release)
    - Result: either "All tests passed" or a short list of failing tests with
      a one-line note each

- Polaris `omega_pr` regression suite
    - Baseline run (develop)
        - Build directory used for -p
        - Work directory used for -w
    - PR branch run
        - Build directory used for -p
        - Work directory used for -w
    - Machine/partition, compiler, and build type (should match between baseline and PR)
    - Result: "All tests passed" or a concise list of diffs/failures with a brief note
    - If there are issues, add the path to the relevant log(s) (e.g., polaris_omega_pr.o$JOBID) and a short excerpt if helpful

- Tests added/modified/impacted
    - Briefly list any new or updated CTest or Polaris tests and what they cover

- Performance PRs
    - Link to the relevant PACE experiment and summarize before/after metrics

Copy-paste template you can use in your PR comment:

```text
## Testing

CTest unit tests
- Machine(s): <machine>[, <partition>]
- Compiler(s): <compiler>[, ...]
- Build type(s): <Debug|Release>[, ...]
- Result: All tests passed | Failures: <name> (<one-line note>), ...

Polaris omega_pr regression suite
- Baseline build (-p): <path/to/develop/build>
- Baseline workdir (-w): <path/to/baseline_omega_pr>
- PR build (-p): <path/to/pr/build>
- PR workdir (-w): <path/to/my_branch_omega_pr>
- Machine/partition: <machine>[, <partition>]
- Compiler/build type: <compiler>, <Debug|Release>
- Result: All tests passed | Diffs/Failures: <brief summary>
- Logs (if applicable): <path/to/polaris_omega_pr.o$JOBID>

Tests added/modified/impacted
- CTest: <name> — <one-line purpose>
- Polaris: <path or name> — <one-line purpose>

Performance (if applicable)
- PACE: <link>
- Before: <metric(s)>
- After: <metric(s)>
- Delta: <summary>
```
