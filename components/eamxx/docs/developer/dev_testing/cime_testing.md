# Full Model Testing

Full model system testing of eamxx is done through CIME test cases
(much like the rest of E3SM).

We offer a number of test suites, including:

* e3sm_scream_v0: Test the full set of V0 (pre-C++) tests
* e3sm_scream_v1: Test the full set of V1 (C++) tests
* e3sm_scream_v1_at: A smaller and quicker set of tests for autotesting
* e3sm_scream_hires: A small number of bigger, longer-running tests to measure performance

Example for running a suite:

```shell
cd $repo/cime/scripts
./create_test e3sm_scream_v1_at --wait
```

Example for running a single test case:

```shell
cd $repo/cime/scripts
./create_test SMS.ne4_ne4.F2010-SCREAMv1 --wait
```

There are many behavioral tweaks you can make to a test case, like
changing the run length, test type, etc. Most of this is not specific
to eamxx and works for any CIME case. This generic stuff
is well-documentated here:
<http://esmci.github.io/cime/versions/master/html/users_guide/testing.html>

When it comes to things specific to eamxx, you have grids, compsets, and
testmods.

Common EAMxx grids are:

* ne4_ne4 (low resolution)
* ne4pg2_ne4pg2 (low resolution with phys grid)
* ne30_ne30 (med resolution)
* ne30pg2_ne30pg2 (med resolution with phys grid)
* ne1024pg2_ne1024pg2 (ultra high with phys grid)

More grid info can be found here:
<https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/933986549/ATM+Grid+Resolution+Summary>

Common EAMxx compsets are:

* F2010-SCREAM-LR: V0 low res compset with eamxx V0 atmosphere
* F2010-SCREAMv1: V1 standard compset with eamxx V1 atmosphere
* FIOP-SCREAMv1-DP: V1 with dpxx (doubly-periodic lateral boundary condition in C++)
* F2010-SCREAMv1-noAero: V1 without aerosol forcing

Full info on supported compsets can be found by looking at this file:
`$scream_repo/components/eamxx/cime_config/config_compsets.xml`

Common EAMxx testmods are:

* small_kernels: Enable smaller-granularity kernels,
  can improve performance on some systems
* scream-output-preset-[1-6]: Our 6 output presets.
  These turn some combination of our three output streams
  (phys_dyn, phys, and diags),
  various remaps, etc.
* bfbhash: Turns on bit-for-bit hash output: <https://acme-climate.atlassian.net/wiki/spaces/NGDNA/pages/3831923056/EAMxx+BFB+hashing>

More info on running EAMxx can be found here:
<https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/3386015745/How+To+Run+EAMxx+SCREAMv1>
