# Continuous Integration and Nightly Testing

## Autotester

EAMxx using GitHub actions and a Sandia product called Autotester 2 (AT2)
to run CI testing on a CPU and GPU machine for every GitHub pull
request.
By default, we run the `e3sm_scream_v1_at` suite and the
standalone EAMxx tests (`test-all-scream`).

## Nightly overview, CDash

Our nightly testing is much more extensive than the CI testing. You
can see our dashboard here under the section [E3SM_EAMXX](https://my.cdash.org/index.php?project=E3SM)

We run a variety of CIME test suites and standalone testing on a number
of platforms.
We even do some performance testing on [Frontier](https://www.olcf.ornl.gov/frontier/).
