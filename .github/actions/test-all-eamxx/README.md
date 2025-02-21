# Composite action to call test-all-eamxx inside a workflow

This action is meant to be used inside a workflow. E.g.,

```yaml
jobs:
  my-testing:
    steps:
      ...
      - name: run-test-all-eamxx
        uses: ./.github/actions/test-all-eamxx
        with:
          build_type: <build-type>
          machine: <machine>
          run_type: <run-type>
```

The input run_type is the only input that this action has to explicitly handle.
As such, this action checks that its value is one of the following.

- nightly: runs tests and then submit to cdash
- at-run: runs tests without submitting to cdash
- bless: runs tests and copy output files to baselines folder

As for build_type and machine, we do not prescribe a list of
valid values, as that will be handled by components/eamxx/scripts/test-all-eamxx.
If their values are not supported, it is up to test-all-scram to handle the error.
As a guideline, however, you may have to ensure that the following exist:

- the file components/eamxx/cmake/machine-files/${machine}.cmake
- the entry ${machine} in the MACHINE_METADATA dict in components/eamxx/scripts/machine_specs.py

 Questions? Contact Luca Bertagna [lbertag@sandia.gov](mailto:lbertag@sandia.gov)
