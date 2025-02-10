# Multi-instance and NBFB testing capabilities

We adopt MVK testing for now.

## MVK testing

In Multivariate Kolmogorov-Smirnov (MVK) test, an ensemble of multi-year runs
with the 2 versions of the code to check that the climate statistics of
important variables are unchanged between ensembles.

For some background, see DOIs:

    - 10.1145/3324989.3325724
    - 10.1145/3468267.3470572

## Multi-instance testing

See following link for background on the multi-instance capability in CIME.

[Link][cime-multi-instance]

[cime-multi-instance]: https://esmci.github.io/cime/versions/master/html/users_guide/multi-instance.html

In SCREAM/EAMxx, it works in exactly the same except that
the namelist file has a different name and location.
Instead of user_nl_eam and atm_in for EAM, in SCREAM,
we do not use user_nl_scream at all (even though it will be there)
and we rely on run/data/scream_input.yaml.

Thus, much like the original implementation in other components,
the user must specify their choice of settings in the corresponding
input files, namely scream_input.yaml_0001, and so on.

## Together

For MVK testing to work, it relies on the multi-instance capability.
We use the default MVK test to instantiate a new MVKxx test in
components/eamxx/cime_config/SystemTests. In the test definition,
we must edit the scream_input.yaml specs.
