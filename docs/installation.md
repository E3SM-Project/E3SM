# Installation Guide

E3SM is not available as a pre-compiled binary.  You install
E3SM by cloning the source code using [git](https://git-scm.com/)
and building the executable on your local platform after making some
choices on model configuration (addressed in the [User Guide](user-guide/index.md)).

It is recommended that you install E3SM on a supported platform.
(See [Hardware/Software Configuration](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/4116447351/Hardware+Software+Configuration))

The E3SM Project can not assist with installation on an un-supported platform but can point you in the right direction.
If you wish to port E3SM to your machine, first check that the target machine has
the software prerequisites detailed in
[Hardware/Software Configuration](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/4116447351/Hardware+Software+Configuration#Software-prerequisites).

E3SM uses the CIME Case Control System (CCS) and a machine must be described to the CCS
in order to build and run E3SM.  See the
[CIME Porting Guide](https://esmci.github.io/cime/versions/master/html/users_guide/porting-cime.html).

Once you are on a supported machine, clone the source code by following the first steps (2-4 and 6)
in the [Development Getting Started Guide](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/1868455/Development+Getting+Started+Guide).

To start configuration cases and building the model, see the User Guide
