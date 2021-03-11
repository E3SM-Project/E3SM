# SCREAM: the Simple Cloud Resolving E3SM Atmosphere Model

SCREAM is a next-gen atmosphere component for the [E3SM project](https://e3sm.org/).
In contrast to previous atmosphere components, which model phenomena on grids
with coarser spatial scales using sub-grid-scale parameterizations, SCREAM
focuses on resolving cloud physics by solving well-established differential
equations using state-of-the-art numerical methods.

E3SM project members can find more information on SCREAM's [Confluence Site](https://acme-climate.atlassian.net/wiki/spaces/NGDNA/overview).
Plenty of documentation is available on the [Documentation page](https://acme-climate.atlassian.net/wiki/spaces/NGDNA/pages/755597313/Documentation),
visible in the left side bar.

SCREAM uses process implementations from legacy Fortran modules from E3SM. You
can find these modules `components/cam/src`. These legacy implementations use
the CAM/E3SMv1 process coupling approach, and are used mostly for comparisons
with the newer C++ implementations. The C++ implementations are exposed via
SCREAM's atmosphere driver.

For information about building SCREAM and running its tests, see
[this page](docs/build.md).

You can also take a brief [tour of SCREAM's source tree](docs/source-tree.md).

SCREAM comes with a python-based set of tools that compliments the CIME tool
suite and also works with standalone SCREAM. See scripts/README.md for more.
