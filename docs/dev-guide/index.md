# Development Guide

Development of E3SM is carefully planned and tied to goals in the proposals that fund
work on E3SM.

## Code development

Most code development occurs within one component and developers should consult the
relevant component's development guide.

More information on code development.

+ [Development Big Picture](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/7997024/Development+Big+Picture)
+ [Getting Started](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/1868455/Development+Getting+Started+Guide)
+ [Development Reference](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/2523163/Development+Reference)

## Shared spherical Earth radius

The shared ``SHR_CONST_REARTH`` is mapped to PhysicalConstantsDictionary's
GRS80 ``mean_radius`` (approximately 6371008.771366666 m), consistent with Omega.
This replaces the legacy 6371220 m shared radius and affects every consumer of
the shared value, not just the coupler. Other shared physical constants and the
shared rotation convention are unchanged by this radius work-up.

E3SM spherical MPAS-SeaIce and MPAS-Ocean runs enforce the common reference.
Existing ocean meshes and SeaIce restarts at legacy radii require coordinated
input-data preparation; changing radius metadata alone is not sufficient.
See the [SeaIce radius and compatibility discussion](../../components/mpas-seaice/docs/dev-guide/index.md#spherical-earth-radius)
for scaling behavior, standalone exceptions, and required coupled validation.

## Creating new input

E3SM's public input data directory contains initial and boundary conditions for many cases.
Creating new initial or boundary conditions typically starts with creating a new grid.

+ [Adding Support for New Grids](adding-grid-support/index.md)

## Testing

+ [Testing on GitHub-hosted runners](gh-hosted-testing.md)
