# Adding Support for New Grids

The purpose of this guide is to outline all the necessary steps for running E3SM on a new grid for the atmosphere and land components. The process is similar for uniform and regionally refined grids, although regionally refined cases will likely require some special considerations which will be noted where appropriate.

If you wish to add a new ocean and sea-ice mesh you will need to use the compass tool to generate the mesh and dynamically adjusted initial condition. This procedure is detailed in a separate tutorial:
<https://mpas-dev.github.io/compass/latest/tutorials/dev_add_rrm.html>

<!-- disable certain linter checks here to allow vertical alignment of links -->
<!-- markdownlint-disable MD039 --> <!-- no-space-in-links -->
<!-- markdownlint-disable MD042 --> <!-- no-empty-links -->

## Step-by-Step Guide

1. [Generate a new grid file                               ](adding-grid-support-step-by-step-guide/generate-new-grid-file.md)
1. [Generate mapping files                                 ](adding-grid-support-step-by-step-guide/generate-mapping-files.md)
1. [Generate domain files                                  ](../../generate_domain_files/index.md)
1. [Generate a topography file                             ](adding-grid-support-step-by-step-guide/generate-topo-file.md)
1. [Generate atmospheric initial condition                 ](adding-grid-support-step-by-step-guide/generate-atm-initial-condition.md)
1. [Generate land input data (*fsurdat*)                   ](adding-grid-support-step-by-step-guide/generate-lnd-input-data.md)
1. [Generate land initial condition (*finidat*)            ](adding-grid-support-step-by-step-guide/generate-lnd-initial-condition.md)
1. [Generate a dry deposition file (*depends on use case*) ](adding-grid-support-step-by-step-guide/generate-dry-deposition.md)
1. [Add new grid configuration to E3SM                     ](adding-grid-support-step-by-step-guide/add-grid-config.md)

## Other Useful Tutorials

1. [Generate a RRM Grid File with SQuadGen                 ](adding-grid-support-step-by-step-guide/generate-RRM-grid-file.md)
