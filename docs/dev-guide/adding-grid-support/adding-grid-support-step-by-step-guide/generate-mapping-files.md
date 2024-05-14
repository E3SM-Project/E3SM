# Generate mapping files

In order to pass data between different components at runtime, a set of mapping files between each component is generated offline.

See [Recommended Mapping Procedures for E3SM Atmosphere Grids](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/178848194/Recommended+Mapping+Procedures+for+E3SM+Atmosphere+Grids) for a discussion of different remap algorithms and when to use each.

TempestRemap and ESMF are the backends that generate the mapping weights, but this is all nicely encapsulated using ncremap. Tempest is the preferred method for creating mapping files. ncremap will call TempestRemap or ESMF depending on the algorithm argument and input file types. If exodus files are provided (i.e. `*.g`) then TempestRemap commands will be used. The ESMF tools are adequate for making atmosphere-only-type component sets for E3SM, but this tool is less conservative than TempestRemap. If you are making grids for a coupled run, then TempestRemap should be used wherever possible. Currently, TempestRemap has trouble with masked grids, such those that are needed for land data generation, so ESMF is still required for certain tasks.
