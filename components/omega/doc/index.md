(omega-home)=
# OMEGA

The Ocean Model for E3SM Global Applications (OMEGA) initially is an eddy-resolving,
global ocean model in the early stages of development by the
[E3SM](https://e3sm.org/) ocean team.  The first release is planned for
Summer 2026.  A non-eddying configuration will be released in early 2027.

The model is written in c++ using the [Kokkos](https://github.com/kokkos)
framework for performance portability.  OMEGA is based on the
[TRSK formulation](https://doi.org/10.1016/j.jcp.2009.08.006) for geophysical
models on unstructured meshes. The first version of OMEGA will primarily be a direct port
of the [MPAS-Ocean](https://e3sm.org/model/e3sm-model-description/v1-description/v1-ocean-sea-ice-land-ice/)
component of E3SM for comparison purposes.

Development is taking place at https://github.com/E3SM-Project/Omega.


```{toctree}
:caption: User's guide
:maxdepth: 2

userGuide/QuickStart
userGuide/OmegaBuild
userGuide/DataTypes
userGuide/MachEnv
userGuide/Config
userGuide/Broadcast
userGuide/Logging
userGuide/Decomp
userGuide/MetaData
userGuide/IO
userGuide/IOField
userGuide/Halo
userGuide/HorzMesh
userGuide/TimeMgr
```

```{toctree}
:caption: Developer's guide
:maxdepth: 2

devGuide/QuickStart
devGuide/CondaEnv
devGuide/Linting
devGuide/Docs
devGuide/BuildDocs
devGuide/DataTypes
devGuide/MachEnv
devGuide/Config
devGuide/Broadcast
devGuide/CMakeBuild
devGuide/Logging
devGuide/Decomp
devGuide/MetaData
devGuide/IO
devGuide/IOField
devGuide/Halo
devGuide/HorzMesh
devGuide/TimeMgr
```

```{toctree}
:caption: Design documents
:maxdepth: 1

design/Broadcast
design/Config
design/DataTypes
design/Decomp
design/Driver
design/Halo
design/HorzMeshClass
design/Logging
design/MachEnv
design/Metadata
design/IO
design/IOStreams
design/Reductions
design/ShallowWaterOmega0
design/State
design/Tendency
design/TimeMgr
design/Timers
design/Tracers

design/OmegaDesignTemplate
```
