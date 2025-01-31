(omega-home)=
# Omega

The Ocean Model for E3SM Global Applications (Omega) initially is an eddy-resolving,
global ocean model in the early stages of development by the
[E3SM](https://e3sm.org/) ocean team.  The first release is planned for
Summer 2026.  A non-eddying configuration will be released in early 2027.

The model is written in c++ using the [Kokkos](https://github.com/kokkos)
framework for performance portability.  Omega is based on the
[TRSK formulation](https://doi.org/10.1016/j.jcp.2009.08.006) for geophysical
models on unstructured meshes. The first version of Omega will primarily be a direct port
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
userGuide/Driver
userGuide/Decomp
userGuide/Dimension
userGuide/Field
userGuide/IO
userGuide/IOStreams
userGuide/Halo
userGuide/HorzMesh
userGuide/HorzOperators
userGuide/AuxiliaryVariables
userGuide/AuxiliaryState
userGuide/TendencyTerms
userGuide/Tendencies
userGuide/OceanState
userGuide/TimeMgr
userGuide/TimeStepping
userGuide/Reductions
userGuide/Tracers
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
devGuide/Driver
devGuide/Broadcast
devGuide/CMakeBuild
devGuide/Logging
devGuide/Decomp
devGuide/Dimension
devGuide/Field
devGuide/IO
devGuide/IOStreams
devGuide/Halo
devGuide/HorzMesh
devGuide/HorzOperators
devGuide/AuxiliaryVariables
devGuide/AuxiliaryState
devGuide/TendencyTerms
devGuide/Tendencies
devGuide/OceanState
devGuide/TimeMgr
devGuide/TimeStepping
devGuide/Reductions
devGuide/Tracers
```

```{toctree}
:caption: Design documents
:maxdepth: 1

design/OmegaV0ShallowWater
design/Broadcast
design/Config
design/DataTypes
design/Decomp
design/Driver
design/Error
design/Halo
design/HorzMeshClass
design/Logging
design/MachEnv
design/Metadata
design/IO
design/IOStreams
design/Reductions
design/State
design/Tendency
design/Tendencies
design/AuxiliaryVariables
design/AuxiliaryState
design/TimeMgr
design/Timers
design/TimeStepping
design/Tracers

design/Template
```
