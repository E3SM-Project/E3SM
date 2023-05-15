(omega-home)=
# OMEGA

The Ocean Model for E3SM Global Applications (OMEGA) is an eddy-resolving, 
global ocean model in the early stages of development by the 
[E3SM](https://e3sm.org/) ocean team.  The first release is planned for
Summer 2026.

The model is written in c++ using the [YAKL](https://github.com/mrnorman/YAKL)
framework for performance portability.  OMEGA is based on the 
[TRSK formulation](https://doi.org/10.1016/j.jcp.2009.08.006) for geophysical 
models on unstructured meshes. The first version of OMEGA will be a direct port
of the [MPAS-Ocean](https://e3sm.org/model/e3sm-model-description/v1-description/v1-ocean-sea-ice-land-ice/)
component of E3SM for comparison purposes.

Development is taking place at https://github.com/E3SM-Project/Omega.


```{toctree}
:caption: User's guide
:maxdepth: 2

userGuide/QuickStart
```

```{toctree}
:caption: Developer's guide
:maxdepth: 2

devGuide/QuickStart
devGuide/CondaEnv
devGuide/Docs
devGuide/BuildDocs
```

```{toctree}
:caption: Design documents
:maxdepth: 1

design/Broadcast
design/Config
design/DataTypes
design/Halo
design/Logging
design/MachEnv
design/TimeMgr

design/OmegaDesignTemplate
```
