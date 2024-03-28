# The E3SM Sea Ice Model (MPAS-seaice)

MPAS-seaice is an unstructured-mesh sea-ice model that uses the Modeling for Prediction Across Scales (MPAS) framework, allowing enhanced horizontal resolution in regions of interest. MPAS-seaice uses many of the methods used in the Los Alamos CICE sea-ice model, but adapted to the Spherical Centroidal Vornoi Tesselation (SCVT) meshes used by the MPAS framework. 

* The [MPAS-seaice User's Guide](user-guide/index.md) explains how to control MPAS-seaice within E3SM.
* The [MPAS-seaice Developer's Guide](dev-guide/index.md) explains MPAS-seaice data structures and how to develop new code.
* The [MPAS-seaice Technical Guide](tech-guide/index.md) explains the science behind MPAS-seaice code.

**Icepack**
-----------

MPAS-seaice incorporates the Icepack software package for sea ice column physics, developed by the [CICE Consortium](https://github.com/cice-consortium), as a submodule. Icepack documentation provides a complete description of the column physics and instructions for using Icepack as a standalone model. The source code for this documentation is maintained in [E3SM's Icepack fork](https://github.com/E3SM-Project/Icepack/) (navigate to the desired branch, then to doc/source/, etc).  This is the documentation associated with the latest Icepack version that has been merged into E3SM, plus any documentation changes made within E3SM itself. This documentation is fully rendered in [E3SM's Icepack readthedocs](https://e3sm-icepack.readthedocs.io/en/latest/).

If needed, documentation for the most recent Icepack release incorporated in E3SM can be found in the CICE Consortium's readthedocs project area:

* Check the [release tags](https://github.com/E3SM-Project/Icepack/tags) to get the release number.
* Choose the release version of the documentation from the [Icepack release table](https://github.com/CICE-Consortium/Icepack/wiki/Icepack-Release-Table).

[Guidance for developing Icepack documentation](https://github.com/CICE-Consortium/About-Us/wiki/Documentation-Workflow-Guide) includes instructions for building the readthedocs documentation yourself.


**MPAS-seaice code structure**
------------------------------

Some MPAS-seaice functionality is sourced from the MPAS Framework:
``E3SM/components/mpas-framework``.  In particular, see ``E3SM/components/mpas-framework/core_seaice``.

Code structure within the ``mpas-seaice/``component-level directory:

| Directories | Function |
| ----------- | -------- |
| ``  bld``         | namelist configuration files |
| ``  cime_config`` | build and configuration scripts |
| ``  docs``        | this documentation |
| ``  driver``      | coupling modules |
| ``  src``         | source code for the model physics and output |
| ``  src/analysis_members`` | source code for model output |
| ``  src/column ``          | source code for the (original) ``column_package`` |
| ``  src/icepack ``         | link to icepack submodule |
| ``  src/model_forward``    | top-level mpas-seaice modules |
| ``  src/shared``           | dynamics and general-purpose modules (e.g. mesh, constants) |
| ``  testing``     | testing scripts |

