# The E3SM Sea Ice Model (MPAS-seaice)

MPAS-seaice is an unstructured-mesh sea-ice model that uses the Modeling for Prediction Across Scales (MPAS) framework, allowing enhanced horizontal resolution in regions of interest. MPAS-seaice uses many of the methods used in the Los Alamos CICE sea-ice model, but adapted to the Spherical Centroidal Vornoi Tesselation (SCVT) meshes used by the MPAS framework. 

* The [MPAS-seaice User's Guide](user-guide/index.md) explains how to control MPAS-seaice within E3SM.
* The [MPAS-seaice Developer's Guide](dev-guide/index.md) explains MPAS-seaice data structures and how to develop new code.
* The [MPAS-seaice Technical Guide](tech-guide/index.md) explains the science behind MPAS-seaice code.

**Icepack**

MPAS-seaice incorporates the Icepack software package for sea ice column physics, developed by the [CICE Consortium](https://github.com/cice-consortium), as a submodule. Icepack documentation provides a complete description of the column physics and instructions for using Icepack as a standalone model. The source code for this documentation is maintained in the Consortium's [Icepack repository](https://github.com/cice-consortium/Icepack), and it is fully rendered at 
[readthedocs](https://cice-consortium-icepack.readthedocs.io/en/main/).

If modifications have been made to the Icepack repository and documentation that have not yet been migrated to E3SM's fork, then the documentation may not match E3SM's version of the Icepack code.  E3SM does not render the readthedocs documentation from within its fork of the Icepack repository, but the source can be viewed directly at [https://github.com/E3SM-Project/Icepack/](https://github.com/E3SM-Project/Icepack/) (navigate to doc/source/ and then likely science_guide/, etc).  This is the documentation associated with the latest Icepack version that has been merged into E3SM, plus any documentation changes made within E3SM itself.

If the source code is difficult to read, the Icepack repository's documentation for the most recent release incorporated in E3SM can be found in readthedocs:

* Check the [release tags](https://github.com/E3SM-Project/Icepack/tags) to get the release number.
* Choose the release version of the documentation from the [Icepack release table](https://github.com/CICE-Consortium/Icepack/wiki/Icepack-Release-Table).

[Guidance for developing Icepack documentation](https://github.com/CICE-Consortium/About-Us/wiki/Documentation-Workflow-Guide) includes instructions for building the readthedocs documentation yourself.