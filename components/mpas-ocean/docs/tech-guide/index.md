# MPAS-Ocean Technical Guide

This Technical Guide describes the governing equations, physics, and numerical discretizations of MPAS-Ocean.

## Guides, Design Documents, and documentation

MPAS-Ocean may be run as either a stand-alone model, or within the E3SM coupled earth system model.
The MPAS-Ocean code for both stand-alone and coupled is housed in the repository [https://github.com/E3SM-Project/E3SM](https://github.com/E3SM-Project/E3SM) within the directory `components/mpas-ocean`. The stand-alone executable may be built within that directory using the make command with the required libraries, as described in Chapter 1 of the [User's Guide](https://zenodo.org/records/11098080).

The [MPAS-Ocean User's Guide](https://zenodo.org/records/11098080) provides a description of the MPAS Framework in Part I, the Governing equations for MPAS-Ocean in Chapter 8, and describes the physics behind each term and parameterization in chapter 11.

All new features are created with design documents. The location of these documents have moved over the years, but can still be found in these locations:

1. [MPAS Documents repository](https://github.com/MPAS-Dev/MPAS-Documents/tree/master/ocean)
2. [E3SM repository MPAS-Ocean docs](https://github.com/E3SM-Project/E3SM/tree/master/components/mpas-ocean/docs/design_docs)
3. Documents for the new Omega model are similar to MPAS-Ocean and may be found in the [Omega documentation](https://docs.e3sm.org/Omega/develop/index.html).

All test cases are housed in the [Compass repository](https://github.com/MPAS-Dev/compass) and, more recently, the [Polaris repository](https://github.com/E3SM-Project/polaris). The corresponding documentation is housed in the [Compass docs](https://mpas-dev.github.io/compass/latest/) and [Polaris docs](http://docs.e3sm.org/polaris/main/) pages.

## Publications

Beyond the documentation, there are many publications that describe the inner workings of MPAS-Ocean:

Ringler, T., Petersen, M., Higdon, R.L., Jacobsen, D., Jones, P.W., Maltrud, M., 2013.
[A multi-resolution approach to global ocean modeling](https://doi.org/10.1016/j.ocemod.2013.04.010). Ocean Modelling 69, 211-232.

Petersen, M.R., D.W. Jacobsen, T.D. Ringler, M.W. Hecht, M.E. Maltrud, [Evaluation of the arbitrary Lagrangian–Eulerian vertical coordinate method in the MPAS-Ocean model](http://dx.doi.org/10.1016/j.ocemod.2014.12.004), Ocean Modelling, Volume 86, February 2015, Pages 93-113, ISSN 1463-5003.

Petersen, M. R., Asay‐Davis, X. S., Berres, A. S., Chen, Q., Feige, N., Hoffman, M. J., et al. (2019). [An evaluation of the ocean and sea ice climate of E3SM using MPAS and interannual CORE‐II forcing](https://doi.org/10.1029/2018MS001373). Journal of Advances in Modeling Earth Systems, 11, 1438– 1458.
