[![E3SM Logo](https://e3sm.org/wp-content/themes/e3sm/assets/images/e3sm-logo.png)](https://e3sm.org)

Energy Exascale Earth System Model (E3SM)
================================================================================
E3SM is a state-of-the-art fully coupled model of the Earth's climate including
important biogeochemical and cryospheric processes. It is intended to address
the most challenging and demanding climate-change research problems and
Department of Energy mission needs while efficiently using DOE Leadership
Computing Facilities.

DOI: [10.11578/E3SM/dc.20180418.36](http://dx.doi.org/10.11578/E3SM/dc.20180418.36)


Please visit the [project website](https://e3sm.org) for further details.

This branch contains code modifications to produce all the  all the long-term (10-year) 
E3SM climate simulations with FV dynamical core and cam4 physics package 
(with different configurations for the macro-physics (condensation) parameterization)
for Wan et al. (2020), Improving time-step convergence in an atmosphere model with
simplified physics: the impacts of closure assumption and process coupling,
Journal of Advances in Modeling Earth Systems, https://doi.org/10.1029/2019MS001982.

Compset and namelist settings to switch between different configurations for large-scale condensation:
---------------------------------------------------------------------------------------------------------

Baseline model (Original spliting) (COMPSET F)
-----------------------------------------------------------------
 ql_incld_opt      = 0

 lc_tend_opt       = 0

 fmin              = 1.0E-4

Revised splitting + original closure (COMPSET F)
-----------------------------------------------------------------
 ql_incld_opt      = 1

 lc_tend_opt       = 0

 fmin              = 1.0E-4

Revised splitting & closure (COMPSET F)
-----------------------------------------------------------------
 ql_incld_opt      = 1

 lc_tend_opt       = 1

 fmin              = 1.0E-4

