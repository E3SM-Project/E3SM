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

This branch contains code modifications to produce all the short-term (1-h -- 12-h) 
convergence test simulations (with different configurations for simple condensation model)
for Wan et al. (2020), Improving time-step convergence in an atmosphere model with 
simplified physics: the impacts of closure assumption and process coupling, 
Journal of Advances in Modeling Earth Systems, https://doi.org/10.1029/2019MS001982.


Compset and namelist settings to active SE dynamical core with/without simple condensation parameterization:
-----------------------------------------------------------------

SE dynamical core only (compset FC5AQUAP)
-----------------------------------------------------------------
 se_ftype           = 1
 
 deep_scheme        = 'off',
 
 shallow_scheme     = 'off',
 
 l_tracer_aero      = .false.
 
 l_vdiff            = .false.
 
 l_rayleigh         = .false.
 
 l_gw_drag          = .false.
 
 l_ac_energy_chk    = .true.
 
 l_bc_energy_fix    = .true.
 
 l_dry_adj          = .false.
 
 l_st_mac           = .false.
 
 l_st_mic           = .false.
 
 l_rad              = .false.


SE dynamical core + simple condensation model (compset FC5AQUAP)
-----------------------------------------------------------------

 se_ftype           = 1
 
 simple_macrop_opt  = 2
 
 deep_scheme        = 'off',
 
 shallow_scheme     = 'off',
 
 reset_init_ql      = .true.
 
 rkz_cldfrc_opt     = 1
 
 rkz_zsmall_opt     = 0
 
 rkz_lmt5_opt       = 0
 
 l_tracer_aero      = .false.
 
 l_vdiff            = .false.
 
 l_rayleigh         = .false.
 
 l_gw_drag          = .false.
 
 l_ac_energy_chk    = .true.
 
 l_bc_energy_fix    = .true.
 
 l_dry_adj          = .false.
 
 l_st_mac           = .true.
 
 l_st_mic           = .false.
 
 l_rad              = .false.
 
 l_rkz_lmt_2        = .false.
 
 l_rkz_lmt_3        = .false.
 
 l_rkz_lmt_4        = .true.
 
 l_rkz_lmt_5        = .false.

Switch in namelist between different configurations of simple condesnation model:
-----------------------------------------------------------------

(A1+B1+C2)

rkz_term_A_opt    = 1

rkz_term_B_opt    = 1

rkz_term_C_opt    = 2


(A1+B3+C2)

rkz_term_A_opt    = 1

rkz_term_B_opt    = 3

rkz_term_C_opt    = 2

Switch in namelist between different calculations of in-cloud liquid water
-----------------------------------------------------------------------------------------

(default)

rkz_term_C_ql_opt = 17


(revised splitting)

rkz_term_C_ql_opt = 19

Switch in namelist between different safeguard value of cloud fraction (fmin)
-----------------------------------------------------------------

(default)

rkz_term_C_fmin   = 1e-3


(test setup)

rkz_term_C_fmin   = 1e-2

rkz_term_C_fmin   = 1e-5

rkz_term_C_fmin   = 1e-8

rkz_term_C_fmin   = 1e-12


