!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
module parameter_indices

! Description:
!   Since f90/95 lacks enumeration, we're stuck numbering each
!   parameter by hand like this.

!   Adding new parameters is relatively simple.  First, the
!   parameter should be added in the common block of the parameters
!   module so it can be used in other parts of the code. Each
!   variable needs a unique number in this module, and nparams must
!   be incremented for the new variable.  Next, the params_list
!   variable in module parameters should have new variable added to
!   it.  The subroutines pack_parameters and uppack_parameters will
!   need to have the variable added to their list, but the order
!   doesn't actually matter, since the i variables in here determine
!   where in the params vector the number is placed.
!   Finally, the namelists clubb_params_nl and initspread will need to
!   have the parameter added to them.
!-------------------------------------------------------------------------------

  implicit none

  private ! Default Scope

  integer, parameter, public ::  & 
    nparams = 82 ! Total tunable parameters

!***************************************************************
!                    ***** IMPORTANT *****
! If you change the order of these parameters, you will need to
! change the order of params_list as well or the tuner will
! break!
!                    ***** IMPORTANT *****
!***************************************************************

  integer, parameter, public :: & 
    iC1               =  1, & 
    iC1b              =  2, & 
    iC1c              =  3, & 
    iC2               =  4, & 
    iC2b              =  5, & 
    iC2c              =  6, & 
    iC2rt             =  7, & 
    iC2thl            =  8, & 
    iC2rtthl          =  9, & 
    iC4               = 10, & 
    iC5               = 11, & 
    iC6rt             = 12, & 
    iC6rtb            = 13, & 
    iC6rtc            = 14, & 
    iC6thl            = 15, & 
    iC6thlb           = 16, & 
    iC6thlc           = 17, & 
    iC7               = 18, & 
    iC7b              = 19, & 
    iC7c              = 20, & 
    iC8               = 21, & 
    iC8b              = 22, & 
    iC10              = 23, & 
    iC11              = 24, & 
    iC11b             = 25, & 
    iC11c             = 26, & 
    iC12              = 27, & 
    iC13              = 28, & 
    iC14              = 29, &
    iC15              = 30, &
    iC_wp2_splat      = 31

  integer, parameter, public :: &
    iC6rt_Lscale0     = 32, &
    iC6thl_Lscale0    = 33, &
    iC7_Lscale0       = 34, &
    iwpxp_L_thresh    = 35

  integer, parameter, public :: & 
    ic_K              = 36, & 
    ic_K1             = 37, & 
    inu1              = 38, & 
    ic_K2             = 39, & 
    inu2              = 40, & 
    ic_K6             = 41, & 
    inu6              = 42, & 
    ic_K8             = 43, & 
    inu8              = 44, & 
    ic_K9             = 45, & 
    inu9              = 46, & 
    inu10             = 47, &
    ic_K_hm           = 48, & 
    ic_K_hmb          = 49, & 
    iK_hm_min_coef    = 50, & 
    inu_hm            = 51 

  integer, parameter, public :: &
    islope_coef_spread_DG_means_w = 52, &
    ipdf_component_stdev_factor_w = 53, &
    icoef_spread_DG_means_rt      = 54, &
    icoef_spread_DG_means_thl     = 55, &
    igamma_coef                   = 56, & 
    igamma_coefb                  = 57, & 
    igamma_coefc                  = 58, & 
    imu                           = 59, & 
    ibeta                         = 60, & 
    ilmin_coef                    = 61, &
    iomicron                      = 62, &
    izeta_vrnce_rat               = 63, &
    iupsilon_precip_frac_rat      = 64, &
    ilambda0_stability_coef       = 65, &
    imult_coef                    = 66, &
    itaumin                       = 67, &
    itaumax                       = 68, &
    iLscale_mu_coef               = 69, &
    iLscale_pert_coef             = 70, &
    ialpha_corr                   = 71, &
    iSkw_denom_coef               = 72, &
    ic_K10                        = 73, &
    ic_K10h                       = 74, &
    ithlp2_rad_coef               = 75, &
    ithlp2_rad_cloud_frac_thresh  = 76, &
    iup2_vp2_factor               = 77, &
    iSkw_max_mag                  = 78, &
    iC_invrs_tau_bkgnd            = 79, &
    iC_invrs_tau_sfc              = 80, &
    iC_invrs_tau_shear            = 81, &
    iC_invrs_tau_N2               = 82


end module parameter_indices
!-----------------------------------------------------------------------
