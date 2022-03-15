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
    nparams = 97 ! Total tunable parameters

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
    iC2rt             =  4, & 
    iC2thl            =  5, & 
    iC2rtthl          =  6, & 
    iC4               =  7, & 
    iC_uu_shr         =  8, &
    iC_uu_buoy        =  9, & 
    iC6rt             = 10, & 
    iC6rtb            = 11, & 
    iC6rtc            = 12, & 
    iC6thl            = 13, & 
    iC6thlb           = 14, & 
    iC6thlc           = 15, & 
    iC7               = 16, & 
    iC7b              = 17, & 
    iC7c              = 18, & 
    iC8               = 19, & 
    iC8b              = 20, & 
    iC10              = 21, & 
    iC11              = 22, & 
    iC11b             = 23, & 
    iC11c             = 24, & 
    iC12              = 25, & 
    iC13              = 26, & 
    iC14              = 27, &
    iC_wp2_pr_dfsn    = 28, &
    iC_wp3_pr_tp      = 29, &
    iC_wp3_pr_turb    = 30, &
    iC_wp3_pr_dfsn    = 31, &
    iC_wp2_splat      = 32

  integer, parameter, public :: &
    iC6rt_Lscale0     = 33, &
    iC6thl_Lscale0    = 34, &
    iC7_Lscale0       = 35, &
    iwpxp_L_thresh    = 36

  integer, parameter, public :: & 
    ic_K              = 37, & 
    ic_K1             = 38, & 
    inu1              = 39, & 
    ic_K2             = 40, & 
    inu2              = 41, & 
    ic_K6             = 42, & 
    inu6              = 43, & 
    ic_K8             = 44, & 
    inu8              = 45, & 
    ic_K9             = 46, & 
    inu9              = 47, & 
    inu10             = 48, &
    ic_K_hm           = 49, & 
    ic_K_hmb          = 50, & 
    iK_hm_min_coef    = 51, & 
    inu_hm            = 52 

  integer, parameter, public :: &
    islope_coef_spread_DG_means_w = 53, &
    ipdf_component_stdev_factor_w = 54, &
    icoef_spread_DG_means_rt      = 55, &
    icoef_spread_DG_means_thl     = 56, &
    igamma_coef                   = 57, & 
    igamma_coefb                  = 58, & 
    igamma_coefc                  = 59, & 
    imu                           = 60, & 
    ibeta                         = 61, & 
    ilmin_coef                    = 62, &
    iomicron                      = 63, &
    izeta_vrnce_rat               = 64, &
    iupsilon_precip_frac_rat      = 65, &
    ilambda0_stability_coef       = 66, &
    imult_coef                    = 67, &
    itaumin                       = 68, &
    itaumax                       = 69, &
    iLscale_mu_coef               = 70, &
    iLscale_pert_coef             = 71, &
    ialpha_corr                   = 72, &
    iSkw_denom_coef               = 73, &
    ic_K10                        = 74, &
    ic_K10h                       = 75, &
    ithlp2_rad_coef               = 76, &
    ithlp2_rad_cloud_frac_thresh  = 77, &
    iup2_sfc_coef                 = 78, &
    iSkw_max_mag                  = 79, &
    iC_invrs_tau_bkgnd            = 80, &
    iC_invrs_tau_sfc              = 81, &
    iC_invrs_tau_shear            = 82, &
    iC_invrs_tau_N2               = 83, &
    iC_invrs_tau_N2_wp2           = 84, &
    iC_invrs_tau_N2_xp2           = 85, &
    iC_invrs_tau_N2_wpxp          = 86, &
    iC_invrs_tau_N2_clear_wp3     = 87, &
    iC_invrs_tau_wpxp_Ri          = 88, &
    iC_invrs_tau_wpxp_N2_thresh   = 89, &
    ixp3_coef_base                = 90, &
    ixp3_coef_slope               = 91, &
    ialtitude_threshold           = 92, &
    irtp2_clip_coef               = 93, &
    iCx_min                       = 94, &
    iCx_max                       = 95, &
    iRichardson_num_min           = 96, &
    iRichardson_num_max           = 97

end module parameter_indices
!-----------------------------------------------------------------------
