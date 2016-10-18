!-----------------------------------------------------------------------
! $Id: parameters_tunable.F90 8220 2016-07-21 18:48:32Z raut@uwm.edu $
!===============================================================================
module parameters_tunable

  ! Description:
  !   This module contains tunable model parameters.  The purpose of the module is to make it
  !   easier for the clubb_tuner code to use the params vector without "knowing" any information
  !   about the individual parameters contained in the vector itself.  It makes it easier to add
  !   new parameters to be tuned for, but does not make the CLUBB_core code itself any simpler.
  !   The parameters within the vector do not need to be the same variables used in the rest of
  !   CLUBB_core (see for e.g. nu1_vert_res_dep or lmin_coef).
  !   The parameters in the params vector only need to be those parameters for which we're not
  !   sure the correct value and we'd like to tune for.
  !
  ! References:
  !   None
  ! 
  ! Notes:
  !   To make it easier to verify of code correctness, please keep the omp threadprivate
  !   directives just after the variable declaration.  All parameters in this
  !   module should be declared threadprivate because of the CLUBB tuner.
  !-----------------------------------------------------------------------

  use parameter_indices, only: nparams ! Variable(s)

  use grid_class, only: gr ! Variable(s)

  use clubb_precision, only: &
      core_rknd ! Variable(s)

  implicit none

  ! Default to private
  private

  public :: setup_parameters, read_parameters, read_param_spread, &
            get_parameters, adj_low_res_nu, cleanup_nu, clubb_param_readnl

  ! NOTE: In CLUBB standalone, as well as some host models, the hardcoded
  !       default values of some or all of the parameters below have no effect,
  !       as the values are simply read in using a namelist or set in host model
  !       specific code.

  ! The parameters below have the same meaning as those without prefix 'clubb_'
  ! They can be specified via namelist to ease parameter tuning
  ! If adding more parameters for tuning via namelist, need to insert blocks
  ! accordingly in subroutines read_parameters and clubb_paramm_readnl
  ! (including updating list of nml variables)

  real( kind = core_rknd ) ::      &
    clubb_C1,                      &
    clubb_C2rt,                    &
    clubb_C2thl,                   &
    clubb_C2rtthl,                 &
    clubb_C6rt,                    &
    clubb_C6rtb,                   &
    clubb_C7,                      &
    clubb_C7b,                     &
    clubb_C8,                      &
    clubb_C11,                     &
    clubb_C11b,                    &
    clubb_C14,                     &
    clubb_beta,                    &
    clubb_gamma_coef,              &
    clubb_gamma_coefb,             &
    clubb_mu,                      &
    clubb_nu1,                     &
    clubb_c_K10,                   &
    clubb_wpxp_L_thresh
    
  ! Model constant parameters
  real( kind = core_rknd ), public :: & 
    C1      = 1.000000_core_rknd,    & ! Low Skewness in C1 Skw. Function    [-]
    C1b     = 1.000000_core_rknd,    & ! High Skewness in C1 Skw. Function   [-]
    C1c     = 1.000000_core_rknd,    & ! Degree of Slope of C1 Skw. Function [-]
    C2      = 1.300000_core_rknd,    & ! Low Skewness in C2 Skw. Function    [-]
    C2rt    = 1.000000_core_rknd,    & ! C2 coef. for the rtp2_dp1 term      [-]
    C2thl   = 1.000000_core_rknd,    & ! C2 coef. for the thlp2_dp1 term     [-]
    C2rtthl = 2.000000_core_rknd,    & ! C2 coef. for the rtpthlp_dp1 term   [-]
    C2b     = 1.300000_core_rknd,    & ! High Skewness in C2 Skw. Function   [-]
    C2c     = 5.000000_core_rknd,    & ! Degree of Slope of C2 Skw. Function [-]
    C4      = 5.200000_core_rknd,    & ! Used only when l_tke_aniso is true  [-]
    C5      = 0.300000_core_rknd,    & ! Coef. in pressure terms: w'^2 eqn   [-]
    C6rt    = 4.000000_core_rknd,    & ! Low Skewness in C6rt Skw. Function  [-]
    C6rtb   = 6.000000_core_rknd,    & ! High Skewness in C6rt Skw. Function [-]
    C6rtc   = 1.000000_core_rknd,    & ! Degree of Slope of C6rt Skw. Fnct.  [-]
    C6thl   = 4.000000_core_rknd,    & ! Low Skewness in C6thl Skw. Function [-]
    C6thlb  = 6.000000_core_rknd,    & ! High Skewness in C6thl Skw. Fnct.   [-]
    C6thlc  = 1.000000_core_rknd,    & ! Degree of Slope of C6thl Skw. Fnct. [-]
    C7      = 0.500000_core_rknd,    & ! Low Skewness in C7 Skw. Function    [-]
    C7b     = 0.800000_core_rknd,    & ! High Skewness in C7 Skw. Function   [-]
    C7c     = 0.500000_core_rknd,    & ! Degree of Slope of C7 Skw. Function [-]
    C8      = 3.000000_core_rknd,    & ! Coef. #1 in C8 Skewness Equation    [-]
    C8b     = 0.000000_core_rknd,    & ! Coef. #2 in C8 Skewness Equation    [-]
    C10     = 3.300000_core_rknd,    & ! Currently Not Used in the Model     [-]
    C11     = 0.80000_core_rknd,     & ! Low Skewness in C11 Skw. Function   [-]
    C11b    = 0.350000_core_rknd,    & ! High Skewness in C11 Skw. Function  [-]
    C11c    = 0.500000_core_rknd,    & ! Degree of Slope of C11 Skw. Fnct.   [-]
    C12     = 1.000000_core_rknd,    & ! Constant in w'^3 Crank-Nich. diff.  [-]
    C13     = 0.100000_core_rknd,    & ! Not currently used in model         [-]
    C14     = 1.000000_core_rknd,    & ! Constant for u'^2 and v'^2 terms    [-]
    C15     = 0.4_core_rknd            ! Coefficient for the wp3_bp2 term    [-]
!$omp threadprivate(C1, C1b, C1c, C2, C2b, C2c, &
!$omp   C2rt, C2thl, C2rtthl, C4, C5, C6rt, C6rtb, C6rtc, &
!$omp   C6thl, C6thlb, C6thlc, &
!$omp   C7, C7b, C7c, C8, C8b, C10, C11, C11b, C11c, C12, &
!$omp   C13, C14, C15)

  real( kind = core_rknd ), public ::    &
    C6rt_Lscale0  = 14.0_core_rknd,      & ! Damp C6rt as a fnct. of Lscale  [-]
    C6thl_Lscale0 = 14.0_core_rknd,      & ! Damp C6thl as a fnct. of Lscale [-]
    C7_Lscale0    = 0.8500000_core_rknd, & ! Damp C7 as a fnct. of Lscale    [-]
    wpxp_L_thresh = 60.0_core_rknd         ! Lscale threshold: damp C6 & C7  [m]
!$omp threadprivate(C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh)

  ! Note: DD 1987 is Duynkerke & Driedonks (1987).
  real( kind = core_rknd ), public :: & 
    c_K         = 0.200000_core_rknd, & ! Constant C_mu^(1/4) in DD 1987 [m^2/s]
    c_K1        = 0.750000_core_rknd, & ! Coef. of Eddy Diffusion: wp2   [m^2/s]
    c_K2        = 0.125000_core_rknd, & ! Coef. of Eddy Diffusion: xp2   [m^2/s]
    c_K6        = 0.375000_core_rknd, & ! Coef. of Eddy Diffusion: wpxp  [m^2/s]
    c_K8        = 1.250000_core_rknd, & ! Coef. of Eddy Diffusion: wp3   [m^2/s]
    c_K9        = 0.250000_core_rknd, & ! Coef. of Eddy Diff.: up2/vp2   [m^2/s]
    c_K_hm      = 0.750000_core_rknd, & ! Coef. of Eddy Diffusion: hmm   [m^2/s]
    c_K_hmb     = 0.10000_core_rknd,  & ! Coef. of Non-Local Factor, Eddy Diffusion: hmm   [m^2/s]
    K_hm_min_coef = 0.10000_core_rknd,& ! Min. of Non-Local Factor, Eddy Diffusion: hmm   [m^2/s]
    gamma_coef  = 0.320000_core_rknd, & ! Low Skw.: gamma coef. Skw. Fnct.   [-]
    gamma_coefb = 0.320000_core_rknd, & ! High Skw.: gamma coef. Skw. Fnct.  [-]
    gamma_coefc = 5.000000_core_rknd, & ! Deg. Slope: gamma coef. Skw. Fnct. [-]
    mu          = 1.000E-3_core_rknd, & ! Fract entrain rate per unit alt  [1/m]
    mult_coef   = 0.500000_core_rknd, & ! Coef. applied to log(avg dz/thresh)[-]
    taumin      = 90.00000_core_rknd, & ! Min. allow. value: time-scale tau  [s]
    taumax      = 3600.000_core_rknd, & ! Max. allow. value: time-scale tau  [s]
    lmin        = 20.00000_core_rknd    ! Min. value for the length scale    [m]
!$omp threadprivate(c_K, c_K1, c_K2, c_K6, &
!$omp   c_K8, c_K9, c_K_hm, c_K_hmb, K_hm_min_coef, gamma_coef, gamma_coefb, gamma_coefc, &
!$omp   mu, mult_coef, taumin, taumax, lmin)

  real( kind = core_rknd ), public :: &
    Lscale_mu_coef   = 2.0_core_rknd, & ! Coef perturb mu: av calc Lscale    [-]
    Lscale_pert_coef = 0.1_core_rknd    ! Coef pert thlm/rtm: av calc Lscale [-]
!$omp threadprivate(Lscale_mu_coef, Lscale_pert_coef)

  real( kind = core_rknd ), public :: &
    alpha_corr = 0.15_core_rknd   ! Coef. for the corr. diagnosis algorithm  [-]

!$omp threadprivate(alpha_corr)

  real( kind = core_rknd ), private :: & 
    nu1   = 20.00000_core_rknd, & ! Bg. Coef. Eddy Diffusion: wp2        [m^2/s]
    nu2   = 5.000000_core_rknd, & ! Bg. Coef. Eddy Diffusion: xp2        [m^2/s]
    nu6   = 5.000000_core_rknd, & ! Bg. Coef. Eddy Diffusion: wpxp       [m^2/s]
    nu8   = 20.00000_core_rknd, & ! Bg. Coef. Eddy Diffusion: wp3        [m^2/s]
    nu9   = 20.00000_core_rknd, & ! Bg. Coef. Eddy Diffusion: up2/vp2    [m^2/s]
    nu10  = 0.000000_core_rknd, & ! Bg. Coef. Eddy Diffusion: edsclrm    [m^2/s]
    nu_hm = 1.500000_core_rknd    ! Bg. Coef. Eddy Diffusion: hmm        [m^2/s]
!$omp threadprivate(nu1, nu2, nu6, nu8, nu9, nu10, nu_hm)


  real( kind = core_rknd ), public, allocatable, dimension(:) :: & 
    nu1_vert_res_dep,   & ! Background Coef. of Eddy Diffusion: wp2      [m^2/s]
    nu2_vert_res_dep,   & ! Background Coef. of Eddy Diffusion: xp2      [m^2/s]
    nu6_vert_res_dep,   & ! Background Coef. of Eddy Diffusion: wpxp     [m^2/s]
    nu8_vert_res_dep,   & ! Background Coef. of Eddy Diffusion: wp3      [m^2/s]
    nu9_vert_res_dep,   & ! Background Coef. of Eddy Diffusion: up2/vp2  [m^2/s]
    nu10_vert_res_dep,  & ! Background Coef. of Eddy Diffusion: edsclrm  [m^2/s]
    nu_hm_vert_res_dep    ! Background Coef. of Eddy Diffusion: hydromet [m^2/s]

!$omp threadprivate(nu1_vert_res_dep, nu2_vert_res_dep, nu6_vert_res_dep, &
!$omp   nu8_vert_res_dep, nu9_vert_res_dep, nu10_vert_res_dep, nu_hm_vert_res_dep)

  ! Vince Larson added a constant to set plume widths for theta_l and rt
  ! beta should vary between 0 and 3.

  real( kind = core_rknd ), public :: &
    beta = 2.400000_core_rknd    ! Beta coefficient     [-]

!$omp threadprivate(beta)

  real( kind = core_rknd ), private :: &
    lmin_coef = 0.500000_core_rknd   ! Coefficient of lmin    [-]

!$omp threadprivate(lmin_coef)

  ! Brian Griffin added a parameter for hydrometeors, omicron, to increase the
  ! standard deviation of each component and decrease the spread between the
  ! component means as the value of omicron inreases.  Valid value are
  ! 0 < omicron <= 1.
  ! A second parameter for hydrometeors, zeta, increases the standard deviation
  ! of component 1 at the expense of the standard deviation of component 2 when
  ! the value of zeta > 0 (and increasingly so as zeta increases).  Valid values
  ! are zeta > -1.
  real( kind = core_rknd ), public :: &
    omicron        = 0.8_core_rknd, & ! Hydromet width/spread-of-means param [-]
    zeta_vrnce_rat = 0.0_core_rknd    ! Ratio sigma^2/mu^2 comp. 1 / comp. 2 [-]

!$omp threadprivate( omicron, zeta_vrnce_rat )

  real( kind = core_rknd ), public :: &
    ! ratio mixt_frac*precip_frac_1/precip_frac (precip_frac_calc_type=2)    [-]
    upsilon_precip_frac_rat = 0.9_core_rknd, &
    ! Intensity of stability correction applied to C1 and C6 [-]
    lambda0_stability_coef = 0.03_core_rknd

!$omp threadprivate( upsilon_precip_frac_rat, lambda0_stability_coef)

  ! Factor to decrease sensitivity in the denominator of Skw calculation
  real( kind = core_rknd ), public :: &
    Skw_denom_coef = 4.0_core_rknd

!$omp threadprivate( Skw_denom_coef )

  ! Momentum coefficient of Kh_zm
  real( kind = core_rknd ), public :: &
    c_K10 = 1.0_core_rknd

!$omp threadprivate( c_K10 )

  ! Thermodynamic coefficient of Kh_zm
  real( kind = core_rknd ), public :: &
    c_K10h = 1.0_core_rknd

!$omp threadprivate( c_K10h )

  real( kind = core_rknd ), public :: &
    thlp2_rad_coef = 1.0_core_rknd, &            ! Coefficient of thlp2_rad                   [-]
    thlp2_rad_cloud_frac_thresh = 0.1_core_rknd ! Minimum cloud fraction for computation
                                                 ! of thlp2_rad                               [-]

!$omp threadprivate( thlp2_rad_coef, thlp2_rad_cloud_frac_thresh )

  real( kind = core_rknd ), public :: &
    up2_vp2_factor = 2.0_core_rknd               ! Coefficients of up2 and vp2    [-]

!$omp threadprivate( up2_vp2_factor )

  ! used in adj_low_res_nu. If .true., avg_deltaz = deltaz
#ifdef GFDL
  logical, public :: l_prescribed_avg_deltaz = .true.
#else
  logical, public :: l_prescribed_avg_deltaz = .false.
#endif

!$omp threadprivate(l_prescribed_avg_deltaz)

  ! Since we lack a devious way to do this just once, this namelist
  ! must be changed as well when a new parameter is added.
  namelist /clubb_params_nl/  & 
    C1, C1b, C1c, C2, C2b, C2c,  & 
    C2rt, C2thl, C2rtthl, C4, C5, & 
    C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, & 
    C7, C7b, C7c, C8, C8b, C10, C11, C11b, C11c, & 
    C12, C13, C14, C15, C6rt_Lscale0, C6thl_Lscale0, &
    C7_Lscale0, wpxp_L_thresh, c_K, c_K1, nu1, c_K2, nu2, & 
    c_K6, nu6, c_K8, nu8, c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, &
    nu_hm, beta, gamma_coef, gamma_coefb, gamma_coefc, lmin_coef, &
    omicron, zeta_vrnce_rat, upsilon_precip_frac_rat, &
    lambda0_stability_coef, mult_coef, taumin, taumax, mu, Lscale_mu_coef, &
    Lscale_pert_coef, alpha_corr, Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef, &
    thlp2_rad_cloud_frac_thresh, up2_vp2_factor

  ! These are referenced together often enough that it made sense to
  ! make a list of them.  Note that lmin_coef is the input parameter,
  ! while the actual lmin model constant is computed from this.
  !***************************************************************
  !                    ***** IMPORTANT *****
  ! If you change the order of the parameters in the parameter_indices,
  ! you will need to change the order of this list as well or the
  ! tuner will break!
  !                    ***** IMPORTANT *****
  !***************************************************************
  character(len=27), dimension(nparams), parameter, public ::  & 
  params_list = & 
     (/"C1                         ", "C1b                        ", &
       "C1c                        ", "C2                         ", &
       "C2b                        ", "C2c                        ", &
       "C2rt                       ", "C2thl                      ", &
       "C2rtthl                    ", "C4                         ", &
       "C5                         ", "C6rt                       ", &
       "C6rtb                      ", "C6rtc                      ", &
       "C6thl                      ", "C6thlb                     ", &
       "C6thlc                     ", "C7                         ", &
       "C7b                        ", "C7c                        ", &
       "C8                         ", "C8b                        ", &
       "C10                        ", "C11                        ", &
       "C11b                       ", "C11c                       ", &
       "C12                        ", "C13                        ", &
       "C14                        ", "C15                        ", &
       "C6rt_Lscale0               ", "C6thl_Lscale0              ", &
       "C7_Lscale0                 ", "wpxp_L_thresh              ", &
       "c_K                        ", "c_K1                       ", &
       "nu1                        ", "c_K2                       ", &
       "nu2                        ", "c_K6                       ", &
       "nu6                        ", "c_K8                       ", &
       "nu8                        ", "c_K9                       ", &
       "nu9                        ", "nu10                       ", &
       "c_K_hm                     ", "c_K_hmb                    ", &
       "K_hm_min_coef              ", "nu_hm                      ", &
       "gamma_coef                 ", "gamma_coefb                ", &
       "gamma_coefc                ", "mu                         ", &
       "beta                       ", "lmin_coef                  ", &
       "omicron                    ", "zeta_vrnce_rat             ", &
       "upsilon_precip_frac_rat    ", "lambda0_stability_coef     ", &
       "mult_coef                  ", "taumin                     ", &
       "taumax                     ", "Lscale_mu_coef             ", &
       "Lscale_pert_coef           ", "alpha_corr                 ", &
       "Skw_denom_coef             ", "c_K10                      ", &
       "c_K10h                     ", "thlp2_rad_coef             ", &
       "thlp2_rad_cloud_frac_thresh", "up2_vp2_factor             "  /)

  real( kind = core_rknd ), parameter, private :: &
    init_value = -999._core_rknd ! Initial value for the parameters, used to detect missing values

  contains

  !=============================================================================
  subroutine clubb_param_readnl(filename)

    ! Description:
    ! Read clubb tunable parameters from namelist to override preset values
    ! To be called by clubb_readnl in clubb_intr.F90
    ! Author: Wuyin Lin

    !-----------------------------------------------------------------------
    use spmd_utils,      only: masterproc
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use cam_abortutils,  only: endrun
    use mpishorthand

    implicit none

    ! Input variables

    character(len=*), intent(in) :: filename

    namelist /clubb_param_nl/      &
    clubb_C1,                      &
    clubb_C2rt,                    &
    clubb_C2thl,                   &
    clubb_C2rtthl,                 &
    clubb_C6rt,                    &
    clubb_C6rtb,                   &
    clubb_C7,                      &
    clubb_C7b,                     &
    clubb_C8,                      &
    clubb_C11,                     &
    clubb_C11b,                    &
    clubb_C14,                     &
    clubb_beta,                    &
    clubb_gamma_coef,              &
    clubb_gamma_coefb,             &
    clubb_mu,                      &
    clubb_nu1,                     &
    clubb_c_K10,                   &
    clubb_wpxp_L_thresh

    integer :: read_status
    integer :: iunit

    ! ---- Begin Code ----

    ! If clubb_tunables_nl present, read in to replace preset values
    ! This is made available for tuning 
     
    clubb_C1 = init_value
    clubb_C2rt = init_value
    clubb_C2thl = init_value
    clubb_C2rtthl = init_value
    clubb_C6rt = init_value
    clubb_C6rtb = init_value
    clubb_C7 = init_value
    clubb_C7b = init_value
    clubb_C8 = init_value
    clubb_C11 = init_value
    clubb_C11b = init_value
    clubb_C14 = init_value
    clubb_beta = init_value
    clubb_gamma_coef = init_value
    clubb_gamma_coefb = init_value
    clubb_mu = init_value
    clubb_nu1 = init_value
    clubb_c_K10 = init_value
    clubb_wpxp_L_thresh = init_value

    if (masterproc) then
      iunit = getunit()
      open( iunit, file=trim(filename), status='old' )
      call find_group_name(iunit, 'clubb_param_nl', status=read_status)
      if (read_status == 0) then
         read(unit=iunit, nml=clubb_param_nl,iostat=read_status)
         if (read_status /= 0) then
            call endrun('clubb_param_readnl:  error reading namelist')
         end if
      endif
      close(unit=iunit)
      call freeunit(iunit)
    end if
#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(clubb_C1,         1, mpir8,  0, mpicom)
   call mpibcast(clubb_C2rt,       1, mpir8,  0, mpicom)
   call mpibcast(clubb_C2thl,      1, mpir8,  0, mpicom)
   call mpibcast(clubb_C2rtthl,    1, mpir8,  0, mpicom)
   call mpibcast(clubb_C6rt,       1, mpir8,  0, mpicom)
   call mpibcast(clubb_C6rtb,      1, mpir8,  0, mpicom)
   call mpibcast(clubb_C7,         1, mpir8,  0, mpicom)
   call mpibcast(clubb_C7b,        1, mpir8,  0, mpicom)
   call mpibcast(clubb_C8,         1, mpir8,  0, mpicom)
   call mpibcast(clubb_C11,        1, mpir8,  0, mpicom)
   call mpibcast(clubb_C11b,       1, mpir8,  0, mpicom)
   call mpibcast(clubb_C14,        1, mpir8,  0, mpicom)
   call mpibcast(clubb_beta,       1, mpir8,  0, mpicom)
   call mpibcast(clubb_gamma_coef, 1, mpir8,  0, mpicom)
   call mpibcast(clubb_gamma_coefb,1, mpir8,  0, mpicom)
   call mpibcast(clubb_mu,         1, mpir8,  0, mpicom)
   call mpibcast(clubb_nu1,        1, mpir8,  0, mpicom)
   call mpibcast(clubb_c_K10,      1, mpir8,  0, mpicom)
   call mpibcast(clubb_wpxp_L_thresh, 1, mpir8,  0, mpicom)
#endif


  end subroutine clubb_param_readnl
  !=============================================================================
  subroutine setup_parameters & 
            ( deltaz, params, nzmax, &
              grid_type, momentum_heights, thermodynamic_heights, &
              err_code )

    ! Description:
    ! Subroutine to setup model parameters

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use constants_clubb, only:  & 
      fstderr ! Variable(s)

    use error_code, only:  & 
      clubb_var_out_of_bounds,  & ! Variable(s)
      clubb_no_error

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none


    ! Constant Parameters
    real( kind = core_rknd ), parameter :: &
      lmin_deltaz = 40.0_core_rknd ! Fixed value for minimum value for the length scale.

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      deltaz  ! Change per height level        [m]

    real( kind = core_rknd ), intent(in), dimension(nparams) :: & 
      params  ! Tuneable model parameters      [-]

    ! Grid definition
    integer, intent(in) :: nzmax  ! Vertical grid levels            [#]

    ! If CLUBB is running on its own, this option determines
    ! if it is using:
    ! 1) an evenly-spaced grid,
    ! 2) a stretched (unevenly-spaced) grid entered on the
    !    thermodynamic grid levels (with momentum levels set
    !    halfway between thermodynamic levels), or
    ! 3) a stretched (unevenly-spaced) grid entered on the
    !    momentum grid levels (with thermodynamic levels set
    !    halfway between momentum levels).
    integer, intent(in) :: grid_type

    ! If the CLUBB parameterization is implemented in a host model,
    ! it needs to use the host model's momentum level altitudes
    ! and thermodynamic level altitudes.
    ! If the CLUBB model is running by itself, but is using a
    ! stretched grid entered on thermodynamic levels (grid_type = 2),
    ! it needs to use the thermodynamic level altitudes as input.
    ! If the CLUBB model is running by itself, but is using a
    ! stretched grid entered on momentum levels (grid_type = 3),
    ! it needs to use the momentum level altitudes as input.
    real( kind = core_rknd ), intent(in), dimension(nzmax) :: &
      momentum_heights,      & ! Momentum level altitudes (input)      [m]
      thermodynamic_heights    ! Thermodynamic level altitudes (input) [m]

    ! Output Variables
    integer, intent(out) ::  &
      err_code ! Error condition

    !-------------------- Begin code --------------------

    call unpack_parameters( params, & 
                            C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl, &
                            C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, &
                            C7, C7b, C7c, C8, C8b, C10, & 
                            C11, C11b, C11c, C12, C13, C14, C15, & 
                            C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
                            c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, & 
                            c_K8, nu8, c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, &
                            nu_hm, gamma_coef, gamma_coefb, gamma_coefc, & 
                            mu, beta, lmin_coef, omicron, zeta_vrnce_rat, &
                            upsilon_precip_frac_rat, lambda0_stability_coef, &
                            mult_coef, taumin, taumax, Lscale_mu_coef, Lscale_pert_coef, &
                            alpha_corr, Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef, &
                            thlp2_rad_cloud_frac_thresh, up2_vp2_factor )


    ! It was decided after some experimentation, that the best
    ! way to produce grid independent results is to set lmin to be
    ! some fixed value. -dschanen 21 May 2007
    !lmin = lmin_coef * deltaz  ! Old
    lmin = lmin_coef * lmin_deltaz ! New fixed value

    ! ### Adjust Constant Diffusivity Coefficients Based On Grid Spacing ###
    call adj_low_res_nu &
           ( nzmax, grid_type, deltaz,  & ! Intent(in)
             momentum_heights, thermodynamic_heights )   ! Intent(in)

    ! Sanity check
    ! Initialize err_code to clubb_no_error.  Only overwrite it if a variable
    ! out-of-bounds error is found.
    err_code = clubb_no_error

    if ( beta < 0.0_core_rknd .or. beta > 3.0_core_rknd ) then

       ! Constraints on beta
       write(fstderr,*) "beta = ", beta
       write(fstderr,*) "beta cannot be < 0 or > 3"
       err_code = clubb_var_out_of_bounds

    endif ! beta < 0 or beta > 3

    if ( omicron <= 0.0_core_rknd .or. omicron > 1.0_core_rknd ) then

       ! Constraints on omicron
       write(fstderr,*) "omicron = ", omicron
       write(fstderr,*) "omicron cannot be <= 0 or > 1"
       err_code = clubb_var_out_of_bounds

    endif ! omicron <= 0 or omicron > 1

    if ( zeta_vrnce_rat <= -1.0_core_rknd ) then

       ! Constraints on zeta_vrnce_rat
       write(fstderr,*) "zeta_vrnce_rat = ", zeta_vrnce_rat
       write(fstderr,*) "zeta_vrnce_rat cannot be <= -1"
       err_code = clubb_var_out_of_bounds

    endif ! zeta_vrnce_rat <= -1

    if ( upsilon_precip_frac_rat < 0.0_core_rknd &
         .or. upsilon_precip_frac_rat > 1.0_core_rknd ) then

       ! Constraints on upsilon_precip_frac_rat
       write(fstderr,*) "upsilon_precip_frac_rat = ", upsilon_precip_frac_rat
       write(fstderr,*) "upsilon_precip_frac_rat cannot be < 0 or > 1"
       err_code = clubb_var_out_of_bounds

    endif ! upsilon_precip_frac_rat < 0 or upsilon_precip_frac_rat > 1

    if ( mu < 0.0_core_rknd ) then

       ! Constraints on entrainment rate, mu.
       write(fstderr,*) "mu = ", mu
       write(fstderr,*) "mu cannot be < 0"
       err_code = clubb_var_out_of_bounds

    endif ! mu < 0.0

    if ( lmin < 4.0_core_rknd ) then

       ! Constraints on mixing length
       write(fstderr,*) "lmin = ", lmin
       write(fstderr,*) "lmin is < 4.0_core_rknd"
       err_code = clubb_var_out_of_bounds

    endif ! lmin < 4.0

!    write(*,nml=clubb_params_nl) ! %% debug


    return

  end subroutine setup_parameters

  !=============================================================================
  subroutine adj_low_res_nu &
               ( nzmax, grid_type, deltaz, & ! Intent(in)
                 momentum_heights, thermodynamic_heights )  ! Intent(in)

    ! Description:
    !   Adjust the values of background eddy diffusivity based on
    !   vertical grid spacing.
    !   This code was made into a public subroutine so that it may be
    !   called multiple times per model run in scenarios where grid
    !   altitudes, and hence average grid spacing, change through space
    !   and/or time.  This occurs, for example, when CLUBB is
    !   implemented in WRF.  --ldgrant Jul 2010
    !----------------------------------------------------------------------

    use constants_clubb, only: &
      fstderr ! Constant(s)

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Constant Parameters

    ! Flag for adjusting the values of the constant background eddy diffusivity
    ! coefficients based on the average vertical grid spacing.  If this flag is
    ! turned off, the values of the various nu coefficients will remain as they
    ! are declared in the tunable_parameters.in file.
    logical, parameter :: l_adj_low_res_nu = .true.

    ! The size of the average vertical grid spacing that serves as a threshold
    ! for when to increase the size of the background eddy diffusivity
    ! coefficients (nus) by a certain factor above what the background
    ! coefficients are specified to be in tunable_parameters.in.  At any average
    ! grid spacing at or below this value, the values of the background
    ! diffusivities remain the same.  However, at any average vertical grid
    ! spacing above this value, the values of the background eddy diffusivities
    ! are increased.  Traditionally, the threshold grid spacing has been set to
    ! 40.0 meters.  This is only relevant if l_adj_low_res_nu is turned on.
    real( kind = core_rknd ), parameter :: &
      grid_spacing_thresh = 40.0_core_rknd  ! grid spacing threshold  [m]

    ! Input Variables

    ! Grid definition
    integer, intent(in) :: nzmax  ! Vertical grid levels            [#]

    ! If CLUBB is running on it's own, this option determines
    ! if it is using:
    ! 1) an evenly-spaced grid,
    ! 2) a stretched (unevenly-spaced) grid entered on the
    !    thermodynamic grid levels (with momentum levels set
    !    halfway between thermodynamic levels), or
    ! 3) a stretched (unevenly-spaced) grid entered on the
    !    momentum grid levels (with thermodynamic levels set
    !    halfway between momentum levels).
    integer, intent(in) :: grid_type

    real( kind = core_rknd ), intent(in) ::  & 
      deltaz  ! Change per height level        [m]

    ! If the CLUBB parameterization is implemented in a host model,
    ! it needs to use the host model's momentum level altitudes
    ! and thermodynamic level altitudes.
    ! If the CLUBB model is running by itself, but is using a
    ! stretched grid entered on thermodynamic levels (grid_type = 2),
    ! it needs to use the thermodynamic level altitudes as input.
    ! If the CLUBB model is running by itself, but is using a
    ! stretched grid entered on momentum levels (grid_type = 3),
    ! it needs to use the momentum level altitudes as input.
    real( kind = core_rknd ), intent(in), dimension(nzmax) :: &
      momentum_heights,      & ! Momentum level altitudes (input)      [m]
      thermodynamic_heights    ! Thermodynamic level altitudes (input) [m]

    ! Local Variables
    real( kind = core_rknd ) :: avg_deltaz  ! Average grid box height   [m]

    ! The factor by which to multiply the coefficients of background eddy
    ! diffusivity if the grid spacing threshold is exceeded and l_adj_low_res_nu
    ! is turned on.
    real( kind = core_rknd ),dimension(gr%nz) :: &
      mult_factor_zt, &  ! Uses gr%dzt for nu values on zt levels
      mult_factor_zm     ! Uses gr%dzm for nu values on zm levels

    ! Flag to enable nu values that are a function of grid spacing
    logical, parameter :: l_nu_grid_dependent = .false.

    integer :: k  ! Loop variable

    !--------------- Begin code -------------------------

    if ( .not. allocated( nu1_vert_res_dep ) ) then
      allocate( nu1_vert_res_dep(1:gr%nz) )
    end if
    if ( .not. allocated( nu2_vert_res_dep ) ) then
      allocate( nu2_vert_res_dep(1:gr%nz) )
    end if
    if ( .not. allocated( nu6_vert_res_dep ) ) then
      allocate( nu6_vert_res_dep(1:gr%nz) )
    end if
    if ( .not. allocated( nu8_vert_res_dep ) ) then
      allocate( nu8_vert_res_dep(1:gr%nz) )
    end if
    if ( .not. allocated( nu9_vert_res_dep ) ) then
      allocate( nu9_vert_res_dep(1:gr%nz) )
    end if
    if ( .not. allocated( nu10_vert_res_dep ) ) then
      allocate( nu10_vert_res_dep(1:gr%nz) )
    end if
    if ( .not. allocated( nu_hm_vert_res_dep ) ) then
      allocate( nu_hm_vert_res_dep(1:gr%nz) )
    end if

    ! Flag for adjusting the values of the constant diffusivity coefficients
    ! based on the grid spacing.  If this flag is turned off, the values of the
    ! various nu coefficients will remain as they are declared in the
    ! parameters.in file.
    if ( l_adj_low_res_nu ) then

      ! ### Adjust Constant Diffusivity Coefficients Based On Grid Spacing ###

      ! All of the background coefficients of eddy diffusivity, as well as the
      ! constant coefficient for 4th-order hyper-diffusion, must be adjusted
      ! based on the size of the grid spacing.  For a case that uses an
      ! evenly-spaced grid, the adjustment is based on the constant grid
      ! spacing deltaz.  For a case that uses a stretched grid, the adjustment
      ! is based on avg_deltaz, which is the average grid spacing over the
      ! vertical domain.
 
      if ( l_prescribed_avg_deltaz ) then
        
        avg_deltaz = deltaz

      else if ( grid_type == 3 ) then

        ! CLUBB is implemented in a host model, or is using grid_type = 3

        ! Find the average deltaz over the grid based on momentum level
        ! inputs.

        avg_deltaz  &
           = ( momentum_heights(nzmax) - momentum_heights(1) )  &
             / real( nzmax - 1, kind = core_rknd )

      else if ( grid_type == 1 ) then

        ! Evenly-spaced grid.

        avg_deltaz = deltaz

      else if ( grid_type == 2 ) then

        ! Stretched (unevenly-spaced) grid:  stretched thermodynamic level
        ! input.

        ! Find the average deltaz over the stretched grid based on
        ! thermodynamic level inputs.

        avg_deltaz  &
          = ( thermodynamic_heights(nzmax) - thermodynamic_heights(1) )  &
             / real( nzmax - 1, kind = core_rknd )
      else
        ! Eric Raut added to remove compiler warning. (Obviously, this value is not used)
        avg_deltaz = 0.0_core_rknd
        write(fstderr,*) "Invalid grid_type:", grid_type
        stop "Fatal error"

      end if ! grid_type

      ! The nu's are chosen for deltaz <= 40 m. Looks like they must
      ! be adjusted for larger grid spacings (Vince Larson)
      if( .not. l_nu_grid_dependent ) then
        ! Use a constant mult_factor so nu does not depend on grid spacing
        if( avg_deltaz > grid_spacing_thresh ) then
          mult_factor_zt = 1.0_core_rknd + mult_coef * log( avg_deltaz / grid_spacing_thresh )
          mult_factor_zm = mult_factor_zt
        else
          mult_factor_zt = 1.0_core_rknd
          mult_factor_zm = 1.0_core_rknd
        end if
      else  ! l_nu_grid_dependent = .true.
        ! mult_factor will vary to create nu values that vary with grid spacing
        do k = 1, gr%nz
          if( gr%dzm(k) > grid_spacing_thresh ) then
            mult_factor_zm(k) = 1.0_core_rknd + mult_coef * log( gr%dzm(k) / grid_spacing_thresh )
          else
            mult_factor_zm(k) = 1.0_core_rknd
          end if

          if( gr%dzt(k) > grid_spacing_thresh ) then
            mult_factor_zt(k) = 1.0_core_rknd + mult_coef * log( gr%dzt(k) / grid_spacing_thresh )
          else
            mult_factor_zt(k) = 1.0_core_rknd
          end if
        end do
      end if ! l_nu_grid_dependent

      !mult_factor = 1.0_core_rknd + mult_coef * log( avg_deltaz / grid_spacing_thresh )
      nu1_vert_res_dep   =  nu1 * mult_factor_zm
      nu2_vert_res_dep   =  nu2 * mult_factor_zm
      nu6_vert_res_dep   =  nu6 * mult_factor_zm
      nu8_vert_res_dep   =  nu8 * mult_factor_zt
      nu9_vert_res_dep   =  nu9 * mult_factor_zm
      nu10_vert_res_dep  =  nu10 * mult_factor_zt !We're unsure of the grid
      nu_hm_vert_res_dep =  nu_hm * mult_factor_zt

    else ! nu values are not adjusted

      nu1_vert_res_dep   =  nu1
      nu2_vert_res_dep   =  nu2
      nu6_vert_res_dep   =  nu6
      nu8_vert_res_dep   =  nu8
      nu9_vert_res_dep   =  nu9
      nu10_vert_res_dep  = nu10
      nu_hm_vert_res_dep = nu_hm

    end if  ! l_adj_low_res_nu

    return
  end subroutine adj_low_res_nu

  !=============================================================================
  subroutine read_parameters( iunit, filename, params )

    ! Description:
    ! Read a namelist containing the model parameters

    ! References:
    ! None
    !-----------------------------------------------------------------------
    use constants_clubb, only: fstderr ! Constant

    implicit none

    ! Input variables
    integer, intent(in) :: iunit

    character(len=*), intent(in) :: filename

    ! Output variables
    real( kind = core_rknd ), intent(out), dimension(nparams) :: params

    ! Local variables
    integer :: i

    logical :: l_error

    ! ---- Begin Code ----

    ! If the filename is empty, assume we're using a `working' set of
    ! parameters that are set statically here (handy for host models).
    ! Read the namelist
    if ( filename /= "" ) then
      ! Read the namelist
      open(unit=iunit, file=filename, status='old', action='read')

      read(unit=iunit, nml=clubb_params_nl)

      close(unit=iunit)

    end if

    if (clubb_C1 /= init_value) then
       C1 = clubb_C1
       C1b = C1
    end if
    ! if clubb_C2thl and clubb_C2rtthl not specified, continue to use C2thl=C2rt, C2rtthl = 1.3*C2rt
    ! to preserve existing compsets that have assumed so and only vary C2rt
    if (clubb_C2rt /= init_value) then
       C2rt = clubb_C2rt
       if (clubb_C2thl == init_value) C2thl = C2rt
       if (clubb_C2rtthl == init_value) C2rtthl = C2rt*1.3_core_rknd
    end if
    ! Allows C2thl and C2rtthl to vary separately
    if (clubb_C2thl /= init_value) C2thl = clubb_C2thl
    if (clubb_C2rtthl /= init_value) C2rtthl = clubb_C2rtthl
    if (clubb_C6rt /= init_value) then
       C6rt = clubb_C6rt
       C6thl = C6rt
    end if
    if (clubb_C6rtb /= init_value) C6rtb = clubb_C6rtb
    if (clubb_C7 /= init_value) C7 = clubb_C7
    if (clubb_C7b /= init_value) C7b = clubb_C7b
    if (clubb_C8 /= init_value) C8 = clubb_C8
    if (clubb_C11 /= init_value) C11 = clubb_C11
    if (clubb_C11b /= init_value) C11b = clubb_C11b
    if (clubb_C14 /= init_value) C14 = clubb_C14
    if (clubb_beta /= init_value) beta = clubb_beta
    ! if clubb_gamma_coefb not specified, continue to use gamma_coefb=gamma_coef
    ! to preserve existing compsets that have assumed so  and only vary gamma_coef
    if (clubb_gamma_coef /= init_value) then
       gamma_coef = clubb_gamma_coef
       if (clubb_gamma_coefb == init_value) gamma_coefb = gamma_coef
    end if
    ! Allows gamma_coefb to vary separately
    if (clubb_gamma_coefb /= init_value) gamma_coefb = clubb_gamma_coefb
    if (clubb_mu /= init_value) mu = clubb_mu
    if (clubb_nu1 /= init_value) nu1 = clubb_nu1
    if (clubb_c_K10 /= init_value) c_K10 = clubb_c_K10
    if (clubb_wpxp_L_thresh /= init_value)wpxp_L_thresh = clubb_wpxp_L_thresh

    ! Put the variables in the output array
    call pack_parameters( C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl, &
                          C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, &
                          C7, C7b, C7c, C8, C8b, C10, &
                          C11, C11b, C11c, C12, C13, C14, C15, &
                          C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
                          c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6,  &
                          c_K8, nu8, c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, &
                          nu_hm, gamma_coef, gamma_coefb, gamma_coefc, &
                          mu, beta, lmin_coef, omicron, zeta_vrnce_rat, &
                          upsilon_precip_frac_rat, lambda0_stability_coef, &
                          mult_coef, taumin, taumax, Lscale_mu_coef, Lscale_pert_coef, &
                          alpha_corr, Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef, &
                          thlp2_rad_cloud_frac_thresh, up2_vp2_factor, params )

    l_error = .false.

    do i = 1, nparams
      if ( params(i) == init_value ) then
        write(fstderr,*) "Tuning parameter "//trim( params_list(i) )// &
          " was missing from "//trim( filename )
        l_error = .true.
      end if
    end do

    if ( l_error ) stop "Fatal error."

    return

  end subroutine read_parameters

  !=============================================================================
  subroutine read_param_spread & 
           ( iunit, filename, nindex, param_spread, ndim )

    ! Description:
    ! Read a namelist containing the amount to vary model parameters.
    ! Used by the downhill simplex / simulated annealing algorithm.

    ! References:
    ! None
    !-----------------------------------------------------------------------
    use constants_clubb, only: fstderr ! Constant

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input variables
    integer, intent(in) :: iunit

    character(len=*), intent(in) :: filename

    ! Output variables

    ! An array of array indices (i.e. which elements of the array `params'
    ! are contained within the simplex and the spread variable)
    integer, intent(out), dimension(nparams) :: nindex

    real( kind = core_rknd ), intent(out), dimension(nparams) ::  & 
      param_spread  ! Amount to vary the parameter in the initial simplex

    integer, intent(out) :: &
        ndim  ! Number of variables, e.g. rcm, to be tuned. Dimension of the init simplex


    ! Local variables
    integer :: i

    logical :: l_error

    ! Amount to change each parameter for the initial simplex
    ! This MUST be changed to match the clubb_params_nl namelist if parameters are added!
    namelist /initspread/  & 
      C1, C1b, C1c, C2, C2b, C2c,  & 
      C2rt, C2thl, C2rtthl, C4, C5, & 
      C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, & 
      C7, C7b, C7c, C8, C8b, C10, C11, C11b, C11c, & 
      C12, C13, C14, C15, C6rt_Lscale0, C6thl_Lscale0, &
      C7_Lscale0, wpxp_L_thresh, c_K, c_K1, nu1, c_K2, nu2,  & 
      c_K6, nu6, c_K8, nu8, c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, &
      nu_hm, beta, gamma_coef, gamma_coefb, gamma_coefc, & 
      lmin_coef, omicron, zeta_vrnce_rat, upsilon_precip_frac_rat, &
      lambda0_stability_coef, mult_coef, taumin, taumax, mu, &
      Lscale_mu_coef, Lscale_pert_coef, alpha_corr, Skw_denom_coef, c_K10, c_K10h, &
      thlp2_rad_coef, thlp2_rad_cloud_frac_thresh, up2_vp2_factor

    ! Initialize values to -999.
    call init_parameters_999( )

    ! Read the namelist
    open(unit=iunit, file=filename, status='old', action='read')

    read(unit=iunit, nml=initspread)

    close(unit=iunit)

    ! Put the variables in the output array
    call pack_parameters( C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl, &
                          C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, &
                          C7, C7b, C7c, C8, C8b, C10, &
                          C11, C11b, C11c, C12, C13, C14, C15, &
                          C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
                          c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6,  &
                          c_K8, nu8, c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, &
                          nu_hm, gamma_coef, gamma_coefb, gamma_coefc, &
                          mu, beta, lmin_coef, omicron, zeta_vrnce_rat, &
                          upsilon_precip_frac_rat, lambda0_stability_coef, &
                          mult_coef, taumin, taumax, Lscale_mu_coef, Lscale_pert_coef, &
                          alpha_corr, Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef, &
                          thlp2_rad_cloud_frac_thresh, up2_vp2_factor, param_spread )

    l_error = .false.

    do i = 1, nparams
      if ( param_spread(i) == init_value ) then
        write(fstderr,*) "A spread parameter "//trim( params_list(i) )// &
          " was missing from "//trim( filename )
        l_error = .true.
      end if
    end do

    if ( l_error ) stop "Fatal error."

    ! Initialize to zero
    nindex(1:nparams) = 0
    ndim = 0

    ! Determine how many variables are being changed
    do i = 1, nparams, 1

      if ( param_spread(i) /= 0.0_core_rknd ) then
        ndim = ndim + 1   ! Increase the total
        nindex(ndim) = i  ! Set the next array index
      endif

    enddo

    return

  end subroutine read_param_spread

  !=============================================================================
  subroutine pack_parameters &
             ( C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl, &
               C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, &
               C7, C7b, C7c, C8, C8b, C10, &
               C11, C11b, C11c, C12, C13, C14, C15, &
               C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
               c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6,  &
               c_K8, nu8, c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, &
               nu_hm, gamma_coef, gamma_coefb, gamma_coefc, &
               mu, beta, lmin_coef, omicron, zeta_vrnce_rat, &
               upsilon_precip_frac_rat, lambda0_stability_coef, &
               mult_coef, taumin, taumax, Lscale_mu_coef, Lscale_pert_coef, &
               alpha_corr, Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef, &
               thlp2_rad_cloud_frac_thresh, up2_vp2_factor, params )

    ! Description:
    ! Takes the list of scalar variables and puts them into a 1D vector.
    ! It is here for the purpose of keeping the code generalized
    ! when new variables are added.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use parameter_indices, only: & 
      iC1,  & ! Variable(s)
      iC1b, & 
      iC1c, & 
      iC2, & 
      iC2b, & 
      iC2c, & 
      iC2rt, & 
      iC2thl, & 
      iC2rtthl, & 
      iC4, & 
      iC5, & 
      iC6rt, & 
      iC6rtb, & 
      iC6rtc, & 
      iC6thl, & 
      iC6thlb, & 
      iC6thlc, & 
      iC7, & 
      iC7b, & 
      iC7c, & 
      iC8, & 
      iC8b, & 
      iC10, & 
      iC11, & 
      iC11b, & 
      iC11c, & 
      iC12, & 
      iC13, & 
      iC14, &
      iC15

    use parameter_indices, only: &
      iC6rt_Lscale0, &
      iC6thl_Lscale0, &
      iC7_Lscale0, &
      iwpxp_L_thresh

    use parameter_indices, only: & 
      ic_K,  & 
      ic_K1, & 
      inu1, & 
      ic_K2, & 
      inu2, & 
      ic_K6, & 
      inu6, & 
      ic_K8, & 
      inu8, & 
      ic_K9, & 
      inu9, & 
      inu10, &
      ic_K_hm, & 
      ic_K_hmb, & 
      iK_hm_min_coef, &
      inu_hm, & 
      igamma_coef, & 
      igamma_coefb, & 
      igamma_coefc, & 
      imu, & 
      ibeta, & 
      ilmin_coef, &
      iomicron, &
      izeta_vrnce_rat, &
      iupsilon_precip_frac_rat, &
      ilambda0_stability_coef, &
      imult_coef, &
      itaumin, & 
      itaumax, & 
      iLscale_mu_coef, &
      iLscale_pert_coef, &
      ialpha_corr, &
      iSkw_denom_coef, &
      ic_K10, &
      ic_K10h, &
      ithlp2_rad_coef, &
      ithlp2_rad_cloud_frac_thresh, &
      iup2_vp2_factor, &
      nparams

    implicit none

    ! Input variables
    real( kind = core_rknd ), intent(in) :: & 
      C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl, & 
      C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, & 
      C7, C7b, C7c, C8, C8b, C10, & 
      C11, C11b, C11c, C12, C13, C14, C15, & 
      C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
      c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, c_K8, nu8,  & 
      c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm, gamma_coef, &
      gamma_coefb, gamma_coefc, mu, beta, lmin_coef, &
      omicron, zeta_vrnce_rat, upsilon_precip_frac_rat, &
      lambda0_stability_coef, mult_coef, taumin, taumax, Lscale_mu_coef, &
      Lscale_pert_coef, alpha_corr, Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef, &
      thlp2_rad_cloud_frac_thresh, up2_vp2_factor

    ! Output variables
    real( kind = core_rknd ), intent(out), dimension(nparams) :: params

    params(iC1)      = C1
    params(iC1b)     = C1b
    params(iC1c)     = C1c
    params(iC2)      = C2
    params(iC2b)     = C2b
    params(iC2c)     = C2c
    params(iC2rt)    = C2rt
    params(iC2thl)   = C2thl
    params(iC2rtthl) = C2rtthl
    params(iC4)      = C4
    params(iC5)      = C5
    params(iC6rt)    = C6rt
    params(iC6rtb)   = C6rtb
    params(iC6rtc)   = C6rtc
    params(iC6thl)   = C6thl
    params(iC6thlb)  = C6thlb
    params(iC6thlc)  = C6thlc
    params(iC7)      = C7
    params(iC7b)     = C7b
    params(iC7c)     = C7c
    params(iC8)      = C8
    params(iC8b)     = C8b
    params(iC10)     = C10
    params(iC11)     = C11
    params(iC11b)    = C11b
    params(iC11c)    = C11c
    params(iC12)     = C12
    params(iC13)     = C13
    params(iC14)     = C14
    params(iC15)     = C15

    params(iC6rt_Lscale0)       = C6rt_Lscale0
    params(iC6thl_Lscale0)      = C6thl_Lscale0
    params(iC7_Lscale0)         = C7_Lscale0
    params(iwpxp_L_thresh)    = wpxp_L_thresh

    params(ic_K)       = c_K
    params(ic_K1)      = c_K1
    params(inu1)       = nu1
    params(ic_K2)      = c_K2
    params(inu2)       = nu2
    params(ic_K6)      = c_K6
    params(inu6)       = nu6
    params(ic_K8)      = c_K8
    params(inu8)       = nu8
    params(ic_K9)      = c_K9
    params(inu9)       = nu9
    params(inu10)      = nu10
    params(ic_K_hm)    = c_K_hm
    params(ic_K_hmb)   = c_K_hmb
    params(iK_hm_min_coef)   = K_hm_min_coef
    params(inu_hm)     = nu_hm

    params(igamma_coef)  = gamma_coef
    params(igamma_coefb) = gamma_coefb
    params(igamma_coefc) = gamma_coefc

    params(imu) = mu

    params(ibeta) = beta

    params(ilmin_coef) = lmin_coef

    params(iomicron) = omicron
    params(izeta_vrnce_rat) = zeta_vrnce_rat

    params(iupsilon_precip_frac_rat) = upsilon_precip_frac_rat
    params(ilambda0_stability_coef) = lambda0_stability_coef
    params(imult_coef) = mult_coef

    params(itaumin) = taumin
    params(itaumax) = taumax

    params(iLscale_mu_coef) = Lscale_mu_coef
    params(iLscale_pert_coef) = Lscale_pert_coef
    params(ialpha_corr) = alpha_corr
    params(iSkw_denom_coef) = Skw_denom_coef
    params(ic_K10) = c_K10
    params(ic_K10h) = c_K10h
    params(ithlp2_rad_coef) = thlp2_rad_coef
    params(ithlp2_rad_cloud_frac_thresh) = thlp2_rad_cloud_frac_thresh
    params(iup2_vp2_factor) = up2_vp2_factor

    return
  end subroutine pack_parameters

  !=============================================================================
  subroutine unpack_parameters & 
             ( params, & 
               C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl, & 
               C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, & 
               C7, C7b, C7c, C8, C8b, C10, & 
               C11, C11b, C11c, C12, C13, C14, C15, & 
               C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
               c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, & 
               c_K8, nu8, c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, &
               nu_hm, gamma_coef, gamma_coefb, gamma_coefc, & 
               mu, beta, lmin_coef, omicron, zeta_vrnce_rat, &
               upsilon_precip_frac_rat, lambda0_stability_coef, &
               mult_coef, taumin, taumax, Lscale_mu_coef, Lscale_pert_coef, &
               alpha_corr, Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef, &
               thlp2_rad_cloud_frac_thresh, up2_vp2_factor )

    ! Description:
    ! Takes the 1D vector and returns the list of scalar variables.
    ! Here for the purposes of keeping the code generalized
    ! when new variables are added.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use parameter_indices, only: & 
      iC1,  & ! Variable(s)
      iC1b, & 
      iC1c, & 
      iC2, & 
      iC2b, & 
      iC2c, & 
      iC2rt, & 
      iC2thl, & 
      iC2rtthl, & 
      iC4, & 
      iC5, & 
      iC6rt, & 
      iC6rtb, & 
      iC6rtc, & 
      iC6thl, & 
      iC6thlb, & 
      iC6thlc, & 
      iC7, & 
      iC7b, & 
      iC7c, & 
      iC8, & 
      iC8b, & 
      iC10, & 
      iC11, & 
      iC11b, & 
      iC11c, & 
      iC12, & 
      iC13, & 
      iC14, &
      iC15

    use parameter_indices, only: &
      iC6rt_Lscale0, &
      iC6thl_Lscale0, &
      iC7_Lscale0, &
      iwpxp_L_thresh

    use parameter_indices, only: & 
      ic_K,  & 
      ic_K1, & 
      inu1, & 
      ic_K2, & 
      inu2, & 
      ic_K6, & 
      inu6, & 
      ic_K8, & 
      inu8, & 
      ic_K9, & 
      inu9, & 
      inu10, &
      ic_K_hm, & 
      ic_K_hmb, & 
      iK_hm_min_coef, & 
      inu_hm, & 
      igamma_coef, & 
      igamma_coefb, & 
      igamma_coefc, & 
      imu, & 
      ibeta, & 
      ilmin_coef, &
      iomicron, &
      izeta_vrnce_rat, &
      iupsilon_precip_frac_rat, &
      ilambda0_stability_coef, &
      imult_coef, &
      itaumin, & 
      itaumax, & 
      iLscale_mu_coef, &
      iLscale_pert_coef, &
      ialpha_corr, &
      iSkw_denom_coef, &
      ic_K10, &
      ic_K10h, & 
      ithlp2_rad_coef, &
      ithlp2_rad_cloud_frac_thresh, &
      iup2_vp2_factor, &
      nparams

    implicit none

    ! Input variables
    real( kind = core_rknd ), intent(in), dimension(nparams) :: params

    ! Output variables
    real( kind = core_rknd ), intent(out) :: & 
      C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl, & 
      C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, & 
      C7, C7b, C7c, C8, C8b, C10, & 
      C11, C11b, C11c, C12, C13, C14, C15, & 
      C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
      c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, & 
      c_K8, nu8, c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm, & 
      gamma_coef, gamma_coefb, gamma_coefc, & 
      mu, beta, lmin_coef, omicron, zeta_vrnce_rat, upsilon_precip_frac_rat, &
      lambda0_stability_coef, mult_coef, taumin, &
      taumax, Lscale_mu_coef, Lscale_pert_coef, alpha_corr, Skw_denom_coef, c_K10, &
      c_K10h, thlp2_rad_coef, thlp2_rad_cloud_frac_thresh, up2_vp2_factor

    C1      = params(iC1)
    C1b     = params(iC1b)
    C1c     = params(iC1c)
    C2      = params(iC2)
    C2b     = params(iC2b)
    C2c     = params(iC2c)
    C2rt    = params(iC2rt)
    C2thl   = params(iC2thl)
    C2rtthl = params(iC2rtthl)
    C4      = params(iC4)
    C5      = params(iC5)
    C6rt    = params(iC6rt)
    C6rtb   = params(iC6rtb)
    C6rtc   = params(iC6rtc)
    C6thl   = params(iC6thl)
    C6thlb  = params(iC6thlb)
    C6thlc  = params(iC6thlc)
    C7      = params(iC7)
    C7b     = params(iC7b)
    C7c     = params(iC7c)
    C8      = params(iC8)
    C8b     = params(iC8b)
    C10     = params(iC10)
    C11     = params(iC11)
    C11b    = params(iC11b)
    C11c    = params(iC11c)
    C12     = params(iC12)
    C13     = params(iC13)
    C14     = params(iC14)
    C15     = params(iC15)

    C6rt_Lscale0       = params(iC6rt_Lscale0)
    C6thl_Lscale0      = params(iC6thl_Lscale0)
    C7_Lscale0         = params(iC7_Lscale0)
    wpxp_L_thresh    = params(iwpxp_L_thresh)

    c_K       = params(ic_K)
    c_K1      = params(ic_K1)
    nu1       = params(inu1)
    c_K2      = params(ic_K2)
    nu2       = params(inu2)
    c_K6      = params(ic_K6)
    nu6       = params(inu6)
    c_K8      = params(ic_K8)
    nu8       = params(inu8)
    c_K9      = params(ic_K9)
    nu9       = params(inu9)
    nu10      = params(inu10)
    c_K_hm    = params(ic_K_hm)
    c_K_hmb   = params(ic_K_hmb)
    K_hm_min_coef   = params(iK_hm_min_coef)
    nu_hm     = params(inu_hm)

    gamma_coef  = params(igamma_coef)
    gamma_coefb = params(igamma_coefb)
    gamma_coefc = params(igamma_coefc)

    mu = params(imu)

    beta = params(ibeta)

    lmin_coef = params(ilmin_coef)

    omicron = params(iomicron)
    zeta_vrnce_rat = params(izeta_vrnce_rat)

    upsilon_precip_frac_rat = params(iupsilon_precip_frac_rat)
    lambda0_stability_coef = params(ilambda0_stability_coef)
    mult_coef = params(imult_coef)

    taumin = params(itaumin)
    taumax = params(itaumax)

    Lscale_mu_coef = params(iLscale_mu_coef)
    Lscale_pert_coef = params(iLscale_pert_coef)
    alpha_corr = params(ialpha_corr)
    Skw_denom_coef = params(iSkw_denom_coef)
    c_K10 = params(ic_K10)
    c_K10h = params(ic_K10h)

    thlp2_rad_coef = params(ithlp2_rad_coef)
    thlp2_rad_cloud_frac_thresh = params(ithlp2_rad_cloud_frac_thresh)
    up2_vp2_factor = params(iup2_vp2_factor)

    return
  end subroutine unpack_parameters

  !=============================================================================
  subroutine get_parameters( params )

    ! Description:
    ! Return an array of all tunable parameters

    ! References:
    ! None
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(out), dimension(nparams) :: params

    call pack_parameters( C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl, &
                          C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, &
                          C7, C7b, C7c, C8, C8b, C10, &
                          C11, C11b, C11c, C12, C13, C14, C15, &
                          C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
                          c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6,  &
                          c_K8, nu8, c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, &
                          nu_hm, gamma_coef, gamma_coefb, gamma_coefc, &
                          mu, beta, lmin_coef, omicron, zeta_vrnce_rat, &
                          upsilon_precip_frac_rat, lambda0_stability_coef, &
                          mult_coef, taumin, taumax, Lscale_mu_coef, Lscale_pert_coef, &
                          alpha_corr, Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef, &
                          thlp2_rad_cloud_frac_thresh, up2_vp2_factor, params )

    return

  end subroutine get_parameters

  !=============================================================================
  subroutine init_parameters_999( )

    ! Description:
    ! Set all tunable parameters to NaN

    ! References:
    ! None
    !-----------------------------------------------------------------------

    implicit none

    ! --- Begin Code ---

    C1                          = init_value
    C1b                         = init_value
    C1c                         = init_value
    C2rt                        = init_value
    C2thl                       = init_value
    C2rtthl                     = init_value
    C2                          = init_value
    C2b                         = init_value
    C2c                         = init_value
    C4                          = init_value
    C5                          = init_value
    C6rt                        = init_value
    C6rtb                       = init_value
    C6rtc                       = init_value
    C6thl                       = init_value
    C6thlb                      = init_value
    C6thlc                      = init_value
    C7                          = init_value
    C7b                         = init_value
    C7c                         = init_value
    C8                          = init_value
    C8b                         = init_value
    C10                         = init_value
    C11                         = init_value
    C11b                        = init_value
    C11c                        = init_value
    C12                         = init_value
    C13                         = init_value
    C14                         = init_value
    C15                         = init_value
    C6rt_Lscale0                = init_value
    C6thl_Lscale0               = init_value
    C7_Lscale0                  = init_value
    wpxp_L_thresh               = init_value
    c_K                         = init_value
    c_K1                        = init_value
    nu1                         = init_value
    c_K2                        = init_value
    nu2                         = init_value
    c_K6                        = init_value
    nu6                         = init_value
    c_K8                        = init_value
    nu8                         = init_value
    c_K9                        = init_value
    nu9                         = init_value
    nu10                        = init_value
    c_K_hm                      = init_value
    c_K_hmb                     = init_value
    K_hm_min_coef               = init_value
    nu_hm                       = init_value
    beta                        = init_value
    gamma_coef                  = init_value
    gamma_coefb                 = init_value
    gamma_coefc                 = init_value
    mult_coef                   = init_value
    taumin                      = init_value
    taumax                      = init_value
    lmin_coef                   = init_value
    omicron                     = init_value
    zeta_vrnce_rat              = init_value
    upsilon_precip_frac_rat     = init_value
    lambda0_stability_coef      = init_value
    mu                          = init_value
    Lscale_mu_coef              = init_value
    Lscale_pert_coef            = init_value
    alpha_corr                  = init_value
    Skw_denom_coef              = init_value
    c_K10                       = init_value
    c_K10h                      = init_value
    thlp2_rad_coef              = init_value
    thlp2_rad_cloud_frac_thresh = init_value
    up2_vp2_factor              = init_value

    return

  end subroutine init_parameters_999

  !=============================================================================
  subroutine cleanup_nu( )

    ! Description:
    !  De-allocates memory used for the nu arrays
    !
    ! References:
    !  None
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
      fstderr  ! Constant

    implicit none

    ! Local Variable(s)
    integer :: ierr

    ! ----- Begin Code -----

    deallocate( nu1_vert_res_dep, nu2_vert_res_dep, nu6_vert_res_dep,  &
                nu8_vert_res_dep, nu9_vert_res_dep, nu10_vert_res_dep, &
                nu_hm_vert_res_dep, stat = ierr )

    if ( ierr /= 0 ) then
      write(fstderr,*) "Deallocation of vertically depedent nu arrays failed."
    end if

    return

  end subroutine cleanup_nu

!===============================================================================

end module parameters_tunable
