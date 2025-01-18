module clubb_intr

  !----------------------------------------------------------------------------------------------------- !
  ! Module to interface CAM with Cloud Layers Unified by Bi-normals (CLUBB), developed                   !
  !    by the University of Wisconsin Milwaukee Group (UWM).                                             !
  !                                                                                                      !
  ! CLUBB replaces the exisiting turbulence, shallow convection, and macrophysics in CAM5                !
  !                                                                                                      !
  ! Lastly, a implicit diffusion solver is called, and tendencies retrieved by                           !
  ! differencing the diffused and initial states.                                                        !
  !                                                                                                      !
  ! Calling sequence:                                                                                    !
  !                                                                                                      !
  !---------------------------Code history-------------------------------------------------------------- !
  ! Authors:  P. Bogenschutz, C. Craig, A. Gettelman                                                     !
  !                                                                                                      !
  !----------------------------------------------------------------------------------------------------- !
  !  2020-01  O. Guba Correct energy density function
  !-----------------------------------------------------------------------------------------------------
  use shr_kind_mod,  only: r8=>shr_kind_r8
  use shr_log_mod ,  only: errMsg => shr_log_errMsg
  use ppgrid,        only: pver, pverp
  use phys_control,  only: phys_getopts, use_od_ss, use_od_fd
  use od_common,     only: od_ls_ncleff, od_bl_ncd, od_ss_sncleff
  use physconst,     only: rair, cpair, gravit, latvap, latice, zvir, rh2o, karman, &
                           tms_orocnst, tms_z0fac, pi
  use cam_logfile,   only: iulog
  use spmd_utils,    only: masterproc
  use constituents,  only: pcnst, cnst_add
  use pbl_utils,     only: calc_ustar, calc_obklen
  use perf_mod,      only: t_startf, t_stopf
  use mpishorthand
  use cam_history_support, only: fillvalue
#ifdef CLUBB_SGS
  use clubb_api_module, only: pdf_parameter, clubb_fatal_error, fstderr
  use clubb_precision,  only: core_rknd
#else
  use shr_kind_mod,     only: core_rknd=>shr_kind_r8
#endif

  use zm_conv,      only: zm_microp

  implicit none

  private
  save

  ! ----------------- !
  ! Public interfaces !
  ! ----------------- !

  public :: clubb_ini_cam, clubb_register_cam, clubb_tend_cam, &
#ifdef CLUBB_SGS
            ! This utilizes CLUBB specific variables in its interface
            stats_init_clubb, &
#endif
            stats_end_timestep_clubb, &
            clubb_surface, &
            clubb_readnl, &
            clubb_init_cnst, &
            clubb_implements_cnst

#ifdef CLUBB_SGS
  ! Both of these utilize CLUBB specific variables in their interface
  private :: stats_zero, stats_avg
#endif

  logical, public :: do_cldcool

  ! ------------ !
  ! Private data !
  ! ------------ !

  integer, parameter :: &
      grid_type    = 3, &               ! The 2 option specifies stretched thermodynamic levels
      hydromet_dim = 0                  ! The hydromet array in SAM-CLUBB is currently 0 elements

  real(core_rknd), dimension(0) :: &
      sclr_tol = 1.e-8_core_rknd        ! Total water in kg/kg

  character(len=6), parameter :: &
      saturation_equation = "flatau"    ! Flatau polynomial approximation for SVP

  real(core_rknd), parameter :: &
      theta0   = 300._core_rknd, &      ! Reference temperature                     [K]
      ts_nudge = 86400._core_rknd, &    ! Time scale for u/v nudging (not used)     [s]
      p0_clubb = 100000._core_rknd

  real(core_rknd), parameter :: &
      host_dx = 100000._core_rknd, &    ! Host model deltax [m]
      host_dy = 100000._core_rknd       ! Host model deltay [m]

  integer, parameter :: &
    sclr_dim = 0                        ! Higher-order scalars, set to zero

  real(r8), parameter :: &
    wp3_const = 1._r8                   ! Constant to add to wp3 when moments are advected

  real(r8), parameter :: &
    wpthlp_const = 10.0_r8              ! Constant to add to wpthlp when moments are advected

  real(r8), parameter :: &
    wprtp_const = 0.01_r8               ! Constant to add to wprtp when moments are advected

  real(r8), parameter :: &
    rtpthlp_const = 0.01_r8             ! Constant to add to rtpthlp when moments are advected

  real(r8), parameter :: unset_r8 = huge(1.0_r8)

  !PMA
  real(r8), parameter :: qsmall = 1.e-18_r8 ! qsmall used in MG

  real(r8), parameter ::  rad_to_deg = 180.0_r8/pi !converts radians to degrees

  real(r8) :: clubb_timestep = unset_r8  ! Default CLUBB timestep, unless overwriten by namelist
  real(r8) :: clubb_rnevap_effic = unset_r8

  !namelist variables
  real(r8) :: clubb_liq_deep = unset_r8
  real(r8) :: clubb_liq_sh   = unset_r8
  real(r8) :: clubb_ice_deep = unset_r8
  real(r8) :: clubb_ice_sh   = unset_r8
  real(r8) :: clubb_tk1      = unset_r8
  real(r8) :: clubb_tk2      = unset_r8

!  Constant parameters
  logical, parameter, private :: &
    l_uv_nudge       = .false.,       &  ! Use u/v nudging (not used)
    l_input_fields   = .false.,       &  ! Always false for EAM-CLUBB.
    l_implemented    = .true.,        &  ! Implemented in a host model (always true)
    l_host_applies_sfc_fluxes = .false.  ! Whether the host model applies the surface fluxes

  logical            :: do_tms
  logical            :: linearize_pbl_winds
  logical            :: lq(pcnst)
  logical            :: lq2(pcnst)
  logical            :: prog_modal_aero
  logical            :: do_rainturb
  logical            :: do_expldiff
  logical            :: clubb_do_adv
  logical            :: clubb_do_deep
  logical            :: micro_do_icesupersat
  logical            :: history_budget
  logical            :: use_sgv !PMA This flag controls tuning for tpert and gustiness
  integer            :: history_budget_histfile_num
  integer            :: edsclr_dim       ! Number of scalars to transport in CLUBB
  integer            :: offset

!  define physics buffer indicies here
  integer :: &
    wp2_idx, &          ! vertical velocity variances
    wp3_idx, &          ! third moment of vertical velocity
    wpthlp_idx, &       ! turbulent flux of thetal
    wprtp_idx, &        ! turbulent flux of total water
    rtpthlp_idx, &      ! covariance of thetal and rt
    rtp2_idx, &         ! variance of total water
    thlp2_idx, &        ! variance of thetal
    up2_idx, &          ! variance of east-west wind
    vp2_idx, &          ! variance of north-south wind
    upwp_idx, &         ! east-west momentum flux
    vpwp_idx, &         ! north-south momentum flux
    um_pert_idx, &      ! perturbed east-west momentum flux
    vm_pert_idx, &      ! perturbed north-south momentum flux
    upwp_pert_idx, &    ! perturbed east-west momentum flux
    vpwp_pert_idx, &    ! perturbed north-south momentum flux
    thlm_idx, &         ! mean thetal
    rtm_idx, &          ! mean total water mixing ratio
    um_idx, &           ! mean of east-west wind
    vm_idx, &           ! mean of north-south wind
    cld_idx, &          ! Cloud fraction
    concld_idx, &       ! Convective cloud fraction
    ast_idx, &          ! Stratiform cloud fraction
    alst_idx, &         ! Liquid stratiform cloud fraction
    aist_idx, &         ! Ice stratiform cloud fraction
    qlst_idx, &         ! Physical in-cloud LWC
    qist_idx, &         ! Physical in-cloud IWC
    dp_frac_idx, &      ! deep convection cloud fraction
    sh_frac_idx, &      ! shallow convection cloud fraction
    rel_idx, &          ! Rel
    kvh_idx, &          ! CLUBB eddy diffusivity on thermo levels
    kvm_idx, &          ! CLUBB eddy diffusivity on mom levels
    pblh_idx, &         ! PBL pbuf
    icwmrdp_idx, &      ! In cloud mixing ratio for deep convection
    tke_idx, &          ! turbulent kinetic energy
    tpert_idx, &        ! temperature perturbation from PBL
    fice_idx, &         ! fice_idx index in physics buffer
    cmeliq_idx, &       ! cmeliq_idx index in physics buffer
    relvar_idx, &       ! relative cloud water variance
    accre_enhan_idx, &  ! optional accretion enhancement factor for MG
    naai_idx, &         ! ice number concentration
    prer_evap_idx, &    ! rain evaporation rate
    qrl_idx, &          ! longwave cooling rate
    radf_idx

 integer :: &          ! newly added pbuf fields for CLUBB
    wpthvp_idx, &       ! < w'th_v' >
    wp2thvp_idx, &      ! < w'^2 th_v' >
    rtpthvp_idx, &      ! < r_t'th_v' >
    thlpthvp_idx, &     ! < th_l'th_v' >
    rcm_idx, &          ! Cloud water mixing ratio
    cloud_frac_idx !, & ! Cloud fraction

  integer :: &          ! added pbuf fields for clubb to have restart bfb when ipdf_call_placement=2
    pdf_zm_w_1_idx, &
    pdf_zm_w_2_idx, &
    pdf_zm_varnce_w_1_idx, &
    pdf_zm_varnce_w_2_idx, &
    pdf_zm_mixt_frac_idx

  integer :: &          !PMA adds pbuf fields for ZM gustiness
    prec_dp_idx, &
    snow_dp_idx, &
    vmag_gust_idx, &
    wsresp_idx, &
    tau_est_idx

  integer, public :: &
    ixthlp2 = 0, &
    ixwpthlp = 0, &
    ixwprtp = 0, &
    ixwp2 = 0, &
    ixwp3 = 0, &
    ixrtpthlp = 0, &
    ixrtp2 = 0, &
    ixup2 = 0, &
    ixvp2 = 0

  integer :: cmfmc_sh_idx = 0

! ZM convective microphysics(zm_microp)
  integer :: &
    dlfzm_idx  = -1,    & ! index of ZM detrainment of convective cloud water mixing ratio.
    difzm_idx  = -1,    & ! index of ZM detrainment of convective cloud ice mixing ratio.
    dsfzm_idx  = -1,    & ! index of ZM detrainment of convective snow mixing ratio.
    dnlfzm_idx = -1,    & ! index of ZM detrainment of convective cloud water num concen.
    dnifzm_idx = -1,    & ! index of ZM detrainment of convective cloud ice num concen.
    dnsfzm_idx = -1       ! index of ZM detrainment of convective snow num concen.

  real(r8) :: dp1 !set in namelist; assigned in cloud_fraction.F90
  !  Output arrays for CLUBB statistics
  real(r8), allocatable, dimension(:,:,:) :: out_zt, out_zm, out_radzt, out_radzm, out_sfc

  character(len=16)  :: eddy_scheme      ! Default set in phys_control.F90
  character(len=16)  :: deep_scheme      ! Default set in phys_control.F90

  integer, parameter :: ncnst=9
  character(len=8)   :: cnst_names(ncnst)
  logical            :: do_cnst=.false.

#ifdef CLUBB_SGS
  type(pdf_parameter), target, allocatable :: pdf_params_chnk(:,:)    ! PDF parameters (thermo. levs.) [units vary]
  type(pdf_parameter), target, allocatable :: pdf_params_zm_chnk(:,:) ! PDF parameters on momentum levs. [units vary]
#endif

  logical :: liqcf_fix = .FALSE.  ! HW for liquid cloud fraction fix
  logical :: relvar_fix = .FALSE. !PMA for relvar fix

  real(r8) :: micro_mg_accre_enhan_fac = huge(1.0_r8) !Accretion enhancement factor from namelist

  contains

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine clubb_register_cam( )
!-------------------------------------------------------------------------------
! Description:
!   Register the constituents and fields in the physics buffer
! Author: P. Bogenschutz, C. Craig, A. Gettelman
!
!-------------------------------------------------------------------------------
#ifdef CLUBB_SGS

    !------------------------------------------------ !
    ! Register physics buffer fields and constituents !
    !------------------------------------------------ !

    !  Add CLUBB fields to pbuf
    use physics_buffer,  only: pbuf_add_field, dtype_r8, dyn_time_lvls
    use ppgrid,          only: pver, pverp, pcols

    call phys_getopts( eddy_scheme_out                 = eddy_scheme, &
                       deep_scheme_out                 = deep_scheme, &
                       do_tms_out                      = do_tms,      &
                       linearize_pbl_winds_out         = linearize_pbl_winds,      &
                       history_budget_out              = history_budget, &
                       history_budget_histfile_num_out = history_budget_histfile_num, &
                       micro_do_icesupersat_out        = micro_do_icesupersat, &
                       micro_mg_accre_enhan_fac_out    = micro_mg_accre_enhan_fac)

    if (clubb_do_adv) then
       cnst_names =(/'THLP2  ','RTP2   ','RTPTHLP','WPTHLP ','WPRTP  ','WP2    ','WP3    ','UP2    ','VP2    '/)
       do_cnst=.true.
       !  If CLUBB moments are advected, do not output them automatically which is typically done.  Some moments
       !    need a constant added to them before they are advected, thus this would corrupt the output.
       !    Users should refer to the "XXXX_CLUBB" (THLP2_CLUBB for instance) output variables for these moments
       call cnst_add(trim(cnst_names(1)),0._r8,0._r8,0._r8,ixthlp2,longname='second moment vertical velocity',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(2)),0._r8,0._r8,0._r8,ixrtp2,longname='second moment rtp',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(3)),0._r8,0._r8,-999999._r8,ixrtpthlp,longname='covariance rtp thlp',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(4)),0._r8,0._r8,-999999._r8,ixwpthlp,longname='CLUBB heat flux',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(5)),0._r8,0._r8,-999999._r8,ixwprtp,longname='CLUBB moisture flux',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(6)),0._r8,0._r8,0._r8,ixwp2,longname='CLUBB wp2',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(7)),0._r8,0._r8,-999999._r8,ixwp3,longname='CLUBB 3rd moment vert velocity',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(8)),0._r8,0._r8,0._r8,ixup2,longname='CLUBB 2nd moment u wind',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(9)),0._r8,0._r8,0._r8,ixvp2,longname='CLUBB 2nd moment v wind',cam_outfld=.false.)
    end if

    !  put pbuf_add calls here (see macrop_driver.F90 for sample) use indicies defined at top
    call pbuf_add_field('pblh',       'global', dtype_r8, (/pcols/),                    pblh_idx)
    call pbuf_add_field('tke',        'global', dtype_r8, (/pcols, pverp/),             tke_idx)
    call pbuf_add_field('kvh',        'global', dtype_r8, (/pcols, pverp/),             kvh_idx)
    call pbuf_add_field('kvm',        'global', dtype_r8, (/pcols, pverp/),             kvm_idx)
    call pbuf_add_field('tpert',      'global', dtype_r8, (/pcols/),                    tpert_idx)
    call pbuf_add_field('AST',        'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    ast_idx)
    call pbuf_add_field('AIST',       'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    aist_idx)
    call pbuf_add_field('ALST',       'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    alst_idx)
    call pbuf_add_field('QIST',       'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    qist_idx)
    call pbuf_add_field('QLST',       'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    qlst_idx)
    call pbuf_add_field('CONCLD',     'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    concld_idx)
    call pbuf_add_field('CLD',        'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    cld_idx)
    call pbuf_add_field('FICE',       'physpkg',dtype_r8, (/pcols,pver/),               fice_idx)
    call pbuf_add_field('RAD_CLUBB',  'global', dtype_r8, (/pcols,pver/),               radf_idx)
    call pbuf_add_field('CMELIQ',     'physpkg',dtype_r8, (/pcols,pver/),                  cmeliq_idx)

    call pbuf_add_field('WP2_nadv',        'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), wp2_idx)
    call pbuf_add_field('WP3_nadv',        'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), wp3_idx)
    call pbuf_add_field('WPTHLP_nadv',     'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), wpthlp_idx)
    call pbuf_add_field('WPRTP_nadv',      'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), wprtp_idx)
    call pbuf_add_field('RTPTHLP_nadv',    'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), rtpthlp_idx)
    call pbuf_add_field('RTP2_nadv',       'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), rtp2_idx)
    call pbuf_add_field('THLP2_nadv',      'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), thlp2_idx)
    call pbuf_add_field('UP2_nadv',        'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), up2_idx)
    call pbuf_add_field('VP2_nadv',        'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), vp2_idx)

    call pbuf_add_field('UPWP',       'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), upwp_idx)
    call pbuf_add_field('VPWP',       'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), vpwp_idx)
    if (linearize_pbl_winds) then
       call pbuf_add_field('UM_PERT',    'physpkg', dtype_r8, (/pcols,pverp/), um_pert_idx)
       call pbuf_add_field('VM_PERT',    'physpkg', dtype_r8, (/pcols,pverp/), vm_pert_idx)
       call pbuf_add_field('UPWP_PERT',  'physpkg', dtype_r8, (/pcols,pverp/), upwp_pert_idx)
       call pbuf_add_field('VPWP_PERT',  'physpkg', dtype_r8, (/pcols,pverp/), vpwp_pert_idx)
    end if
    call pbuf_add_field('THLM',       'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), thlm_idx)
    call pbuf_add_field('RTM',        'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), rtm_idx)
    call pbuf_add_field('UM',         'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), um_idx)
    call pbuf_add_field('VM',         'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), vm_idx)

    call pbuf_add_field('WPTHVP',        'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), wpthvp_idx)
    call pbuf_add_field('WP2THVP',       'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), wp2thvp_idx)
    call pbuf_add_field('RTPTHVP',       'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), rtpthvp_idx)
    call pbuf_add_field('THLPTHVP',      'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), thlpthvp_idx)
    call pbuf_add_field('RCM',           'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), rcm_idx)
    call pbuf_add_field('CLOUD_FRAC',    'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), cloud_frac_idx)

    call pbuf_add_field('pdf_zm_w_1',    'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), pdf_zm_w_1_idx)
    call pbuf_add_field('pdf_zm_w_2',    'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), pdf_zm_w_2_idx)
    call pbuf_add_field('pdf_zm_var_w_1', 'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), pdf_zm_varnce_w_1_idx)
    call pbuf_add_field('pdf_zm_var_w_2', 'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), pdf_zm_varnce_w_2_idx)
    call pbuf_add_field('pdf_zm_mixt_frac',  'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), pdf_zm_mixt_frac_idx)

    call pbuf_add_field('vmag_gust',       'global', dtype_r8, (/pcols/),      vmag_gust_idx) !PMA total gustiness

    if (linearize_pbl_winds) then
       call pbuf_add_field('wsresp',          'global', dtype_r8, (/pcols/),      wsresp_idx)
       call pbuf_add_field('tau_est',         'global', dtype_r8, (/pcols/),      tau_est_idx)
    end if
#endif

  end subroutine clubb_register_cam
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

function clubb_implements_cnst(name)

  !----------------------------------------------------------------------------- !
  !                                                                              !
  ! Return true if specified constituent is implemented by this package          !
  !                                                                              !
  !----------------------------------------------------------------------------- !

   character(len=*), intent(in) :: name      ! constituent name
   logical :: clubb_implements_cnst     ! return value

   !-----------------------------------------------------------------------

   clubb_implements_cnst = (do_cnst .and. any(name == cnst_names))

end function clubb_implements_cnst


  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

subroutine clubb_init_cnst(name, q, gcid)

#ifdef CLUBB_SGS
    use constants_clubb,        only: w_tol_sqd, rt_tol, thl_tol
#endif

   !----------------------------------------------------------------------- !
   !                                                                        !
   ! Initialize the state if clubb_do_adv                                   !
   !                                                                        !
   !----------------------------------------------------------------------- !

   character(len=*), intent(in)  :: name     ! constituent name
   real(r8),         intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
   integer,          intent(in)  :: gcid(:)  ! global column id
   !-----------------------------------------------------------------------

#ifdef CLUBB_SGS
   if (clubb_do_adv) then
      if (trim(name) == trim(cnst_names(1))) q = real(thl_tol**2, kind = r8)
      if (trim(name) == trim(cnst_names(2))) q = real(rt_tol**2, kind = r8)
      if (trim(name) == trim(cnst_names(3))) q = 0.0_r8
      if (trim(name) == trim(cnst_names(4))) q = 0.0_r8
      if (trim(name) == trim(cnst_names(5))) q = 0.0_r8
      if (trim(name) == trim(cnst_names(6))) q = real(w_tol_sqd, kind = r8)
      if (trim(name) == trim(cnst_names(7))) q = 0.0_r8
      if (trim(name) == trim(cnst_names(8))) q = real(w_tol_sqd, kind = r8)
      if (trim(name) == trim(cnst_names(9))) q = real(w_tol_sqd, kind = r8)
   end if
#endif

end subroutine clubb_init_cnst


  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine clubb_readnl(nlfile)

#ifdef CLUBB_SGS
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use cam_abortutils,  only: endrun
    use stats_variables, only: l_stats, l_output_rad_files
    use mpishorthand
    use model_flags,     only: l_diffuse_rtm_and_thlm, l_stability_correct_Kh_N2_zm, &
                               l_vert_avg_closure, l_trapezoidal_rule_zt, &
                               l_trapezoidal_rule_zm, l_call_pdf_closure_twice,&
                               ipdf_call_placement

    use parameters_tunable, only: clubb_param_readnl
#endif

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

#ifdef CLUBB_SGS
    logical :: clubb_history, clubb_rad_history, clubb_cloudtop_cooling, clubb_rainevap_turb, &
               clubb_stabcorrect, clubb_expldiff ! Stats enabled (T/F)
    logical :: clubb_use_sgv !PMA This flag controls tuning for tpert and gustiness
    logical :: clubb_vert_avg_closure !XZheng This flag sets four clubb config flags for pdf_closure and the trapezoidal rule to  compute the varibles that are output from high order closure
    integer :: clubb_ipdf_call_placement  !XZheng This flag sets options for the placement of the call to CLUBB's PDF.

    integer :: iunit, read_status

    namelist /clubb_his_nl/ clubb_history, clubb_rad_history
    namelist /clubbpbl_diff_nl/ clubb_cloudtop_cooling, clubb_rainevap_turb, clubb_expldiff, &
                                clubb_do_adv, clubb_do_deep, clubb_timestep, clubb_stabcorrect, &
                                clubb_rnevap_effic, clubb_liq_deep, clubb_liq_sh, clubb_ice_deep, &
                                clubb_ice_sh, clubb_tk1, clubb_tk2, relvar_fix, clubb_use_sgv, &
                                clubb_vert_avg_closure, clubb_ipdf_call_placement


    !----- Begin Code -----

    !  Determine if we want clubb_history to be output
    clubb_history      = .false.   ! Initialize to false
    l_stats            = .false.   ! Initialize to false
    l_output_rad_files = .false.   ! Initialize to false
    do_cldcool         = .false.   ! Initialize to false
    do_rainturb        = .false.   ! Initialize to false
    do_expldiff        = .false.   ! Initialize to false
    relvar_fix         = .false.   ! Initialize to false
    clubb_do_adv       = .false.   ! Initialize to false
    clubb_do_deep      = .false.   ! Initialize to false
    use_sgv            = .false.
    clubb_vert_avg_closure = .true.
    clubb_ipdf_call_placement = -999


    !  Read namelist to determine if CLUBB history should be called
    if (masterproc) then
      iunit = getunit()
      open( iunit, file=trim(nlfile), status='old' )

      call find_group_name(iunit, 'clubb_his_nl', status=read_status)
      if (read_status == 0) then
         read(unit=iunit, nml=clubb_his_nl, iostat=read_status)
         if (read_status /= 0) then
            call endrun('clubb_readnl:  error reading namelist'//errmsg(__FILE__,__LINE__))
         end if
      end if

      call find_group_name(iunit, 'clubbpbl_diff_nl', status=read_status)
      if (read_status == 0) then
         read(unit=iunit, nml=clubbpbl_diff_nl, iostat=read_status)
         if (read_status /= 0) then
            call endrun('clubb_readnl:  error reading namelist'//errmsg(__FILE__,__LINE__))
         end if
      end if

      close(unit=iunit)
      call freeunit(iunit)
    end if

#ifdef SPMD
! Broadcast namelist variables
      call mpibcast(clubb_history,            1,   mpilog,   0, mpicom)
      call mpibcast(clubb_rad_history,        1,   mpilog,   0, mpicom)
      call mpibcast(clubb_cloudtop_cooling,   1,   mpilog,   0, mpicom)
      call mpibcast(clubb_rainevap_turb,      1,   mpilog,   0, mpicom)
      call mpibcast(clubb_expldiff,           1,   mpilog,   0, mpicom)
      call mpibcast(clubb_do_adv,             1,   mpilog,   0, mpicom)
      call mpibcast(clubb_do_deep,            1,   mpilog,   0, mpicom)
      call mpibcast(clubb_timestep,           1,   mpir8,   0, mpicom)
      call mpibcast(clubb_stabcorrect,        1,   mpilog,   0, mpicom)
      call mpibcast(clubb_rnevap_effic,       1,   mpir8,   0, mpicom)
      call mpibcast(clubb_liq_deep,           1,   mpir8,   0, mpicom)
      call mpibcast(clubb_liq_sh,             1,   mpir8,   0, mpicom)
      call mpibcast(clubb_ice_deep,           1,   mpir8,   0, mpicom)
      call mpibcast(clubb_ice_sh,             1,   mpir8,   0, mpicom)
      call mpibcast(clubb_tk1,                1,   mpir8,   0, mpicom)
      call mpibcast(clubb_tk2,                1,   mpir8,   0, mpicom)
      call mpibcast(relvar_fix,               1,   mpilog,  0, mpicom)
      call mpibcast(clubb_use_sgv,            1,   mpilog,   0, mpicom)
      call mpibcast(clubb_vert_avg_closure,   1,   mpilog,   0, mpicom)
      call mpibcast(clubb_ipdf_call_placement,   1,   mpiint,   0, mpicom)
#endif

    !  Overwrite defaults if they are true
    if (clubb_ipdf_call_placement > 0) ipdf_call_placement = clubb_ipdf_call_placement
    if (clubb_history) l_stats = .true.
    if (clubb_rad_history) l_output_rad_files = .true.
    if (clubb_cloudtop_cooling) do_cldcool = .true.
    if (clubb_rainevap_turb) do_rainturb = .true.
    if (clubb_expldiff) do_expldiff = .true.
    if (clubb_use_sgv) use_sgv =.true.
    if (clubb_stabcorrect .and. clubb_expldiff)  then
      call endrun('clubb_readnl: clubb_stabcorrect and clubb_expldiff may not both be set to true at the same time'//errmsg(__FILE__,__LINE__))
    end if

    if (clubb_stabcorrect) then
      l_diffuse_rtm_and_thlm       = .true.   ! CLUBB flag set to true
      l_stability_correct_Kh_N2_zm = .true.   ! CLUBB flag set to true
    endif

    if (clubb_vert_avg_closure) then
      l_vert_avg_closure       = .true.   ! CLUBB flag set to true
      l_trapezoidal_rule_zt    = .true.   ! CLUBB flag set to true
      l_trapezoidal_rule_zm    = .true.   ! CLUBB flag set to true
      l_call_pdf_closure_twice = .true.   ! CLUBB flag set to true
    else
      l_vert_avg_closure       = .false.   ! CLUBB flag set to false
      l_trapezoidal_rule_zt    = .false.   ! CLUBB flag set to false
      l_trapezoidal_rule_zm    = .false.   ! CLUBB flag set to false
      l_call_pdf_closure_twice = .false.   ! CLUBB flag set to false
    endif

    ! read tunable parameters from namelist, handlings of masterproc vs others
    ! are done within clubb_param_readnl
    call clubb_param_readnl(nlfile)
#endif
  end subroutine clubb_readnl

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine clubb_ini_cam(pbuf2d, dp1_in)
!-------------------------------------------------------------------------------
! Description:
!   Initialize UWM CLUBB.
! Author: Cheryl Craig March 2011
! Modifications: Pete Bogenschutz 2011 March and onward
! Origin: Based heavily on UWM clubb_init.F90
! References:
!   None
!-------------------------------------------------------------------------------



#ifdef CLUBB_SGS

    !  From CAM libraries
    use physics_types,          only: physics_state, physics_ptend
    use cam_history,            only: addfld, horiz_only, add_default
    use ppgrid,                 only: pver, pverp, pcols, begchunk, endchunk
    use ref_pres,               only: pref_mid
    use hb_diff,                only: init_hb_diff
    use trb_mtn_stress,         only: init_tms
    use rad_constituents,       only: rad_cnst_get_info, rad_cnst_get_mode_num_idx, rad_cnst_get_mam_mmr_idx

    !  From the CLUBB libraries

    !  From the CLUBB libraries
    use clubb_api_module, only: &
         setup_clubb_core_api, &
         time_precision, &
         core_rknd, &
         set_clubb_debug_level_api, &
         nparams, &
         read_parameters_api, &
         l_stats, &
         l_stats_samp, &
         l_grads, &
         stats_zt, &
         stats_zm, &
         stats_sfc, &
         stats_rad_zt, &
         stats_rad_zm, &
         w_tol_sqd, &
         rt_tol, &
         l_do_expldiff_rtm_thlm, &
         init_pdf_params_api
    use stats_variables,           only: l_output_rad_files

    use units,                     only: getunit, freeunit
    use error_messages,            only: handle_errmsg
    use time_manager,              only: is_first_step
    use constants_clubb,           only: thl_tol


    !  These are only needed if we're using a passive scalar
    use array_index,            only: iisclr_rt, iisclr_thl, iisclr_CO2, &    ! [kg/kg]/[K]/[1e6 mol/mol]
                                      iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2 ! "    "
    use constituents,           only: cnst_get_ind
    use phys_control,           only: phys_getopts

    use parameters_tunable,     only: params_list
    use cam_abortutils,         only: endrun

#endif

    use physics_buffer,         only: pbuf_get_index, pbuf_set_field, physics_buffer_desc
    implicit none
    !  Input Variables
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    real(r8) :: dp1_in

#ifdef CLUBB_SGS

    real(kind=time_precision) :: dum1, dum2, dum3

    real(core_rknd), dimension(nparams)  :: clubb_params    ! These adjustable CLUBB parameters (C1, C2 ...)

    logical :: clubb_history, clubb_rad_history, clubb_cloudtop_cooling, clubb_rainevap_turb, clubb_expldiff ! Stats enabled (T/F)

    ! The similar name to clubb_history is unfortunate...
    logical :: history_amwg, history_clubb

    character(len=128) :: errstring             ! error status for CLUBB init

    integer :: err_code, iunit                  ! Code for when CLUBB fails
    integer :: i, j, k, l, idx_chunk, idx_pcols ! Indices
    integer :: read_status                      ! Length of a string
    integer :: ntop_eddy                        ! Top    interface level to which eddy vertical diffusion is applied ( = 1 )
    integer :: nbot_eddy                        ! Bottom interface level to which eddy vertical diffusion is applied ( = pver )
    integer :: nmodes, nspec, pmam_ncnst, m
    integer :: ixnumliq
    integer :: lptr

    real(core_rknd)  :: zt_g(pverp)                        ! Height dummy array
    real(core_rknd)  :: zi_g(pverp)                        ! Height dummy array


    !----- Begin Code -----
    !$OMP PARALLEL
    l_do_expldiff_rtm_thlm = do_expldiff
    !$OMP END PARALLEL

    allocate( &
       pdf_params_chnk(pcols,begchunk:endchunk),   &
       pdf_params_zm_chnk(pcols,begchunk:endchunk) )

    do idx_chunk = begchunk, endchunk
        do idx_pcols = 1, pcols
            call init_pdf_params_api( pverp, pdf_params_chnk(idx_pcols,idx_chunk) )
            call init_pdf_params_api( pverp, pdf_params_zm_chnk(idx_pcols,idx_chunk) )
        end do
    end do

    ! ----------------------------------------------------------------- !
    ! Determine how many constituents CLUBB will transport.  Note that
    ! CLUBB does not transport aerosol consituents.  Therefore, need to
    ! determine how many aerosols constituents there are and subtract that
    ! off of pcnst (the total consituents)
    ! ----------------------------------------------------------------- !

    call phys_getopts(prog_modal_aero_out=prog_modal_aero, &
                      history_amwg_out=history_amwg, &
                      history_clubb_out=history_clubb,&
                      liqcf_fix_out   = liqcf_fix)

    !  Select variables to apply tendencies back to CAM

    ! Initialize all consituents to true to start
    lq(1:pcnst) = .true.
    edsclr_dim  = pcnst

    if (prog_modal_aero) then
       ! Turn off modal aerosols and decrement edsclr_dim accordingly
       call rad_cnst_get_info(0, nmodes=nmodes)

       do m = 1, nmodes
          call rad_cnst_get_mode_num_idx(m, lptr)
          lq(lptr)=.false.
          edsclr_dim = edsclr_dim-1

          call rad_cnst_get_info(0, m, nspec=nspec)
          do l = 1, nspec
             call rad_cnst_get_mam_mmr_idx(m, l, lptr)
             lq(lptr)=.false.
             edsclr_dim = edsclr_dim-1
          end do
       end do

       !  In addition, if running with MAM, droplet number is transported
       !  in dropmixnuc, therefore we do NOT want CLUBB to apply transport
       !  tendencies to avoid double counted.  Else, we apply tendencies.
       call cnst_get_ind('NUMLIQ',ixnumliq)
       lq(ixnumliq) = .false.
       edsclr_dim = edsclr_dim-1
    endif

    ! ----------------------------------------------------------------- !
    ! Set the debug level.  Level 2 has additional computational expense since
    ! it checks the array variables in CLUBB for invalid values.
    ! ----------------------------------------------------------------- !
    call set_clubb_debug_level_api( 0 )

    ! ----------------------------------------------------------------- !
    ! use pbuf_get_fld_idx to get existing physics buffer fields from other
    ! physics packages (e.g. tke)
    ! ----------------------------------------------------------------- !


    !  Defaults
    l_stats_samp = .false.
    l_grads = .false.

    !  Overwrite defaults if needbe
    if (l_stats) l_stats_samp = .true.

    !  Define physics buffers indexes
    cld_idx     = pbuf_get_index('CLD')         ! Cloud fraction
    concld_idx  = pbuf_get_index('CONCLD')      ! Convective cloud cover
    ast_idx     = pbuf_get_index('AST')         ! Stratiform cloud fraction
    alst_idx    = pbuf_get_index('ALST')        ! Liquid stratiform cloud fraction
    aist_idx    = pbuf_get_index('AIST')        ! Ice stratiform cloud fraction
    qlst_idx    = pbuf_get_index('QLST')        ! Physical in-stratus LWC
    qist_idx    = pbuf_get_index('QIST')        ! Physical in-stratus IWC
    dp_frac_idx = pbuf_get_index('DP_FRAC')     ! Deep convection cloud fraction
    icwmrdp_idx = pbuf_get_index('ICWMRDP')     ! In-cloud deep convective mixing ratio
    sh_frac_idx = pbuf_get_index('SH_FRAC')     ! Shallow convection cloud fraction
    relvar_idx  = pbuf_get_index('RELVAR')      ! Relative cloud water variance
    accre_enhan_idx = pbuf_get_index('ACCRE_ENHAN') ! accretion enhancement for MG
    prer_evap_idx   = pbuf_get_index('PRER_EVAP')
    qrl_idx         = pbuf_get_index('QRL')
    cmfmc_sh_idx    = pbuf_get_index('CMFMC_SH')

    prec_dp_idx = pbuf_get_index('PREC_DP') !PMA ZM precip for gustiness
    snow_dp_idx = pbuf_get_index('SNOW_DP') !PMA ZM snow for gustiness
    vmag_gust_idx = pbuf_get_index('vmag_gust') !PMA ZM snow for gustiness

    iisclr_rt  = -1
    iisclr_thl = -1
    iisclr_CO2 = -1

    iiedsclr_rt  = -1
    iiedsclr_thl = -1
    iiedsclr_CO2 = -1

    if (zm_microp) then
        dlfzm_idx = pbuf_get_index('DLFZM')
        difzm_idx = pbuf_get_index('DIFZM')
        dsfzm_idx = pbuf_get_index('DSFZM')
        dnlfzm_idx = pbuf_get_index('DNLFZM')
        dnifzm_idx = pbuf_get_index('DNIFZM')
        dnsfzm_idx = pbuf_get_index('DNSFZM')
    end if

    ! ----------------------------------------------------------------- !
    ! Define number of tracers for CLUBB to diffuse
    ! ----------------------------------------------------------------- !

    if (do_expldiff) then
       offset = 2 ! diffuse temperature and moisture explicitly
       edsclr_dim = edsclr_dim + offset
    endif

    ! ----------------------------------------------------------------- !
    ! Setup CLUBB core
    ! ----------------------------------------------------------------- !

    !  Read in parameters for CLUBB.  Pack the default and updated (via nml)
    !  tunable parameters into clubb_params
!$OMP PARALLEL
    call read_parameters_api( -99, "", clubb_params )
!$OMP END PARALLEL

    ! Print the list of CLUBB parameters, if multi-threaded, it may print by each thread
    if (masterproc) then
       write(iulog,*)'CLUBB tunable parameters: total ',nparams
       write(iulog,*)'--------------------------------------------------'
       do i = 1, nparams
          write(iulog,*) params_list(i), " = ", clubb_params(i)
       enddo
    endif


    !  Fill in dummy arrays for height.  Note that these are overwrote
    !  at every CLUBB step to physical values.
    do k=1,pverp
       zt_g(k) = ((k-1)*1000._core_rknd)-500._core_rknd  !  this is dummy garbage
       zi_g(k) = (k-1)*1000._core_rknd                   !  this is dummy garbage
    enddo

    !  Set up CLUBB core.  Note that some of these inputs are overwrote
    !  when clubb_tend_cam is called.  The reason is that heights can change
    !  at each time step, which is why dummy arrays are read in here for heights
    !  as they are immediately overwrote.
!$OMP PARALLEL
    call setup_clubb_core_api     &
         ( pverp, theta0, ts_nudge, &                                 ! In
           hydromet_dim,  sclr_dim, &                                 ! In
           sclr_tol, edsclr_dim, clubb_params, &                      ! In
           l_host_applies_sfc_fluxes, &                               ! In
           l_uv_nudge, saturation_equation, l_input_fields,  &        ! In
           l_implemented, grid_type, zi_g(2), zi_g(1), zi_g(pverp), & ! In
           zi_g(1:pverp), zt_g(1:pverp), zi_g(1), &
           err_code )
!$OMP END PARALLEL

    ! ----------------------------------------------------------------- !
    ! Set-up HB diffusion.  Only initialized to diagnose PBL depth      !
    ! ----------------------------------------------------------------- !

    ! Initialize eddy diffusivity module

    ntop_eddy = 1    ! if >1, must be <= nbot_molec
    nbot_eddy = pver ! currently always pver

    call init_hb_diff( gravit, cpair, ntop_eddy, nbot_eddy, pref_mid, karman, eddy_scheme )

    ! ----------------------------------------------------------------- !
    ! Initialize turbulent mountain stress module                       !
    ! ------------------------------------------------------------------!

    if ( do_tms) then
       call init_tms( r8, tms_orocnst, tms_z0fac, karman, gravit, rair, errstring)
       call handle_errmsg(errstring, subname="init_tms")

       call addfld( 'TAUTMSX' ,  horiz_only,  'A','N/m2',  'Zonal      turbulent mountain surface stress' )
       call addfld( 'TAUTMSY' ,  horiz_only,  'A','N/m2',  'Meridional turbulent mountain surface stress' )
       if (history_amwg) then
          call add_default( 'TAUTMSX ', 1, ' ' )
          call add_default( 'TAUTMSY ', 1, ' ' )
       end if
       if (masterproc) then
          write(iulog,*)'Using turbulent mountain stress module'
          write(iulog,*)'  tms_orocnst = ',tms_orocnst
          write(iulog,*)'  tms_z0fac = ',tms_z0fac
       end if
    endif

    ! ----------------------------------------------------------------- !
    ! Add output fields for the history files
    ! ----------------------------------------------------------------- !

    if (clubb_do_deep) then
       call addfld ('MU_CLUBB',horiz_only,'A','1/m','CLUBB value of entrainment')
    endif

    !  These are default CLUBB output.  Not the higher order history budgets
    call addfld ('RHO_CLUBB',    (/ 'ilev' /), 'A',        'kg/m3', 'Air Density')
    call addfld ('UP2_CLUBB',    (/ 'ilev' /), 'A',        'm2/s2', 'Zonal Velocity Variance')
    call addfld ('VP2_CLUBB',    (/ 'ilev' /), 'A',        'm2/s2', 'Meridional Velocity Variance')
    call addfld ('WP2_CLUBB',    (/ 'ilev' /), 'A',        'm2/s2', 'Vertical Velocity Variance')
    call addfld ('UPWP_CLUBB',    (/ 'ilev' /), 'A',       'm2/s2', 'Zonal Momentum Flux')
    call addfld ('VPWP_CLUBB',    (/ 'ilev' /), 'A',       'm2/s2', 'Meridional Momentum Flux')
    call addfld ('WP3_CLUBB',    (/ 'ilev' /), 'A',        'm3/s3', 'Third Moment Vertical Velocity')
    call addfld ('WPTHLP_CLUBB',     (/ 'ilev' /), 'A',     'W/m2', 'Heat Flux')
    call addfld ('WPRTP_CLUBB',     (/ 'ilev' /), 'A',      'W/m2', 'Moisture Flux')
    call addfld ('RTP2_CLUBB', (/ 'ilev' /), 'A',       'g^2/kg^2', 'Moisture Variance')
    call addfld ('THLP2_CLUBB',      (/ 'ilev' /), 'A',      'K^2', 'Temperature Variance')
    call addfld ('RTPTHLP_CLUBB',   (/ 'ilev' /), 'A',    'K g/kg', 'Temp. Moist. Covariance')
    call addfld ('RCM_CLUBB',     (/ 'ilev' /), 'A',        'g/kg', 'Cloud Water Mixing Ratio')
    call addfld ('WPRCP_CLUBB',     (/ 'ilev' /), 'A',      'W/m2', 'Liquid Water Flux')
    call addfld ('CLOUDFRAC_CLUBB', (/ 'lev' /),  'A',  '1', 'Cloud Fraction')
    call addfld ('RCMINLAYER_CLUBB',     (/ 'ilev' /), 'A', 'g/kg', 'Cloud Water in Layer')
    call addfld ('CLOUDCOVER_CLUBB', (/ 'ilev' /), 'A', '1', 'Cloud Cover')
    call addfld ('WPTHVP_CLUBB',     (/ 'lev' /),  'A',     'W/m2', 'Buoyancy Flux')
    call addfld ('RVMTEND_CLUBB',  (/ 'lev' /),  'A',    'g/kg /s', 'Water vapor tendency')
    call addfld ('TTEND_CLUBB',      (/ 'lev' /),  'A',      'k/s', 'Temperature tendency')
    call addfld ('RCMTEND_CLUBB',  (/ 'lev' /),  'A',    'g/kg /s', 'Cloud Liquid Water Tendency')
    call addfld ('RIMTEND_CLUBB',  (/ 'lev' /),  'A',    'g/kg /s', 'Cloud Ice Tendency')
    call addfld ('UTEND_CLUBB',   (/ 'lev' /),  'A',      'm/s /s', 'U-wind Tendency')
    call addfld ('VTEND_CLUBB',   (/ 'lev' /),  'A',      'm/s /s', 'V-wind Tendency')
    call addfld ('ZT_CLUBB',        (/ 'ilev' /), 'A',         'm', 'Thermodynamic Heights')
    call addfld ('ZM_CLUBB',        (/ 'ilev' /), 'A',         'm', 'Momentum Heights')
    call addfld ('UM_CLUBB',      (/ 'ilev' /), 'A',         'm/s', 'Zonal Wind')
    call addfld ('VM_CLUBB',      (/ 'ilev' /), 'A',         'm/s', 'Meridional Wind')
    call addfld ('THETAL',        (/ 'lev' /),  'A',           'K', 'Liquid Water Potential Temperature')
    call addfld ('PBLH',        horiz_only,     'A',             'm', 'PBL height')
    call addfld ('QT',    (/ 'lev' /),  'A',               'kg/kg', 'Total water mixing ratio')
    call addfld ('SL',     (/ 'lev' /),  'A',               'J/kg', 'Liquid water static energy')
    call addfld ('CLDST', (/ 'lev' /),  'A',            'fraction', 'Stratus cloud fraction')
    call addfld ('ZMDLF',  (/ 'lev' /),  'A',            'kg/kg/s', 'Detrained liquid water from ZM convection')
    call addfld ('TTENDICE',     (/ 'lev' /),  'A',         'K/s', 'T tendency from Ice Saturation Adjustment')
    call addfld ('QVTENDICE', (/ 'lev' /),  'A',        'kg/kg/s', 'Q tendency from Ice Saturation Adjustment')
    call addfld ('QITENDICE', (/ 'lev' /),  'A',        'kg/kg/s', 'CLDICE tendency from Ice Saturation Adjustment')
    call addfld ('NITENDICE', (/ 'lev' /),  'A',        'kg/kg/s', 'NUMICE tendency from Ice Saturation Adjustment')
    call addfld ('DPDLFLIQ', (/ 'lev' /),  'A',        'kg/kg/s', 'Detrained liquid water from deep convection')
    call addfld ('DPDLFICE', (/ 'lev' /),  'A',        'kg/kg/s', 'Detrained ice from deep convection')
    call addfld ('DPDLFT', (/ 'lev' /),  'A',        'K/s', 'T-tendency due to deep convective detrainment')
    call addfld ('RELVAR', (/ 'lev' /),  'A',        '-', 'Relative cloud water variance')
    call addfld ('RELVARC', (/ 'lev' /),  'A',        '-', 'Relative cloud water variance', flag_xyfill=.true.,fill_value=fillvalue)
    call addfld ('CONCLD', (/ 'lev' /),  'A',        'fraction', 'Convective cloud cover')
    call addfld ('CMELIQ', (/ 'lev' /),  'A',        'kg/kg/s', 'Rate of cond-evap of liq within the cloud')
!PMA gustiness output fields
    call addfld ('VMAGGUST',       horiz_only,     'A',             '-', 'Total gustiness enhancement')
    call addfld ('VMAGDP',        horiz_only,     'A',             '-', 'ZM gustiness enhancement')
    call addfld ('VMAGCL',        horiz_only,     'A',             '-', 'CLUBB gustiness enhancement')
    call addfld ('TPERTBLT',        horiz_only,     'A',             'K', 'perturbation temperature at PBL top')
    !
    if (use_od_fd) then
       !added for turbulent orographic form drag (TOFD) output
       call addfld ('DTAUX3_FD',(/'lev'/),'A','m/s2','U tendency - fd orographic drag')
       call addfld ('DTAUY3_FD',(/'lev'/),'A','m/s2','V tendency - fd orographic drag')
       call addfld ('DUSFC_FD',horiz_only,'A','N/m2','fd zonal oro surface stress')
       call addfld ('DVSFC_FD',horiz_only,'A','N/m2','fd merio oro surface stress')
       call add_default('DTAUX3_FD', 1,  ' ')
       call add_default('DTAUY3_FD', 1,  ' ')
       call add_default('DUSFC_FD',  1,  ' ')
       call add_default('DVSFC_FD',  1,  ' ')
       if (masterproc) then
          write(iulog,*)'Using turbulent orographic form drag scheme (TOFD)'
       end if
       if (use_od_fd.and.do_tms) then
          call endrun("clubb_intr: Both TMS and TOFD are turned on, please turn one off&
          &by setting use_od_fd or do_tms as .false.")
       end if
    end if
    !  Initialize statistics, below are dummy variables
    dum1 = 300._r8
    dum2 = 1200._r8
    dum3 = 300._r8

    if (l_stats) then

       call stats_init_clubb( .true., dum1, dum2, &
                         pverp, pverp, pverp, dum3 )

       allocate(out_zt(pcols,pverp,stats_zt%num_output_fields))
       allocate(out_zm(pcols,pverp,stats_zm%num_output_fields))
       allocate(out_sfc(pcols,1,stats_sfc%num_output_fields))

       allocate(out_radzt(pcols,pverp,stats_rad_zt%num_output_fields))
       allocate(out_radzm(pcols,pverp,stats_rad_zm%num_output_fields))

    endif

    ! ----------------------------------------------------------------- !
    ! Make all of this output default, this is not CLUBB history
    ! ----------------------------------------------------------------- !
    if (clubb_do_adv .or. history_clubb) then
       call add_default('WP2_CLUBB',        1, ' ')
       call add_default('WP3_CLUBB',        1, ' ')
       call add_default('WPTHLP_CLUBB',     1, ' ')
       call add_default('WPRTP_CLUBB',      1, ' ')
       call add_default('RTP2_CLUBB',       1, ' ')
       call add_default('THLP2_CLUBB',      1, ' ')
       call add_default('RTPTHLP_CLUBB',    1, ' ')
       call add_default('UP2_CLUBB',        1, ' ')
       call add_default('VP2_CLUBB',        1, ' ')
    end if

    if (history_clubb) then

       if (clubb_do_deep) then
          call add_default('MU_CLUBB',         1, ' ')
       endif

       call add_default('RELVAR',           1, ' ')
       call add_default('RHO_CLUBB',        1, ' ')
       call add_default('UPWP_CLUBB',       1, ' ')
       call add_default('VPWP_CLUBB',       1, ' ')
       call add_default('RCM_CLUBB',        1, ' ')
       call add_default('WPRCP_CLUBB',      1, ' ')
       call add_default('CLOUDFRAC_CLUBB',  1, ' ')
       call add_default('RCMINLAYER_CLUBB', 1, ' ')
       call add_default('CLOUDCOVER_CLUBB', 1, ' ')
       call add_default('WPTHVP_CLUBB',     1, ' ')
       call add_default('RVMTEND_CLUBB',    1, ' ')
       call add_default('TTEND_CLUBB',      1, ' ')
       call add_default('RCMTEND_CLUBB',    1, ' ')
       call add_default('RIMTEND_CLUBB',    1, ' ')
       call add_default('UTEND_CLUBB',      1, ' ')
       call add_default('VTEND_CLUBB',      1, ' ')
       call add_default('ZT_CLUBB',         1, ' ')
       call add_default('ZM_CLUBB',         1, ' ')
       call add_default('UM_CLUBB',         1, ' ')
       call add_default('VM_CLUBB',         1, ' ')
       call add_default('SL',               1, ' ')
       call add_default('QT',               1, ' ')
       call add_default('CONCLD',           1, ' ')
    else
       call add_default('CLOUDFRAC_CLUBB',  1, ' ')
       call add_default('CONCLD',           1, ' ')
    end if

    if (history_amwg) then
       call add_default('PBLH',             1, ' ')
    end if

    if (history_budget) then
       call add_default('DPDLFLIQ',         history_budget_histfile_num, ' ')
       call add_default('DPDLFICE',         history_budget_histfile_num, ' ')
       call add_default('DPDLFT',           history_budget_histfile_num, ' ')
       call add_default('TTEND_CLUBB',      history_budget_histfile_num, ' ')
       call add_default('RCMTEND_CLUBB',    history_budget_histfile_num, ' ')
       call add_default('RIMTEND_CLUBB',    history_budget_histfile_num, ' ')
       call add_default('RVMTEND_CLUBB',    history_budget_histfile_num, ' ')
       call add_default('UTEND_CLUBB',      history_budget_histfile_num, ' ')
       call add_default('VTEND_CLUBB',      history_budget_histfile_num, ' ')
    endif


    ! --------------- !
    ! First step?     !
    ! Initialization  !
    ! --------------- !

    !  Is this the first time step?  If so then initialize CLUBB variables as follows
    if (is_first_step()) then

       call pbuf_set_field(pbuf2d, wp2_idx,     real(w_tol_sqd, kind = r8))
       call pbuf_set_field(pbuf2d, wp3_idx,     0.0_r8)
       call pbuf_set_field(pbuf2d, wpthlp_idx,  0.0_r8)
       call pbuf_set_field(pbuf2d, wprtp_idx,   0.0_r8)
       call pbuf_set_field(pbuf2d, rtpthlp_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, rtp2_idx,    real(rt_tol**2, kind = r8))
       call pbuf_set_field(pbuf2d, thlp2_idx,   real(thl_tol**2, kind = r8))
       call pbuf_set_field(pbuf2d, up2_idx,     real(w_tol_sqd, kind = r8))
       call pbuf_set_field(pbuf2d, vp2_idx,     real(w_tol_sqd, kind = r8))

       call pbuf_set_field(pbuf2d, upwp_idx,    0.0_r8)
       call pbuf_set_field(pbuf2d, vpwp_idx,    0.0_r8)
       if (linearize_pbl_winds) then
          call pbuf_set_field(pbuf2d, um_pert_idx, 0.0_r8)
          call pbuf_set_field(pbuf2d, vm_pert_idx, 0.0_r8)
          call pbuf_set_field(pbuf2d, upwp_pert_idx, 0.0_r8)
          call pbuf_set_field(pbuf2d, vpwp_pert_idx, 0.0_r8)
       end if
       call pbuf_set_field(pbuf2d, tke_idx,     0.0_r8)
       call pbuf_set_field(pbuf2d, kvh_idx,     0.0_r8)
       call pbuf_set_field(pbuf2d, fice_idx,    0.0_r8)
       call pbuf_set_field(pbuf2d, radf_idx,    0.0_r8)
       call pbuf_set_field(pbuf2d, qrl_idx,     0.0_r8)

       call pbuf_set_field(pbuf2d, wpthvp_idx,     0.0_r8)
       call pbuf_set_field(pbuf2d, wp2thvp_idx,    0.0_r8)
       call pbuf_set_field(pbuf2d, rtpthvp_idx,    0.0_r8)
       call pbuf_set_field(pbuf2d, thlpthvp_idx,   0.0_r8)
       call pbuf_set_field(pbuf2d, rcm_idx,        0.0_r8)
       call pbuf_set_field(pbuf2d, cloud_frac_idx, 0.0_r8)

       call pbuf_set_field(pbuf2d, pdf_zm_w_1_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, pdf_zm_w_2_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, pdf_zm_varnce_w_1_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, pdf_zm_varnce_w_2_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, pdf_zm_mixt_frac_idx, 0.0_r8)

       call pbuf_set_field(pbuf2d, vmag_gust_idx,    1.0_r8)

       if (linearize_pbl_winds) then
          call pbuf_set_field(pbuf2d, wsresp_idx,    0.0_r8)
          call pbuf_set_field(pbuf2d, tau_est_idx,   0.0_r8)
       end if
    endif

    ! --------------- !
    ! End             !
    ! Initialization  !
    ! --------------- !

#endif
    dp1 = dp1_in !set via namelist, assigned in cloud_fraction.F90
    end subroutine clubb_ini_cam


  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

   subroutine clubb_tend_cam( &
                              state,   ptend_all,   pbuf,  diag,   hdtime, &
                              cmfmc,   cam_in,   cam_out,  sgh30,          &
                              macmic_it, cld_macmic_num_steps,dlf, det_s, det_ice, alst_o)

!-------------------------------------------------------------------------------
! Description: Provide tendencies of shallow convection, turbulence, and
!              macrophysics from CLUBB to CAM
!
! Author: Cheryl Craig, March 2011
! Modifications: Pete Bogenschutz, March 2011 and onward
! Origin: Based heavily on UWM clubb_init.F90
! References:
!   None
!-------------------------------------------------------------------------------

   use physics_types,  only: physics_state, physics_ptend, &
                             physics_state_copy, physics_ptend_init, &
                             physics_ptend_sum, set_dry_to_wet

   use physics_update_mod, only: physics_update

   use physics_buffer, only: pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field, &
                             pbuf_set_field, physics_buffer_desc

   use ppgrid,         only: pver, pverp, pcols
   use constituents,   only: cnst_get_ind, cnst_type
   use co2_cycle,      only: co2_cycle_set_cnst_type
   use camsrfexch,     only: cam_in_t, cam_out_t
   use ref_pres,       only: top_lev => trop_cloud_top_lev
   use time_manager,   only: is_first_step, is_first_restart_step, get_nstep
   use cam_abortutils, only: endrun
   use wv_saturation,  only: qsat
   use micro_mg_cam,   only: micro_mg_version

   use conditional_diag,      only: cnd_diag_t
   use conditional_diag_main, only: cnd_diag_checkpoint

#ifdef CLUBB_SGS
   use hb_diff,                   only: pblintd
   use iop_data_mod,              only: single_column
   use phys_grid,                 only: get_gcol_p
   use cldfrc2m,                  only: aist_vector
   use cam_history,               only: outfld
   use trb_mtn_stress,            only: compute_tms
   use macrop_driver,             only: ice_macro_tend

   use parameters_tunable,        only: mu
   use clubb_api_module, only: &
        cleanup_clubb_core_api, &
        nparams, &
        read_parameters_api, &
        setup_parameters_api, &
        setup_grid_heights_api, &
        w_tol_sqd, &
        rt_tol, &
        thl_tol, &
        l_stats, &
        stats_tsamp, &
        stats_tout, &
        stats_zt, &
        stats_sfc, &
        stats_zm, &
        stats_rad_zt, &
        stats_rad_zm, &
        l_output_rad_files, &
        pdf_parameter, &
        stats_begin_timestep_api, &
        advance_clubb_core_api, &
        calculate_thlp2_rad_api, &
        update_xp2_mc_api, &
        zt2zm_api, zm2zt_api
   use model_flags, only: ipdf_call_placement
   use advance_clubb_core_module, only: ipdf_post_advance_fields
#endif
   use od_common,          only: grid_size, oro_drag_interface
   use hycoef,             only: etamid
   use physconst,          only: rh2o,pi,rearth,r_universal
   implicit none

   ! --------------- !
   ! Input Auguments !
   ! --------------- !

   type(physics_state), intent(in)    :: state                    ! Physics state variables                 [vary]
   type(cam_in_t),      intent(in)    :: cam_in
   real(r8),            intent(in)    :: hdtime                   ! Host model timestep                     [s]
   real(r8),            intent(in)    :: dlf(pcols,pver)          ! Detraining cld H20 from deep convection [kg/ks/s]
   real(r8),            intent(in)    :: cmfmc(pcols,pverp)       ! convective mass flux--m sub c           [kg/m2/s]
   real(r8),            intent(in)    :: sgh30(pcols)             ! std deviation of orography              [m]
   integer,             intent(in)    :: cld_macmic_num_steps     ! number of mac-mic iterations
   integer,             intent(in)    :: macmic_it                ! number of mac-mic iterations
   type(cam_out_t),     intent(in)    :: cam_out

   ! ---------------------- !
   ! Input-Output Auguments !
   ! ---------------------- !

   type(physics_buffer_desc), pointer :: pbuf(:)
   type(cnd_diag_t),    intent(inout) :: diag                     !  conditionally sampled fields

   ! ---------------------- !
   ! Output Auguments !
   ! ---------------------- !

   type(physics_ptend), intent(out)   :: ptend_all                 ! package tendencies

   ! These two variables are needed for energy check
   real(r8),            intent(out)   :: det_s(pcols)              ! Integral of detrained static energy from ice
   real(r8),            intent(out)   :: det_ice(pcols)            ! Integral of detrained ice for energy check

   real(r8), intent(out) :: alst_o(pcols,pver)  ! H. Wang: for old liquid status fraction

   ! --------------- !
   ! Local Variables !
   ! --------------- !

#ifdef CLUBB_SGS

   type(physics_state) :: state1                ! Local copy of state variable
   type(physics_ptend) :: ptend_loc             ! Local tendency from processes, added up to return as ptend_all

   integer :: i, j, k, t, ixind, nadv
   integer :: ixcldice, ixcldliq, ixnumliq, ixnumice, ixq
   integer :: itim_old
   integer :: ncol, lchnk                       ! # of columns, and chunk identifier
   integer :: err_code                          ! Diagnostic, for if some calculation goes amiss.
   integer :: icnt, clubbtop

   real(r8) :: frac_limit, ic_limit

  !=====================================================================================
  ! The variables defined as core_rknd is required by the advance_clubb_core_api()
  ! subroutine, the changes here is to change the precision in the CLUBB calculation
  !=====================================================================================
   real(core_rknd) :: dtime                            ! CLUBB time step                              [s]
   real(core_rknd) :: edsclr_in(pverp,edsclr_dim)      ! Scalars to be diffused through CLUBB         [units vary]
   real(core_rknd) :: wp2_in(pverp)                    ! vertical velocity variance (CLUBB)           [m^2/s^2]
   real(core_rknd) :: wp3_in(pverp)                    ! third moment vertical velocity               [m^3/s^3]
   real(core_rknd) :: wpthlp_in(pverp)                 ! turbulent flux of thetal                     [K m/s]
   real(core_rknd) :: wprtp_in(pverp)                  ! turbulent flux of total water                [kg/kg m/s]
   real(core_rknd) :: rtpthlp_in(pverp)                ! covariance of thetal and qt                  [kg/kg K]
   real(core_rknd) :: rtp2_in(pverp)                   ! total water variance                         [kg^2/k^2]
   real(core_rknd) :: thlp2_in(pverp)                  ! thetal variance                              [K^2]
   real(core_rknd) :: wp2thvp_inout(pverp)                ! thermodynamic levels (< w'^2 th_v' >)        [m^2/s^2 K]
   real(core_rknd) :: wpthvp_inout(pverp)                 ! momentum levels (< w' th_v' > )              [kg/kg K]
   real(core_rknd) :: rtpthvp_inout(pverp)                ! momentum levels (< r_t' th_v' > )            [kg/kg K]
   real(core_rknd) :: rtp3_in(pverp)                   ! thermodynamic levels (r_t'^3 )               [(kg/kg)^3]
   real(core_rknd) :: thlp3_in(pverp)                  ! thermodynamic levels (th_l'^3)               [K^3]
   real(core_rknd) :: thlpthvp_inout(pverp)               ! momentum levels (< th_l' th_v' >)            [K^2]
   real(core_rknd) :: up2_in(pverp)                    ! meridional wind variance                     [m^2/s^2]
   real(core_rknd) :: vp2_in(pverp)                    ! zonal wind variance                          [m^2/s^2]
   real(core_rknd) :: upwp_in(pverp)                   ! meridional wind flux                         [m^2/s^2]
   real(core_rknd) :: vpwp_in(pverp)                   ! zonal wind flux                              [m^2/s^2]
   real(core_rknd) :: thlm_in(pverp)                   ! liquid water potential temperature (thetal)  [K]
   real(core_rknd) :: rtm_in(pverp)                    ! total water mixing ratio                     [kg/kg]
   real(core_rknd) :: rvm_in(pverp)                    ! water vapor mixing ratio                     [kg/kg]
   real(core_rknd) :: um_in(pverp)                     ! meridional wind                              [m/s]
   real(core_rknd) :: vm_in(pverp)                     ! zonal wind                                   [m/s]
   real(core_rknd) :: rho_in(pverp)                    ! mid-point density                            [kg/m^3]
   real(core_rknd) :: pre_in(pverp)                    ! input for precip evaporation
   real(core_rknd) :: rtp2_mc_out(pverp)               ! total water tendency from rain evap
   real(core_rknd) :: thlp2_mc_out(pverp)              ! thetal tendency from rain evap
   real(core_rknd) :: wprtp_mc_out(pverp)
   real(core_rknd) :: wpthlp_mc_out(pverp)
   real(core_rknd) :: rtpthlp_mc_out(pverp)
   real(core_rknd) :: rcm_inout(pverp)                 ! CLUBB output of liquid water mixing ratio     [kg/kg]
   real(core_rknd) :: rcm_out_zm(pverp)
   real(core_rknd) :: wprcp_out(pverp)                 ! CLUBB output of flux of liquid water          [kg/kg m/s]
   real(core_rknd) :: cloud_frac_inout(pverp)            ! CLUBB output of cloud fraction                [fraction]
   real(core_rknd) :: rcm_in_layer_out(pverp)          ! CLUBB output of in-cloud liq. wat. mix. ratio [kg/kg]
   real(core_rknd) :: cloud_cover_out(pverp)           ! CLUBB output of in-cloud cloud fraction       [fraction]
   real(core_rknd) :: thlprcp_out(pverp)
   real(core_rknd) :: rho_ds_zm(pverp)                 ! Dry, static density on momentum levels        [kg/m^3]
   real(core_rknd) :: rho_ds_zt(pverp)                 ! Dry, static density on thermodynamic levels   [kg/m^3]
   real(core_rknd) :: invrs_rho_ds_zm(pverp)           ! Inv. dry, static density on momentum levels   [m^3/kg]
   real(core_rknd) :: invrs_rho_ds_zt(pverp)           ! Inv. dry, static density on thermo. levels    [m^3/kg]
   real(core_rknd) :: thv_ds_zm(pverp)                 ! Dry, base-state theta_v on momentum levels    [K]
   real(core_rknd) :: thv_ds_zt(pverp)                 ! Dry, base-state theta_v on thermo. levels     [K]
   real(core_rknd) :: rfrzm(pverp)
   real(core_rknd) :: radf(pverp)
   real(core_rknd) :: wprtp_forcing(pverp)
   real(core_rknd) :: wpthlp_forcing(pverp)
   real(core_rknd) :: rtp2_forcing(pverp)
   real(core_rknd) :: thlp2_forcing(pverp)
   real(core_rknd) :: rtpthlp_forcing(pverp)
   real(core_rknd) :: ice_supersat_frac(pverp)
   real(core_rknd) :: zt_g(pverp)                      ! Thermodynamic grid of CLUBB                   [m]
   real(core_rknd) :: zi_g(pverp)                      ! Momentum grid of CLUBB                        [m]
   real(core_rknd) :: fcor                             ! Coriolis forcing                              [s^-1]
   real(core_rknd) :: sfc_elevation                    ! Elevation of ground                           [m AMSL]
   real(core_rknd) :: thlm_forcing(pverp)              ! theta_l forcing (thermodynamic levels)        [K/s]
   real(core_rknd) :: rtm_forcing(pverp)               ! r_t forcing (thermodynamic levels)            [(kg/kg)/s]
   real(core_rknd) :: um_forcing(pverp)                ! u wind forcing (thermodynamic levels)         [m/s/s]
   real(core_rknd) :: vm_forcing(pverp)                ! v wind forcing (thermodynamic levels)         [m/s/s]
   real(core_rknd) :: wm_zm(pverp)                     ! w mean wind component on momentum levels      [m/s]
   real(core_rknd) :: wm_zt(pverp)                     ! w mean wind component on thermo. levels       [m/s]
   real(core_rknd) :: p_in_Pa(pverp)                   ! Air pressure (thermodynamic levels)           [Pa]
   real(core_rknd) :: rho_zt(pverp)                    ! Air density on thermo levels                  [kt/m^3]
   real(core_rknd) :: rho_zm(pverp)                    ! Air density on momentum levels                [kg/m^3]
   real(core_rknd) :: exner(pverp)                     ! Exner function (thermodynamic levels)         [-]
   real(core_rknd) :: wpthlp_sfc                       ! w' theta_l' at surface                        [(m K)/s]
   real(core_rknd) :: wprtp_sfc                        ! w' r_t' at surface                            [(kg m)/( kg s)]
   real(core_rknd) :: upwp_sfc                         ! u'w' at surface                               [m^2/s^2]
   real(core_rknd) :: vpwp_sfc                         ! v'w' at surface                               [m^2/s^2]
   real(core_rknd) :: sclrpthvp_inout(pverp,sclr_dim)     ! momentum levels (< sclr' th_v' >)             [units vary]
   real(core_rknd) :: sclrm_forcing(pverp,sclr_dim)    ! Passive scalar forcing                        [{units vary}/s]
   real(core_rknd) :: wpsclrp_sfc(sclr_dim)            ! Scalar flux at surface                        [{units vary} m/s]
   real(core_rknd) :: edsclrm_forcing(pverp,edsclr_dim)! Eddy passive scalar forcing                   [{units vary}/s]
   real(core_rknd) :: wpedsclrp_sfc(edsclr_dim)        ! Eddy-scalar flux at surface                   [{units vary} m/s]
   real(core_rknd) :: sclrm(pverp,sclr_dim)            ! Passive scalar mean (thermo. levels)          [units vary]
   real(core_rknd) :: wpsclrp(pverp,sclr_dim)          ! w'sclr' (momentum levels)                     [{units vary} m/s]
   real(core_rknd) :: sclrp2(pverp,sclr_dim)           ! sclr'^2 (momentum levels)                     [{units vary}^2]
   real(core_rknd) :: sclrprtp(pverp,sclr_dim)         ! sclr'rt' (momentum levels)                    [{units vary} (kg/kg)]
   real(core_rknd) :: sclrpthlp(pverp,sclr_dim)        ! sclr'thlp' (momentum levels)                  [{units vary} (K)]
   real(core_rknd) :: hydromet(pverp,hydromet_dim)
   real(core_rknd) :: wphydrometp(pverp,hydromet_dim)
   real(core_rknd) :: wp2hmp(pverp,hydromet_dim)
   real(core_rknd) :: rtphmp_zt(pverp,hydromet_dim)
   real(core_rknd) :: thlphmp_zt (pverp,hydromet_dim)
   real(core_rknd) :: C_10                             ! transfer coefficient                          [-]
   real(core_rknd) :: khzm_out(pverp)                  ! eddy diffusivity on momentum grids            [m^2/s]
   real(core_rknd) :: khzt_out(pverp)                  ! eddy diffusivity on thermo grids              [m^2/s]
   real(core_rknd) :: qclvar_out(pverp)                ! cloud water variance                          [kg^2/kg^2]
   real(core_rknd) :: varmu2
   real(core_rknd) :: qrl_clubb(pverp)
   real(core_rknd) :: qrl_zm(pverp)
   real(core_rknd) :: thlp2_rad_out(pverp)

   real(core_rknd), dimension(nparams)  :: clubb_params ! These adjustable CLUBB parameters (C1, C2 ...)
   real(core_rknd), dimension(sclr_dim) :: sclr_tol     ! Tolerance on passive scalar       [units vary]

   real(core_rknd) :: dum_core_rknd                    ! dummy variable  [units vary]
   real(core_rknd) :: hdtime_core_rknd                  ! host model time step in core_rknd

   real(core_rknd), pointer :: upwp_sfc_pert    ! u'w' at surface                               [m^2/s^2]
   real(core_rknd), pointer :: vpwp_sfc_pert    ! v'w' at surface                               [m^2/s^2]
   ! Pointers to temporary copies of particular columns of um_pert/upwp_pert fields
   real(core_rknd), pointer, dimension(:) :: um_pert_col ! Pointer to a particular column of um
   real(core_rknd), pointer, dimension(:) :: vm_pert_col
   real(core_rknd), pointer, dimension(:) :: upwp_pert_col
   real(core_rknd), pointer, dimension(:) :: vpwp_pert_col

   !===========================================================================================================================
   ! End of defining the variables for the change of precision in the CLUBB calculation
   !===========================================================================================================================

   real(r8) :: apply_const
   real(r8) :: qclvar(pcols,pverp)              ! cloud water variance                          [kg^2/kg^2]
   real(r8) :: newfice(pcols,pver)              ! fraction of ice in cloud at CLUBB start       [-]
   real(r8) :: bflx22                           ! Variable for buoyancy flux for pbl            [K m/s]
   real(r8) :: invrs_hdtime                     ! Preculate 1/hdtime to reduce divide operations
   real(r8) :: invrs_gravit                     ! Preculate 1/gravit to reduce divide operations
   real(r8) :: ubar                             ! surface wind                                  [m/s]
   real(r8) :: ustar                            ! surface stress                                [m/s]
   real(r8) :: z0                               ! roughness height                              [m]
   real(r8) :: zo                               ! roughness height                              [m]
   real(r8) :: dz_g(pver)                       ! thickness of layer                            [m]
   real(r8) :: minqn                            ! minimum total cloud liquid + ice threshold    [kg/kg]
   real(r8) :: tempqn                           ! temporary total cloud liquid + ice            [kg/kg]
   real(r8) :: cldthresh                        ! threshold to determin cloud fraction          [kg/kg]
   real(r8) :: relvarmax,relvarmin
   real(r8) :: qmin
   real(r8) :: varmu(pcols)
   real(r8) :: zt_out(pcols,pverp)              ! output for the thermo CLUBB grid              [m]
   real(r8) :: zi_out(pcols,pverp)              ! output for momentum CLUBB grid                [m]

   real(r8) :: pdf_zm_w_1_inout(pverp)          ! work array for pdf_params_zm%w_1
   real(r8) :: pdf_zm_w_2_inout(pverp)          ! work array for pdf_params_zm%w_2
   real(r8) :: pdf_zm_varnce_w_1_inout(pverp)   ! work array for pdf_params_zm%varnce_w_1
   real(r8) :: pdf_zm_varnce_w_2_inout(pverp)   ! work array for pdf_params_zm%varnce_w_2
   real(r8) :: pdf_zm_mixt_frac_inout(pverp)    ! work array for pdf_params_zm%mixt_frac

   ! Variables below are needed to compute energy integrals for conservation
   real(r8) :: ke_a(pcols), ke_b(pcols), te_a(pcols), te_b(pcols)
   real(r8) :: wv_a(pcols), wv_b(pcols), wl_b(pcols), wl_a(pcols)
   real(r8) :: se_dis, se_a(pcols), se_b(pcols), clubb_s(pver), enthalpy

   real(r8) :: exner_clubb(pcols,pverp)         ! Exner function consistent with CLUBB          [-]
   real(r8) :: wpthlp_output(pcols,pverp)       ! Heat flux output variable                     [W/m2]
   real(r8) :: wprtp_output(pcols,pverp)        ! Total water flux output variable              [W/m2]
   real(r8) :: wp3_output(pcols,pverp)          ! wp3 output                                    [m^3/s^3]
   real(r8) :: rtpthlp_output(pcols,pverp)      ! rtpthlp ouptut                                [K kg/kg]
   real(r8) :: qt_output(pcols,pver)            ! Total water mixing ratio for output           [kg/kg]
   real(r8) :: thetal_output(pcols,pver)        ! Liquid water potential temperature output     [K]
   real(r8) :: sl_output(pcols,pver)            ! Liquid water static energy                    [J/kg]
   real(r8) :: ustar2(pcols)                    ! Surface stress for PBL height                 [m2/s2]
   real(r8) :: rho(pcols,pverp)                 ! Midpoint density in CAM                       [kg/m^3]
   real(r8) :: thv(pcols,pver)                  ! virtual potential temperature                 [K]
   real(r8) :: edsclr_out(pverp,edsclr_dim)     ! Scalars to be diffused through CLUBB          [units vary]
   real(r8) :: sclrpthvp(pcols,pverp,sclr_dim)  ! sclr'th_v' (momentum levels)                  [{units vary} K]
   real(r8) :: rcm_in_layer(pcols,pverp)        ! CLUBB in-cloud liquid water mixing ratio      [kg/kg]
   real(r8) :: cloud_cover(pcols,pverp)         ! CLUBB in-cloud cloud fraction                 [fraction]
   real(r8) :: wprcp(pcols,pverp)               ! CLUBB liquid water flux                       [m/s kg/kg]
   real(r8) :: wpthvp_diag(pcols,pverp)              ! CLUBB buoyancy flux                           [W/m^2]
   real(r8) :: rvm(pcols,pverp)
   real(r8) :: dlf2(pcols,pver)                 ! Detraining cld H20 from shallow convection    [kg/kg/day]
   real(r8) :: eps                              ! Rv/Rd                                         [-]
   real(r8) :: dum1                             ! dummy variable                                [units vary]
   real(r8) :: obklen(pcols)                    ! Obukov length                                 [m]
   real(r8) :: kbfs(pcols)                      ! Kinematic Surface heat flux                   [K m/s]
   real(r8) :: th(pcols,pver)                   ! potential temperature                         [K]
   real(r8) :: dummy2(pcols)                    ! dummy variable                                [units vary]
   real(r8) :: dummy3(pcols)                    ! dummy variable                                [units vary]
   real(r8) :: kinheat(pcols)                   ! Kinematic Surface heat flux                   [K m/s]
   real(r8) :: ksrftms(pcols)                   ! Turbulent mountain stress surface drag        [kg/s/m2]
   real(r8) :: tautmsx(pcols)                   ! U component of turbulent mountain stress      [N/m2]
   real(r8) :: tautmsy(pcols)                   ! V component of turbulent mountain stress      [N/m2]
   real(r8) :: rrho                             ! Inverse of air density                        [1/kg/m^3]
   real(r8) :: kinwat(pcols)                    ! Kinematic water vapor flux                    [m/s]
   real(r8) :: latsub

   integer  :: ktop(pcols,pver)
   integer  :: ncvfin(pcols)
   real(r8) :: chs(pcols,pverp)
   real(r8) :: lwp_CL(pver)
   real(r8) :: opt_depth_CL(pver)
   real(r8) :: radinvfrac_CL(pver)
   real(r8) :: radf_CL(pver)
   real(r8) :: radf_out(pver)
   real(r8) :: es(pcols,pver)
   real(r8) :: qs(pcols,pver)
   real(r8) :: gam(pcols,pver)
   real(r8) :: tmp_array(state%ncol,pverp)
   real(r8) :: bfact, orgparam, delpavg
   character(len=6) :: choice_radf

   integer                               :: time_elapsed                ! time keep track of stats          [s]
   character(len=200)                    :: temp1, sub                  ! Strings needed for CLUBB output
   logical                               :: l_Lscale_plume_centered, l_use_ice_latent
   character(len=3), dimension(pcnst)    :: cnst_type_loc               ! local override option for constituents cnst_type


   ! --------------- !
   ! Pointers        !
   ! --------------- !

   real(r8), pointer, dimension(:,:) :: wp2      ! vertical velocity variance                   [m^2/s^2]
   real(r8), pointer, dimension(:,:) :: wp3      ! third moment of vertical velocity            [m^3/s^3]
   real(r8), pointer, dimension(:,:) :: wpthlp   ! turbulent flux of thetal                     [m/s K]
   real(r8), pointer, dimension(:,:) :: wprtp    ! turbulent flux of moisture                   [m/s kg/kg]
   real(r8), pointer, dimension(:,:) :: rtpthlp  ! covariance of thetal and qt                  [kg/kg K]
   real(r8), pointer, dimension(:,:) :: rtp2     ! moisture variance                            [kg^2/kg^2]
   real(r8), pointer, dimension(:,:) :: thlp2    ! temperature variance                         [K^2]
   real(r8), pointer, dimension(:,:) :: up2      ! east-west wind variance                      [m^2/s^2]
   real(r8), pointer, dimension(:,:) :: vp2      ! north-south wind variance                    [m^2/s^2]

   real(r8), pointer, dimension(:,:) :: wpthvp     ! < w'th_v' > (momentum levels)                [m/s K]
   real(r8), pointer, dimension(:,:) :: wp2thvp    ! < w'^2 th_v' > (thermodynamic levels)        [m^2/s^2 K]
   real(r8), pointer, dimension(:,:) :: rtpthvp    ! < r_t'th_v' > (momentum levels)              [kg/kg K]
   real(r8), pointer, dimension(:,:) :: thlpthvp   ! < th_l'th_v' > (momentum levels)             [K^2]
   real(r8), pointer, dimension(:,:) :: rcm        ! CLUBB cloud water mixing ratio               [kg/kg]
   real(r8), pointer, dimension(:,:) :: cloud_frac ! CLUBB cloud fraction                       [-]
   real(r8), pointer, dimension(:,:) :: pdf_zm_w_1(:,:)        !work pointer for pdf_params_zm
   real(r8), pointer, dimension(:,:) :: pdf_zm_w_2(:,:)        !work pointer for pdf_params_zm
   real(r8), pointer, dimension(:,:) :: pdf_zm_varnce_w_1(:,:) !work pointer for pdf_params_zm
   real(r8), pointer, dimension(:,:) :: pdf_zm_varnce_w_2(:,:) !work pointer for pdf_params_zm
   real(r8), pointer, dimension(:,:) :: pdf_zm_mixt_frac(:,:)  !work pointer for pdf_params_zm

   real(r8), pointer, dimension(:,:) :: upwp     ! east-west momentum flux                      [m^2/s^2]
   real(r8), pointer, dimension(:,:) :: vpwp     ! north-south momentum flux                    [m^2/s^2]
   real(r8), pointer, dimension(:,:) :: um_pert  ! perturbed meridional wind                    [m/s]
   real(r8), pointer, dimension(:,:) :: vm_pert  ! perturbed zonal wind                         [m/s]
   real(r8), pointer, dimension(:,:) :: upwp_pert! perturbed meridional wind flux               [m^2/s^2]
   real(r8), pointer, dimension(:,:) :: vpwp_pert! perturbed zonal wind flux                    [m^2/s^2]
   real(r8), pointer, dimension(:,:) :: thlm     ! mean temperature                             [K]
   real(r8), pointer, dimension(:,:) :: rtm      ! mean moisture mixing ratio                   [kg/kg]
   real(r8), pointer, dimension(:,:) :: um       ! mean east-west wind                          [m/s]
   real(r8), pointer, dimension(:,:) :: vm       ! mean north-south wind                        [m/s]
   real(r8), pointer, dimension(:,:) :: cld      ! cloud fraction                               [fraction]
   real(r8), pointer, dimension(:,:) :: concld   ! convective cloud fraction                    [fraction]
   real(r8), pointer, dimension(:,:) :: ast      ! stratiform cloud fraction                    [fraction]
   real(r8), pointer, dimension(:,:) :: alst     ! liquid stratiform cloud fraction             [fraction]
   real(r8), pointer, dimension(:,:) :: aist     ! ice stratiform cloud fraction                [fraction]
   real(r8), pointer, dimension(:,:) :: qlst     ! Physical in-stratus LWC                      [kg/kg]
   real(r8), pointer, dimension(:,:) :: qist     ! Physical in-stratus IWC                      [kg/kg]
   real(r8), pointer, dimension(:,:) :: deepcu   ! deep convection cloud fraction               [fraction]
   real(r8), pointer, dimension(:,:) :: shalcu   ! shallow convection cloud fraction            [fraction]
   real(r8), pointer, dimension(:,:) :: khzt     ! eddy diffusivity on thermo levels            [m^2/s]
   real(r8), pointer, dimension(:,:) :: khzm     ! eddy diffusivity on momentum levels          [m^2/s]
   real(r8), pointer, dimension(:) :: pblh     ! planetary boundary layer height                [m]
   real(r8), pointer, dimension(:,:) :: tke      ! turbulent kinetic energy                     [m^2/s^2]
   real(r8), pointer, dimension(:,:) :: dp_icwmr ! deep convection in cloud mixing ratio        [kg/kg]
   real(r8), pointer, dimension(:,:) :: relvar   ! relative cloud water variance                [-]
   real(r8), pointer, dimension(:,:) :: accre_enhan ! accretion enhancement factor              [-]
   real(r8), pointer, dimension(:,:) :: cmeliq
   real(r8), pointer, dimension(:,:) :: cmfmc_sh ! Shallow convective mass flux--m subc (pcols,pverp) [kg/m2/s/]

   type(pdf_parameter), pointer :: pdf_params    ! PDF parameters (thermo. levs.) [units vary]
   type(pdf_parameter), pointer :: pdf_params_zm ! PDF parameters on momentum levs. [units vary]

   real(r8), pointer, dimension(:,:) :: naai
   real(r8), pointer, dimension(:,:) :: prer_evap
   real(r8), pointer, dimension(:,:) :: qrl
   real(r8), pointer, dimension(:,:) :: radf_clubb

   ! ZM convective microphysics(zm_microp)
   real(r8), pointer :: dlfzm(:,:)  ! ZM detrainment of convective cloud water mixing ratio.
   real(r8), pointer :: difzm(:,:)  ! ZM detrainment of convective cloud ice mixing ratio.
   real(r8), pointer :: dsfzm(:,:)  ! ZM detrainment of convective snow mixing ratio.
   real(r8), pointer :: dnlfzm(:,:) ! ZM detrainment of convective cloud water num concen.
   real(r8), pointer :: dnifzm(:,:) ! ZM detrainment of convective cloud ice num concen.
   real(r8), pointer :: dnsfzm(:,:) ! ZM detrainment of convective snow num concen.

!PMA
   real(r8)  relvarc(pcols,pver)
   real(r8)  stend(pcols,pver)
   real(r8)  qvtend(pcols,pver)
   real(r8)  qitend(pcols,pver)
   real(r8)  initend(pcols,pver)
   logical            :: lqice(pcnst)
   integer :: ktopi(pcols)

   integer :: ixorg

   character(len=2) :: char_macmic_it

   intrinsic :: selected_real_kind, max

!PMA adds gustiness and tpert
   real(r8), pointer :: prec_dp(:)                 ! total precipitation from ZM convection
   real(r8), pointer :: snow_dp(:)                 ! snow precipitation from ZM convection
   real(r8), pointer :: vmag_gust(:)
   real(r8), pointer :: wsresp(:)
   real(r8), pointer :: tau_est(:)
   real(r8), pointer :: tpert(:)

   real(r8) :: ugust  ! function: gustiness as a function of convective rainfall
   real(r8) :: gfac
   real(r8) :: gprec
   real(r8) :: prec_gust(pcols)
   real(r8) :: vmag_gust_dp(pcols),vmag_gust_cl(pcols)
   real(r8) :: gust_fac(pcols)
   real(r8) :: umb(pcols), vmb(pcols),up2b(pcols),vp2b(pcols)
   real(r8),parameter :: gust_facl = 1.2_r8 !gust fac for land
   real(r8),parameter :: gust_faco = 0.9_r8 !gust fac for ocean
   real(r8),parameter :: gust_facc = 1.5_r8 !gust fac for clubb

   real(r8) :: sfc_v_diff_tau(pcols) ! Response to tau perturbation, m/s
   real(r8), parameter :: pert_tau = 0.1_r8 ! tau perturbation, Pa

   !variables for turbulent orographic form drag (TOFD) interface
   real(r8) :: dtaux3_fd(pcols,pver)
   real(r8) :: dtauy3_fd(pcols,pver)
   real(r8) :: dusfc_fd(pcols)
   real(r8) :: dvsfc_fd(pcols)
   logical  :: gwd_ls,gwd_bl,gwd_ss,gwd_fd
   real(r8) :: dummy_nm(pcols,pver)
   real(r8) :: dummy_utgw(pcols,pver)
   real(r8) :: dummy_vtgw(pcols,pver)
   real(r8) :: dummy_ttgw(pcols,pver)
   real(r8) :: dummx_ls(pcols,pver)
   real(r8) :: dummx_bl(pcols,pver)
   real(r8) :: dummx_ss(pcols,pver)
   real(r8) :: dummy_ls(pcols,pver)
   real(r8) :: dummy_bl(pcols,pver)
   real(r8) :: dummy_ss(pcols,pver)
   real(r8) :: dummx3_ls(pcols,pver)
   real(r8) :: dummx3_bl(pcols,pver)
   real(r8) :: dummx3_ss(pcols,pver)
   real(r8) :: dummy3_ls(pcols,pver)
   real(r8) :: dummy3_bl(pcols,pver)
   real(r8) :: dummy3_ss(pcols,pver)

   real(r8) :: inv_exner_clubb_surf


! ZM gustiness equation below from Redelsperger et al. (2000)
! numbers are coefficients of the empirical equation

   ugust(gprec) = log(1._R8+57801.6_R8*gprec-3.55332096e7_R8*(gprec**2.0_R8))

#endif
   det_s(:)   = 0.0_r8
   det_ice(:) = 0.0_r8

   if (macmic_it > 99) then
      call endrun('clubb_tend_cam: macmic_it > 99. Revise checkpoint name for cnd_diag_checkpoint.')
   end if

#ifdef CLUBB_SGS

   !-----------------------------------------------------------------------------------------------!
   !-----------------------------------------------------------------------------------------------!
   !-----------------------------------------------------------------------------------------------!
   !       MAIN COMPUTATION BEGINS HERE                                                            !
   !-----------------------------------------------------------------------------------------------!
   !-----------------------------------------------------------------------------------------------!
   !-----------------------------------------------------------------------------------------------!

   write(char_macmic_it,'(i2.2)') macmic_it

   call t_startf('clubb_tend_cam_init')
   invrs_hdtime = 1._r8 / hdtime
   invrs_gravit = 1._r8 / gravit
   frac_limit = 0.01_r8
   ic_limit   = 1.e-12_r8

   if (clubb_do_adv) then
     apply_const = 1._r8  ! Initialize to one, only if CLUBB's moments are advected
   else
     apply_const = 0._r8  ! Never want this if CLUBB's moments are not advected
   endif

   !  Define forcings from CAM to CLUBB as zero for momentum and thermo,
   !  forcings already applied through CAM
   thlm_forcing(1:pverp) = 0._core_rknd
   rtm_forcing(1:pverp)  = 0._core_rknd
   um_forcing(1:pverp)   = 0._core_rknd
   vm_forcing(1:pverp)   = 0._core_rknd

   wprtp_forcing(1:pverp)   = 0._core_rknd
   wpthlp_forcing(1:pverp)  = 0._core_rknd
   rtp2_forcing(1:pverp)    = 0._core_rknd
   thlp2_forcing(1:pverp)   = 0._core_rknd
   rtpthlp_forcing(1:pverp) = 0._core_rknd

   ! rtp3_in and thlp3_in are not currently used in CLUBB's default code.
   rtp3_in(:)  = 0.0_core_rknd
   thlp3_in(:) = 0.0_core_rknd

   !  Define surface sources for transported variables for diffusion, will
   !  be zero as these tendencies are done in clubb_surface
   do ixind=1,edsclr_dim
      wpedsclrp_sfc(ixind) = 0._core_rknd
   enddo

   ice_supersat_frac(1:pverp) = 0._core_rknd

   !  Higher order scalar inputs, set to zero
   wpsclrp_sfc(:)      = 0._core_rknd
   hydromet(:,:)       = 0._core_rknd
   wphydrometp(:,:)    = 0._core_rknd
   wp2hmp(:,:)         = 0._core_rknd
   rtphmp_zt(:,:)      = 0._core_rknd
   thlphmp_zt(:,:)     = 0._core_rknd

   !  Initialize forcings for transported scalars to zero
   sclrm_forcing(:,:)   = 0._core_rknd
   edsclrm_forcing(:,:) = 0._core_rknd

   !  Determine Coriolis force at given latitude.  This is never used
   !  when CLUBB is implemented in a host model, therefore just set
   !  to zero.
   fcor = 0._core_rknd

 !  Get indicees for cloud and ice mass and cloud and ice number

   call cnst_get_ind('Q',ixq)
   call cnst_get_ind('CLDLIQ',ixcldliq)
   call cnst_get_ind('CLDICE',ixcldice)
   call cnst_get_ind('NUMLIQ',ixnumliq)
   call cnst_get_ind('NUMICE',ixnumice)

 !  Initialize physics tendency arrays, copy the state to state1 array to use in this routine

   if (.not. micro_do_icesupersat) then
     call physics_ptend_init(ptend_loc,state%psetcols, 'clubb_ice1', ls=.true., lu=.true., lv=.true., lq=lq)
   endif

   call physics_state_copy(state,state1)

   ! constituents are all treated as wet mmr by clubb
   ! don't convert co2 tracers to wet mixing ratios
   cnst_type_loc(:) = cnst_type(:)
   call co2_cycle_set_cnst_type(cnst_type_loc, 'wet')
   call set_dry_to_wet(state1, cnst_type_loc)

   if (micro_do_icesupersat) then
     naai_idx      = pbuf_get_index('NAAI')
     call pbuf_get_field(pbuf, naai_idx, naai)
     call physics_ptend_init(ptend_all, state%psetcols, 'clubb_ice2')
   endif

   !  Determine number of columns and which chunk computation is to be performed on

   ncol = state%ncol
   lchnk = state%lchnk

   !  Determine time step of physics buffer

   itim_old = pbuf_old_tim_idx()

   !  Establish associations between pointers and physics buffer fields

   call pbuf_get_field(pbuf, wp2_idx,     wp2,     start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, wp3_idx,     wp3,     start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, wpthlp_idx,  wpthlp,  start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, wprtp_idx,   wprtp,   start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, rtpthlp_idx, rtpthlp, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, rtp2_idx,    rtp2,    start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, thlp2_idx,   thlp2,   start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, up2_idx,     up2,     start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, vp2_idx,     vp2,     start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))

   call pbuf_get_field(pbuf, upwp_idx,    upwp,    start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, vpwp_idx,    vpwp,    start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   if (linearize_pbl_winds) then
      call pbuf_get_field(pbuf, um_pert_idx, um_pert, start=(/1,1/),          kount=(/pcols,pverp/))
      call pbuf_get_field(pbuf, vm_pert_idx, vm_pert, start=(/1,1/),          kount=(/pcols,pverp/))
      call pbuf_get_field(pbuf, upwp_pert_idx,upwp_pert,start=(/1,1/),        kount=(/pcols,pverp/))
      call pbuf_get_field(pbuf, vpwp_pert_idx,vpwp_pert,start=(/1,1/),        kount=(/pcols,pverp/))
   else
      nullify(um_pert)
      nullify(vm_pert)
      nullify(upwp_pert)
      nullify(vpwp_pert)
   end if
   call pbuf_get_field(pbuf, thlm_idx,    thlm,    start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, rtm_idx,     rtm,     start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, um_idx,      um,      start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, vm_idx,      vm,      start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))

   call pbuf_get_field(pbuf, wpthvp_idx,     wpthvp,     start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, wp2thvp_idx,    wp2thvp,    start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, rtpthvp_idx,    rtpthvp,    start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, thlpthvp_idx,   thlpthvp,   start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, rcm_idx,        rcm,        start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, cloud_frac_idx, cloud_frac, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))

   call pbuf_get_field(pbuf, pdf_zm_w_1_idx, pdf_zm_w_1, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, pdf_zm_w_2_idx, pdf_zm_w_2, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, pdf_zm_varnce_w_1_idx, pdf_zm_varnce_w_1, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, pdf_zm_varnce_w_2_idx, pdf_zm_varnce_w_2, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, pdf_zm_mixt_frac_idx, pdf_zm_mixt_frac, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))

   call pbuf_get_field(pbuf, tke_idx,     tke)
   call pbuf_get_field(pbuf, qrl_idx,     qrl)
   call pbuf_get_field(pbuf, radf_idx,    radf_clubb)

   call pbuf_get_field(pbuf, cld_idx,     cld,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, concld_idx,  concld,  start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, ast_idx,     ast,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, alst_idx,    alst,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, aist_idx,    aist,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, qlst_idx,    qlst,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, qist_idx,    qist,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   call pbuf_get_field(pbuf, prer_evap_idx, prer_evap)
   call pbuf_get_field(pbuf, accre_enhan_idx, accre_enhan)
   call pbuf_get_field(pbuf, cmeliq_idx,  cmeliq)
   call pbuf_get_field(pbuf, relvar_idx,  relvar)
   call pbuf_get_field(pbuf, dp_frac_idx, deepcu)
   call pbuf_get_field(pbuf, sh_frac_idx, shalcu)
   call pbuf_get_field(pbuf, kvm_idx,     khzt)
   call pbuf_get_field(pbuf, kvh_idx,     khzm)
   call pbuf_get_field(pbuf, pblh_idx,    pblh)
   call pbuf_get_field(pbuf, icwmrdp_idx, dp_icwmr)
   call pbuf_get_field(pbuf, cmfmc_sh_idx, cmfmc_sh)

!PMA adds fields for gustiness
   call pbuf_get_field(pbuf, tpert_idx,   tpert)
   call pbuf_get_field(pbuf, prec_dp_idx, prec_dp)
   call pbuf_get_field(pbuf, snow_dp_idx, snow_dp)
   call pbuf_get_field(pbuf, vmag_gust_idx, vmag_gust)

   if (linearize_pbl_winds) then
      call pbuf_get_field(pbuf, wsresp_idx, wsresp)
      call pbuf_get_field(pbuf, tau_est_idx, tau_est)
   end if

   ! Intialize the apply_const variable (note special logic is due to eularian backstepping)
   if (clubb_do_adv .and. (is_first_step() .or. all(wpthlp(1:ncol,1:pver) .eq. 0._r8))) then
      apply_const = 0._r8  ! On first time through do not remove constant
                           !  from moments since it has not been added yet
   endif

   if (micro_do_icesupersat) then

     ! -------------------------------------- !
     ! Ice Saturation Adjustment Computation  !
     ! -------------------------------------- !

     lq2(:)  = .FALSE.
     lq2(1)  = .TRUE.
     lq2(ixcldice) = .TRUE.
     lq2(ixnumice) = .TRUE.

     latsub = latvap + latice

     call physics_ptend_init(ptend_loc, state%psetcols, 'iceadj', ls=.true., lq=lq2 )

     stend(:ncol,:)=0._r8
     qvtend(:ncol,:)=0._r8
     qitend(:ncol,:)=0._r8
     initend(:ncol,:)=0._r8

     call t_startf('ice_macro_tend')
     call ice_macro_tend(naai(:ncol,top_lev:pver),state1%t(:ncol,top_lev:pver), &
        state1%pmid(:ncol,top_lev:pver),state1%q(:ncol,top_lev:pver,1),state1%q(:ncol,top_lev:pver,ixcldice),&
        state1%q(:ncol,top_lev:pver,ixnumice),latsub,hdtime,&
        stend(:ncol,top_lev:pver),qvtend(:ncol,top_lev:pver),qitend(:ncol,top_lev:pver),&
        initend(:ncol,top_lev:pver))
     call t_stopf('ice_macro_tend')

     ! update local copy of state with the tendencies
     ptend_loc%q(:ncol,top_lev:pver,1)=qvtend(:ncol,top_lev:pver)
     ptend_loc%q(:ncol,top_lev:pver,ixcldice)=qitend(:ncol,top_lev:pver)
     ptend_loc%q(:ncol,top_lev:pver,ixnumice)=initend(:ncol,top_lev:pver)
     ptend_loc%s(:ncol,top_lev:pver)=stend(:ncol,top_lev:pver)

    ! Add the ice tendency to the output tendency
     call physics_ptend_sum(ptend_loc, ptend_all, ncol)

    ! ptend_loc is reset to zero by this call
     call physics_update(state1, ptend_loc, hdtime)

    !Write output for tendencies:
    !        oufld: QVTENDICE,QITENDICE,NITENDICE
     call outfld( 'TTENDICE',  stend/cpair, pcols, lchnk )
     call outfld( 'QVTENDICE', qvtend, pcols, lchnk )
     call outfld( 'QITENDICE', qitend, pcols, lchnk )
     call outfld( 'NITENDICE', initend, pcols, lchnk )

   endif

   call cnd_diag_checkpoint(diag, 'ICEMAC'//char_macmic_it, state1, pbuf, cam_in, cam_out)

   !  Determine CLUBB time step and make it sub-step friendly
   !  For now we want CLUBB time step to be 5 min since that is
   !  what has been scientifically validated.  However, there are certain
   !  instances when a 5 min time step will not be possible (based on
   !  host model time step or on macro-micro sub-stepping

   dtime = clubb_timestep
   hdtime_core_rknd = real(hdtime, kind = core_rknd)

   !  Now check to see if dtime is greater than the host model
   !    (or sub stepped) time step.  If it is, then simply
   !    set it equal to the host (or sub step) time step.
   !    This section is mostly to deal with small host model
   !    time steps (or small sub-steps)

   if (dtime .gt. hdtime_core_rknd) then
     dtime = hdtime_core_rknd
   endif

   !  Now check to see if CLUBB time step divides evenly into
   !    the host model time step.  If not, force it to divide evenly.
   !    We also want it to be 5 minutes or less.  This section is
   !    mainly for host model time steps that are not evenly divisible
   !    by 5 minutes

   if (mod(hdtime_core_rknd,dtime) .ne. 0) then
     dtime = hdtime_core_rknd/2._core_rknd
     do while (dtime .gt. 300._core_rknd)
       dtime = dtime/2._core_rknd
     end do
   endif

   !  If resulting host model time step and CLUBB time step do not divide evenly
   !    into each other, have model throw a fit.

   if (mod(hdtime_core_rknd,dtime) .ne. 0) then
     call endrun('clubb_tend_cam:  CLUBB time step and HOST time step NOT compatible'//errmsg(__FILE__,__LINE__))
   endif

   !  determine number of timesteps CLUBB core should be advanced,
   !  host time step divided by CLUBB time step
   nadv = max(hdtime_core_rknd/dtime,1._core_rknd)

   minqn = 0._r8
   newfice(:,:) = 0._r8
   where(state1%q(:ncol,:pver,3) .gt. minqn) &
       newfice(:ncol,:pver) = state1%q(:ncol,:pver,3)/(state1%q(:ncol,:pver,2)+state1%q(:ncol,:pver,3))

   !  Compute exner function consistent with CLUBB's definition, which uses a constant
   !  surface pressure.  CAM's exner (in state does not).  Therefore, for consistent
   !  treatment with CLUBB code, anytime exner is needed to treat CLUBB variables
   !  (such as thlm), use "exner_clubb" other wise use the exner in state

   do k=1,pver
     do i=1,ncol
       exner_clubb(i,k) = (real(p0_clubb, kind = r8 )/state1%pmid(i,k))**(rair/cpair)
     enddo
   enddo

   !  At each CLUBB call, initialize mean momentum  and thermo CLUBB state
   !  from the CAM state

   rvm = 0._r8
   do k=1,pver   ! loop over levels
     do i=1,ncol ! loop over columns

       rtm(i,k)     = state1%q(i,k,ixq)+state1%q(i,k,ixcldliq)
       rvm(i,k)     = state1%q(i,k,ixq)
       um(i,k)      = state1%u(i,k)
       vm(i,k)      = state1%v(i,k)

#define NEWTHETAL
#ifndef NEWTHETAL
       thlm(i,k)    = state1%t(i,k)*exner_clubb(i,k)-(latvap/cpair)*state1%q(i,k,ixcldliq)
#else
!NCAR
!       thlm(i,k) = ( state1%t(i,k) &
!                     - (latvap/cpairv(i,k,lchnk))*state1%q(i,k,ixcldliq) ) &
!                   * inv_exner_clubb(i,k)

       thlm(i,k) = ( state1%t(i,k) &
                     - (latvap/cpair)*state1%q(i,k,ixcldliq) ) &
                   * exner_clubb(i,k)
#endif

       if (clubb_do_adv) then
          if (macmic_it .eq. 1) then

            !  Note that some of the moments below can be positive or negative.
            !    Remove a constant that was added to prevent dynamics from clipping
            !    them to prevent dynamics from making them positive.
            thlp2(i,k)   = state1%q(i,k,ixthlp2)
            rtp2(i,k)    = state1%q(i,k,ixrtp2)
            rtpthlp(i,k) = state1%q(i,k,ixrtpthlp) - (rtpthlp_const*apply_const)
            wpthlp(i,k)  = state1%q(i,k,ixwpthlp) - (wpthlp_const*apply_const)
            wprtp(i,k)   = state1%q(i,k,ixwprtp) - (wprtp_const*apply_const)
            wp2(i,k)     = state1%q(i,k,ixwp2)
            wp3(i,k)     = state1%q(i,k,ixwp3) - (wp3_const*apply_const)
            up2(i,k)     = state1%q(i,k,ixup2)
            vp2(i,k)     = state1%q(i,k,ixvp2)
          endif
       endif

     enddo
   enddo

   if (clubb_do_adv) then
     ! If not last step of macmic loop then set apply_const back to
     !   zero to prevent output from being corrupted.
     if (macmic_it .eq. cld_macmic_num_steps) then
       apply_const = 1._r8
     else
       apply_const = 0._r8
     endif
   endif

   rtm(1:ncol,pverp)  = rtm(1:ncol,pver)
   um(1:ncol,pverp)   = state1%u(1:ncol,pver)
   vm(1:ncol,pverp)   = state1%v(1:ncol,pver)
   thlm(1:ncol,pverp) = thlm(1:ncol,pver)

   if (clubb_do_adv) then
      thlp2(1:ncol,pverp)=thlp2(1:ncol,pver)
      rtp2(1:ncol,pverp)=rtp2(1:ncol,pver)
      rtpthlp(1:ncol,pverp)=rtpthlp(1:ncol,pver)
      wpthlp(1:ncol,pverp)=wpthlp(1:ncol,pver)
      wprtp(1:ncol,pverp)=wprtp(1:ncol,pver)
      wp2(1:ncol,pverp)=wp2(1:ncol,pver)
      wp3(1:ncol,pverp)=wp3(1:ncol,pver)
      up2(1:ncol,pverp)=up2(1:ncol,pver)
      vp2(1:ncol,pverp)=vp2(1:ncol,pver)
   endif

   ! Compute integrals of static energy, kinetic energy, water vapor, and liquid water
   ! for the computation of total energy before CLUBB is called.  This is for an
   ! effort to conserve energy since liquid water potential temperature (which CLUBB
   ! conserves) and static energy (which CAM conserves) are not exactly equal.
   se_b = 0._r8
   ke_b = 0._r8
   wv_b = 0._r8
   wl_b = 0._r8
   do k=1,pver
     do i=1,ncol
       ! use s=c_pT+g*z, total energy needs term c_pT but not gz
       se_b(i) = se_b(i) + (state1%s(i,k) - gravit*state1%zm(i,k) - state1%phis(i)) &
                         *  state1%pdel(i,k)*invrs_gravit
       ke_b(i) = ke_b(i) + 0.5_r8*(um(i,k)**2+vm(i,k)**2)*state1%pdel(i,k)*invrs_gravit
       wv_b(i) = wv_b(i) + state1%q(i,k,ixq)*state1%pdel(i,k)*invrs_gravit
       wl_b(i) = wl_b(i) + state1%q(i,k,ixcldliq)*state1%pdel(i,k)*invrs_gravit
     enddo
   enddo

   !  Compute virtual potential temperature, which is needed for CLUBB
   do k=1,pver
     do i=1,ncol
       thv(i,k) = state1%t(i,k)*exner_clubb(i,k)*(1._r8+zvir*state1%q(i,k,ixq)&
                  -state1%q(i,k,ixcldliq))
     enddo
   enddo
   call t_stopf('clubb_tend_cam_init')

   ! ------------------------------------------------- !
   ! Begin module to compute turbulent mountain stress !
   ! ------------------------------------------------- !
    if ( do_tms) then
       call t_startf('compute_tms')
       call compute_tms( pcols,        pver,      ncol,                   &
                     state1%u,     state1%v,  state1%t,  state1%pmid, &
                     state1%exner, state1%zm, sgh30,     ksrftms,     &
                     tautmsx,      tautmsy,   cam_in%landfrac )
       call t_stopf('compute_tms')
    endif

    if (use_od_fd) then
        gwd_ls    =.false.
        gwd_bl    =.false.
        gwd_ss    =.false.
        gwd_fd    =use_od_fd
        dummy_nm  =0.0_r8
        dummy_utgw=0.0_r8
        dummy_vtgw=0.0_r8
        dummy_ttgw=0.0_r8
        !sgh30 as the input for turbulent orographic form drag (TOFD) instead of sgh
	call oro_drag_interface(state,cam_in,sgh30,pbuf,hdtime,dummy_nm,&
                                gwd_ls,gwd_bl,gwd_ss,gwd_fd,&
                                od_ls_ncleff,od_bl_ncd,od_ss_sncleff,&
                                dummy_utgw,dummy_vtgw,dummy_ttgw,& 
                                dtaux3_ls=dummx3_ls,dtauy3_ls=dummy3_ls,&
                                dtaux3_bl=dummx3_bl,dtauy3_bl=dummy3_bl,&
                                dtaux3_ss=dummx3_ss,dtauy3_ss=dummy3_ss,&
                                dtaux3_fd=dtaux3_fd,dtauy3_fd=dtauy3_fd,&
                                dusfc_ls=dummx_ls,dvsfc_ls=dummy_ls,&
                                dusfc_bl=dummx_bl,dvsfc_bl=dummy_bl,&
                                dusfc_ss=dummx_ss,dvsfc_ss=dummy_ss,&
                                dusfc_fd=dusfc_fd,dvsfc_fd=dvsfc_fd)

        call outfld ('DTAUX3_FD', dtaux3_fd,  pcols, lchnk)
        call outfld ('DTAUY3_FD', dtauy3_fd,  pcols, lchnk)
        call outfld ('DUSFC_FD', dusfc_fd,  pcols, lchnk)
        call outfld ('DVSFC_FD', dvsfc_fd,  pcols, lchnk)
   endif

   if (micro_do_icesupersat) then
     call physics_ptend_init(ptend_loc,state%psetcols, 'clubb_ice3', ls=.true., lu=.true., lv=.true., lq=lq)
   endif

   ! ------------------------------------------------- !
   ! End module to compute turbulent mountain stress   !
   ! ------------------------------------------------- !

   call t_startf('adv_clubb_core_col_loop')
   !  Loop over all columns in lchnk to advance CLUBB core
   do i=1,ncol   ! loop over columns

      !  Set time_elapsed to host model time step, this is for
      !  CLUBB's budget stats
      time_elapsed = hdtime_core_rknd

      !  Define the CLUBB momentum grid (in height, units of m)
      do k=1,pverp
         dum1 = state1%zi(i,pverp-k+1)-state1%zi(i,pver+1)
         zi_g(k) = real(dum1, kind = core_rknd)
      enddo

      !  Define the CLUBB thermodynamic grid (in units of m)
      do k=1,pver
         dum1 = state1%zm(i,pver-k+1)-state1%zi(i,pver+1)
         zt_g(k+1) = real(dum1, kind = core_rknd)
         dz_g(k) = state1%zi(i,k)-state1%zi(i,k+1)  ! compute thickness
      enddo

      !  Thermodynamic ghost point is below surface
      zt_g(1) = -1._core_rknd*zt_g(2)

      !  Set the elevation of the surface
      sfc_elevation = real(state1%zi(i,pver+1), kind = core_rknd)

      !  Compute thermodynamic stuff needed for CLUBB on thermo levels.
      !  Inputs for the momentum levels are set below setup_clubb core
      do k=1,pver
         p_in_Pa(k+1)         = real(state1%pmid(i,pver-k+1), kind = core_rknd)             ! Pressure profile
         exner(k+1)           = 1._core_rknd/real(exner_clubb(i,pver-k+1), kind = core_rknd)
         rho(i,k+1)           = invrs_gravit*state1%pdel(i,pver-k+1)/dz_g(pver-k+1)
         rho_ds_zt(k+1)       = real(rho(i,k+1), kind = core_rknd)
         invrs_rho_ds_zt(k+1) = 1._core_rknd/(rho_ds_zt(k+1))                               ! Inverse ds rho at thermo
         thv_ds_zt(k+1)       = real(thv(i,pver-k+1), kind = core_rknd)                     ! thetav on thermo
         rfrzm(k+1)           = real(state1%q(i,pver-k+1,ixcldice), kind = core_rknd)
         radf(k+1)            = real(radf_clubb(i,pver-k+1), kind = core_rknd)
         dum1                 = qrl(i,pver-k+1)/(cpair*state1%pdel(i,pver-k+1))
         qrl_clubb(k+1)       = real(dum1, kind = core_rknd)
      enddo

      !  Below computes the same stuff for the ghost point.  May or may
      !  not be needed, just to be safe to avoid NaN's
      rho_ds_zt(1)       = rho_ds_zt(2)
      invrs_rho_ds_zt(1) = invrs_rho_ds_zt(2)
      rho(i,1)           = rho(i,2)     !rho_ds_zt(2)
      thv_ds_zt(1)       = thv_ds_zt(2)
      rho_zt(:)          = rho_ds_zt(:) !rho(i,:)
      p_in_Pa(1)         = p_in_Pa(2)
      exner(1)           = exner(2)
      rfrzm(1)           = rfrzm(2)
      radf(1)            = radf(2)
      qrl_clubb(1)       = qrl_clubb(2)

      !  Compute mean w wind on thermo grid, convert from omega to w
      wm_zt(1) = 0._core_rknd
      do k=1,pver
        dum1 = -1._r8*state1%omega(i,pver-k+1)*real(invrs_rho_ds_zt(k+1), kind = r8)*invrs_gravit
        wm_zt(k+1) = real(dum1, kind = core_rknd)
      enddo

      !  Set stats output and increment equal to CLUBB and host dt
      stats_tsamp = dtime
      stats_tout  = hdtime_core_rknd

      !  Heights need to be set at each timestep.  Therefore, recall
      !  setup_grid and setup_parameters for this.

      !  Read in parameters for CLUBB.  Pack the default and updated (via nml)
      !  tunable parameters into clubb_params
      call read_parameters_api( -99, "", clubb_params )

      !  Set-up CLUBB core at each CLUBB call because heights can change
      call setup_grid_heights_api(l_implemented, grid_type, zi_g(2), &
         zi_g(1), zi_g, zt_g)

      call setup_parameters_api(zi_g(2), clubb_params, pverp, grid_type, &
        zi_g, zt_g, err_code)

      !  Compute some inputs from the thermodynamic grid
      !  to the momentum grid
      rho_ds_zm       = zt2zm_api(rho_ds_zt)
      rho_zm          = zt2zm_api(rho_zt)
      invrs_rho_ds_zm = zt2zm_api(invrs_rho_ds_zt)
      thv_ds_zm       = zt2zm_api(thv_ds_zt)
      wm_zm           = zt2zm_api(wm_zt)

      !  Surface fluxes provided by host model
      wpthlp_sfc = real(cam_in%shf(i), kind = core_rknd)/(real(cpair, kind = core_rknd)*rho_ds_zm(1)) ! Sensible heat flux
#if 1
      inv_exner_clubb_surf = 1._r8/((state1%pmid(i,pver)/p0_clubb)**(rair/cpair)) !phl Option 2
      wpthlp_sfc = wpthlp_sfc*inv_exner_clubb_surf
#endif
#if 0
      inv_exner_clubb_surf = 1._r8/((state1%pint(i,pverp)/p0_clubb)**(rair/cpair)) !Peter B option
      wpthlp_sfc = wpthlp_sfc*inv_exner_clubb_surf
#endif

      wprtp_sfc  = real(cam_in%cflx(i,1), kind = core_rknd)/rho_ds_zm(1)                              ! Latent heat flux
      upwp_sfc   = real(cam_in%wsx(i), kind = core_rknd)/rho_ds_zm(1)                                 ! Surface meridional momentum flux
      vpwp_sfc   = real(cam_in%wsy(i), kind = core_rknd)/rho_ds_zm(1)                                 ! Surface zonal momentum flux

      ! ------------------------------------------------- !
      ! Apply TMS                                         !
      ! ------------------------------------------------- !
       if ( do_tms ) then
         dum_core_rknd = real((ksrftms(i)*state1%u(i,pver)), kind = core_rknd)
         upwp_sfc      = upwp_sfc-(dum_core_rknd/rho_ds_zm(1))
         dum_core_rknd = real((ksrftms(i)*state1%v(i,pver)), kind = core_rknd)
         vpwp_sfc      = vpwp_sfc-(dum_core_rknd/rho_ds_zm(1))
       endif
      ! ------------------------------------------------- !
      ! Apply TOFD
      ! ------------------------------------------------- !
      ! tendency is flipped already
      if (use_od_fd) then
        um_forcing(2:pverp)=dtaux3_fd(i,pver:1:-1)
        vm_forcing(2:pverp)=dtauy3_fd(i,pver:1:-1)
      endif
      !  Need to flip arrays around for CLUBB core
      do k=1,pverp
         um_in(k)      = real(um(i,pverp-k+1), kind = core_rknd)
         vm_in(k)      = real(vm(i,pverp-k+1), kind = core_rknd)
         upwp_in(k)    = real(upwp(i,pverp-k+1), kind = core_rknd)
         vpwp_in(k)    = real(vpwp(i,pverp-k+1), kind = core_rknd)
         up2_in(k)     = real(up2(i,pverp-k+1), kind = core_rknd)
         vp2_in(k)     = real(vp2(i,pverp-k+1), kind = core_rknd)
         wp2_in(k)     = real(wp2(i,pverp-k+1), kind = core_rknd)
         wp3_in(k)     = real(wp3(i,pverp-k+1), kind = core_rknd)
         rtp2_in(k)    = real(rtp2(i,pverp-k+1), kind = core_rknd)
         thlp2_in(k)   = real(thlp2(i,pverp-k+1), kind = core_rknd)
         thlm_in(k)    = real(thlm(i,pverp-k+1), kind = core_rknd)
         rtm_in(k)     = real(rtm(i,pverp-k+1), kind = core_rknd)
         rvm_in(k)     = real(rvm(i,pverp-k+1), kind = core_rknd)
         wprtp_in(k)   = real(wprtp(i,pverp-k+1), kind = core_rknd)
         wpthlp_in(k)  = real(wpthlp(i,pverp-k+1), kind = core_rknd)
         rtpthlp_in(k) = real(rtpthlp(i,pverp-k+1), kind = core_rknd)

         wpthvp_inout(k)     = real(wpthvp(i,pverp-k+1), kind = core_rknd)
         wp2thvp_inout(k)    = real(wp2thvp(i,pverp-k+1), kind = core_rknd)
         rtpthvp_inout(k)    = real(rtpthvp(i,pverp-k+1), kind = core_rknd)
         thlpthvp_inout(k)   = real(thlpthvp(i,pverp-k+1), kind = core_rknd)
         rcm_inout(k)        = real(rcm(i,pverp-k+1), kind = core_rknd)
         cloud_frac_inout(k) = real(cloud_frac(i,pverp-k+1), kind = core_rknd)

         ! also flip the arrays for the pdf_params_zm variables to have
         ! consistent vertical orientation for all variables in restart file
         ! also need to flip back after calling advance_clubb_core

         pdf_zm_w_1_inout(k) = pdf_zm_w_1(i,pverp-k+1)
         pdf_zm_w_2_inout(k) = pdf_zm_w_2(i,pverp-k+1)
         pdf_zm_varnce_w_1_inout(k) = pdf_zm_varnce_w_1(i,pverp-k+1)
         pdf_zm_varnce_w_2_inout(k) = pdf_zm_varnce_w_2(i,pverp-k+1)
         pdf_zm_mixt_frac_inout(k) =  pdf_zm_mixt_frac(i,pverp-k+1)

         !  Higher order scalar inouts, set to zero
         sclrpthvp_inout(k,:)= 0._core_rknd
         sclrm(k,:)          = 0._core_rknd
         wpsclrp(k,:)        = 0._core_rknd
         sclrp2(k,:)         = 0._core_rknd
         sclrprtp(k,:)       = 0._core_rknd
         sclrpthlp(k,:)      = 0._core_rknd

         if (k .ne. 1) then
            pre_in(k)    = real(prer_evap(i,pverp-k+1), kind = core_rknd)
         endif

      enddo

      if (linearize_pbl_winds) then
         ! Each host model time step, reset the perturbed variables to be equal to
         ! the unperturbed values.
         if (macmic_it == 1) then
            um_pert(i,:) = um_in
            vm_pert(i,:) = vm_in
            upwp_pert(i,:) = upwp_in
            vpwp_pert(i,:) = vpwp_in
         end if
         allocate(um_pert_col(pverp))
         allocate(vm_pert_col(pverp))
         allocate(upwp_pert_col(pverp))
         allocate(vpwp_pert_col(pverp))
         um_pert_col = um_pert(i,:)
         vm_pert_col = vm_pert(i,:)
         upwp_pert_col = upwp_pert(i,:)
         vpwp_pert_col = vpwp_pert(i,:)

         allocate(upwp_sfc_pert)
         allocate(vpwp_sfc_pert)
         ! Prefer to perturb wind/stress in the direction of the existing stress.
         ! However, if there's no existing surface stress, just perturb zonal
         ! wind/stress.
         if (abs(cam_in%wsx(i)) < 1.e-12 .and. abs(cam_in%wsy(i)) < 1.e-12) then
            upwp_sfc_pert = upwp_sfc + pert_tau / rho_ds_zm(1)
            vpwp_sfc_pert = vpwp_sfc
         else
            upwp_sfc_pert = upwp_sfc + cam_in%wsx(i) * &
                 (pert_tau / (rho_ds_zm(1) * hypot(cam_in%wsx(i), cam_in%wsy(i))))
            vpwp_sfc_pert = vpwp_sfc + cam_in%wsy(i) * &
                 (pert_tau / (rho_ds_zm(1) * hypot(cam_in%wsx(i), cam_in%wsy(i))))
         end if
      else
         nullify(upwp_sfc_pert)
         nullify(vpwp_sfc_pert)
         nullify(um_pert_col)
         nullify(vm_pert_col)
         nullify(upwp_pert_col)
         nullify(vpwp_pert_col)
      end if

      pre_in(1) = pre_in(2)

      !  Initialize these to prevent crashing behavior
      edsclr_in(:,:)      = 0._core_rknd
      edsclr_out(:,:)     = 0._r8


      if (clubb_do_adv) then
        if (macmic_it .eq. 1) then
          wp2_in=zt2zm_api(wp2_in)
          wpthlp_in=zt2zm_api(wpthlp_in)
          wprtp_in=zt2zm_api(wprtp_in)
          up2_in=zt2zm_api(up2_in)
          vp2_in=zt2zm_api(vp2_in)
          thlp2_in=zt2zm_api(thlp2_in)
          rtp2_in=zt2zm_api(rtp2_in)
          rtpthlp_in=zt2zm_api(rtpthlp_in)

          do k=1,pverp
            thlp2_in(k)=max(thl_tol**2,thlp2_in(k))
            rtp2_in(k)=max(rt_tol**2,rtp2_in(k))
            wp2_in(k)=max(w_tol_sqd,wp2_in(k))
            up2_in(k)=max(w_tol_sqd,up2_in(k))
            vp2_in(k)=max(w_tol_sqd,vp2_in(k))
          enddo
        endif
      endif

      !  Do the same for tracers
      icnt=0
      do ixind=1,pcnst
         if (lq(ixind))  then
            icnt=icnt+1
            do k=1,pver
               edsclr_in(k+1,icnt) = real(state1%q(i,pver-k+1,ixind), kind = core_rknd)
            enddo
            edsclr_in(1,icnt) = edsclr_in(2,icnt)
         end if
      enddo

      if (do_expldiff) then
        do k=1,pver
          edsclr_in(k+1,icnt+1) = real(thlm(i,pver-k+1), kind = core_rknd)
          edsclr_in(k+1,icnt+2) = real(rtm(i,pver-k+1), kind = core_rknd)
        enddo

        edsclr_in(1,icnt+1) = edsclr_in(2,icnt+1)
        edsclr_in(1,icnt+2) = edsclr_in(2,icnt+2)
      endif

      rho_in(:) = real(rho(i,:), kind = core_rknd)

      ! --------------------------------------------------------- !
      ! Compute cloud-top radiative cooling contribution to CLUBB !
      ! --------------------------------------------------------- !

      ! Sandbox version of code to take into account meso organization

      if (clubb_do_deep) then
         orgparam = 0._r8
         delpavg = 0._r8

         do k = 1, pver
           if (abs(prer_evap(i,k)) .gt. 0._r8) then
             orgparam = orgparam + (abs(prer_evap(i,k)) * 1000._r8 * 1000._r8 * 2._r8 ) * state1%pdel(i,k)
             delpavg = delpavg + state1%pdel(i,k)
           endif
         enddo

         if (delpavg .gt. 0._r8) then
           orgparam = orgparam/delpavg
         endif

         ! Now compute new entrainment rate based on organization
         varmu(i) = mu / (1._r8 + orgparam * 100._r8)
         varmu2   = real(varmu(i), kind = core_rknd)

      endif

      ! --------------------------------------------------------- !
      ! End cloud-top radiative cooling contribution to CLUBB     !
      ! --------------------------------------------------------- !

      pdf_params    => pdf_params_chnk(i,lchnk)
      pdf_params_zm => pdf_params_zm_chnk(i,lchnk)

      if ( is_first_restart_step() .and. ipdf_call_placement .eq. ipdf_post_advance_fields ) then
         ! assign the values read back from restart file
         ! This is necessary when ipdf_call_placement = 2
         pdf_params_zm%w_1 = pdf_zm_w_1_inout
         pdf_params_zm%w_2 = pdf_zm_w_2_inout
         pdf_params_zm%varnce_w_1 = pdf_zm_varnce_w_1_inout
         pdf_params_zm%varnce_w_2 = pdf_zm_varnce_w_2_inout
         pdf_params_zm%mixt_frac = pdf_zm_mixt_frac_inout
      end if

      call t_startf('adv_clubb_core_ts_loop')
      do t=1,nadv    ! do needed number of "sub" timesteps for each CAM step

         !  Increment the statistics then being stats timestep
         if (l_stats) then
            time_elapsed = time_elapsed+dtime
            call stats_begin_timestep_api(time_elapsed, 1, 1)
         endif

         !  Advance CLUBB CORE one timestep in the future
         call t_startf('advance_clubb_core')
         !Balli- to do: check whether initent-ins and intent-inouts are actually what they say

         call advance_clubb_core_api &
              ( l_implemented, dtime, fcor, sfc_elevation, hydromet_dim, & ! intent(in)
              thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &         ! intent(in)
              sclrm_forcing, edsclrm_forcing, wprtp_forcing, &             ! intent(in)
              wpthlp_forcing, rtp2_forcing, thlp2_forcing, &               ! intent(in)
              rtpthlp_forcing, wm_zm, wm_zt, &                             ! intent(in)
              wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &                 ! intent(in)
              wpsclrp_sfc, wpedsclrp_sfc, &                                ! intent(in)
              p_in_Pa, rho_zm, rho_in, exner, &                            ! intent(in)
              rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &                     ! intent(in)
              invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, hydromet, &           ! intent(in)
              rfrzm, radf, &                                               ! intent(in)
#ifdef CLUBBND_CAM
              varmu2, &                                                    ! intent(in)
#endif
              wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt, &                ! intent(in)
              host_dx, host_dy, &                                          ! intent(in)
              um_in, vm_in, upwp_in, &                                     ! intent(inout)
              vpwp_in, up2_in, vp2_in, &                                   ! intent(inout)
              thlm_in, rtm_in, wprtp_in, wpthlp_in, &                      ! intent(inout)
              wp2_in, wp3_in, rtp2_in, &                                   ! intent(inout)
              rtp3_in, thlp2_in, thlp3_in, rtpthlp_in, &                   ! intent(inout)
              sclrm,   &                                                   ! intent(inout)
              sclrp2, sclrprtp, sclrpthlp, &                               ! intent(inout)
              wpsclrp, edsclr_in, err_code, &                              ! intent(inout)
              rcm_inout, cloud_frac_inout, &                               ! intent(inout)
              wpthvp_inout, wp2thvp_inout, rtpthvp_inout, thlpthvp_inout, & ! intent(inout)
              sclrpthvp_inout, &                                            ! intent(inout)
              pdf_params, pdf_params_zm, &                                 ! intent(inout)
              khzm_out, khzt_out, qclvar_out, thlprcp_out, &               ! intent(out)
              wprcp_out, ice_supersat_frac, &                              ! intent(out)
              rcm_in_layer_out, cloud_cover_out, &                         ! intent(out)
              upwp_sfc_pert, vpwp_sfc_pert, &                              ! intent(in)
              um_pert_col, vm_pert_col, upwp_pert_col, vpwp_pert_col)      ! intent(inout)
         call t_stopf('advance_clubb_core')

         if ( err_code == clubb_fatal_error ) then
            write(fstderr,*) "Fatal error in CLUBB: at timestep ", get_nstep(), &
                 "LAT (radians): ", state1%lat(i), &
                 "LON (radians): ", state1%lon(i), &
                 "LAT (degrees): ", rad_to_deg*state1%lat(i), &
                 "LON (degrees): ", rad_to_deg*state1%lon(i), &
                 "Global Column Number: ", get_gcol_p(lchnk,i)
            call endrun('clubb_tend_cam:  Fatal error in CLUBB library'//errmsg(__FILE__,__LINE__))
         end if

         pdf_zm_w_1_inout = pdf_params_zm%w_1
         pdf_zm_w_2_inout = pdf_params_zm%w_2
         pdf_zm_varnce_w_1_inout = pdf_params_zm%varnce_w_1
         pdf_zm_varnce_w_2_inout = pdf_params_zm%varnce_w_2
         pdf_zm_mixt_frac_inout = pdf_params_zm%mixt_frac

         if (do_rainturb) then
            rvm_in = rtm_in - rcm_inout
            call update_xp2_mc_api(pverp, dtime, cloud_frac_inout, &
            rcm_inout, rvm_in, thlm_in, wm_zt, exner, pre_in, pdf_params, &
            rtp2_mc_out, thlp2_mc_out, &
            wprtp_mc_out, wpthlp_mc_out, &
            rtpthlp_mc_out)

            if (clubb_do_deep) then
               dum_core_rknd = 1._core_rknd
            else
               dum_core_rknd = (1._core_rknd - real(cam_in%landfrac(i), kind = core_rknd))
            end if

            ! update turbulent moments based on rain evaporation
            rtp2_in  = rtp2_in + real(clubb_rnevap_effic, kind = core_rknd) * dum_core_rknd * rtp2_mc_out * dtime
            thlp2_in = thlp2_in + real(clubb_rnevap_effic, kind = core_rknd) * dum_core_rknd * thlp2_mc_out * dtime
            if (.not. clubb_do_deep) then
               wprtp_in = wprtp_in + real(clubb_rnevap_effic, kind = core_rknd) * dum_core_rknd * wprtp_mc_out * dtime
               wpthlp_in = wpthlp_in + real(clubb_rnevap_effic, kind = core_rknd) * dum_core_rknd * wpthlp_mc_out * dtime
            endif
!                     rtpthlp_in = rtpthlp_in + rtpthlp_mc_out * dtime

         endif

         if (do_cldcool) then

            rcm_out_zm = zt2zm_api(rcm_inout)
            qrl_zm = zt2zm_api(qrl_clubb)
            thlp2_rad_out(:) = 0._r8
            call calculate_thlp2_rad_api(pverp, rcm_out_zm, thlprcp_out, qrl_zm, thlp2_rad_out)
            thlp2_in = thlp2_in + thlp2_rad_out * dtime
            thlp2_in = max(thl_tol**2,thlp2_in)
          endif

          !  Check to see if stats should be output, here stats are read into
          !  output arrays to make them conformable to CAM output
          if (l_stats) call stats_end_timestep_clubb(lchnk,i,out_zt,out_zm,&
                                                     out_radzt,out_radzm,out_sfc)

      enddo  ! end time loop
      call t_stopf('adv_clubb_core_ts_loop')

      if (clubb_do_adv) then
         if (macmic_it .eq. cld_macmic_num_steps) then
            wp2_in=zm2zt_api(wp2_in)
            wpthlp_in=zm2zt_api(wpthlp_in)
            wprtp_in=zm2zt_api(wprtp_in)
            up2_in=zm2zt_api(up2_in)
            vp2_in=zm2zt_api(vp2_in)
            thlp2_in=zm2zt_api(thlp2_in)
            rtp2_in=zm2zt_api(rtp2_in)
            rtpthlp_in=zm2zt_api(rtpthlp_in)

            do k=1,pverp
               thlp2_in(k)=max(thl_tol**2,thlp2_in(k))
               rtp2_in(k)=max(rt_tol**2,rtp2_in(k))
               wp2_in(k)=max(w_tol_sqd,wp2_in(k))
               up2_in(k)=max(w_tol_sqd,up2_in(k))
               vp2_in(k)=max(w_tol_sqd,vp2_in(k))
            enddo
         endif
      endif

      if (linearize_pbl_winds) then
         ! Copy column variables back to pbuf arrays.
         um_pert(i,:) = um_pert_col
         vm_pert(i,:) = vm_pert_col
         upwp_pert(i,:) = upwp_pert_col
         vpwp_pert(i,:) = vpwp_pert_col
         deallocate(um_pert_col)
         deallocate(vm_pert_col)
         deallocate(upwp_pert_col)
         deallocate(vpwp_pert_col)
         deallocate(upwp_sfc_pert)
         deallocate(vpwp_sfc_pert)
         if (abs(cam_in%wsx(i)) < 1.e-12 .and. abs(cam_in%wsy(i)) < 1.e-12) then
            sfc_v_diff_tau(i) = um_pert(i,2) - um_in(2)
         else
            sfc_v_diff_tau(i) = ((um_pert(i,2) - um_in(2))*cam_in%wsx(i) &
                 + (vm_pert(i,2) - vm_in(2))*cam_in%wsy(i)) &
                 / hypot(cam_in%wsx(i), cam_in%wsy(i))
         end if
      end if

      !  Arrays need to be "flipped" to CAM grid
      do k=1,pverp

          um(i,k)           = real(um_in(pverp-k+1), kind = r8)
          vm(i,k)           = real(vm_in(pverp-k+1), kind = r8)
          upwp(i,k)         = real(upwp_in(pverp-k+1), kind = r8)
          vpwp(i,k)         = real(vpwp_in(pverp-k+1), kind = r8)
          wpthvp(i,k)       = real(wpthvp_inout(pverp-k+1), kind = r8)
          wp2thvp(i,k)      = real(wp2thvp_inout(pverp-k+1), kind = r8)
          rtpthvp(i,k)      = real(rtpthvp_inout(pverp-k+1), kind = r8)
          thlpthvp(i,k)     = real(thlpthvp_inout(pverp-k+1), kind = r8)
          up2(i,k)          = real(up2_in(pverp-k+1), kind = r8)
          vp2(i,k)          = real(vp2_in(pverp-k+1), kind = r8)
          thlm(i,k)         = real(thlm_in(pverp-k+1), kind = r8)
          rtm(i,k)          = real(rtm_in(pverp-k+1), kind = r8)
          wprtp(i,k)        = real(wprtp_in(pverp-k+1), kind = r8)
          wpthlp(i,k)       = real(wpthlp_in(pverp-k+1), kind = r8)
          wp2(i,k)          = real(wp2_in(pverp-k+1), kind = r8)
          wp3(i,k)          = real(wp3_in(pverp-k+1), kind = r8)
          rtp2(i,k)         = real(rtp2_in(pverp-k+1), kind = r8)
          thlp2(i,k)        = real(thlp2_in(pverp-k+1), kind = r8)
          rtpthlp(i,k)      = real(rtpthlp_in(pverp-k+1), kind = r8)
          rcm(i,k)          = real(rcm_inout(pverp-k+1), kind = r8)
          wprcp(i,k)        = real(wprcp_out(pverp-k+1), kind = r8)
          cloud_frac(i,k)   = min(real(cloud_frac_inout(pverp-k+1), kind = r8),1._r8)
          rcm_in_layer(i,k) = real(rcm_in_layer_out(pverp-k+1), kind = r8)
          cloud_cover(i,k)  = min(real(cloud_cover_out(pverp-k+1), kind = r8),1._r8)
          zt_out(i,k)       = real(zt_g(pverp-k+1), kind = r8)
          zi_out(i,k)       = real(zi_g(pverp-k+1), kind = r8)
          khzm(i,k)         = real(khzm_out(pverp-k+1), kind = r8)
          khzt(i,k)         = real(khzt_out(pverp-k+1), kind = r8)
          qclvar(i,k)       = min(1._r8,real(qclvar_out(pverp-k+1), kind = r8))
          sclrpthvp(i,k,:)  = real(sclrpthvp_inout(pverp-k+1,:), kind = r8)
          pdf_zm_w_1(i,k) = pdf_zm_w_1_inout(pverp-k+1)
          pdf_zm_w_2(i,k) = pdf_zm_w_2_inout(pverp-k+1)
          pdf_zm_varnce_w_1(i,k) = pdf_zm_varnce_w_1_inout(pverp-k+1)
          pdf_zm_varnce_w_2(i,k) = pdf_zm_varnce_w_2_inout(pverp-k+1)
          pdf_zm_mixt_frac(i,k) =  pdf_zm_mixt_frac_inout(pverp-k+1)

          do ixind=1,edsclr_dim
              edsclr_out(k,ixind) = real(edsclr_in(pverp-k+1,ixind), kind = r8)
          enddo

      enddo

      !  Fill up arrays needed for McICA.  Note we do not want the ghost point,
      !   thus why the second loop is needed.

      zi_out(i,1) = 0._r8

      ! Compute integrals for static energy, kinetic energy, water vapor, and liquid water
      ! after CLUBB is called.  This is for energy conservation purposes.
      se_a = 0._r8
      ke_a = 0._r8
      wv_a = 0._r8
      wl_a = 0._r8
      do k=1,pver
#ifdef NEWTHETAL
         enthalpy = cpair*thlm(i,k)/exner_clubb(i,k) + latvap*rcm(i,k)
#else
         enthalpy = cpair*((thlm(i,k)+(latvap/cpair)*rcm(i,k))/exner_clubb(i,k))
#endif
         clubb_s(k) = enthalpy + gravit*state1%zm(i,k)+state1%phis(i)
!         se_a(i) = se_a(i) + clubb_s(k)*state1%pdel(i,k)*invrs_gravit
         se_a(i) = se_a(i) + enthalpy * state1%pdel(i,k)*invrs_gravit
         ke_a(i) = ke_a(i) + 0.5_r8*(um(i,k)**2+vm(i,k)**2)*state1%pdel(i,k)*invrs_gravit
         wv_a(i) = wv_a(i) + (rtm(i,k)-rcm(i,k))*state1%pdel(i,k)*invrs_gravit
         wl_a(i) = wl_a(i) + (rcm(i,k))*state1%pdel(i,k)*invrs_gravit
      enddo

      ! Based on these integrals, compute the total energy before and after CLUBB call
      ! TE as in Williamson2015, E= \int_{whole domain} (K+c_p*T) +
      ! \int_{surface} p_s\phi_s (up to water forms), but we ignore surface term
      ! under assumption that CLUBB does not change surface pressure
      do k=1,pver
         te_a(i) = se_a(i) + ke_a(i) + (latvap+latice)*wv_a(i)+latice*wl_a(i)
         te_b(i) = se_b(i) + ke_b(i) + (latvap+latice)*wv_b(i)+latice*wl_b(i)
      enddo

      ! Take into account the surface fluxes of heat and moisture
      te_b(i) = te_b(i)+(cam_in%shf(i)+(cam_in%cflx(i,1))*(latvap+latice))*hdtime

      ! Limit the energy fixer to find highest layer where CLUBB is active
      ! Find first level where wp2 is higher than lowest threshold
      clubbtop = 1
      do while (wp2(i,clubbtop) .eq. w_tol_sqd .and. clubbtop .lt. pver-1)
         clubbtop = clubbtop + 1
      enddo

      ! Compute the disbalance of total energy, over depth where CLUBB is active
      se_dis = (te_a(i) - te_b(i))/(state1%pint(i,pverp)-state1%pint(i,clubbtop))

      ! Apply this fixer throughout the column evenly, but only at layers where
      ! CLUBB is active.
      do k=clubbtop,pver
         clubb_s(k) = clubb_s(k) - se_dis*gravit
      enddo

      !  Now compute the tendencies of CLUBB to CAM, note that pverp is the ghost point
      !  for all variables and therefore is never called in this loop
      do k=1,pver

         ptend_loc%u(i,k)   = (um(i,k)-state1%u(i,k))*invrs_hdtime                  ! east-west wind
         ptend_loc%v(i,k)   = (vm(i,k)-state1%v(i,k))*invrs_hdtime                  ! north-south wind
         ptend_loc%q(i,k,ixq) = (rtm(i,k)-rcm(i,k)-state1%q(i,k,ixq))*invrs_hdtime  ! water vapor
         ptend_loc%q(i,k,ixcldliq) = (rcm(i,k)-state1%q(i,k,ixcldliq))*invrs_hdtime ! Tendency of liquid water
         ptend_loc%s(i,k)   = (clubb_s(k)-state1%s(i,k))*invrs_hdtime               ! Tendency of static energy

         if (clubb_do_adv) then
            if (macmic_it .eq. cld_macmic_num_steps) then

               ! Here add a constant to moments which can be either positive or
               !  negative.  This is to prevent clipping when dynamics tries to
               !  make all constituents positive
               wp3(i,k) = wp3(i,k) + wp3_const
               rtpthlp(i,k) = rtpthlp(i,k) + rtpthlp_const
               wpthlp(i,k) = wpthlp(i,k) + wpthlp_const
               wprtp(i,k) = wprtp(i,k) + wprtp_const

               ptend_loc%q(i,k,ixthlp2)=(thlp2(i,k)-state1%q(i,k,ixthlp2))*invrs_hdtime         ! THLP Variance
               ptend_loc%q(i,k,ixrtp2)=(rtp2(i,k)-state1%q(i,k,ixrtp2))*invrs_hdtime            ! RTP Variance
               ptend_loc%q(i,k,ixrtpthlp)=(rtpthlp(i,k)-state1%q(i,k,ixrtpthlp))*invrs_hdtime   ! RTP THLP covariance
               ptend_loc%q(i,k,ixwpthlp)=(wpthlp(i,k)-state1%q(i,k,ixwpthlp))*invrs_hdtime      ! WPTHLP
               ptend_loc%q(i,k,ixwprtp)=(wprtp(i,k)-state1%q(i,k,ixwprtp))*invrs_hdtime         ! WPRTP
               ptend_loc%q(i,k,ixwp2)=(wp2(i,k)-state1%q(i,k,ixwp2))*invrs_hdtime               ! WP2
               ptend_loc%q(i,k,ixwp3)=(wp3(i,k)-state1%q(i,k,ixwp3))*invrs_hdtime               ! WP3
               ptend_loc%q(i,k,ixup2)=(up2(i,k)-state1%q(i,k,ixup2))*invrs_hdtime               ! UP2
               ptend_loc%q(i,k,ixvp2)=(vp2(i,k)-state1%q(i,k,ixvp2))*invrs_hdtime               ! VP2
            else
               ptend_loc%q(i,k,ixthlp2)=0._r8
               ptend_loc%q(i,k,ixrtp2)=0._r8
               ptend_loc%q(i,k,ixrtpthlp)=0._r8
               ptend_loc%q(i,k,ixwpthlp)=0._r8
               ptend_loc%q(i,k,ixwprtp)=0._r8
               ptend_loc%q(i,k,ixwp2)=0._r8
               ptend_loc%q(i,k,ixwp3)=0._r8
               ptend_loc%q(i,k,ixup2)=0._r8
               ptend_loc%q(i,k,ixvp2)=0._r8
            endif

         endif

         !  Apply tendencies to ice mixing ratio, liquid and ice number, and aerosol constituents.
         !  Loading up this array doesn't mean the tendencies are applied.
         ! edsclr_out is compressed with just the constituents being used, ptend and state are not compressed

         icnt=0
         do ixind=1,pcnst
            if (lq(ixind)) then
               icnt=icnt+1
               if ((ixind /= ixq)       .and. (ixind /= ixcldliq) .and.&
                   (ixind /= ixthlp2)   .and. (ixind /= ixrtp2)   .and.&
                   (ixind /= ixrtpthlp) .and. (ixind /= ixwpthlp) .and.&
                   (ixind /= ixwprtp)   .and. (ixind /= ixwp2)    .and.&
                   (ixind /= ixwp3)     .and. (ixind /= ixup2)    .and. (ixind /= ixvp2) ) then
                       ptend_loc%q(i,k,ixind) = (edsclr_out(k,icnt)-state1%q(i,k,ixind))*invrs_hdtime ! transported constituents
               end if
            end if
         enddo

      enddo

   enddo  ! end column loop
   call t_stopf('adv_clubb_core_col_loop')


   call outfld('fixerCLUBB', te_a(:ncol)-te_b(:ncol), ncol, lchnk )


   ! Add constant to ghost point so that output is not corrupted
   if (clubb_do_adv) then
      if (macmic_it .eq. cld_macmic_num_steps) then
         wp3(:,pverp) = wp3(:,pverp) + wp3_const
         rtpthlp(:,pverp) = rtpthlp(:,pverp) + rtpthlp_const
         wpthlp(:,pverp) = wpthlp(:,pverp) + wpthlp_const
         wprtp(:,pverp) = wprtp(:,pverp) + wprtp_const
      endif
   endif

   cmeliq(:,:) = ptend_loc%q(:,:,ixcldliq)

   ! ------------------------------------------------- !
   ! End column computation of CLUBB, begin to apply   !
   ! and compute output, etc                           !
   ! ------------------------------------------------- !

   !  Output CLUBB tendencies
   call outfld( 'RVMTEND_CLUBB', ptend_loc%q(:,:,ixq), pcols, lchnk)
   call outfld( 'RCMTEND_CLUBB', ptend_loc%q(:,:,ixcldliq), pcols, lchnk)
   call outfld( 'RIMTEND_CLUBB', ptend_loc%q(:,:,ixcldice), pcols, lchnk)
   call outfld( 'TTEND_CLUBB',   ptend_loc%s/cpair,pcols, lchnk)
   call outfld( 'UTEND_CLUBB',   ptend_loc%u,pcols, lchnk)
   call outfld( 'VTEND_CLUBB',   ptend_loc%v,pcols, lchnk)

   if (clubb_do_deep) call outfld( 'MU_CLUBB',      varmu      ,pcols, lchnk)

   call outfld( 'CMELIQ',        cmeliq, pcols, lchnk)

   !  Update physics tendencies
   if (.not. micro_do_icesupersat) then
      call physics_ptend_init(ptend_all, state%psetcols, 'clubb_ice4')
   endif
   call physics_ptend_sum(ptend_loc,ptend_all,ncol)
   call physics_update(state1,ptend_loc,hdtime)

   call cnd_diag_checkpoint(diag, 'CLUBB'//char_macmic_it, state1, pbuf, cam_in, cam_out)

   ! ------------------------------------------------------------ !
   ! ------------------------------------------------------------ !
   ! ------------------------------------------------------------ !
   ! The rest of the code deals with diagnosing variables         !
   ! for microphysics/radiation computation and macrophysics      !
   ! ------------------------------------------------------------ !
   ! ------------------------------------------------------------ !
   ! ------------------------------------------------------------ !
   call t_startf('clubb_tend_cam_diag')

   ! --------------------------------------------------------------------------------- !
   !  COMPUTE THE ICE CLOUD DETRAINMENT                                                !
   !  Detrainment of convective condensate into the environment or stratiform cloud    !
   ! --------------------------------------------------------------------------------- !

   !  Initialize the shallow convective detrainment rate, will always be zero
   dlf2(:,:) = 0.0_r8

   lqice(:)        = .false.
   lqice(ixcldliq) = .true.
   lqice(ixcldice) = .true.
   lqice(ixnumliq) = .true.
   lqice(ixnumice) = .true.

   call physics_ptend_init(ptend_loc,state%psetcols, 'clubb_det', ls=.true., lq=lqice)

   if (zm_microp) then
      call pbuf_get_field(pbuf, dlfzm_idx, dlfzm)
      call pbuf_get_field(pbuf, difzm_idx, difzm)
      call pbuf_get_field(pbuf, dsfzm_idx, dsfzm)
      call pbuf_get_field(pbuf, dnlfzm_idx, dnlfzm)
      call pbuf_get_field(pbuf, dnifzm_idx, dnifzm)
      call pbuf_get_field(pbuf, dnsfzm_idx, dnsfzm)
   end if

   call t_startf('ice_cloud_detrain_diag')
   do k=1,pver
      do i=1,ncol
         if( state1%t(i,k) > clubb_tk1 ) then
            dum1 = 0.0_r8
         elseif ( state1%t(i,k) < clubb_tk2 ) then
            dum1 = 1.0_r8
         else
            !Note: Denominator is changed from 30.0_r8 to (clubb_tk1 - clubb_tk2),
            !(clubb_tk1 - clubb_tk2) is also 30.0 but it introduced a non-bfb change
            dum1 = ( clubb_tk1 - state1%t(i,k) ) /(clubb_tk1 - clubb_tk2)
         endif

         if (zm_microp) then
           ptend_loc%q(i,k,ixcldliq) = dlfzm(i,k) + dlf2(i,k) * ( 1._r8 - dum1 )
           ptend_loc%q(i,k,ixcldice) = difzm(i,k) + dsfzm(i,k) +  dlf2(i,k) * dum1
           ptend_loc%q(i,k,ixnumliq) = dnlfzm(i,k) + 3._r8 * ( dlf2(i,k) * ( 1._r8 - dum1 ) )   &
                                                   / (4._r8*pi*clubb_liq_sh**3._r8*997._r8)      ! Shallow Convection
           ptend_loc%q(i,k,ixnumice) = dnifzm(i,k) + dnsfzm(i,k) + 3._r8 * ( dlf2(i,k) * dum1 ) &
                                                   / (4._r8*pi*clubb_ice_sh**3._r8*500._r8)      ! Shallow Convection
           ptend_loc%s(i,k)          = dlf2(i,k) * dum1 * latice
         else

           ptend_loc%q(i,k,ixcldliq) = dlf(i,k) * ( 1._r8 - dum1 )
           ptend_loc%q(i,k,ixcldice) = dlf(i,k) * dum1
           ptend_loc%q(i,k,ixnumliq) = 3._r8 * ( max(0._r8, ( dlf(i,k) - dlf2(i,k) )) * ( 1._r8 - dum1 ) ) &
                                       / (4._r8*3.14_r8* clubb_liq_deep**3._r8*997._r8) + & ! Deep    Convection
                                       3._r8 * (                         dlf2(i,k)    * ( 1._r8 - dum1 ) ) &
                                       / (4._r8*3.14_r8*clubb_liq_sh**3._r8*997._r8)     ! Shallow Convection
           ptend_loc%q(i,k,ixnumice) = 3._r8 * ( max(0._r8, ( dlf(i,k) - dlf2(i,k) )) *  dum1 ) &
                                       / (4._r8*3.14_r8*clubb_ice_deep**3._r8*500._r8) + & ! Deep    Convection
                                       3._r8 * (                         dlf2(i,k)    *  dum1 ) &
                                       / (4._r8*3.14_r8*clubb_ice_sh**3._r8*500._r8)     ! Shallow Convection
           ptend_loc%s(i,k)          = dlf(i,k) * dum1 * latice

         end if 

         ! Only rliq is saved from deep convection, which is the reserved liquid.  We need to keep
         !   track of the integrals of ice and static energy that is effected from conversion to ice
         !   so that the energy checker doesn't complain.
         det_s(i)                  = det_s(i) + ptend_loc%s(i,k)*state1%pdel(i,k)*invrs_gravit
         det_ice(i)                = det_ice(i) - ptend_loc%q(i,k,ixcldice)*state1%pdel(i,k)*invrs_gravit

      enddo
   enddo

   det_ice(:ncol) = det_ice(:ncol)/1000._r8  ! divide by density of water
   call t_stopf('ice_cloud_detrain_diag')

   call outfld( 'DPDLFLIQ', ptend_loc%q(:,:,ixcldliq), pcols, lchnk)
   call outfld( 'DPDLFICE', ptend_loc%q(:,:,ixcldice), pcols, lchnk)
   call outfld( 'DPDLFT',   ptend_loc%s(:,:)/cpair, pcols, lchnk)

   call physics_ptend_sum(ptend_loc,ptend_all,ncol)
   call physics_update(state1,ptend_loc,hdtime)
   call cnd_diag_checkpoint(diag, 'CUDET'//char_macmic_it, state1, pbuf, cam_in, cam_out)

   ! ptend_all now has all accumulated tendencies.  Convert the tendencies for the
   ! dry constituents to dry air basis.
   do ixind = 1, pcnst
      if (lq(ixind) .and. cnst_type(ixind).eq.'dry') then
         do k = 1, pver
            do i = 1, ncol
               ptend_all%q(i,k,ixind) = ptend_all%q(i,k,ixind)*state1%pdel(i,k)/state1%pdeldry(i,k)
            end do
         end do
      end if
   end do

   ! ------------------------------------------------- !
   ! Diagnose relative cloud water variance            !
   ! ------------------------------------------------- !

   if (deep_scheme .eq. 'CLUBB_SGS') then
      relvarmax = 2.0_r8
   else
      relvarmax = 10.0_r8
   endif

   relvar(:,:) = relvarmax  ! default
!
!PMA c20161114: The lower bound of 0.7 is the mean of scattered Cu in Barker et al (1996).
!     With the new formulation the lower bound and is rarely reached.
!
   relvarmin   = 0.7_r8

!PMA c20161114: Xue Zheng identified the issue with small relvar: the original
!               code uses grid mean variance and water content instead of in-cloud
!               quantities.
!               Following equation A7 in Guo et al (2014), relvar is now  calculated
!               using in-cloud variance and in-cloud total water instead of grid
!               mean. This effectively reduces autoconversion rate especially
!               for thin clouds.
!
!

   relvarc(:ncol,:pver)=fillvalue

   if (deep_scheme .ne. 'CLUBB_SGS') then
      if (relvar_fix) then
         where (rcm(:ncol,:pver) > qsmall .and. qclvar(:ncol,:pver) /= 0._r8)  &
              relvar(:ncol,:pver) = min(relvarmax,max(relvarmin,rcm(:ncol,:pver)**2/max(qsmall,  &
              cloud_frac(:ncol,:pver)*qclvar(:ncol,:pver)-  &
              (1._r8-cloud_frac(:ncol,:pver))*rcm(:ncol,:pver)**2)))
              relvarc(:ncol,:pver) = min(relvarmax,max(relvarmin,rcm(:ncol,:pver)**2/max(qsmall,  &
              cloud_frac(:ncol,:pver)*qclvar(:ncol,:pver)-  &
              (1._r8-cloud_frac(:ncol,:pver))*rcm(:ncol,:pver)**2)))
      else

         where (rcm(:ncol,:pver) /= 0 .and. qclvar(:ncol,:pver) /= 0) &
              relvar(:ncol,:pver) = min(relvarmax,max(0.001_r8,rcm(:ncol,:pver)**2/qclvar(:ncol,:pver)))
      endif
   endif

   ! ------------------------------------------------- !
   ! Optional Accretion enhancement factor             !
   ! ------------------------------------------------- !

     accre_enhan(:ncol,:pver) = micro_mg_accre_enhan_fac !default is 1._r8


   ! ------------------------------------------------- !
   ! Diagnose some output variables                    !
   ! ------------------------------------------------- !

   !  density
   rho(:ncol,1:pver) = state1%pmid(:ncol,1:pver)/(rair*state1%t(:ncol,1:pver))
   rho(:ncol,pverp)  = state1%ps(:ncol)/(rair*state1%t(:ncol,pver))

   eps = rair/rh2o
   wpthvp_diag(:,:) = 0.0_r8
   do k=1,pver
      do i=1,ncol
         !  buoyancy flux
         wpthvp_diag(i,k) = (wpthlp(i,k)-(apply_const*wpthlp_const))+((1._r8-eps)/eps)*real(theta0, kind = r8)* &
                       (wprtp(i,k)-(apply_const*wprtp_const))+((latvap/cpair)* &
                       state1%exner(i,k)-(1._r8/eps)*real(theta0, kind = r8))*wprcp(i,k)

         !  total water mixing ratio
         qt_output(i,k) = state1%q(i,k,ixq)+state1%q(i,k,ixcldliq)+state1%q(i,k,ixcldice)
         !  liquid water potential temperature
         thetal_output(i,k) = (state1%t(i,k)*state1%exner(i,k))-(latvap/cpair)*state1%q(i,k,ixcldliq)
         !  liquid water static energy
         sl_output(i,k) = cpair*state1%t(i,k)+gravit*state1%zm(i,k)-latvap*state1%q(i,k,ixcldliq)
      enddo
   enddo

   do k=1,pverp
      do i=1,ncol
         wpthlp_output(i,k)  = (wpthlp(i,k)-(apply_const*wpthlp_const))*rho(i,k)*cpair !  liquid water potential temperature flux
         wprtp_output(i,k)   = (wprtp(i,k)-(apply_const*wprtp_const))*rho(i,k)*latvap  !  total water mixig ratio flux
         rtpthlp_output(i,k) = rtpthlp(i,k)-(apply_const*rtpthlp_const)                !  rtpthlp output
         wp3_output(i,k)     = wp3(i,k) - (apply_const*wp3_const)                      !  wp3 output
         tke(i,k)            = 0.5_r8*(up2(i,k)+vp2(i,k)+wp2(i,k))                     !  turbulent kinetic energy
      enddo
   enddo

   ! --------------------------------------------------------------------------------- !
   !  Diagnose some quantities that are computed in macrop_tend here.                  !
   !  These are inputs required for the microphysics calculation.                      !
   !                                                                                   !
   !  FIRST PART COMPUTES THE STRATIFORM CLOUD FRACTION FROM CLUBB CLOUD FRACTION      !
   ! --------------------------------------------------------------------------------- !

  ! HW: set alst to alst_o before getting updated
  if(liqcf_fix) then
      if(.not.is_first_step()) alst_o(:ncol,:pver) = alst(:ncol,:pver)
   endif

   !  initialize variables
   alst(:,:) = 0.0_r8
   qlst(:,:) = 0.0_r8

   do k=1,pver
      do i=1,ncol
         alst(i,k) = cloud_frac(i,k)
         qlst(i,k) = rcm(i,k)/max(0.01_r8,alst(i,k))  ! Incloud stratus condensate mixing ratio
      enddo
   enddo

   ! HW
   if(liqcf_fix) then
      if(is_first_step()) alst_o(:ncol,:pver) = alst(:ncol,:pver)
   endif
   !HW

   ! --------------------------------------------------------------------------------- !
   !  THIS PART COMPUTES CONVECTIVE AND DEEP CONVECTIVE CLOUD FRACTION                 !
   ! --------------------------------------------------------------------------------- !

   deepcu(:,pver) = 0.0_r8
   shalcu(:,pver) = 0.0_r8

   do k=1,pver-1
      do i=1,ncol
         !  diagnose the deep convective cloud fraction, as done in macrophysics based on the
         !  deep convective mass flux, read in from pbuf.  Since shallow convection is never
         !  called, the shallow convective mass flux will ALWAYS be zero, ensuring that this cloud
         !  fraction is purely from deep convection scheme.
         deepcu(i,k) = max(0.0_r8,min(dp1*log(1.0_r8+500.0_r8*(cmfmc(i,k+1)-cmfmc_sh(i,k+1))),0.6_r8))
         shalcu(i,k) = 0._r8

         if (deepcu(i,k) <= frac_limit .or. dp_icwmr(i,k) < ic_limit) then
            deepcu(i,k) = 0._r8
         endif

         !  using the deep convective cloud fraction, and CLUBB cloud fraction (variable
         !  "cloud_frac"), compute the convective cloud fraction.  This follows the formulation
         !  found in macrophysics code.  Assumes that convective cloud is all nonstratiform cloud
         !  from CLUBB plus the deep convective cloud fraction
         concld(i,k) = min(cloud_frac(i,k)-alst(i,k)+deepcu(i,k),0.80_r8)
      enddo
   enddo

   ! --------------------------------------------------------------------------------- !
   !  COMPUTE THE ICE CLOUD FRACTION PORTION                                           !
   !  use the aist_vector function to compute the ice cloud fraction                   !
   ! --------------------------------------------------------------------------------- !

   call t_startf('ice_cloud_frac_diag')
   do k=1,pver
      call aist_vector(state1%q(:,k,ixq),state1%t(:,k),state1%pmid(:,k),state1%q(:,k,ixcldice), &
           state1%q(:,k,ixnumice),cam_in%landfrac(:),cam_in%snowhland(:),aist(:,k),ncol)
   enddo
   call t_stopf('ice_cloud_frac_diag')

   ! --------------------------------------------------------------------------------- !
   !  THIS PART COMPUTES THE LIQUID STRATUS FRACTION                                   !
   !                                                                                   !
   !  For now leave the computation of ice stratus fraction from macrop_driver intact  !
   !  because CLUBB does nothing with ice.  Here I simply overwrite the liquid stratus !
   !  fraction that was coded in macrop_driver                                         !
   ! --------------------------------------------------------------------------------- !

   !  Recompute net stratus fraction using maximum over-lapping assumption, as done
   !  in macrophysics code, using alst computed above and aist read in from physics buffer

   cldthresh=1.e-18_r8

   do k=1,pver
      do i=1,ncol

         ast(i,k) = max(alst(i,k),aist(i,k))

         qist(i,k) = state1%q(i,k,ixcldice)/max(0.01_r8,aist(i,k))
      enddo
   enddo

   !  Probably need to add deepcu cloud fraction to the cloud fraction array, else would just
   !  be outputting the shallow convective cloud fraction

   do k=1,pver
      do i=1,ncol
         cloud_frac(i,k) = min(ast(i,k)+deepcu(i,k),1.0_r8)
      enddo
   enddo

   ! --------------------------------------------------------------------------------- !
   !  DIAGNOSE THE PBL DEPTH                                                           !
   !  this is needed for aerosol code                                                  !
   ! --------------------------------------------------------------------------------- !

   call t_startf('pbl_depth_diag')
   do i=1,ncol
      do k=1,pver
         th(i,k) = state1%t(i,k)*state1%exner(i,k)
         if (use_sgv) then
           thv(i,k) = th(i,k)*(1.0_r8+zvir*state1%q(i,k,ixq) &
                    - state1%q(i,k,ixcldliq))  !PMA corrects thv formula
         else
           thv(i,k) = th(i,k)*(1.0_r8+zvir*state1%q(i,k,ixq))
         end if
      enddo
   enddo

   ! diagnose surface friction and obukhov length (inputs to diagnose PBL depth)
   kbfs = 0._r8
   do i=1,ncol
      rrho = invrs_gravit*(state1%pdel(i,pver)/dz_g(pver))
      call calc_ustar( state1%t(i,pver), state1%pmid(i,pver), cam_in%wsx(i), cam_in%wsy(i), &
                       rrho, ustar2(i) )
      call calc_obklen( th(i,pver), thv(i,pver), cam_in%cflx(i,1), cam_in%shf(i), rrho, ustar2(i), &
                        kinheat(i), kinwat(i), kbfs(i), obklen(i) )
   enddo

   dummy2(:) = 0._r8
   dummy3(:) = 0._r8

   !  Compute PBL depth according to Holtslag-Boville Scheme
   call pblintd(ncol, thv, state1%zm, state1%u, state1%v, &
                ustar2, obklen, kbfs, pblh, dummy2, &
                state1%zi, cloud_frac(:,1:pver), 1._r8-cam_in%landfrac, dummy3)
   call t_stopf('pbl_depth_diag')

   !  Output the PBL depth
   call outfld('PBLH', pblh, pcols, lchnk)

   ! Assign the first pver levels of cloud_frac back to cld
   cld(:,1:pver) = cloud_frac(:,1:pver)

!PMA adds gustiness and tpert

    vmag_gust(:)    = 0._r8
    vmag_gust_dp(:) = 0._r8
    vmag_gust_cl(:) = 0._r8
    ktopi(:)        = pver

    if (use_sgv) then
       do i=1,ncol
           up2b(i)          = up2(i,pver)
           vp2b(i)          = vp2(i,pver)
           umb(i)           = state1%u(i,pver)
           vmb(i)           = state1%v(i,pver)
           prec_gust(i)     = max(0._r8,prec_dp(i)-snow_dp(i))*1.e3_r8
           if (cam_in%landfrac(i).gt.0.95_r8) then
             gust_fac(i)   = gust_facl
           else
             gust_fac(i)   = gust_faco
           endif
           vmag_gust_dp(i) = ugust(min(prec_gust(i),6.94444e-4_r8)) ! Limit for the ZM gustiness equation set in Redelsperger et al. (2000)
           vmag_gust_dp(i) = max(0._r8, vmag_gust_dp(i) )
           vmag_gust_cl(i) = sqrt(max(0._r8,up2b(i)+vp2b(i)))
           vmag_gust(i)    = sqrt(gust_facc * vmag_gust_cl(i)**2 &
                + gust_fac(i) * vmag_gust_dp(i)**2)
          do k=1,pver
             if (state1%zi(i,k)>pblh(i).and.state1%zi(i,k+1)<=pblh(i)) then
                ktopi(i) = k
                exit
             end if
          end do
          tpert(i) = min(2._r8,(sqrt(thlp2(i,ktopi(i)))+(latvap/cpair)*state1%q(i,ktopi(i),ixcldliq)) &
                    /max(state1%exner(i,ktopi(i)),1.e-3_r8)) !proxy for tpert
       end do
    end if

   if (linearize_pbl_winds .and. macmic_it == cld_macmic_num_steps) then
      do i = 1, ncol
         wsresp(i) = sfc_v_diff_tau(i) / pert_tau
         ! Estimated tau in balance with wind is the tau we just used.
         if (cam_in%wsx(i) == 0._r8 .or. cam_in%wsy(i) == 0._r8) then
            ! Work around an odd FPE issue with intel compiler.
            tau_est(i) = abs(cam_in%wsx(i)) + abs(cam_in%wsy(i))
         else
            tau_est(i) = hypot(cam_in%wsx(i), cam_in%wsy(i))
         end if
      end do
   end if

   call outfld('VMAGGUST', vmag_gust, pcols, lchnk)
   call outfld('VMAGDP', vmag_gust_dp, pcols, lchnk)
   call outfld('VMAGCL', vmag_gust_cl, pcols, lchnk)
   call outfld('TPERTBLT', tpert, pcols, lchnk)

   call t_stopf('clubb_tend_cam_diag')

   ! --------------------------------------------------------------------------------- !
   !  END CLOUD FRACTION DIAGNOSIS, begin to store variables back into buffer          !
   ! --------------------------------------------------------------------------------- !

   !  Output calls of variables goes here
   call outfld( 'RELVAR',           relvar,                  pcols, lchnk )
   call outfld( 'RELVARC',          relvarc,                 pcols, lchnk )
   call outfld( 'RHO_CLUBB',        rho,                     pcols, lchnk )
   call outfld( 'WP2_CLUBB',        wp2,                     pcols, lchnk )
   call outfld( 'UP2_CLUBB',        up2,                     pcols, lchnk )
   call outfld( 'VP2_CLUBB',        vp2,                     pcols, lchnk )
   call outfld( 'WP3_CLUBB',        wp3_output,              pcols, lchnk )
   call outfld( 'UPWP_CLUBB',       upwp,                    pcols, lchnk )
   call outfld( 'VPWP_CLUBB',       vpwp,                    pcols, lchnk )
   call outfld( 'WPTHLP_CLUBB',     wpthlp_output,           pcols, lchnk )
   call outfld( 'WPRTP_CLUBB',      wprtp_output,            pcols, lchnk )
   tmp_array = rtp2(:ncol,:)*1000._r8
   call outfld( 'RTP2_CLUBB',       tmp_array,               ncol,  lchnk )
   call outfld( 'THLP2_CLUBB',      thlp2,                   pcols, lchnk )
   tmp_array = rtpthlp_output(:ncol,:)*1000._r8
   call outfld( 'RTPTHLP_CLUBB',    tmp_array,               ncol,  lchnk )
   tmp_array = rcm(:ncol,:)*1000._r8
   call outfld( 'RCM_CLUBB',        tmp_array,               ncol,  lchnk )
   tmp_array = wprcp(:ncol,:)*latvap
   call outfld( 'WPRCP_CLUBB',      tmp_array,               ncol,  lchnk )
   call outfld( 'CLOUDFRAC_CLUBB',  alst,                    pcols, lchnk )
   tmp_array = rcm_in_layer(:ncol,:)*1000._r8
   call outfld( 'RCMINLAYER_CLUBB', tmp_array,               ncol,  lchnk )
   call outfld( 'CLOUDCOVER_CLUBB', cloud_frac,              pcols, lchnk )
   tmp_array = wpthvp_diag(:ncol,:)*cpair
   call outfld( 'WPTHVP_CLUBB',     tmp_array,               ncol,  lchnk )
   tmp_array = 1._r8*zt_out(:ncol,:)
   call outfld( 'ZT_CLUBB',         tmp_array,               ncol,  lchnk )
   tmp_array = 1._r8*zi_out(:ncol,:)
   call outfld( 'ZM_CLUBB',         tmp_array,               ncol,  lchnk )
   call outfld( 'UM_CLUBB',         um,                      pcols, lchnk )
   call outfld( 'VM_CLUBB',         vm,                      pcols, lchnk )
   call outfld( 'THETAL',           thetal_output,           pcols, lchnk )
   call outfld( 'QT',               qt_output,               pcols, lchnk )
   call outfld( 'SL',               sl_output,               pcols, lchnk )
   call outfld( 'CONCLD',           concld,                  pcols, lchnk )

   !  Output CLUBB history here
   if (l_stats) then

      do i=1,stats_zt%num_output_fields

         temp1 = trim(stats_zt%file%var(i)%name)
         sub   = temp1
         if (len(temp1) .gt. 16) sub = temp1(1:16)

         call outfld(trim(sub), out_zt(:,:,i), pcols, lchnk )
      enddo

      do i=1,stats_zm%num_output_fields

         temp1 = trim(stats_zm%file%var(i)%name)
         sub   = temp1
         if (len(temp1) .gt. 16) sub = temp1(1:16)

         call outfld(trim(sub),out_zm(:,:,i), pcols, lchnk)
      enddo

      if (l_output_rad_files) then
         do i=1,stats_rad_zt%num_output_fields
            call outfld(trim(stats_rad_zt%file%var(i)%name), out_radzt(:,:,i), pcols, lchnk)
         enddo

         do i=1,stats_rad_zm%num_output_fields
            call outfld(trim(stats_rad_zm%file%var(i)%name), out_radzm(:,:,i), pcols, lchnk)
         enddo
      endif

      do i=1,stats_sfc%num_output_fields
         call outfld(trim(stats_sfc%file%var(i)%name), out_sfc(:,:,i), pcols, lchnk)
      enddo

   endif

   call cnd_diag_checkpoint(diag, 'MACDIAG'//char_macmic_it, state1, pbuf, cam_in, cam_out)

   return
#endif
  end subroutine clubb_tend_cam

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

    subroutine clubb_surface (state, cam_in, pbuf, ustar, obklen)

!-------------------------------------------------------------------------------
! Description: Provide the obukhov length and the surface friction velocity
!              for the dry deposition code in routine tphysac.  Since University
!              of Washington Moist Turbulence (UWMT) scheme is not called when
!              CLUBB is turned on the obukov length and ustar are never initialized
!              nor computed (sometimes never updated from NaN).  In addition, surface
!              fluxes are applied to the constituents.
!
! Author: Peter Bogenschutz, August 2011
! Origin: Based heavily on UWMT code (eddy_diff.F90)
! References:
!   None
!-------------------------------------------------------------------------------

    use physics_types,          only: physics_state
    use physconst,              only: zvir,gravit
    use ppgrid,                 only: pver, pcols
    use constituents,           only: cnst_get_ind
    use camsrfexch,             only: cam_in_t
    use hb_diff,                only: pblintd_ri
    use physics_buffer,         only: pbuf_get_index, pbuf_get_field, physics_buffer_desc

    implicit none

    ! --------------- !
    ! Input Auguments !
    ! --------------- !

    type(physics_state),                intent(inout)  :: state ! Physics state variables
    type(cam_in_t),                     intent(in)     :: cam_in
    type(physics_buffer_desc), pointer, intent(in)     :: pbuf(:)

    ! ---------------- !
    ! Output Auguments !
    ! ---------------- !

    real(r8),            intent(out)    :: obklen(pcols)        ! Obukhov length [ m ]
    real(r8),            intent(out)    :: ustar(pcols)         ! Surface friction velocity [ m/s ]

#if defined(CLUBB_SGS) || defined(SHOC_SGS)

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    integer :: i                                                ! indicees
    integer :: k
    integer :: ncol                                             ! # of atmospheric columns

    real(r8) :: th(pcols)                                       ! surface potential temperature
    real(r8) :: thv(pcols)                                      ! surface virtual potential temperature
    real(r8) :: th_lv(pcols,pver)                               ! level potential temperature
    real(r8) :: thv_lv(pcols,pver)                              ! level virtual potential temperature
    real(r8) :: kinheat                                         ! kinematic surface heat flux
    real(r8) :: kinwat                                          ! kinematic surface vapor flux
    real(r8) :: kbfs                                            ! kinematic surface buoyancy flux
    real(r8) :: kbfs_pcol(pcols)                                ! kinematic surface buoyancy flux stored for all pcols
    integer  :: ixq,ixcldliq !PMA fix for thv
    real(r8) :: rrho                                            ! Inverse air density

    integer  :: oro_drag_ribulk_idx                             ! pbuf index of bulk richardson number for oro drag
    real(r8), pointer :: oro_drag_ribulk(:)                     ! pbuf pointer for bulk richardson number


#endif
    obklen(pcols) = 0.0_r8
    ustar(pcols)  = 0.0_r8
#if defined(CLUBB_SGS) || defined(SHOC_SGS)

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !
    call cnst_get_ind('Q',ixq)
    if (use_sgv) then
       call cnst_get_ind('CLDLIQ',ixcldliq)
    endif

    ncol = state%ncol

    ! Compute the surface friction velocity and obukov length

    do i = 1, ncol
       th(i) = state%t(i,pver)*state%exner(i,pver)         ! diagnose potential temperature
       if (use_sgv) then
         thv(i) = th(i)*(1._r8+zvir*state%q(i,pver,ixq) & ! PMA corrects virtual potential temperature formula
                       - state%q(i,pver,ixcldliq))
       else
         thv(i) = th(i)*(1._r8+zvir*state%q(i,pver,ixq))  ! diagnose virtual potential temperature
       end if
    enddo

    do i = 1, ncol
       call calc_ustar( state%t(i,pver), state%pmid(i,pver), cam_in%wsx(i), cam_in%wsy(i), &
                        rrho, ustar(i) )
       call calc_obklen( th(i), thv(i), cam_in%cflx(i,1), cam_in%shf(i), rrho, ustar(i), &
                        kinheat, kinwat, kbfs, obklen(i) )
    enddo

    if (use_od_ss) then
      !add calculation of bulk richardson number here
      !compute the whole level th and thv for diagnose of bulk richardson number
      thv_lv=0.0_r8
      th_lv =0.0_r8

      !use the same virtual potential temperature formula as above (thv) except for all vertical levels
      !used for bulk richardson number below in pblintd_ri
      do i=1,ncol
        do k=1,pver
          th_lv(i,k) = state%t(i,k)*state%exner(i,k)
            if (use_sgv) then
              thv_lv(i,k) = th_lv(i,k)*(1.0_r8+zvir*state%q(i,k,ixq) &
                            - state%q(i,k,ixcldliq))  
            else
              thv_lv(i,k) = th_lv(i,k)*(1.0_r8+zvir*state%q(i,k,ixq))
            end if
        enddo
      enddo

      !recalculate the kbfs stored in kbfs_pcol for bulk richardson number in pblintd_ri
      kbfs_pcol=0.0_r8
      do i=1,ncol
        call calc_ustar( state%t(i,pver), state%pmid(i,pver), cam_in%wsx(i), cam_in%wsy(i), rrho, ustar(i) )
        call calc_obklen( th(i), thv(i), cam_in%cflx(i,1), cam_in%shf(i), rrho, ustar(i), &
                        kinheat, kinwat, kbfs, obklen(i) )
        kbfs_pcol(i)=kbfs
      enddo

      oro_drag_ribulk_idx = pbuf_get_index('oro_drag_ribulk')
      call pbuf_get_field(pbuf, oro_drag_ribulk_idx, oro_drag_ribulk)

      !calculate the bulk richardson number
      call pblintd_ri(ncol, gravit, thv_lv, state%zm, state%u, state%v, &
                      ustar, obklen, kbfs_pcol, oro_drag_ribulk)
    endif

    return

#endif

    end subroutine clubb_surface

#ifdef CLUBB_SGS
! ----------------------------------------------------------------------
!
! DISCLAIMER : this code appears to be correct but has not been
!              very thouroughly tested. If you do notice any
!              anomalous behaviour then please contact Andy and/or
!              Bjorn
!
! Function diag_ustar:  returns value of ustar using the below
! similarity functions and a specified buoyancy flux (bflx) given in
! kinematic units
!
! phi_m (zeta > 0) =  (1 + am * zeta)
! phi_m (zeta < 0) =  (1 - bm * zeta)^(-1/4)
!
! where zeta = z/lmo and lmo = (theta_rev/g*vonk) * (ustar^2/tstar)
!
! Ref: Businger, 1973, Turbulent Transfer in the Atmospheric Surface
! Layer, in Workshop on Micormeteorology, pages 67-100.
!
! Code writen March, 1999 by Bjorn Stevens
!

real(r8) function diag_ustar( z, bflx, wnd, z0 )

use shr_const_mod, only : shr_const_karman, shr_const_pi, shr_const_g

implicit none

real(r8), parameter      :: am   =  4.8_r8   !   "          "         "
real(r8), parameter      :: bm   = 19.3_r8  !   "          "         "

real(r8), parameter      :: grav = shr_const_g
real(r8), parameter      :: vonk = shr_const_karman
real(r8), parameter      :: pi   = shr_const_pi

real(r8), intent (in)    :: z             ! height where u locates
real(r8), intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
real(r8), intent (in)    :: wnd           ! wind speed at z
real(r8), intent (in)    :: z0            ! momentum roughness height


integer :: iterate
real(r8)    :: lnz, klnz, c1, x, psi1, zeta, lmo, ustar

lnz   = log( z / z0 )
klnz  = vonk/lnz
c1    = pi / 2.0_r8 - 3.0_r8*log( 2.0_r8 )

ustar =  wnd*klnz
if (abs(bflx) > 1.e-6_r8) then
   do iterate=1,4

      if (ustar > 1.e-6_r8) then
         lmo   = -ustar**3 / ( vonk * bflx )
         zeta  = z/lmo
         if (zeta > 0._r8) then
            ustar =  vonk*wnd  /(lnz + am*zeta)
         else
            x     = sqrt( sqrt( 1.0_r8 - bm*zeta ) )
            psi1  = 2._r8*log( 1.0_r8+x ) + log( 1.0_r8+x*x ) - 2._r8*atan( x ) + c1
            ustar = wnd*vonk/(lnz - psi1)
         end if

      endif

   end do
end if


diag_ustar = ustar

return


end function diag_ustar
#endif

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

#ifdef CLUBB_SGS

  subroutine stats_init_clubb( l_stats_in, stats_tsamp_in, stats_tout_in, &
                         nnzp, nnrad_zt,nnrad_zm, delt )
    !
    ! Description: Initializes the statistics saving functionality of
    !   the CLUBB model.  This is for purpose of CAM-CLUBB interface.  Here
    !   the traditional stats_init of CLUBB is not called, as it is not compatible
    !   with CAM output.

    !-----------------------------------------------------------------------


    use stats_variables, only: &
      stats_zt,      & ! Variables
      ztscr01, &
      ztscr02, &
      ztscr03, &
      ztscr04, &
      ztscr05, &
      ztscr06, &
      ztscr07, &
      ztscr08, &
      ztscr09, &
      ztscr10, &
      ztscr11, &
      ztscr12, &
      ztscr13, &
      ztscr14, &
      ztscr15, &
      ztscr16, &
      ztscr17, &
      ztscr18, &
      ztscr19, &
      ztscr20, &
      ztscr21

    use stats_variables, only: &
      stats_zm,      &
      zmscr01, &
      zmscr02, &
      zmscr03, &
      zmscr04, &
      zmscr05, &
      zmscr06, &
      zmscr07, &
      zmscr08, &
      zmscr09, &
      zmscr10, &
      zmscr11, &
      zmscr12, &
      zmscr13, &
      zmscr14, &
      zmscr15, &
      zmscr16, &
      zmscr17, &
      stats_rad_zt,  &
      stats_rad_zm,  &
      stats_sfc,     &
      l_stats, &
      l_output_rad_files, &
      stats_tsamp,   &
      stats_tout,    &
      l_stats_samp,  &
      l_stats_last, &
      fname_rad_zt, &
      fname_rad_zm, &
      fname_sfc, &
      l_netcdf, &
      l_grads

    use clubb_precision,        only: time_precision, core_rknd   !
    use stats_zm_module,        only: nvarmax_zm, stats_init_zm !
    use stats_zt_module,        only: nvarmax_zt, stats_init_zt !
    use stats_rad_zt_module,    only: nvarmax_rad_zt, stats_init_rad_zt !
    use stats_rad_zm_module,    only: nvarmax_rad_zm, stats_init_rad_zm !
    use stats_sfc_module,       only: nvarmax_sfc, stats_init_sfc !
    use error_code,             only: clubb_at_least_debug_level !
    use constants_clubb,        only: fstderr, var_length !
    use cam_history,            only: addfld, horiz_only
    use namelist_utils,         only: find_group_name
    use units,                  only: getunit, freeunit
    use cam_abortutils,         only: endrun

    implicit none

    ! Input Variables

    logical, intent(in) :: l_stats_in ! Stats on? T/F

    real(kind=time_precision), intent(in) ::  &
      stats_tsamp_in,  & ! Sampling interval   [s]
      stats_tout_in      ! Output interval     [s]

    integer, intent(in) :: nnzp     ! Grid points in the vertical [count]
    integer, intent(in) :: nnrad_zt ! Grid points in the radiation grid [count]
    integer, intent(in) :: nnrad_zm ! Grid points in the radiation grid [count]

    real(kind=time_precision), intent(in) ::   delt         ! Timestep (dtmain in CLUBB)         [s]


    !  Local Variables

    !  Namelist Variables

    character(len=var_length), dimension(nvarmax_zt)     ::   clubb_vars_zt      ! Variables on the thermodynamic levels
    character(len=var_length), dimension(nvarmax_zm)     ::   clubb_vars_zm      ! Variables on the momentum levels
    character(len=var_length), dimension(nvarmax_rad_zt) ::   clubb_vars_rad_zt  ! Variables on the radiation levels
    character(len=var_length), dimension(nvarmax_rad_zm) ::   clubb_vars_rad_zm  ! Variables on the radiation levels
    character(len=var_length), dimension(nvarmax_sfc)    ::   clubb_vars_sfc     ! Variables at the model surface

    namelist /clubb_stats_nl/ &
      clubb_vars_zt, &
      clubb_vars_zm, &
      clubb_vars_rad_zt, &
      clubb_vars_rad_zm, &
      clubb_vars_sfc

    !  Local Variables

    logical :: l_error

    character(len=200) :: fname, temp1, sub

    integer :: i, ntot, read_status
    integer :: iunit

    !  Initialize
    l_error = .false.

    !  Set stats_variables variables with inputs from calling subroutine
    l_stats = l_stats_in

    stats_tsamp = stats_tsamp_in
    stats_tout  = stats_tout_in

    if ( .not. l_stats ) then
       l_stats_samp  = .false.
       l_stats_last  = .false.
       return
    end if

    !  Initialize namelist variables

    clubb_vars_zt     = ''
    clubb_vars_zm     = ''
    clubb_vars_rad_zt = ''
    clubb_vars_rad_zm = ''
    clubb_vars_sfc    = ''

    !  Read variables to compute from the namelist
    if (masterproc) then
       iunit= getunit()
       open(unit=iunit,file="atm_in",status='old')
       call find_group_name(iunit, 'clubb_stats_nl', status=read_status)
       if (read_status == 0) then
          read(unit=iunit, nml=clubb_stats_nl, iostat=read_status)
          if (read_status /= 0) then
             call endrun('stats_init_clubb:  error reading namelist'//errmsg(__FILE__,__LINE__))
          end if
       end if
       close(unit=iunit)
       call freeunit(iunit)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(clubb_vars_zt,      var_length*nvarmax_zt,       mpichar,   0, mpicom)
    call mpibcast(clubb_vars_zm,      var_length*nvarmax_zm,       mpichar,   0, mpicom)
    call mpibcast(clubb_vars_rad_zt,  var_length*nvarmax_rad_zt,   mpichar,   0, mpicom)
    call mpibcast(clubb_vars_rad_zm,  var_length*nvarmax_rad_zm,   mpichar,   0, mpicom)
    call mpibcast(clubb_vars_sfc,     var_length*nvarmax_sfc,      mpichar,   0, mpicom)
#endif

    !  Hardcode these for use in CAM-CLUBB, don't want either
    l_netcdf = .false.
    l_grads  = .false.

    !  Check sampling and output frequencies

    !  The model time step length, delt (which is dtmain), should multiply
    !  evenly into the statistical sampling time step length, stats_tsamp.
    if ( abs( stats_tsamp/delt - floor(stats_tsamp/delt) ) > 1.e-8_r8 ) then
       l_error = .true.  ! This will cause the run to stop.
       write(fstderr,*) 'Error:  stats_tsamp should be an even multiple of ',  &
                        'delt (which is dtmain).  Check the appropriate ',  &
                        'model.in file.'
       write(fstderr,*) 'stats_tsamp = ', stats_tsamp
       write(fstderr,*) 'delt = ', delt
    endif

    !  Initialize zt (mass points)

    i = 1
    do while ( ichar(clubb_vars_zt(i)(1:1)) /= 0 .and. &
               len_trim(clubb_vars_zt(i))   /= 0 .and. &
               i <= nvarmax_zt )
       i = i + 1
    enddo
    ntot = i - 1
    if ( ntot == nvarmax_zt ) then
       write(fstderr,*) "There are more statistical variables listed in ",  &
                        "clubb_vars_zt than allowed for by nvarmax_zt."
       write(fstderr,*) "Check the number of variables listed for clubb_vars_zt ",  &
                        "in the stats namelist, or change nvarmax_zt."
       write(fstderr,*) "nvarmax_zt = ", nvarmax_zt
       call endrun ("stats_init_clubb:  number of zt statistical variables exceeds limit"//errmsg(__FILE__,__LINE__))
    endif

    stats_zt%num_output_fields = ntot
    stats_zt%kk = nnzp

    allocate( stats_zt%z( stats_zt%kk ) )

    allocate( stats_zt%accum_field_values( 1, 1, stats_zt%kk, stats_zt%num_output_fields ) )
    allocate( stats_zt%accum_num_samples( 1, 1, stats_zt%kk, stats_zt%num_output_fields ) )
    allocate( stats_zt%l_in_update( 1, 1, stats_zt%kk, stats_zt%num_output_fields ) )
    call stats_zero( stats_zt%kk, stats_zt%num_output_fields, stats_zt%accum_field_values, &
                     stats_zt%accum_num_samples, stats_zt%l_in_update )

    allocate( stats_zt%file%var( stats_zt%num_output_fields ) )
    allocate( stats_zt%file%z( stats_zt%kk ) )

    !  Allocate scratch space

    allocate( ztscr01(stats_zt%kk) )
    allocate( ztscr02(stats_zt%kk) )
    allocate( ztscr03(stats_zt%kk) )
    allocate( ztscr04(stats_zt%kk) )
    allocate( ztscr05(stats_zt%kk) )
    allocate( ztscr06(stats_zt%kk) )
    allocate( ztscr07(stats_zt%kk) )
    allocate( ztscr08(stats_zt%kk) )
    allocate( ztscr09(stats_zt%kk) )
    allocate( ztscr10(stats_zt%kk) )
    allocate( ztscr11(stats_zt%kk) )
    allocate( ztscr12(stats_zt%kk) )
    allocate( ztscr13(stats_zt%kk) )
    allocate( ztscr14(stats_zt%kk) )
    allocate( ztscr15(stats_zt%kk) )
    allocate( ztscr16(stats_zt%kk) )
    allocate( ztscr17(stats_zt%kk) )
    allocate( ztscr18(stats_zt%kk) )
    allocate( ztscr19(stats_zt%kk) )
    allocate( ztscr20(stats_zt%kk) )
    allocate( ztscr21(stats_zt%kk) )

    ztscr01 = 0.0_core_rknd
    ztscr02 = 0.0_core_rknd
    ztscr03 = 0.0_core_rknd
    ztscr04 = 0.0_core_rknd
    ztscr05 = 0.0_core_rknd
    ztscr06 = 0.0_core_rknd
    ztscr07 = 0.0_core_rknd
    ztscr08 = 0.0_core_rknd
    ztscr09 = 0.0_core_rknd
    ztscr10 = 0.0_core_rknd
    ztscr11 = 0.0_core_rknd
    ztscr12 = 0.0_core_rknd
    ztscr13 = 0.0_core_rknd
    ztscr14 = 0.0_core_rknd
    ztscr15 = 0.0_core_rknd
    ztscr16 = 0.0_core_rknd
    ztscr17 = 0.0_core_rknd
    ztscr18 = 0.0_core_rknd
    ztscr19 = 0.0_core_rknd
    ztscr20 = 0.0_core_rknd
    ztscr21 = 0.0_core_rknd

    !  Default initialization for array indices for zt

    call stats_init_zt( clubb_vars_zt, l_error )

    !  Initialize zm (momentum points)

    i = 1
    do while ( ichar(clubb_vars_zm(i)(1:1)) /= 0  .and. &
               len_trim(clubb_vars_zm(i)) /= 0    .and. &
               i <= nvarmax_zm )
       i = i + 1
    end do
    ntot = i - 1
    if ( ntot == nvarmax_zm ) then
       write(fstderr,*) "There are more statistical variables listed in ",  &
                        "clubb_vars_zm than allowed for by nvarmax_zm."
       write(fstderr,*) "Check the number of variables listed for clubb_vars_zm ",  &
                        "in the stats namelist, or change nvarmax_zm."
       write(fstderr,*) "nvarmax_zm = ", nvarmax_zm
       call endrun ("stats_init_clubb:  number of zm statistical variables exceeds limit"//errmsg(__FILE__,__LINE__))
    endif

    stats_zm%num_output_fields = ntot
    stats_zm%kk = nnzp

    allocate( stats_zm%z( stats_zm%kk ) )

    allocate( stats_zm%accum_field_values( 1, 1, stats_zm%kk, stats_zm%num_output_fields ) )
    allocate( stats_zm%accum_num_samples( 1, 1, stats_zm%kk, stats_zm%num_output_fields ) )
    allocate( stats_zm%l_in_update( 1, 1, stats_zm%kk, stats_zm%num_output_fields ) )
    call stats_zero( stats_zm%kk, stats_zm%num_output_fields, stats_zm%accum_field_values, &
                     stats_zm%accum_num_samples, stats_zm%l_in_update )

    allocate( stats_zm%file%var( stats_zm%num_output_fields ) )
    allocate( stats_zm%file%z( stats_zm%kk ) )

    !  Allocate scratch space

    allocate( zmscr01(stats_zm%kk) )
    allocate( zmscr02(stats_zm%kk) )
    allocate( zmscr03(stats_zm%kk) )
    allocate( zmscr04(stats_zm%kk) )
    allocate( zmscr05(stats_zm%kk) )
    allocate( zmscr06(stats_zm%kk) )
    allocate( zmscr07(stats_zm%kk) )
    allocate( zmscr08(stats_zm%kk) )
    allocate( zmscr09(stats_zm%kk) )
    allocate( zmscr10(stats_zm%kk) )
    allocate( zmscr11(stats_zm%kk) )
    allocate( zmscr12(stats_zm%kk) )
    allocate( zmscr13(stats_zm%kk) )
    allocate( zmscr14(stats_zm%kk) )
    allocate( zmscr15(stats_zm%kk) )
    allocate( zmscr16(stats_zm%kk) )
    allocate( zmscr17(stats_zm%kk) )

    zmscr01 = 0.0_core_rknd
    zmscr02 = 0.0_core_rknd
    zmscr03 = 0.0_core_rknd
    zmscr04 = 0.0_core_rknd
    zmscr05 = 0.0_core_rknd
    zmscr06 = 0.0_core_rknd
    zmscr07 = 0.0_core_rknd
    zmscr08 = 0.0_core_rknd
    zmscr09 = 0.0_core_rknd
    zmscr10 = 0.0_core_rknd
    zmscr11 = 0.0_core_rknd
    zmscr12 = 0.0_core_rknd
    zmscr13 = 0.0_core_rknd
    zmscr14 = 0.0_core_rknd
    zmscr15 = 0.0_core_rknd
    zmscr16 = 0.0_core_rknd
    zmscr17 = 0.0_core_rknd

    call stats_init_zm( clubb_vars_zm, l_error )

    !  Initialize rad_zt (radiation points)

    if (l_output_rad_files) then

       i = 1
       do while ( ichar(clubb_vars_rad_zt(i)(1:1)) /= 0  .and. &
                  len_trim(clubb_vars_rad_zt(i))   /= 0  .and. &
                  i <= nvarmax_rad_zt )
          i = i + 1
       end do
       ntot = i - 1
       if ( ntot == nvarmax_rad_zt ) then
          write(fstderr,*) "There are more statistical variables listed in ",  &
                           "clubb_vars_rad_zt than allowed for by nvarmax_rad_zt."
          write(fstderr,*) "Check the number of variables listed for clubb_vars_rad_zt ",  &
                           "in the stats namelist, or change nvarmax_rad_zt."
          write(fstderr,*) "nvarmax_rad_zt = ", nvarmax_rad_zt
          call endrun ("stats_init_clubb:  number of rad_zt statistical variables exceeds limit"//errmsg(__FILE__,__LINE__))
       endif

      stats_rad_zt%num_output_fields = ntot
      stats_rad_zt%kk = nnrad_zt

      allocate( stats_rad_zt%z( stats_rad_zt%kk ) )

      allocate( stats_rad_zt%accum_field_values( 1, 1, stats_rad_zt%kk, stats_rad_zt%num_output_fields ) )
      allocate( stats_rad_zt%accum_num_samples( 1, 1, stats_rad_zt%kk, stats_rad_zt%num_output_fields ) )
      allocate( stats_rad_zt%l_in_update( 1, 1, stats_rad_zt%kk, stats_rad_zt%num_output_fields ) )

      call stats_zero( stats_rad_zt%kk, stats_rad_zt%num_output_fields, stats_rad_zt%accum_field_values, &
                     stats_rad_zt%accum_num_samples, stats_rad_zt%l_in_update )

      allocate( stats_rad_zt%file%var( stats_rad_zt%num_output_fields ) )
      allocate( stats_rad_zt%file%z( stats_rad_zt%kk ) )

       fname = trim( fname_rad_zt )

       call stats_init_rad_zt( clubb_vars_rad_zt, l_error )

       !  Initialize rad_zm (radiation points)

       i = 1
       do while ( ichar(clubb_vars_rad_zm(i)(1:1)) /= 0 .and. &
                  len_trim(clubb_vars_rad_zm(i))   /= 0 .and. &
                  i <= nvarmax_rad_zm )
          i = i + 1
       end do
       ntot = i - 1
       if ( ntot == nvarmax_rad_zm ) then
          write(fstderr,*) "There are more statistical variables listed in ",  &
                           "clubb_vars_rad_zm than allowed for by nvarmax_rad_zm."
          write(fstderr,*) "Check the number of variables listed for clubb_vars_rad_zm ",  &
                           "in the stats namelist, or change nvarmax_rad_zm."
          write(fstderr,*) "nvarmax_rad_zm = ", nvarmax_rad_zm
          call endrun ("stats_init_clubb:  number of rad_zm statistical variables exceeds limit"//errmsg(__FILE__,__LINE__))
       endif

       stats_rad_zm%num_output_fields = ntot
       stats_rad_zm%kk = nnrad_zm

       allocate( stats_rad_zm%z( stats_rad_zm%kk ) )

       allocate( stats_rad_zm%accum_field_values( 1, 1, stats_rad_zm%kk, stats_rad_zm%num_output_fields ) )
       allocate( stats_rad_zm%accum_num_samples( 1, 1, stats_rad_zm%kk, stats_rad_zm%num_output_fields ) )
       allocate( stats_rad_zm%l_in_update( 1, 1, stats_rad_zm%kk, stats_rad_zm%num_output_fields ) )

       call stats_zero( stats_rad_zm%kk, stats_rad_zm%num_output_fields, stats_rad_zm%accum_field_values, &
                     stats_rad_zm%accum_num_samples, stats_rad_zm%l_in_update )

       allocate( stats_rad_zm%file%var( stats_rad_zm%num_output_fields ) )
       allocate( stats_rad_zm%file%z( stats_rad_zm%kk ) )

       fname = trim( fname_rad_zm )

       call stats_init_rad_zm( clubb_vars_rad_zm, l_error )
    end if ! l_output_rad_files


    !  Initialize sfc (surface point)

    i = 1
    do while ( ichar(clubb_vars_sfc(i)(1:1)) /= 0 .and. &
               len_trim(clubb_vars_sfc(i))   /= 0 .and. &
               i <= nvarmax_sfc )
       i = i + 1
    end do
    ntot = i - 1
    if ( ntot == nvarmax_sfc ) then
       write(fstderr,*) "There are more statistical variables listed in ",  &
                        "clubb_vars_sfc than allowed for by nvarmax_sfc."
       write(fstderr,*) "Check the number of variables listed for clubb_vars_sfc ",  &
                        "in the stats namelist, or change nvarmax_sfc."
       write(fstderr,*) "nvarmax_sfc = ", nvarmax_sfc
       call endrun ("stats_init_clubb:  number of sfc statistical variables exceeds limit"//errmsg(__FILE__,__LINE__))
    endif

    stats_sfc%num_output_fields = ntot
    stats_sfc%kk = 1

    allocate( stats_sfc%z( stats_sfc%kk ) )

    allocate( stats_sfc%accum_field_values( 1, 1, stats_sfc%kk, stats_sfc%num_output_fields ) )
    allocate( stats_sfc%accum_num_samples( 1, 1, stats_sfc%kk, stats_sfc%num_output_fields ) )
    allocate( stats_sfc%l_in_update( 1, 1, stats_sfc%kk, stats_sfc%num_output_fields ) )

    call stats_zero( stats_sfc%kk, stats_sfc%num_output_fields, stats_sfc%accum_field_values, &
                     stats_sfc%accum_num_samples, stats_sfc%l_in_update )

    allocate( stats_sfc%file%var( stats_sfc%num_output_fields ) )
    allocate( stats_sfc%file%z( stats_sfc%kk ) )

    fname = trim( fname_sfc )

    call stats_init_sfc( clubb_vars_sfc, l_error )

    ! Check for errors

    if ( l_error ) then
       call endrun ('stats_init:  errors found'//errmsg(__FILE__,__LINE__))
    endif

!   Now call add fields
    do i = 1, stats_zt%num_output_fields

      temp1 = trim(stats_zt%file%var(i)%name)
      sub   = temp1
      if (len(temp1) .gt. 16) sub = temp1(1:16)

       call addfld(trim(sub),(/ 'ilev' /),&
            'A',trim(stats_zt%file%var(i)%units),trim(stats_zt%file%var(i)%description))
    enddo

    do i = 1, stats_zm%num_output_fields

      temp1 = trim(stats_zm%file%var(i)%name)
      sub   = temp1
      if (len(temp1) .gt. 16) sub = temp1(1:16)

      call addfld(trim(sub),(/ 'ilev' /),&
           'A',trim(stats_zm%file%var(i)%units),trim(stats_zm%file%var(i)%description))
    enddo

    if (l_output_rad_files) then
      do i = 1, stats_rad_zt%num_output_fields
        call addfld(trim(stats_rad_zt%file%var(i)%name),(/ 'foobar' /),&
           'A',trim(stats_rad_zt%file%var(i)%units),trim(stats_rad_zt%file%var(i)%description))
      enddo

      do i = 1, stats_rad_zm%num_output_fields
        call addfld(trim(stats_rad_zm%file%var(i)%name),(/ 'foobar' /),&
           'A',trim(stats_rad_zm%file%var(i)%units),trim(stats_rad_zm%file%var(i)%description))
      enddo
    endif

    do i = 1, stats_sfc%num_output_fields
      call addfld(trim(stats_sfc%file%var(i)%name),horiz_only,&
           'A',trim(stats_sfc%file%var(i)%units),trim(stats_sfc%file%var(i)%description))
    enddo

    return


  end subroutine stats_init_clubb

#endif

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !


    !-----------------------------------------------------------------------
  subroutine stats_end_timestep_clubb(lchnk,thecol,out_zt,out_zm,out_radzt,out_radzm,out_sfc)

    !     Description: Called when the stats timestep has ended. This subroutine
    !     is responsible for calling statistics to be written to the output
    !     format.
    !-----------------------------------------------------------------------

#ifdef CLUBB_SGS

    use constants_clubb, only: &
        fstderr ! Constant(s)

    use stats_variables, only: &
        stats_zt,  & ! Variable(s)
        stats_zm, &
        stats_rad_zt, &
        stats_rad_zm, &
        stats_sfc, &
        l_stats_last, &
        stats_tsamp, &
        stats_tout, &
        l_output_rad_files

    use error_code, only: &
        clubb_at_least_debug_level ! Procedure(s)

    use cam_history, only: outfld

    use ppgrid,      only: pcols, pverp

    use cam_abortutils,  only: endrun

    implicit none


#endif

    integer :: lchnk
    integer :: thecol

    real(r8), intent(inout) :: out_zt(:,:,:)     ! (pcols,pverp,zt%nn)
    real(r8), intent(inout) :: out_zm(:,:,:)     ! (pcols,pverp,zt%nn)
    real(r8), intent(inout) :: out_radzt(:,:,:)  ! (pcols,pverp,rad_zt%nn)
    real(r8), intent(inout) :: out_radzm(:,:,:)  ! (pcols,pverp,rad_zm%nn)
    real(r8), intent(inout) :: out_sfc(:,:,:)    ! (pcols,1,sfc%nn)

#ifdef CLUBB_SGS
    ! Local Variables

    integer :: i, k
    logical :: l_error

    !  Check if it is time to write to file

    if ( .not. l_stats_last ) return

    !  Initialize
    l_error = .false.

    !  Look for errors by checking the number of sampling points
    !  for each variable in the zt statistics at each vertical level.
    do i = 1, stats_zt%num_output_fields
      do k = 1, stats_zt%kk

        if ( stats_zt%accum_num_samples(1,1,k,i) /= 0 .and.  &
             stats_zt%accum_num_samples(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

          l_error = .true.  ! This will stop the run

          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) 'Possible sampling error for variable ',  &
                             trim(stats_zt%file%var(i)%name), ' in zt ',  &
                             'at k = ', k,  &
                             '; stats_zt%accum_num_samples(',k,',',i,') = ', stats_zt%accum_num_samples(1,1,k,i)
          endif

        endif

      enddo
    enddo

    !  Look for errors by checking the number of sampling points
    !  for each variable in the zm statistics at each vertical level.
    do i = 1, stats_zm%num_output_fields
      do k = 1, stats_zm%kk

        if ( stats_zm%accum_num_samples(1,1,k,i) /= 0 .and.  &
             stats_zm%accum_num_samples(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

          l_error = .true.  ! This will stop the run

          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) 'Possible sampling error for variable ',  &
                             trim(stats_zm%file%var(i)%name), ' in zm ',  &
                             'at k = ', k,  &
                             '; stats_zm%accum_num_samples(',k,',',i,') = ', stats_zm%accum_num_samples(1,1,k,i)
          endif

        endif

      enddo
    enddo

    if (l_output_rad_files) then
      !  Look for errors by checking the number of sampling points
      !  for each variable in the rad_zt statistics at each vertical level.
      do i = 1, stats_rad_zt%num_output_fields
        do k = 1, stats_rad_zt%kk

          if ( stats_rad_zt%accum_num_samples(1,1,k,i) /= 0 .and.  &
               stats_rad_zt%accum_num_samples(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

            l_error = .true.  ! This will stop the run

            if ( clubb_at_least_debug_level( 1 ) ) then
              write(fstderr,*) 'Possible sampling error for variable ',  &
                               trim(stats_rad_zt%file%var(i)%name), ' in rad_zt ',  &
                               'at k = ', k,  &
                               '; stats_rad_zt%accum_num_samples(',k,',',i,') = ', stats_rad_zt%accum_num_samples(1,1,k,i)
            endif

          endif

        enddo
      enddo

      !  Look for errors by checking the number of sampling points
      !  for each variable in the rad_zm statistics at each vertical level.
      do i = 1, stats_rad_zm%num_output_fields
        do k = 1, stats_rad_zm%kk

          if ( stats_rad_zm%accum_num_samples(1,1,k,i) /= 0 .and.  &
               stats_rad_zm%accum_num_samples(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

            l_error = .true.  ! This will stop the run

            if ( clubb_at_least_debug_level( 1 ) ) then
              write(fstderr,*) 'Possible sampling error for variable ',  &
                               trim(stats_rad_zm%file%var(i)%name), ' in rad_zm ',  &
                               'at k = ', k,  &
                               '; stats_rad_zm%accum_num_samples(',k,',',i,') = ', stats_rad_zm%accum_num_samples(1,1,k,i)
            endif

          endif

        enddo
      enddo
    end if ! l_output_rad_files

    !  Look for errors by checking the number of sampling points
    !  for each variable in the sfc statistics at each vertical level.
    do i = 1, stats_sfc%num_output_fields
      do k = 1, stats_sfc%kk

        if ( stats_sfc%accum_num_samples(1,1,k,i) /= 0 .and.  &
             stats_sfc%accum_num_samples(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

          l_error = .true.  ! This will stop the run

          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) 'Possible sampling error for variable ',  &
                             trim(stats_sfc%file%var(i)%name), ' in sfc ',  &
                             'at k = ', k,  &
                             '; stats_sfc%accum_num_samples(',k,',',i,') = ', stats_sfc%accum_num_samples(1,1,k,i)
          endif

        endif

      enddo
    enddo

    !  Stop the run if errors are found.
    if ( l_error ) then
       write(fstderr,*) 'Possible statistical sampling error'
       write(fstderr,*) 'For details, set debug_level to a value of at ',  &
                        'least 1 in the appropriate model.in file.'
       call endrun ('stats_end_timestep:  error(s) found'//errmsg(__FILE__,__LINE__))
    endif

    !  Compute averages
    call stats_avg( stats_zt%kk, stats_zt%num_output_fields, stats_zt%accum_field_values, stats_zt%accum_num_samples )
    call stats_avg( stats_zm%kk, stats_zm%num_output_fields, stats_zm%accum_field_values, stats_zm%accum_num_samples )
    if (l_output_rad_files) then
      call stats_avg( stats_rad_zt%kk, stats_rad_zt%num_output_fields, stats_rad_zt%accum_field_values, &
                      stats_rad_zt%accum_num_samples )
      call stats_avg( stats_rad_zm%kk, stats_rad_zm%num_output_fields, stats_rad_zm%accum_field_values, &
                      stats_rad_zm%accum_num_samples )
    end if
    call stats_avg( stats_sfc%kk, stats_sfc%num_output_fields, stats_sfc%accum_field_values, stats_sfc%accum_num_samples )

   !  Here we are not outputting the data, rather reading the stats into
   !  arrays which are conformable to CAM output.  Also, the data is "flipped"
   !  in the vertical level to be the same as CAM output.
    do i = 1, stats_zt%num_output_fields
      do k = 1, stats_zt%kk
         out_zt(thecol,k,i) = stats_zt%accum_field_values(1,1,stats_zt%kk-k+1,i)
         if(out_zt(thecol,k,i) .ne. out_zt(thecol,k,i)) out_zt(thecol,k,i) = 0.0_r8
      enddo
    enddo

    do i = 1, stats_zm%num_output_fields
      do k = 1, stats_zt%kk
         out_zm(thecol,k,i) = stats_zm%accum_field_values(1,1,stats_zt%kk-k+1,i)
         if(out_zm(thecol,k,i) .ne. out_zm(thecol,k,i)) out_zm(thecol,k,i) = 0.0_r8
      enddo
    enddo

    if (l_output_rad_files) then
      do i = 1, stats_rad_zt%num_output_fields
        do k = 1, stats_rad_zt%kk
          out_radzt(thecol,k,i) = stats_rad_zt%accum_field_values(1,1,stats_zt%kk-k+1,i)
          if(out_radzt(thecol,k,i) .ne. out_radzt(thecol,k,i)) out_radzt(thecol,k,i) = 0.0_r8
        enddo
      enddo

      do i = 1, stats_rad_zm%num_output_fields
        do k = 1, stats_rad_zm%kk
          out_radzm(thecol,k,i) = stats_rad_zm%accum_field_values(1,1,stats_zt%kk-k+1,i)
          if(out_radzm(thecol,k,i) .ne. out_radzm(thecol,k,i)) out_radzm(thecol,k,i) = 0.0_r8
        enddo
      enddo
    endif

    do i = 1, stats_sfc%num_output_fields
      out_sfc(thecol,1,i) = stats_sfc%accum_field_values(1,1,1,i)
      if(out_sfc(thecol,1,i) .ne. out_sfc(thecol,1,i)) out_sfc(thecol,1,i) = 0.0_r8
    enddo

    !  Reset sample fields
    call stats_zero( stats_zt%kk, stats_zt%num_output_fields, stats_zt%accum_field_values, &
                     stats_zt%accum_num_samples, stats_zt%l_in_update )
    call stats_zero( stats_zm%kk, stats_zm%num_output_fields, stats_zm%accum_field_values, &
                     stats_zm%accum_num_samples, stats_zm%l_in_update )
    if (l_output_rad_files) then
      call stats_zero( stats_rad_zt%kk, stats_rad_zt%num_output_fields, stats_rad_zt%accum_field_values, &
                       stats_rad_zt%accum_num_samples, stats_rad_zt%l_in_update )
      call stats_zero( stats_rad_zm%kk, stats_rad_zm%num_output_fields, stats_rad_zm%accum_field_values, &
                       stats_rad_zm%accum_num_samples, stats_rad_zm%l_in_update )
    end if
    call stats_zero( stats_sfc%kk, stats_sfc%num_output_fields, stats_sfc%accum_field_values, &
                     stats_sfc%accum_num_samples, stats_sfc%l_in_update )

    return

#endif

  end subroutine stats_end_timestep_clubb


  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

#ifdef CLUBB_SGS

    !-----------------------------------------------------------------------
  subroutine stats_zero( kk, nn, x, n, l_in_update )

    !     Description:
    !     Initialize stats to zero
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        stat_rknd,   & ! Variable(s)
        stat_nknd


    implicit none

    !  Input
    integer, intent(in) :: kk, nn

    !  Output
    real(kind=stat_rknd),    dimension(1,1,kk,nn), intent(out) :: x
    integer(kind=stat_nknd), dimension(1,1,kk,nn), intent(out) :: n
    logical,                 dimension(1,1,kk,nn), intent(out) :: l_in_update

    !  Zero out arrays

    if ( nn > 0 ) then
       x(:,:,:,:) = 0.0_r8
       n(:,:,:,:) = 0
       l_in_update(:,:,:,:) = .false.
    end if

    return

  end subroutine stats_zero

#endif

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !


#ifdef CLUBB_SGS
    !-----------------------------------------------------------------------
  subroutine stats_avg( kk, nn, x, n )

    !     Description:
    !     Compute the average of stats fields
    !-----------------------------------------------------------------------
    use clubb_precision, only: &
        stat_rknd,   & ! Variable(s)
        stat_nknd

    implicit none

    !  Input
    integer, intent(in) :: nn, kk
    integer(kind=stat_nknd), dimension(1,1,kk,nn), intent(in) :: n

    !  Output
    real(kind=stat_rknd), dimension(1,1,kk,nn), intent(inout)  :: x

    !  Internal

    integer k,m

    !  Compute averages

    do m=1,nn
       do k=1,kk

          if ( n(1,1,k,m) > 0 ) then
             x(1,1,k,m) = x(1,1,k,m) / real( n(1,1,k,m) )
          end if

       end do
    end do

    return

  end subroutine stats_avg

#endif

end module clubb_intr
