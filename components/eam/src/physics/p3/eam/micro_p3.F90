!__________________________________________________________________________________________
! This module contains the Predicted Particle Property (P3) bulk microphysics scheme.      !
!                                                                                          !
! This code was originally written by H. Morrison,  MMM Division, NCAR (Dec 2012).         !
! Modifications were made by J. Milbrandt, RPN, Environment Canada (July 2014).            !
! Peter Caldwell (caldwell19@llnl.gov) further modified this code to remove multiple       !
! ice categories and to clean up/simplify the code for conversion to C++ (9/11/18)         !
! Jacob Shpund (jacob.shpund@pnnl.gov) implemented and further modified/cleaned            !
! the scheme in E3SMv2                                                                     ! 
!                                                                                          !
! Three configurations of the P3 scheme are currently available:                           !
!  1) specified droplet number (i.e. 1-moment cloud water), 1 ice category                 !
!  2) predicted droplet number (i.e. 2-moment cloud water), 1 ice category                 !
!                                                                                          !
!  The  2-moment cloud version is based on a specified aerosol distribution and            !
!  does not include a subgrid-scale vertical velocity for droplet activation. Hence,       !
!  this version should only be used for high-resolution simulations that resolve           !
!  vertical motion driving droplet activation.                                             !
!                                                                                          !
! For details see: Morrison and Milbrandt (2015) [J. Atmos. Sci., 72, 287-311]             !
!                  Milbrandt and Morrison (2016) [J. Atmos. Sci., 73, 975-995]             !
!                                                                                          !
! For questions or bug reports, please contact:                                            !
!    Hugh Morrison   (morrison@ucar.edu), or                                               !
!    Jason Milbrandt (jason.milbrandt@canada.ca)                                           !
!__________________________________________________________________________________________!
!                                                                                          !
! Version:       2.8.2.4 + Peter/Aaron's fixes                                             !
! Last updated:  2018-02-04                                                                !
!__________________________________________________________________________________________!
! Comments from Aaron Donahue:                                                             !
! 1) Need to change the dz coordinate system in sedimentation to be consistent             !
! with E3SM's pressure based coordinate system, i.e. dp.                                   !
! 2) Move all physical constants into a micro_p3_util module and match them to             !
! universal constants in E3SM for consistency.                                            !
! 3) Need to include extra in/out values which correspond with microphysics PBUF           !
! variables and outputs expected in E3SM.                                                  !
!__________________________________________________________________________________________!

! Include bit-for-bit math macros.
! #include "bfb_math.inc"

module micro_p3

   ! get real kind from utils
   use physics_utils, only: rtype,rtype8,btype

   use phys_control,  only: use_hetfrz_classnuc

   ! physical and mathematical constants
   use micro_p3_utils, only: rho_1000mb,rho_600mb,ar,br,f1r,f2r,rho_h2o,kr,kc,aimm,mi0,nccnst,  &
       eci,eri,bcn,cpw,cons1,cons3,cons4,cons5,cons6,cons7,         &
       inv_rho_h2o,inv_dropmass,qsmall,nsmall,cp,g,rd,rv,ep_2,inv_cp,   &
       thrd,sxth,piov6,rho_rimeMin,     &
       rho_rimeMax,inv_rho_rimeMax,max_total_ni,dbrk,nmltratio,clbfact_sub,  &
       clbfact_dep,iparam, isize, densize, rimsize, rcollsize, ice_table_size, collect_table_size, &
       get_latent_heat, T_zerodegc, pi=>pi_e3sm, dnu, &
       T_rainfrz, T_icenuc, T_homogfrz, iulog=>iulog_e3sm, &
       masterproc=>masterproc_e3sm, calculate_incloud_mixingratios, mu_r_constant, &
       lookup_table_1a_dum1_c, rho_h2o, &
       do_Cooper_inP3

   use wv_sat_scream, only:qv_sat

  ! Bit-for-bit math functions.
#ifdef SCREAM_CONFIG_IS_CMAKE
  use physics_share_f2c, only: cxx_pow, cxx_sqrt, cxx_cbrt, cxx_gamma, cxx_log, &
                                 cxx_log10, cxx_exp, cxx_expm1, cxx_tanh
#endif

  implicit none
  save

  public  :: p3_init,p3_main

  ! protected items should be treated as private for everyone except tests

  real(rtype), protected, dimension(densize,rimsize,isize,ice_table_size) :: ice_table_vals   !ice lookup table values

  !ice lookup table values for ice-rain collision/collection
  real(rtype), protected, dimension(densize,rimsize,isize,rcollsize,collect_table_size) :: collect_table_vals

  ! lookup table values for rain shape parameter mu_r
  real(rtype), protected, dimension(150) :: mu_r_table_vals

  ! lookup table values for rain number- and mass-weighted fallspeeds and ventilation parameters
  real(rtype), protected, dimension(300,10) :: vn_table_vals,vm_table_vals,revap_table_vals

!<shanyp 20221218
! set lookup tables as threadprivate variables to see if this change can address the threading issue.
!  !$OMP THREADPRIVATE(ice_table_vals,collect_table_vals,mu_r_table_vals,vn_table_vals,vm_table_vals,revap_table_vals)
! !$OMP THREADPRIVATE(ice_table_vals) 
 !,vn_table_vals,vm_table_vals)
!shanyp 20221218>


  type realptr
     real(rtype), dimension(:), pointer :: p
  end type realptr

contains

  SUBROUTINE p3_init(lookup_file_dir,version_p3)
    !------------------------------------------------------------------------------------------!
    ! This subroutine initializes all physical constants and parameters needed by the P3       !
    ! scheme, including reading in two lookup table files and creating a third.                !
    ! 'P3_INIT' be called at the first model time step, prior to first call to 'P3_MAIN'.      !
    !------------------------------------------------------------------------------------------!

    implicit none

    ! Passed arguments:
    character*(*), intent(in)    :: lookup_file_dir                !directory of the lookup tables
    character(len=16), intent(in) :: version_p3  !version number of P3 package

  END SUBROUTINE p3_init

  SUBROUTINE p3_main(qc,nc,qr,nr,th_atm,qv,dt,qi,qm,ni,bm,                                                                                                               &
       pres,dz,nc_nuceat_tend,nccn_prescribed,ni_activated,frzimm,frzcnt,frzdep,inv_qc_relvar,it,precip_liq_surf,precip_ice_surf,its,ite,kts,kte,diag_eff_radius_qc,     &
       diag_eff_radius_qi,rho_qi,do_predict_nc, do_prescribed_CCN,p3_autocon_coeff,p3_accret_coeff,p3_qc_autocon_expon,p3_nc_autocon_expon,p3_qc_accret_expon,           &
       p3_wbf_coeff,p3_max_mean_rain_size,p3_embryonic_rain_size,                                                                                                        &
       dpres,exner,qv2qi_depos_tend,precip_total_tend,nevapr,qr_evap_tend,precip_liq_flux,precip_ice_flux,rflx,sflx,cflx,cld_frac_r,cld_frac_l,cld_frac_i,               &
       p3_tend_out,mu_c,lamc,liq_ice_exchange,vap_liq_exchange,                                                                                                          &
       vap_ice_exchange,qv_prev,t_prev,col_location,diag_equiv_reflectivity,diag_ze_rain,diag_ze_ice                                                                     &
#ifdef SCREAM_CONFIG_IS_CMAKE
       ,elapsed_s &
#endif
      )

    !----------------------------------------------------------------------------------------!
    !                                                                                        !
    ! This is the main subroutine for the P3 microphysics scheme.  It is called from the     !
    ! wrapper subroutine ('MP_P3_WRAPPER') and is passed i,k slabs of all prognostic         !
    ! variables -- hydrometeor fields, potential temperature, and water vapor mixing ratio.  !
    ! Microphysical process rates are computed first.  These tendencies are then used to     !
    ! computed updated values of the prognostic variables.  The hydrometeor variables are    !
    ! then updated further due to sedimentation.                                             !
    !                                                                                        !
    ! Several diagnostic values are also computed and returned to the wrapper subroutine,    !
    ! including precipitation rates.                                                         !
    !                                                                                        !
    !----------------------------------------------------------------------------------------!

    implicit none

    !----- Input/ouput arguments:  ----------------------------------------------------------!

    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: qc         ! cloud, mass mixing ratio         kg kg-1
    ! note: Nc may be specified or predicted (set by do_predict_nc)
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: nc         ! cloud, number mixing ratio       #  kg-1
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: qr         ! rain, mass mixing ratio          kg kg-1
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: nr         ! rain, number mixing ratio        #  kg-1
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: qi      ! ice, total mass mixing ratio     kg kg-1
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: qm      ! ice, rime mass mixing ratio      kg kg-1
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: ni      ! ice, total number mixing ratio   #  kg-1
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: bm      ! ice, rime volume mixing ratio    m3 kg-1

    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: qv         ! water vapor mixing ratio         kg kg-1
    real(rtype), intent(inout), dimension(its:ite,kts:kte)      :: th_atm         ! potential temperature            K
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: pres       ! pressure                         Pa
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: dz        ! vertical grid spacing            m
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: nc_nuceat_tend      ! IN ccn activated number tendency kg-1 s-1
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: nccn_prescribed
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: ni_activated       ! IN actived ice nuclei concentration  1/kg
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: frzimm,frzcnt,frzdep ! From macrophysics aerop (CNT scheme) [#/cm3]
    real(rtype), intent(in)                                     :: dt         ! model time step                  s

    real(rtype), intent(out),   dimension(its:ite)              :: precip_liq_surf    ! precipitation rate, liquid       m s-1
    real(rtype), intent(out),   dimension(its:ite)              :: precip_ice_surf    ! precipitation rate, solid        m s-1
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: diag_eff_radius_qc  ! effective radius, cloud          m
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: diag_eff_radius_qi  ! effective radius, ice            m
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: rho_qi  ! bulk density of ice              kg m-3
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: mu_c       ! Size distribution shape parameter for radiation
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: lamc       ! Size distribution slope parameter for radiation

    integer, intent(in)                                  :: its,ite    ! array bounds (horizontal)
    integer, intent(in)                                  :: kts,kte    ! array bounds (vertical)
    integer, intent(in)                                  :: it         ! time step counter NOTE: starts at 1 for first time step

    logical(btype), intent(in)                           :: do_predict_nc ! .T. (.F.) for prediction (specification) of Nc

    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: dpres       ! pressure thickness               Pa
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: exner      ! Exner expression

    ! OUTPUT for PBUF variables used by other parameterizations
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: qv2qi_depos_tend    ! qitend due to deposition/sublimation
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: precip_total_tend      ! Total precipitation (rain + snow)
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: nevapr     ! evaporation of total precipitation (rain + snow)
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: qr_evap_tend  ! evaporation of rain
    real(rtype), intent(out),   dimension(its:ite,kts:kte+1)    :: precip_liq_flux       ! grid-box average rain flux (kg m^-2 s^-1) pverp
    real(rtype), intent(out),   dimension(its:ite,kts:kte+1)    :: precip_ice_flux       ! grid-box average ice/snow flux (kg m^-2 s^-1) pverp
    real(rtype), intent(out),   dimension(its:ite,kts:kte+1)    :: rflx       ! grid-box average rain flux (kg m^-2 s^-1) pverp
    real(rtype), intent(out),   dimension(its:ite,kts:kte+1)    :: sflx       ! grid-box average ice/snow flux (kg m^-2 s^-1) pverp
    real(rtype), intent(out),   dimension(its:ite,kts:kte+1)    :: cflx       ! grid-box average cloud droplets flux (kg m^-2 s^-1) pverp
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: liq_ice_exchange ! sum of liq-ice phase change tendenices
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: vap_liq_exchange ! sum of vap-liq phase change tendenices
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: vap_ice_exchange ! sum of vap-ice phase change tendenices

    ! INPUT for prescribed CCN option
    logical(btype), intent(in)                                  :: do_prescribed_CCN

    ! INPUT for p3 tuning parameters
    real(rtype), intent(in)                                     :: p3_autocon_coeff         ! autconversion coefficient
    real(rtype), intent(in)                                     :: p3_accret_coeff          ! accretion coefficient
    real(rtype), intent(in)                                     :: p3_qc_autocon_expon      ! autconversion qc exponent
    real(rtype), intent(in)                                     :: p3_nc_autocon_expon      ! autconversion nc exponent
    real(rtype), intent(in)                                     :: p3_qc_accret_expon       ! accretion qc and qr exponent
    real(rtype), intent(in)                                     :: p3_wbf_coeff             ! WBF coefficient
    real(rtype), intent(in)                                     :: p3_max_mean_rain_size    ! max mean rain size allowed
    real(rtype), intent(in)                                     :: p3_embryonic_rain_size   ! embryonic rain size from autoconversion

    ! INPUT needed for PBUF variables used by other parameterizations

    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: cld_frac_i, cld_frac_l, cld_frac_r ! Ice, Liquid and Rain cloud fraction
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: qv_prev, t_prev                    ! qv and t from previous p3_main call
    ! AaronDonahue, the following variable (p3_tend_out) is a catch-all for passing P3-specific variables outside of p3_main
    ! so that they can be written as ouput.  NOTE TO C++ PORT: This variable is entirely optional and doesn't need to be
    ! included in the port to C++, or can be changed if desired.
    real(rtype), intent(out),   dimension(its:ite,kts:kte,49)   :: p3_tend_out ! micro physics tendencies
    real(rtype), intent(in),    dimension(its:ite,3)            :: col_location
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: inv_qc_relvar
    real(rtype), intent(out),   dimension(its:ite,kts:kte)      :: diag_equiv_reflectivity,diag_ze_rain,diag_ze_ice  ! equivalent reflectivity [dBZ]

#ifdef SCREAM_CONFIG_IS_CMAKE
    real(rtype), optional, intent(out) :: elapsed_s ! duration of main loop in seconds
#endif

    !
    !----- Local variables and parameters:  -------------------------------------------------!
    !

    ! These outputs are no longer provided by p3_main.
    real(rtype), dimension(its:ite,kts:kte) :: diag_vm_qi ! mass-weighted fall speed of ice  m s-1
    real(rtype), dimension(its:ite,kts:kte) :: diag_diam_qi  ! mean diameter of ice             m
    real(rtype), dimension(its:ite,kts:kte) :: pratot   ! accretion of cloud by rain
    real(rtype), dimension(its:ite,kts:kte) :: prctot   ! autoconversion of cloud to rain

    real(rtype), dimension(its:ite,kts:kte) :: mu_r  ! shape parameter of rain
    real(rtype), dimension(its:ite,kts:kte) :: t_atm     ! temperature at the beginning of the microhpysics step [K]

    ! 2D size distribution and fallspeed parameters:

    real(rtype), dimension(its:ite,kts:kte) :: lamr
    real(rtype), dimension(its:ite,kts:kte) :: logn0r

    real(rtype), dimension(its:ite,kts:kte) :: nu
    real(rtype), dimension(its:ite,kts:kte) :: cdist
    real(rtype), dimension(its:ite,kts:kte) :: cdist1
    real(rtype), dimension(its:ite,kts:kte) :: cdistr

    ! Variables needed for in-cloud calculations
    real(rtype), dimension(its:ite,kts:kte) :: inv_cld_frac_i, inv_cld_frac_l, inv_cld_frac_r ! Inverse cloud fractions (1/cld)
    real(rtype), dimension(its:ite,kts:kte) :: qc_incld, qr_incld, qi_incld, qm_incld ! In cloud mass-mixing ratios
    real(rtype), dimension(its:ite,kts:kte) :: nc_incld, nr_incld, ni_incld, bm_incld ! In cloud number concentrations

    real(rtype), dimension(its:ite,kts:kte)      :: inv_dz,inv_rho,ze_ice,ze_rain,prec,rho,       &
         rhofacr,rhofaci,acn,latent_heat_sublim,latent_heat_vapor,latent_heat_fusion,qv_sat_l,qv_sat_i,qv_supersat_i,       &
         tmparr1,inv_exner

    ! -- scalar locals -- !

    real(rtype) :: inv_dt, timeScaleFactor

    integer :: ktop,kbot,kdir,i

    logical(btype) :: is_nucleat_possible, is_hydromet_present

    !--These will be added as namelist parameters in the future
    logical(btype), parameter :: debug_ON     = .true.  !.true. to switch on debugging checks/traps throughout code  TODO: Turn this back off as default once the tlay error is found.
    logical(btype), parameter :: debug_ABORT  = .false.  !.true. will result in forced abort in s/r 'check_values'

    real(rtype),dimension(its:ite,kts:kte) :: qc_old, nc_old, qr_old, nr_old, qi_old, ni_old, qv_old, th_atm_old

#ifdef SCREAM_CONFIG_IS_CMAKE
    integer :: clock_count1, clock_count_rate, clock_count_max, clock_count2, clock_count_diff
#endif

  END SUBROUTINE p3_main

end module micro_p3
