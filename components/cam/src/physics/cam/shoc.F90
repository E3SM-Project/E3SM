!c_doubli--------------------------------------------------------------
! SHOC parameterization
!   SHOC = Simplified Higher Order Closure
!   reference, Bogenschutz and Krueger 2013
!
! PDF-based parameterization for low clouds and turbulence
! email: bogenschutz1@llnl.gov
!--------------------------------------------------------------

! Include bit-for-bit math macros.
#include "bfb_math.inc"

module shoc

  use physics_utils, only: rtype, rtype8, itype, btype
  use scream_abortutils, only: endscreamrun

! Bit-for-bit math functions.
#ifdef SCREAM_CONFIG_IS_CMAKE
  use physics_share_f2c, only: cxx_pow, cxx_sqrt, cxx_cbrt, cxx_gamma, cxx_log, &
                               cxx_log10, cxx_exp
#endif

implicit none
save ! for module variables

public  :: shoc_init, shoc_main

logical :: use_cxx = .true.

real(rtype), parameter, public :: largeneg = -99999999.99_rtype

!=========================================================
! Physical constants used in SHOC
!=========================================================

! These are set in initialization and should be set to
!  to the values used in whatever host model SHOC is
!  implemented in
real(rtype) :: ggr   ! gravity [m/s^2]
real(rtype) :: rgas  ! dry air gas constant [J/kg.K]
real(rtype) :: rv    ! water vapor gas constant [J/kg.K]
real(rtype) :: cp    ! specific heat of dry air [J/kg.K]
real(rtype) :: lcond ! latent heat of vaporization [J/kg]
real(rtype) :: lice  ! latent heat of fusion [J/kg]
real(rtype) :: eps   ! rh2o/rair - 1 [-]
real(rtype) :: vk    ! von karmann constant [-]

!=========================================================
! Private module parameters
!=========================================================

! These are adjustable tunable parameters for SHOC
!  for certain variables and turbulent moments.
!  All are unitless

! temperature variance
real(rtype), parameter :: thl2tune=1.0_rtype
! moisture variance
real(rtype), parameter :: qw2tune=1.0_rtype
! temp moisture covariance
real(rtype), parameter :: qwthl2tune=1.0_rtype
! vertical velocity variance
real(rtype), parameter :: w2tune=1.0_rtype
! third moment of vertical velocity
real(rtype), parameter :: w3clip=1.2_rtype
! mixing length scaling parameter
real(rtype), parameter :: length_fac=0.5_rtype
! coefficient for diag third moment parameters
real(rtype), parameter :: c_diag_3rd_mom = 7.0_rtype

! =========
! Below are options to activate certain features in SHOC

! Allow temperature skewness to be independent of moisture
!  variance
logical(btype), parameter :: dothetal_skew = .false.

! ========
! Below define some parameters for SHOC

! Reference temperature [K]
real(rtype), parameter :: basetemp = 300._rtype
! Reference pressure [Pa]
real(rtype), parameter :: basepres = 100000._rtype
! Lower troposphere pressure [Pa]
real(rtype), parameter :: troppres = 80000._rtype
! Minimum surface friction velocity
real(rtype), parameter :: ustar_min = 0.01_rtype
! PBL max depth in pressure units
real(rtype), parameter :: pblmaxp = 4.e4_rtype

! ========
! Set upper limits for certain SHOC quantities
! Note that these upper limits are quite high
! and they are rarely reached in a stable simulation

! Mixing length [m]
real(rtype), parameter :: maxlen = 20000.0_rtype
! Minimum Mixing length [m]
real(rtype), parameter :: minlen = 20.0_rtype
! Maximum TKE [m2/s2]
real(rtype), parameter :: maxtke = 50.0_rtype
! Minimum TKE [m2/s2]
real(rtype), parameter :: mintke = 0.0004_rtype

!===================
! const parameter for Diagnosis of PBL depth
real(rtype), parameter :: tiny = 1.e-36_rtype     ! lower bound for wind magnitude
real(rtype), parameter :: fac  = 100._rtype       ! ustar parameter in height diagnosis
real(rtype), parameter :: ricr  =  0.3_rtype      ! Critical richardson number

! Maximum number of levels in pbl from surface
integer :: npbl

!==============================================================
! Begin SHOC parameterization code!
contains
!==============================================================

subroutine shoc_init( &
         nlev, gravit, rair, rh2o, cpair, &
         zvir, latvap, latice, karman, &
         pref_mid, nbot_shoc, ntop_shoc)

  implicit none

  ! Purpose:  Initialize constants for SHOC
  ! These should be set to the same values used in
  ! whatever host model SHOC is implemented in

  integer, intent(in)   :: nlev ! number of levels

  real(rtype), intent(in)  :: gravit ! gravity
  real(rtype), intent(in)  :: rair   ! dry air gas constant
  real(rtype), intent(in)  :: rh2o   ! water vapor gas constant
  real(rtype), intent(in)  :: cpair  ! specific heat of dry air
  real(rtype), intent(in)  :: zvir   ! rh2o/rair - 1
  real(rtype), intent(in)  :: latvap ! latent heat of vaporization
  real(rtype), intent(in)  :: latice ! latent heat of fusion
  real(rtype), intent(in)  :: karman ! Von Karman's constant

  real(rtype), intent(in) :: pref_mid(nlev) ! reference pressures at midpoints

  integer, intent(in)   :: nbot_shoc ! Bottom level to which SHOC is applied
  integer, intent(in)   :: ntop_shoc ! Top level to which SHOC is applied

  integer :: k

  ggr = gravit   ! [m/s2]
  rgas = rair    ! [J/kg.K]
  rv = rh2o      ! [J/kg.K]
  cp = cpair     ! [J/kg.K]
  eps = zvir     ! [-]
  lcond = latvap ! [J/kg]
  lice = latice  ! [J/kg]
  vk = karman    ! [-]

   ! Limit pbl height to regions below 400 mb
   ! npbl = max number of levels (from bottom) in pbl

   npbl = 0
   do k=nbot_shoc,ntop_shoc,-1
      if (pref_mid(k) >= pblmaxp) then
         npbl = npbl + 1
      end if
   end do
   npbl = max(npbl,1)

   return

end subroutine shoc_init

!==============================================================
! Main driver for the SHOC scheme
! Host models should call the following routine to call SHOC

subroutine shoc_main ( &
     shcol, nlev, nlevi, dtime, nadv, &   ! Input
     host_dx, host_dy,thv, &              ! Input
     zt_grid,zi_grid,pres,presi,pdel,&    ! Input
     wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, & ! Input
     wtracer_sfc,num_qtracers,w_field, &  ! Input
     exner,phis, &                        ! Input
     host_dse, tke, thetal, qw, &         ! Input/Output
     u_wind, v_wind,qtracers,&            ! Input/Output
     wthv_sec,tkh,tk,&                    ! Input/Output
     shoc_ql,shoc_cldfrac,&               ! Input/Output
     pblh,&                               ! Output
     shoc_mix, isotropy,&                 ! Output (diagnostic)
     w_sec, thl_sec, qw_sec, qwthl_sec,&  ! Output (diagnostic)
     wthl_sec, wqw_sec, wtke_sec,&        ! Output (diagnostic)
     uw_sec, vw_sec, w3,&                 ! Output (diagnostic)
     wqls_sec, brunt, shoc_ql2)           ! Output (diagnostic)

  implicit none

! INPUT VARIABLES
  ! number of SHOC columns in the array
  integer, intent(in) :: shcol
  ! number of levels [-]
  integer, intent(in) :: nlev
  ! number of levels on interface grid [-]
  integer, intent(in) :: nlevi
  ! number of tracers [-]
  integer, intent(in) :: num_qtracers
  ! number of times to loop SHOC
  integer, intent(in) :: nadv

  ! SHOC timestep [s]
  real(rtype), intent(in) :: dtime
  ! grid spacing of host model in x direction [m]
  real(rtype), intent(in) :: host_dx(shcol)
  ! grid spacing of host model in y direction [m]
  real(rtype), intent(in) :: host_dy(shcol)
  ! heights, for thermo grid [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! heights, for interface grid [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)
  ! pressure levels on thermo grid [Pa]
  real(rtype), intent(in) :: pres(shcol,nlev)
  ! pressure levels on interface grid [Pa]
  real(rtype), intent(in) :: presi(shcol,nlevi)
  ! Differences in pressure levels [Pa]
  real(rtype), intent(in) :: pdel(shcol,nlev)
  ! virtual potential temperature [K]
  real(rtype), intent(in) :: thv(shcol,nlev)
  ! large scale vertical velocity [m/s]
  real(rtype), intent(in) :: w_field(shcol,nlev)
  ! Surface sensible heat flux [K m/s]
  real(rtype), intent(in) :: wthl_sfc(shcol)
  ! Surface latent heat flux [kg/kg m/s]
  real(rtype), intent(in) :: wqw_sfc(shcol)
  ! Surface momentum flux (u-direction) [m2/s2]
  real(rtype), intent(in) :: uw_sfc(shcol)
  ! Surface momentum flux (v-direction) [m2/s2]
  real(rtype), intent(in) :: vw_sfc(shcol)
  ! Surface flux for tracers [varies]
  real(rtype), intent(in) :: wtracer_sfc(shcol,num_qtracers)
  ! Exner function [-]
  real(rtype), intent(in) :: exner(shcol,nlev)
  ! Host model surface geopotential height
  real(rtype), intent(in) :: phis(shcol)

! INPUT/OUTPUT VARIABLES
  ! prognostic temp variable of host model
  ! dry static energy [J/kg]
  ! dse = Cp*T + g*z + phis
  real(rtype), intent(inout) :: host_dse(shcol,nlev)
  ! turbulent kinetic energy [m2/s2]
  real(rtype), intent(inout) :: tke(shcol,nlev)
  ! liquid water potential temperature [K]
  real(rtype), intent(inout) :: thetal(shcol,nlev)
  ! total water mixing ratio [kg/kg]
  real(rtype), intent(inout) :: qw(shcol,nlev)
  ! u wind component [m/s]
  real(rtype), intent(inout) :: u_wind(shcol,nlev)
  ! v wind component [m/s]
  real(rtype), intent(inout) :: v_wind(shcol,nlev)
  ! buoyancy flux [K m/s]
  real(rtype), intent(inout) :: wthv_sec(shcol,nlev)
  ! tracers [varies]
  real(rtype), intent(inout) :: qtracers(shcol,nlev,num_qtracers)
  ! eddy coefficient for momentum [m2/s]
  real(rtype), intent(inout) :: tk(shcol,nlev)
  ! eddy coefficent for heat [m2/s]
  real(rtype), intent(inout) :: tkh(shcol,nlev)
  ! Cloud fraction [-]
  real(rtype), intent(inout) :: shoc_cldfrac(shcol,nlev)
  ! cloud liquid mixing ratio [kg/kg]
  real(rtype), intent(inout) :: shoc_ql(shcol,nlev)

  ! OUTPUT VARIABLES

  ! planetary boundary layer depth [m]
  real(rtype), intent(out) :: pblh(shcol)
  ! cloud liquid mixing ratio variance [kg^2/kg^2]
  real(rtype), intent(out) :: shoc_ql2(shcol, nlev)

  ! also output variables, but part of the SHOC diagnostics
  !  to be output to history file by host model (if desired)

  ! Turbulent length scale [m]
  real(rtype) :: shoc_mix(shcol,nlev)
  ! vertical velocity variance [m2/s2]
  real(rtype) :: w_sec(shcol,nlev)
  ! temperature variance [K^2]
  real(rtype) :: thl_sec(shcol,nlevi)
  ! moisture variance [kg2/kg2]
  real(rtype) :: qw_sec(shcol,nlevi)
  ! temp moisture covariance [K kg/kg]
  real(rtype) :: qwthl_sec(shcol,nlevi)
  ! vertical heat flux [K m/s]
  real(rtype) :: wthl_sec(shcol,nlevi)
  ! vertical moisture flux [K m/s]
  real(rtype) :: wqw_sec(shcol,nlevi)
  ! vertical tke flux [m3/s3]
  real(rtype) :: wtke_sec(shcol,nlevi)
  ! vertical zonal momentum flux [m2/s2]
  real(rtype) :: uw_sec(shcol,nlevi)
  ! vertical meridional momentum flux [m2/s2]
  real(rtype) :: vw_sec(shcol,nlevi)
  ! third moment vertical velocity [m3/s3]
  real(rtype) :: w3(shcol,nlevi)
  ! liquid water flux [kg/kg m/s]
  real(rtype) :: wqls_sec(shcol,nlev)
  ! brunt vaisala frequency [s-1]
  real(rtype) :: brunt(shcol,nlev)
  ! return to isotropic timescale [s]
  real(rtype) :: isotropy(shcol,nlev)

  !============================================================================
! LOCAL VARIABLES

  ! time counter
  integer :: t

  ! air density on thermo grid [kg/m3]
  real(rtype) :: rho_zt(shcol,nlev)

  ! Grid difference centereted on thermo grid [m]
  real(rtype) :: dz_zt(shcol,nlev)
  ! Grid difference centereted on interface grid [m]
  real(rtype) :: dz_zi(shcol,nlevi)

  ! Surface friction velocity [m/s]
  real(rtype) :: ustar(shcol)
  ! Monin Obukhov length [m]
  real(rtype) :: obklen(shcol)
  ! Kinematic surface buoyancy flux [m^2/s^3]
  real(rtype) :: kbfs(shcol)

  ! Variables related to energy conservation
  real(rtype) :: se_b(shcol),ke_b(shcol),&
              wv_b(shcol),wl_b(shcol),&
              se_a(shcol),ke_a(shcol),&
              wv_a(shcol),wl_a(shcol)

  ! Compute integrals of static energy, kinetic energy, water vapor, and liquid water
  ! for the computation of total energy before SHOC is called.  This is for an
  ! effort to conserve energy since liquid water potential temperature (which SHOC
  ! conserves) and static energy (which E3SM conserves) are not exactly equal.
  call shoc_energy_integrals(&
     shcol,nlev,host_dse,pdel,&             ! Input
     qw,shoc_ql,u_wind,v_wind,&             ! Input
     se_b,ke_b,wv_b,wl_b)                   ! Input/Output

  do t=1,nadv

    ! Check TKE to make sure values lie within acceptable
    !  bounds after host model performs horizontal advection
    call check_tke(shcol,nlev,&                 ! Input
           tke)                                 ! Input/Output

    ! Define vertical grid arrays needed for
    !   vertical derivatives in SHOC, also
    !   define air density
    call shoc_grid( &
       shcol,nlev,nlevi,&                   ! Input
       zt_grid,zi_grid,pdel,&               ! Input
       dz_zt,dz_zi,rho_zt)          ! Output

    ! Compute the planetary boundary layer height, which is an
    !   input needed for the length scale calculation.

    call shoc_diag_obklen(&
       shcol,uw_sfc,vw_sfc,&                          ! Input
       wthl_sfc,wqw_sfc,thetal(:shcol,nlev),&         ! Input
       shoc_ql(:shcol,nlev),qtracers(:shcol,nlev,1),& ! Input
       ustar,kbfs,obklen)                             ! Output

    call pblintd(&
       shcol,nlev,nlevi,&                   ! Input
       zt_grid,zi_grid,thetal,shoc_ql,&     ! Input
       qtracers(:shcol,:,1),u_wind,v_wind,& ! Input
       ustar,obklen,kbfs,shoc_cldfrac,&     ! Input
       pblh)                                ! Output

    ! Update the turbulent length scale
    call shoc_length(&
       shcol,nlev,nlevi,tke,&               ! Input
       host_dx,host_dy,pblh,&               ! Input
       zt_grid,zi_grid,dz_zt,dz_zi,&        ! Input
       thetal,wthv_sec,thv,&                ! Input
       brunt,shoc_mix)                      ! Output

    ! Advance the SGS TKE equation
    call shoc_tke(&
       shcol,nlev,nlevi,dtime,&             ! Input
       wthv_sec,shoc_mix,&                  ! Input
       dz_zi,dz_zt,pres,&                   ! Input
       u_wind,v_wind,brunt,obklen,&         ! Input
       zt_grid,zi_grid,pblh,&               ! Input
       tke,tk,tkh,&                         ! Input/Output
       isotropy)                            ! Output

    ! Update SHOC prognostic variables here
    !   via implicit diffusion solver
    call update_prognostics_implicit(&      ! Input
       shcol,nlev,nlevi,num_qtracers,&      ! Input
       dtime,dz_zt,dz_zi,rho_zt,&           ! Input
       zt_grid,zi_grid,tk,tkh,&             ! Input
       uw_sfc,vw_sfc,wthl_sfc,wqw_sfc,&     ! Input
       thetal,qw,qtracers,tke,&             ! Input/Output
       u_wind,v_wind)                       ! Input/Output

    ! Diagnose the second order moments
    call diag_second_shoc_moments(&
       shcol,nlev,nlevi, &                    ! Input
       num_qtracers,thetal,qw, &              ! Input
       u_wind,v_wind,qtracers,tke, &          ! Input
       isotropy,tkh,tk,&                      ! Input
       dz_zi,zt_grid,zi_grid,shoc_mix, &      ! Input
       wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, &   ! Input
       wtracer_sfc, &                         ! Input
       thl_sec, qw_sec,wthl_sec,wqw_sec,&     ! Output
       qwthl_sec, uw_sec, vw_sec, wtke_sec, & ! Output
       w_sec)                                 ! Output

    ! Diagnose the third moment of vertical velocity,
    !  needed for the PDF closure
    call diag_third_shoc_moments(&
       shcol,nlev,nlevi,&                   ! Input
       w_sec,thl_sec,qw_sec,qwthl_sec,&     ! Input
       wthl_sec,isotropy,brunt,&            ! Input
       thetal,tke,wthv_sec,&                ! Input
       dz_zt,dz_zi,zt_grid,zi_grid,&        ! Input
       w3)                                  ! Output

    ! Call the PDF to close on SGS cloud and turbulence
    call shoc_assumed_pdf(&
       shcol,nlev,nlevi,&                   ! Input
       thetal,qw,w_field,thl_sec,qw_sec,&   ! Input
       wthl_sec,w_sec,&                     ! Input
       wqw_sec,qwthl_sec,w3,pres,&          ! Input
       zt_grid,zi_grid,&                    ! Input
       shoc_cldfrac,shoc_ql,&               ! Output
       wqls_sec,wthv_sec,shoc_ql2)          ! Output

    ! Check TKE to make sure values lie within acceptable
    !  bounds after vertical advection, etc.
    call check_tke(shcol,nlev,tke)

  enddo ! end time loop

  ! End SHOC parameterization

  ! Use SHOC outputs to update the host model
  !  temperature
  call update_host_dse(&
     shcol,nlev,thetal,&                   ! Input
     shoc_ql,exner,zt_grid,phis,&          ! Input
     host_dse)                             ! Output

  call shoc_energy_integrals(&             ! Input
     shcol,nlev,host_dse,pdel,&            ! Input
     qw,shoc_ql,u_wind,v_wind,&            ! Input
     se_a,ke_a,wv_a,wl_a)                  ! Output

  call shoc_energy_fixer(&
     shcol,nlev,nlevi,dtime,nadv,&         ! Input
     zt_grid,zi_grid,&                     ! Input
     se_b,ke_b,wv_b,wl_b,&                 ! Input
     se_a,ke_a,wv_a,wl_a,&                 ! Input
     wthl_sfc,wqw_sfc,pdel,&               ! Input
     rho_zt,tke,presi,&                    ! Input
     host_dse)                             ! Input/Output

  ! Remaining code is to diagnose certain quantities
  !  related to PBL.  No answer changing subroutines
  !  should be placed at this point onward.

  ! Update PBLH, as other routines outside of SHOC
  !  may require this variable.
  call shoc_diag_obklen(&
     shcol,uw_sfc,vw_sfc,&                          ! Input
     wthl_sfc,wqw_sfc,thetal(:shcol,nlev),&         ! Input
     shoc_ql(:shcol,nlev),qtracers(:shcol,nlev,1),& ! Input
     ustar,kbfs,obklen)                             ! Output

  call pblintd(&
     shcol,nlev,nlevi,&                   ! Input
     zt_grid,zi_grid,thetal,shoc_ql,&     ! Input
     qtracers(:shcol,:,1),u_wind,v_wind,& ! Input
     ustar,obklen,kbfs,shoc_cldfrac,&     ! Input
     pblh)                                ! Output
  return

end subroutine shoc_main

!==============================================================
! Define grid variables needed for the parameterization

subroutine shoc_grid( &
          shcol,nlev,nlevi,&           ! Input
          zt_grid,zi_grid,pdel,&       ! Input
          dz_zt,dz_zi,rho_zt)          ! Output

  ! Purpose of this subroutine is to define the thickness
  !  arrays of each column, to be used for finite differencing
  !  throughout the SHOC parameterization, also define air
  !  density in SHOC

  implicit none

! INPUT VARIABLES
  ! number of columns [-]
  integer, intent(in) :: shcol
  ! number of mid-point levels [-]
  integer, intent(in) :: nlev
  ! number of interface levels [-]
  integer, intent(in) :: nlevi
  ! mid-point grid heights [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! interface grid heights [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)
  ! pressure differences centered on mid-point grid [Pa]
  real(rtype), intent(in) :: pdel(shcol,nlev)

! OUTPUT VARIABLES
  ! thickness (dz) on the thermo grid [m]
  real(rtype), intent(out) :: dz_zt(shcol,nlev)
  ! thickness (dz) on the interface grid [m]
  real(rtype), intent(out) :: dz_zi(shcol,nlevi)
  ! air density on the thermo grid [kg/m3]
  real(rtype), intent(out) :: rho_zt(shcol,nlev)

  ! local variables
  integer :: i, k
  do k=1,nlev
    do i=1,shcol
      ! define thickness of the thermodynamic gridpoints
      dz_zt(i,k) = zi_grid(i,k) - zi_grid(i,k+1)

      ! define thickness of the interface grid points
      if (k .eq. 1) then
        dz_zi(i,k) = 0._rtype ! never used
      else
        dz_zi(i,k) = zt_grid(i,k-1) - zt_grid(i,k)
      endif

      ! Define the air density on the thermo grid
      rho_zt(i,k) = (1._rtype/ggr)*(pdel(i,k)/dz_zt(i,k))

    enddo ! end i loop (column loop)
  enddo ! end k loop (vertical loop)

  ! Set lower condition for dz_zi
  dz_zi(:shcol,nlevi) = zt_grid(:shcol,nlev)

  return

end subroutine shoc_grid

!==============================================================
! Update T, q, tracers, tke, u, and v based on implicit diffusion
! Here we use a backward Euler scheme.

subroutine update_prognostics_implicit( &
         shcol,nlev,nlevi,num_tracer,&    ! Input
         dtime,dz_zt,dz_zi,rho_zt,&       ! Input
         zt_grid,zi_grid,tk,tkh,&         ! Input
         uw_sfc,vw_sfc,wthl_sfc,wqw_sfc,& ! Input
         thetal,qw,tracer,tke,&           ! Input/Output
         u_wind,v_wind)                   ! Input/Output

  implicit none

! INPUT VARIABLES
  ! number of SHOC columns
  integer, intent(in) :: shcol
  ! number of vertical levels
  integer, intent(in) :: nlev
  ! number of interface levels
  integer, intent(in) :: nlevi
  ! number of tracers
  integer, intent(in) :: num_tracer

  ! SHOC timestep [s]
  real(rtype), intent(in) :: dtime
  ! Eddy coefficient for momentum [m2/s]
  real(rtype), intent(in) :: tk(shcol,nlev)
  ! Eddy coefficient for heat [m2/s]
  real(rtype), intent(in) :: tkh(shcol,nlev)
  ! Air density on thermo grid [kg/m3]
  real(rtype), intent(in) :: rho_zt(shcol,nlev)
  ! height thickness centered on thermo grid [m]
  real(rtype), intent(in) :: dz_zt(shcol,nlev)
  ! height thickness centered on interface grid [m]
  real(rtype), intent(in) :: dz_zi(shcol,nlevi)
  ! vertical zonal momentum flux at surface [m3/s3]
  real(rtype), intent(in) :: uw_sfc(shcol)
  ! vertical meridional momentum flux at surface [m3/s3]
  real(rtype), intent(in) :: vw_sfc(shcol)
  ! vertical heat flux at surface [K m/s]
  real(rtype), intent(in) :: wthl_sfc(shcol)
  ! vertical moisture flux at surface [kg/kg m/s]
  real(rtype), intent(in) :: wqw_sfc(shcol)
  ! heights of mid-point [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! heights at interfaces [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)

! IN/OUT VARIABLES
  ! liquid water potential temperature [K]
  real(rtype), intent(inout) :: thetal(shcol,nlev)
  ! total water mixing ratio [kg/kg]
  real(rtype), intent(inout) :: qw(shcol,nlev)
  ! tracers [varies]
  real(rtype), intent(inout) :: tracer(shcol,nlev,num_tracer)
  ! zonal wind [m/s]
  real(rtype), intent(inout) :: u_wind(shcol,nlev)
  ! meridional wind [m/s]
  real(rtype), intent(inout) :: v_wind(shcol,nlev)
  ! turbulent kinetic energy [m2/s2]
  real(rtype), intent(inout) :: tke(shcol,nlev)

! LOCAL VARIABLES
  integer     :: p
  real(rtype) :: rdp_zt(shcol,nlev)
  real(rtype) :: tmpi(shcol,nlevi)
  real(rtype) :: tkh_zi(shcol,nlevi)
  real(rtype) :: tk_zi(shcol,nlevi)
  real(rtype) :: rho_zi(shcol,nlevi)

  real(rtype) :: flux_dummy(shcol)
  real(rtype) :: ksrf(shcol), wtke_flux(shcol)

  real(rtype) :: ca(shcol,nlev) ! superdiagonal for solver
  real(rtype) :: cc(shcol,nlev) ! subdiagonal for solver
  real(rtype) :: denom(shcol,nlev) ! denominator in solver
  real(rtype) :: ze(shcol,nlev)

  ! linearly interpolate tkh, tk, and air density onto the interface grids
  call linear_interp(zt_grid,zi_grid,tkh,tkh_zi,nlev,nlevi,shcol,0._rtype)
  call linear_interp(zt_grid,zi_grid,tk,tk_zi,nlev,nlevi,shcol,0._rtype)
  call linear_interp(zt_grid,zi_grid,rho_zt,rho_zi,nlev,nlevi,shcol,0._rtype)

  tmpi(:,1) = 0._rtype
  ! Define the tmpi variable, which is really dt*(g*rho)**2/dp
  !  at interfaces. Substitue dp = g*rho*dz in the above equation
  call compute_tmpi(nlevi, shcol, dtime, rho_zi, dz_zi, tmpi)

  ! compute 1/dp term, needed in diffusion solver
  call dp_inverse(nlev, nlevi, shcol, rho_zt, dz_zt, rdp_zt)

  ! compute terms needed for the implicit surface stress (ksrf)
  ksrf(1:shcol)      = impli_srf_stress_term(shcol, nlev, nlevi, rho_zi, &
      uw_sfc, vw_sfc, u_wind, v_wind)

  !compute term needed for tke flux calc (wtke_flux)
  wtke_flux(1:shcol) = tke_srf_flux_term(shcol, uw_sfc, vw_sfc)

  ! compute surface fluxes for liq. potential temp, water and tke
  call sfc_fluxes(shcol, nlev, nlevi, dtime, rho_zi, rdp_zt, &
       wthl_sfc, wqw_sfc, wtke_flux, &
       thetal, qw, tke)

  ! Call decomp for momentum variables
  call vd_shoc_decomp(shcol,nlev,nlevi,tk_zi,tmpi,rdp_zt,dtime,&
     ksrf,ca,cc,denom,ze)

  ! march u_wind one step forward using implicit solver
  call vd_shoc_solve(shcol,nlev,nlevi,ca,cc,denom,ze,u_wind)

  ! march v_wind one step forward using implicit solver
  call vd_shoc_solve(shcol,nlev,nlevi,ca,cc,denom,ze,v_wind)

! Call decomp for thermo variables
  flux_dummy(:) = 0._rtype ! fluxes applied explicitly, so zero fluxes out
                           ! for implicit solver decomposition
  call vd_shoc_decomp(shcol,nlev,nlevi,tkh_zi,tmpi,rdp_zt,dtime,&
     flux_dummy,ca,cc,denom,ze)

  ! march temperature one step forward using implicit solver
  call vd_shoc_solve(shcol,nlev,nlevi,ca,cc,denom,ze,thetal)

  ! march total water one step forward using implicit solver
  call vd_shoc_solve(shcol,nlev,nlevi,ca,cc,denom,ze,qw)

  ! march tke one step forward using implicit solver
  call vd_shoc_solve(shcol,nlev,nlevi,ca,cc,denom,ze,tke)

  ! march tracers one step forward using implicit solver
  do p=1,num_tracer
    call vd_shoc_solve(shcol,nlev,nlevi,ca,cc,denom,ze,tracer(:shcol,:nlev,p))
  enddo

  return

end subroutine update_prognostics_implicit

subroutine compute_tmpi(nlevi, shcol, dtime, rho_zi, dz_zi, tmpi)

  !intent-ins
  integer,     intent(in) :: nlevi, shcol
  !time step [s]
  real(rtype), intent(in) :: dtime
  !air density at interfaces [kg/m3]
  real(rtype), intent(in) :: rho_zi(shcol,nlevi)
  !height thickness at interfaces [m]
  real(rtype), intent(in) :: dz_zi(shcol,nlevi)

  !intent-out
  real(rtype), intent(out) :: tmpi(shcol,nlevi)

  !local vars
  integer :: i, k

  tmpi(:,1) = 0._rtype
  ! eqn: tmpi = dt*(g*rho)**2/dp, where dp = g*rho*dz, therefore tmpi = dt*g*rho/dz
  do k = 2, nlevi
    do i = 1, shcol
       tmpi(i,k) = dtime * (ggr*rho_zi(i,k)) / dz_zi(i,k)
    enddo
  enddo

end subroutine compute_tmpi

subroutine dp_inverse(nlev, nlevi, shcol, rho_zt, dz_zt, rdp_zt)

  !intent-ins
  integer,     intent(in) :: nlev, nlevi, shcol
  ! Air density on thermo grid [kg/m3]
  real(rtype), intent(in) :: rho_zt(shcol,nlev)
  ! height thickness centered on thermo grid [m]
  real(rtype), intent(in) :: dz_zt(shcol,nlev)

  !intent-out
  real(rtype), intent(out) :: rdp_zt(shcol,nlev)

  !local vars
  integer :: i, k

  do k = 1, nlev
    do i = 1, shcol
      rdp_zt(i,k) = 1._rtype/(ggr*rho_zt(i,k)*dz_zt(i,k))
    enddo
  enddo

end subroutine dp_inverse

pure function impli_srf_stress_term(shcol, nlev, nlevi, rho_zi, uw_sfc, &
     vw_sfc, u_wind, v_wind) result (ksrf)

  !intent-ins
  integer,     intent(in) :: shcol, nlev, nlevi

  !air density at interfaces [kg/m3]
  real(rtype), intent(in) :: rho_zi(shcol,nlevi)
  !vertical zonal momentum flux at surface [m3/s3]
  real(rtype), intent(in) :: uw_sfc(shcol)
  !vertical meridional momentum flux at surface [m3/s3]
  real(rtype), intent(in) :: vw_sfc(shcol)
  !zonal wind [m/s]
  real(rtype), intent(in) :: u_wind(shcol,nlev)
  !meridional wind [m/s]
  real(rtype), intent(in) :: v_wind(shcol,nlev)

  !function return value
  real(rtype) :: ksrf(shcol)

  !local vars
  integer :: i

  real(rtype) :: taux, tauy !stresses (N/m2)
  real(rtype) :: ws         !wind speed (m/s)
  real(rtype) :: rho, tau, uw, vw, rho_zi_srf(shcol)

  real(rtype), parameter :: wsmin    = 1._rtype    ! Minimum wind speed for ksrfturb computation [ m/s ]
  real(rtype), parameter :: ksrfmin  = 1.e-4_rtype ! Minimum surface drag coefficient  [ kg/s/m^2 ]

  !store surface values of rho in a 1d array
  rho_zi_srf(1:shcol) = rho_zi(1:shcol,nlevi)

  do i = 1, shcol
     rho          = rho_zi_srf(i)
     uw           = uw_sfc(i)
     vw           = vw_sfc(i)

     taux         = rho*uw ! stress in N/m2
     tauy         = rho*vw ! stress in N/m2
     ! compute the wind speed
     ws           = max(sqrt(u_wind(i,nlev)**2._rtype + v_wind(i,nlev)**2._rtype),wsmin)
     tau          = sqrt( taux**2._rtype + tauy**2._rtype )
     ksrf(i)      = max(tau/ws, ksrfmin)
  enddo

  return
end function impli_srf_stress_term

pure function tke_srf_flux_term(shcol, uw_sfc, vw_sfc) result(wtke_flux)

  !intent-ins
  integer,     intent(in) :: shcol

  !vertical zonal momentum flux at surface [m3/s3]
  real(rtype), intent(in) :: uw_sfc(shcol)
  !vertical meridional momentum flux at surface [m3/s3]
  real(rtype), intent(in) :: vw_sfc(shcol)

  !function return value
  real(rtype) :: wtke_flux(shcol)

  !local vars
  integer :: i

  real(rtype) :: ustar, uw, vw

  real(rtype), parameter :: ustarmin = 0.01_rtype  ! Minimum ustar

  do i = 1, shcol
     uw           = uw_sfc(i)
     vw           = vw_sfc(i)
     ustar        = max(sqrt(sqrt(uw**2._rtype + vw**2._rtype)),ustarmin)
     wtke_flux(i) = ustar**3
  enddo

  return
end function tke_srf_flux_term


subroutine sfc_fluxes(shcol, nlev, nlevi, dtime, rho_zi, rdp_zt, &
     wthl_sfc, wqw_sfc, wtke_flux, thetal, qw, tke)

  implicit none

  !intent-ins
  integer,     intent(in) :: shcol, nlev, nlevi
  !time step [s]
  real(rtype), intent(in) :: dtime
  !air density at interfaces [kg/m3]
  real(rtype), intent(in) :: rho_zi(shcol,nlevi)
  !inverse of dp
  real(rtype), intent(in) :: rdp_zt(shcol,nlev)
  !vertical heat flux at surface [K m/s]
  real(rtype), intent(in) :: wthl_sfc(shcol)
  !vertical moisture flux at surface [kg/kg m/s]
  real(rtype), intent(in) :: wqw_sfc(shcol)
  !vertical tke flux at surface [m3/s3]
  real(rtype), intent(in) :: wtke_flux(shcol)

  !intent-inouts
  !liquid water potential temperature [K]
  real(rtype), intent(inout) :: thetal(shcol,nlev)
  !total water mixing ratio [kg/kg]
  real(rtype), intent(inout) :: qw(shcol,nlev)
  !turbulent kinetic energy [m2/s2]
  real(rtype), intent(inout) :: tke(shcol,nlev)

  !local variables
  integer :: i
  real(rtype) :: cmnfac

  ! Apply the surface fluxes explicitly for temperature and moisture
  do i = 1, shcol
     cmnfac       =  dtime * (ggr * rho_zi(i,nlevi) * rdp_zt(i,nlev)) !a common factor for the following 3 equations

     thetal(i,nlev) = thetal(i,nlev) + cmnfac * wthl_sfc(i)
     qw(i,nlev)     = qw(i,nlev)     + cmnfac * wqw_sfc(i)
     tke(i,nlev)    = tke(i,nlev)    + cmnfac * wtke_flux(i)
  enddo

end subroutine sfc_fluxes


!=======================================================
! SHOC Diagnose the second order moments,
!  main routine

subroutine diag_second_shoc_moments(&
         shcol,nlev,nlevi, &                    ! Input
         num_tracer,thetal,qw, &                ! Input
         u_wind,v_wind,tracer,tke, &            ! Input
         isotropy,tkh,tk,&                      ! Input
         dz_zi,zt_grid,zi_grid,shoc_mix, &      ! Input
         wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, &   ! Input
         wtracer_sfc, &                         ! Input
         thl_sec,qw_sec,wthl_sec,wqw_sec,&      ! Output
         qwthl_sec, uw_sec, vw_sec, wtke_sec, & ! Output
         w_sec)                                 ! Output

  ! This is the main routine to compute the second
  !   order moments in SHOC.

  implicit none

! INPUT VARIABLES
  ! number of SHOC columns
  integer, intent(in) :: shcol
  ! number of midpoint levels
  integer, intent(in) :: nlev
  ! number of interface levels
  integer, intent(in) :: nlevi
  ! number of tracers
  integer, intent(in) :: num_tracer

  ! liquid water potential temperature [K]
  real(rtype), intent(in) :: thetal(shcol,nlev)
  ! total water mixing ratio [kg/kg]
  real(rtype), intent(in) :: qw(shcol,nlev)
  ! zonal wind component [m/s]
  real(rtype), intent(in) :: u_wind(shcol,nlev)
  ! meridional wind component [m/s]
  real(rtype), intent(in) :: v_wind(shcol,nlev)
  ! turbulent kinetic energy [m2/s2]
  real(rtype), intent(in) :: tke(shcol,nlev)
  ! return to isotropy timescale [s]
  real(rtype), intent(in) :: isotropy(shcol,nlev)
  ! eddy coefficient for heat [m2/s]
  real(rtype), intent(in) :: tkh(shcol,nlev)
  ! eddy coefficient for momentum [m2/s]
  real(rtype), intent(in) :: tk(shcol,nlev)
  ! tracers [varies]
  real(rtype), intent(in) :: tracer(shcol,nlev,num_tracer) ! tracers
  ! heights of mid-point grid [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! heights of interface grid [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)
  ! thickness centered on interface grid [m]
  real(rtype), intent(in) :: dz_zi(shcol,nlevi)
  ! Mixing length [m]
  real(rtype), intent(in) :: shoc_mix(shcol,nlev)
  ! Surface sensible heat flux [K m/s]
  real(rtype), intent(in) :: wthl_sfc(shcol)
  ! Surface latent heat flux [kg/kg m/s]
  real(rtype), intent(in) :: wqw_sfc(shcol)
  ! Surface momentum flux (u-direction) [m2/s2]
  real(rtype), intent(in) :: uw_sfc(shcol)
  ! Surface momentum flux (v-direction) [m2/s2]
  real(rtype), intent(in) :: vw_sfc(shcol)
  ! Tracer flux [varies m/s]
  real(rtype), intent(in) :: wtracer_sfc(shcol,num_tracer)

! OUTPUT VARIABLES
  ! second order liquid wat. potential temp. [K^2]
  real(rtype), intent(out) :: thl_sec(shcol,nlevi)
  ! second order total water mixing rat. [kg^2/kg^2]
  real(rtype), intent(out) :: qw_sec(shcol,nlevi)
  ! covariance of temp and moisture [K kg/kg]
  real(rtype), intent(out) :: qwthl_sec(shcol,nlevi)
  ! vertical flux of heat [K m/s]
  real(rtype), intent(out) :: wthl_sec(shcol,nlevi)
  ! vertical flux of total water [kg/kg m/s]
  real(rtype), intent(out) :: wqw_sec(shcol,nlevi)
  ! vertical flux of zonal wind [m2/s2]
  real(rtype), intent(out) :: uw_sec(shcol,nlevi)
  ! vertical flux of meridional wind [m2/s2]
  real(rtype), intent(out) :: vw_sec(shcol,nlevi)
  ! vertical flux of tke [m3/s3]
  real(rtype), intent(out) :: wtke_sec(shcol,nlevi)
  ! second order vertical velocity [m2/s2]
  real(rtype), intent(out) :: w_sec(shcol,nlev)

! LOCAL VARIABLES
  real(rtype) :: wstar(shcol)
  real(rtype) :: ustar2(shcol)

  ! Calculate surface properties needed for lower
  !  boundary conditions
  call diag_second_moments_srf(&
     shcol,       &                         ! Input
     wthl_sfc, uw_sfc, vw_sfc, &            ! Input
     ustar2,wstar)                          ! Output

  ! Diagnose the second order moments flux,
  !  for the lower boundary
  call diag_second_moments_lbycond(&
     shcol, num_tracer,&                             ! Input
     wthl_sfc, wqw_sfc, uw_sfc, vw_sfc,&             ! Input
     wtracer_sfc,ustar2,wstar,&                      ! Input
     wthl_sec(:shcol,nlevi),wqw_sec(:shcol,nlevi),&  ! Output
     uw_sec(:shcol,nlevi), vw_sec(:shcol,nlevi),&    ! Output
     wtke_sec(:shcol,nlevi), thl_sec(:shcol,nlevi),& ! Output
     qw_sec(:shcol,nlevi), qwthl_sec(:shcol,nlevi))  ! Output

  ! Diagnose the second order moments,
  !  for points away from boundaries.  this is
  !  the main computation for the second moments
  call diag_second_moments(&
     shcol,nlev,nlevi, &                    ! Input
     num_tracer,thetal,qw, &                ! Input
     u_wind,v_wind,tracer,tke, &            ! Input
     isotropy,tkh,tk,&                      ! Input
     dz_zi,zt_grid,zi_grid,shoc_mix, &      ! Input
     thl_sec, qw_sec,wthl_sec,wqw_sec,&     ! Input/Output
     qwthl_sec, uw_sec, vw_sec, wtke_sec, & ! Input/Output
     w_sec)                                 ! Output

  ! Diagnose the second order moments,
  !  calculate the upper boundary conditions
  call diag_second_moments_ubycond(&
     shcol,                              &  ! Input
     thl_sec(:shcol,1), qw_sec(:shcol,1),&  ! Output
     wthl_sec(:shcol,1),wqw_sec(:shcol,1),& ! Output
     qwthl_sec(:shcol,1), uw_sec(:shcol,1),&! Output
     vw_sec(:shcol,1), wtke_sec(:shcol,1))  ! Output

  return
end subroutine diag_second_shoc_moments

!==============================================================
! SHOC Diagnose the second order moments,
!  lower boundary conditions

subroutine diag_second_moments_srf(&
         shcol,       &                         ! Input
         wthl_sfc, uw_sfc, vw_sfc, &            ! Input
         ustar2,wstar)                          ! Output

  ! Purpose of this subroutine is to diagnose surface
  !  properties needed for the the lower
  !  boundary condition for the second order moments needed
  !  for the SHOC parameterization.
#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: shoc_diag_second_moments_srf_f
#endif

  implicit none

! INPUT VARIABLES
  ! number of SHOC columns
  integer, intent(in) :: shcol

  ! Surface sensible heat flux [K m/s]
  real(rtype), intent(in) :: wthl_sfc(shcol)
  ! Surface momentum flux (u-direction) [m2/s2]
  real(rtype), intent(in) :: uw_sfc(shcol)
  ! Surface momentum flux (v-direction) [m2/s2]
  real(rtype), intent(in) :: vw_sfc(shcol)

! OUTPUT VARIABLES
  ! Surface friction velocity [m4/s4]
  real(rtype), intent(out) :: ustar2(shcol)
  ! Surface convective velocity [m/s]
  real(rtype), intent(out) :: wstar(shcol)

! LOCAL VARIABLES
  integer :: i, p

  ! Constants to parameterize surface variances
  real(rtype), parameter :: z_const = 1.0_rtype

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
      call shoc_diag_second_moments_srf_f(shcol,wthl_sfc, uw_sfc, vw_sfc, &            ! Input
            ustar2,wstar)                          ! Output
      return
   endif
#endif

  ! apply the surface conditions to diagnose turbulent
  !  moments at the surface
  do i=1,shcol

    ! Parameterize thermodyanmics variances via Andre et al. 1978
    ustar2(i) = bfb_sqrt(uw_sfc(i) * uw_sfc(i) + vw_sfc(i) * vw_sfc(i))
    if (wthl_sfc(i) > 0._rtype) then
      wstar(i) = bfb_pow((1._rtype/basetemp * ggr * wthl_sfc(i) * z_const), (1._rtype/3._rtype))
    else
      wstar(i) = 0._rtype
    endif

  enddo ! end i loop (column loop)
  return
end subroutine diag_second_moments_srf

!==============================================================
! SHOC Diagnose the second order moments flux,
!  lower boundary conditions

subroutine diag_second_moments_lbycond(&
         shcol,num_tracer,&                           ! Input
         wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, &         ! Input
         wtracer_sfc,ustar2,wstar,&                   ! Input
         wthl_sec,wqw_sec,&                           ! Output
         uw_sec, vw_sec, wtke_sec,&                   ! Output
         thl_sec,qw_sec,qwthl_sec)                    ! Output

  ! Purpose of this subroutine is to diagnose the lower
  !  boundary condition for the second order moments needed
  !  for the SHOC parameterization.
  ! The thermodymnamic, tracer, and momentum fluxes are set
  !  to the surface fluxes for the host model, while the
  !  thermodynamic variances and covariances are computed
  !  according to that of Andre et al. 1978.

  implicit none

! INPUT VARIABLES
  ! number of SHOC columns
  integer, intent(in) :: shcol
  ! number of tracers
  integer, intent(in) :: num_tracer

  ! Surface sensible heat flux [K m/s]
  real(rtype), intent(in) :: wthl_sfc(shcol)
  ! Surface latent heat flux [kg/kg m/s]
  real(rtype), intent(in) :: wqw_sfc(shcol)
  ! Surface momentum flux (u-direction) [m2/s2]
  real(rtype), intent(in) :: uw_sfc(shcol)
  ! Surface momentum flux (v-direction) [m2/s2]
  real(rtype), intent(in) :: vw_sfc(shcol)
  ! Tracer flux [varies m/s]
  real(rtype), intent(in) :: wtracer_sfc(shcol,num_tracer)
  ! Surface friction velocity squared [m4/s4]
  real(rtype), intent(in) :: ustar2(shcol)
  ! Surface convective velocity scale [m/s]
  real(rtype), intent(in) :: wstar(shcol)

! OUTPUT VARIABLES
  ! vertical flux of heat [K m/s]
  real(rtype), intent(out) :: wthl_sec(shcol)
  ! vertical flux of total water [kg/kg m/s]
  real(rtype), intent(out) :: wqw_sec(shcol)
  ! vertical flux of zonal wind [m2/s2]
  real(rtype), intent(out) :: uw_sec(shcol)
  ! vertical flux of meridional wind [m2/s2]
  real(rtype), intent(out) :: vw_sec(shcol)
  ! vertical flux of tke [m3/s3]
  real(rtype), intent(out) :: wtke_sec(shcol)
  ! second order liquid wat. potential temp. [K^2]
  real(rtype), intent(out) :: thl_sec(shcol)
  ! second order total water mixing rat. [kg^2/kg^2]
  real(rtype), intent(out) :: qw_sec(shcol)
  ! covariance of temp and moisture [K kg/kg]
  real(rtype), intent(out) :: qwthl_sec(shcol)

! LOCAL VARIABLES
  integer :: i, p
  real(rtype) :: uf

  ! Constants to parameterize surface variances
  real(rtype), parameter :: a_const = 1.8_rtype
  real(rtype), parameter :: ufmin = 0.01_rtype

  ! apply the surface conditions to diagnose turbulent
  !  moments at the surface
  do i=1,shcol

    uf = sqrt(ustar2(i) + 0.3_rtype * wstar(i) * wstar(i))
    uf = max(ufmin,uf)

    ! Diagnose thermodynamics variances and covariances
    thl_sec(i) = 0.4_rtype * a_const * (wthl_sfc(i)/uf)**2
    qw_sec(i) = 0.4_rtype * a_const * (wqw_sfc(i)/uf)**2
    qwthl_sec(i) = 0.2_rtype * a_const * (wthl_sfc(i)/uf) * &
                         (wqw_sfc(i)/uf)

    ! Vertical fluxes of heat and moisture, simply
    !  use the surface fluxes given by host model
    wthl_sec(i) = wthl_sfc(i)
    wqw_sec(i) = wqw_sfc(i)
    uw_sec(i) = uw_sfc(i)
    vw_sec(i) = vw_sfc(i)
    wtke_sec(i) = max(sqrt(ustar2(i)),0.01_rtype)**3

  enddo ! end i loop (column loop)
  return
end subroutine diag_second_moments_lbycond

subroutine diag_second_moments(&
         shcol,nlev,nlevi, &                    ! Input
         num_tracer,thetal,qw, &                ! Input
         u_wind,v_wind,tracer,tke, &            ! Input
         isotropy,tkh,tk,&                      ! Input
         dz_zi,zt_grid,zi_grid,shoc_mix, &      ! Input
         thl_sec,qw_sec,wthl_sec,wqw_sec,&      ! Input/Output
         qwthl_sec, uw_sec, vw_sec, wtke_sec, & ! Input/Output
         w_sec)                                 ! Output

  ! Purpose of this subroutine is to diagnose the second
  !  order moments needed for the SHOC parameterization.
  !  Namely these are variances of thetal, qw, and vertical
  !  velocity.  In addition the vertical fluxes of thetal, qw,
  !  u, v, TKE, and tracers are computed here as well as the
  !  correlation of qw and thetal.

  implicit none

! INPUT VARIABLES
  ! number of SHOC columns
  integer, intent(in) :: shcol
  ! number of midpoint levels
  integer, intent(in) :: nlev
  ! number of interface levels
  integer, intent(in) :: nlevi
  ! number of tracers
  integer, intent(in) :: num_tracer

  ! liquid water potential temperature [K]
  real(rtype), intent(in) :: thetal(shcol,nlev)
  ! total water mixing ratio [kg/kg]
  real(rtype), intent(in) :: qw(shcol,nlev)
  ! zonal wind component [m/s]
  real(rtype), intent(in) :: u_wind(shcol,nlev)
  ! meridional wind component [m/s]
  real(rtype), intent(in) :: v_wind(shcol,nlev)
  ! turbulent kinetic energy [m2/s2]
  real(rtype), intent(in) :: tke(shcol,nlev)
  ! return to isotropy timescale [s]
  real(rtype), intent(in) :: isotropy(shcol,nlev)
  ! eddy coefficient for heat [m2/s]
  real(rtype), intent(in) :: tkh(shcol,nlev)
  ! eddy coefficient for momentum [m2/s]
  real(rtype), intent(in) :: tk(shcol,nlev)
  ! tracers [varies]
  real(rtype), intent(in) :: tracer(shcol,nlev,num_tracer) ! tracers
  ! heights of mid-point grid [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! heights of interface grid [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)
  ! thickness centered on interface grid [m]
  real(rtype), intent(in) :: dz_zi(shcol,nlevi)
  ! Mixing length [m]
  real(rtype), intent(in) :: shoc_mix(shcol,nlev)

! INPUT/OUTPUT VARIABLES
  ! second order liquid wat. potential temp. [K^2]
  real(rtype), intent(inout) :: thl_sec(shcol,nlevi)
  ! second order total water mixing rat. [kg^2/kg^2]
  real(rtype), intent(inout) :: qw_sec(shcol,nlevi)
  ! covariance of temp and moisture [K kg/kg]
  real(rtype), intent(inout) :: qwthl_sec(shcol,nlevi)
  ! vertical flux of heat [K m/s]
  real(rtype), intent(inout) :: wthl_sec(shcol,nlevi)
  ! vertical flux of total water [kg/kg m/s]
  real(rtype), intent(inout) :: wqw_sec(shcol,nlevi)
  ! vertical flux of zonal wind [m2/s2]
  real(rtype), intent(inout) :: uw_sec(shcol,nlevi)
  ! vertical flux of meridional wind [m2/s2]
  real(rtype), intent(inout) :: vw_sec(shcol,nlevi)
  ! vertical flux of tke [m3/s3]
  real(rtype), intent(inout) :: wtke_sec(shcol,nlevi)

! OUTPUT VARIABLES
  ! second order vertical velocity [m2/s2]
  real(rtype), intent(out) :: w_sec(shcol,nlev)

  ! LOCAL VARIABLES
  integer :: p
  real(rtype) :: isotropy_zi(shcol,nlevi)
  real(rtype) :: tkh_zi(shcol,nlevi)
  real(rtype) :: tk_zi(shcol,nlevi)

  ! Interpolate some variables from the midpoint grid to the interface grid
  call linear_interp(zt_grid,zi_grid,isotropy,isotropy_zi,nlev,nlevi,shcol,0._rtype)
  call linear_interp(zt_grid,zi_grid,tkh,tkh_zi,nlev,nlevi,shcol,0._rtype)
  call linear_interp(zt_grid,zi_grid,tk,tk_zi,nlev,nlevi,shcol,0._rtype)

  ! Vertical velocity variance is assumed to be propotional
  !  to the TKE
  w_sec = w2tune*(2._rtype/3._rtype)*tke

  ! Calculate the temperature variance
  call calc_shoc_varorcovar(&
         shcol,nlev,nlevi,thl2tune,&              ! Input
         isotropy_zi,tkh_zi,dz_zi,thetal,thetal,& ! Input
         thl_sec)                                 ! Input/Output

  ! Calculate the moisture variance
  call calc_shoc_varorcovar(&
         shcol,nlev,nlevi,qw2tune,&               ! Input
         isotropy_zi,tkh_zi,dz_zi,qw,qw,&         ! Input
         qw_sec)                                  ! Input/Output

  ! Calculate the temperature and moisture covariance
  call calc_shoc_varorcovar(&
         shcol,nlev,nlevi,qwthl2tune,&            ! Input
         isotropy_zi,tkh_zi,dz_zi,thetal,qw,&     ! Input
         qwthl_sec)                               ! Input/Output

  ! Calculate vertical flux for heat
  call calc_shoc_vertflux(&
         shcol,nlev,nlevi,tkh_zi,dz_zi,thetal,&   ! Input
         wthl_sec)                                ! Input/Output

  ! Calculate vertical flux for moisture
  call calc_shoc_vertflux(&
         shcol,nlev,nlevi,tkh_zi,dz_zi,qw,&       ! Input
         wqw_sec)                                 ! Input/Output

  ! Calculate vertical flux for TKE
  call calc_shoc_vertflux(&
         shcol,nlev,nlevi,tkh_zi,dz_zi,tke,&      ! Input
         wtke_sec)                                ! Input/Output

  ! Calculate vertical flux for momentum (zonal wind)
  call calc_shoc_vertflux(&
         shcol,nlev,nlevi,tk_zi,dz_zi,u_wind,&    ! Input
         uw_sec)                                  ! Input/Output

  ! Calculate vertical flux for momentum (meridional wind)
  call calc_shoc_vertflux(&
         shcol,nlev,nlevi,tk_zi,dz_zi,v_wind,&    ! Input
         vw_sec)                                  ! Input/Output

  return
end subroutine diag_second_moments

subroutine calc_shoc_varorcovar(&
         shcol,nlev,nlevi,tunefac,&                ! Input
         isotropy_zi,tkh_zi,dz_zi,invar1,invar2,&  ! Input
         varorcovar)                               ! Input/Output

  ! Compute either the variance or covariance
  !  (depending on if invar1 is the same as invar2)
  !  for a given set of inputs

#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: calc_shoc_varorcovar_f
#endif

  implicit none

! INPUT VARIABLES
  ! number of SHOC columns
  integer, intent(in) :: shcol
  ! number of midpoint levels
  integer, intent(in) :: nlev
  ! number of interface levels
  integer, intent(in) :: nlevi
  ! tuning factor for (co)variance []
  real(rtype), intent(in) :: tunefac
  ! Return to isotopic timescale [s]
  real(rtype), intent(in) :: isotropy_zi(shcol,nlevi)
  ! Eddy diffusivity for heat [ms-2]
  real(rtype), intent(in) :: tkh_zi(shcol,nlevi)
  ! delta z centerend on zi grid [m]
  real(rtype), intent(in) :: dz_zi(shcol,nlevi)
  ! Input variable 1 [units vary]
  real(rtype), intent(in) :: invar1(shcol,nlev)
  ! Input variable 2 [units vary]
  real(rtype), intent(in) :: invar2(shcol,nlev)

! INPUT/OUTPUT VARIABLES
  ! variance or covariance [units vary]
  real(rtype), intent(inout) :: varorcovar(shcol,nlevi)

! LOCAL VARIABLES
  integer :: i, k, kt
  real(rtype) :: sm, grid_dz2

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
      call calc_shoc_varorcovar_f(shcol,nlev,nlevi,tunefac,isotropy_zi,tkh_zi,dz_zi,invar1,invar2,&  ! Input
           varorcovar)                              ! Input/Output)
      return
   endif
#endif
   
  do k=2,nlev

    kt=k-1 ! define upper grid point indicee
    do i=1,shcol

      grid_dz2=bfb_square(1._rtype/dz_zi(i,k)) ! vertical grid diff squared
      sm=isotropy_zi(i,k)*tkh_zi(i,k) ! coefficient for variances

      ! Compute the variance or covariance
      varorcovar(i,k)=tunefac*sm*grid_dz2*(invar1(i,kt)-invar1(i,k))*&
        (invar2(i,kt)-invar2(i,k))

    enddo
  enddo

  return
end subroutine calc_shoc_varorcovar

subroutine calc_shoc_vertflux(&
         shcol,nlev,nlevi,tkh_zi,dz_zi,invar,&  ! Input
         vertflux)                              ! Input/Output

  ! Compute either the vertical flux via
  !  downgradient diffusion for a given set of
  !  input variables

#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: calc_shoc_vertflux_f
#endif

  implicit none

! INPUT VARIABLES
  ! number of SHOC columns
  integer, intent(in) :: shcol
  ! number of midpoint levels
  integer, intent(in) :: nlev
  ! number of interface levels
  integer, intent(in) :: nlevi
  ! Eddy diffusivity for heat [ms-2]
  real(rtype), intent(in) :: tkh_zi(shcol,nlevi)
  ! delta z centerend on zi grid [m]
  real(rtype), intent(in) :: dz_zi(shcol,nlevi)
  ! Input variable [units vary]
  real(rtype), intent(in) :: invar(shcol,nlev)

! INPUT/OUTPUT VARIABLES
  real(rtype), intent(out) :: vertflux(shcol,nlevi)

! LOCAL VARIABLES
  integer :: i, k, kt
  real(rtype) :: grid_dz

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
      call calc_shoc_vertflux_f(shcol,nlev,nlevi,tkh_zi,dz_zi,invar,&  ! Input
           vertflux)                              ! Input/Output)
      return
   endif
#endif

  do k=2,nlev

    kt=k-1 ! define upper grid point indicee
    do i=1,shcol

      grid_dz=1._rtype/dz_zi(i,k) ! vertical grid diff squared

      ! Compute the vertical flux via downgradient diffusion
      vertflux(i,k)=-1._rtype*tkh_zi(i,k)*grid_dz*&
        (invar(i,kt)-invar(i,k))

    enddo
  enddo

  return
end subroutine calc_shoc_vertflux

subroutine diag_second_moments_ubycond(&
         shcol, &                               ! Input
         thl_sec, qw_sec,&                      ! Output
         wthl_sec,wqw_sec,&                     ! Output
         qwthl_sec, uw_sec, vw_sec, wtke_sec)   ! Output

  ! Purpose of this subroutine is to diagnose the upper
  !  boundary condition for the second order moments
  !  needed for the SHOC parameterization.  Currently
  !  set all to zero.

#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: shoc_diag_second_moments_ubycond_f
#endif
  implicit none

  ! INPUT VARIABLES
  ! number of SHOC columns
  integer, intent(in) :: shcol

  ! OUTPUT VARIABLES
  ! second order liquid wat. potential temp. [K^2]
  real(rtype), intent(out) :: thl_sec(shcol)
  ! second order total water mixing rat. [kg^2/kg^2]
  real(rtype), intent(out) :: qw_sec(shcol)
  ! covariance of temp and moisture [K kg/kg]
  real(rtype), intent(out) :: qwthl_sec(shcol)
  ! vertical flux of heat [K m/s]
  real(rtype), intent(out) :: wthl_sec(shcol)
  ! vertical flux of total water [kg/kg m/s]
  real(rtype), intent(out) :: wqw_sec(shcol)
  ! vertical flux of zonal wind [m2/s2]
  real(rtype), intent(out) :: uw_sec(shcol)
  ! vertical flux of meridional wind [m2/s2]
  real(rtype), intent(out) :: vw_sec(shcol)
  ! vertical flux of tke [m3/s3]
  real(rtype), intent(out) :: wtke_sec(shcol)

  ! LOCAL VARIABLES
  integer :: i

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
       call shoc_diag_second_moments_ubycond_f(&
                             shcol, &                    ! Input
                             thl_sec, qw_sec,&                      ! Output
                             wthl_sec,wqw_sec,&                     ! Output
                             qwthl_sec, uw_sec, vw_sec, wtke_sec)   ! Output
      return
   endif
#endif

  ! apply the upper boundary condition
  do i=1,shcol
    wthl_sec(i) = 0._rtype
    wqw_sec(i) = 0._rtype
    uw_sec(i) = 0._rtype
    vw_sec(i) = 0._rtype
    wtke_sec(i) = 0._rtype

    thl_sec(i) = 0._rtype
    qw_sec(i) = 0._rtype
    qwthl_sec(i) = 0._rtype
  enddo ! end i loop (column loop)
  return
end subroutine diag_second_moments_ubycond

!==============================================================
! SHOC Diagnose the third order moment of vertical velocity

subroutine diag_third_shoc_moments(&
          shcol,nlev,nlevi, &                 ! Input
          w_sec, thl_sec, qw_sec, qwthl_sec,& ! Input
          wthl_sec, isotropy, brunt,&         ! Input
          thetal,tke,wthv_sec,&               ! Input
          dz_zt, dz_zi, zt_grid, zi_grid,&    ! Input
          w3)                                 ! Output

  ! Purpose of this subroutine is to diagnose the third
  !  order moment of the vertical velocity, needed
  !  for the skewness calculation in the PDF.
  !  This calculation follows that of Canuto et al. (2001)

  implicit none

! INPUT VARIABLES
  ! number of SHOC columns
  integer, intent(in) :: shcol
  ! number of midpoint levels
  integer, intent(in) :: nlev
  ! number of interface levels
  integer, intent(in) :: nlevi

  ! second order vertical velocity [m2/s2]
  real(rtype), intent(in) :: w_sec(shcol,nlev)
  ! second order liquid wat. potential temperature [K^2]
  real(rtype), intent(in) :: thl_sec(shcol,nlevi)
  ! second order total water mixing ratio [kg2/kg2]
  real(rtype), intent(in) :: qw_sec(shcol,nlevi)
  ! covariance of temp and moisture [K kg/kg]
  real(rtype), intent(in) :: qwthl_sec(shcol,nlevi)
  ! vertical flux of heat [K m/s]
  real(rtype), intent(in) :: wthl_sec(shcol,nlevi)
  ! return to isotropy timescale [s]
  real(rtype), intent(in) :: isotropy(shcol,nlev)
  ! brunt vaisallia frequency [s]
  real(rtype), intent(in) :: brunt(shcol,nlev)
  ! liquid water potential temperature [K]
  real(rtype), intent(in) :: thetal(shcol,nlev)
  ! turbulent kinetic energy [m2/s2]
  real(rtype), intent(in) :: tke(shcol,nlev)
  ! buoyancy flux [K m/s]
  real(rtype), intent(in) :: wthv_sec(shcol,nlev)
  ! thickness centered on thermodynamic grid [m]
  real(rtype), intent(in) :: dz_zt(shcol,nlev)
  ! thickness centered on interface grid [m]
  real(rtype), intent(in) :: dz_zi(shcol,nlevi)
  ! heights of thermodynamics points [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! heights of interface points [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)

! OUTPUT VARIABLES
  ! third moment of vertical velocity
  real(rtype), intent(out) :: w3(shcol,nlevi)

! LOCAL VARIABLES
  real(rtype) :: w_sec_zi(shcol,nlevi)    ! second order vertical velocity
  real(rtype) :: isotropy_zi(shcol,nlevi)
  real(rtype) :: brunt_zi(shcol,nlevi)
  real(rtype) :: thetal_zi(shcol,nlevi)
  real(rtype) :: wthv_sec_zi(shcol,nlevi)

  ! Interpolate variables onto the interface levels
  call linear_interp(zt_grid,zi_grid,isotropy,isotropy_zi,nlev,nlevi,shcol,0._rtype)
  call linear_interp(zt_grid,zi_grid,brunt,brunt_zi,nlev,nlevi,shcol,largeneg)
  call linear_interp(zt_grid,zi_grid,w_sec,w_sec_zi,nlev,nlevi,shcol,(2._rtype/3._rtype)*mintke)
  call linear_interp(zt_grid,zi_grid,thetal,thetal_zi,nlev,nlevi,shcol,0._rtype)
  call linear_interp(zt_grid,zi_grid,wthv_sec,wthv_sec_zi,nlev,nlevi,shcol,largeneg)

  !Diagnose the third moment of the vertical-velocity
  call compute_diag_third_shoc_moment(&
          shcol,nlev,nlevi, &                 ! Input
          w_sec,thl_sec, qw_sec, qwthl_sec,&  ! Input
          wthl_sec, tke, dz_zt, dz_zi,&       ! Input
          zt_grid,zi_grid, isotropy_zi,&      ! Input
          brunt_zi,w_sec_zi,thetal_zi,&       ! Input
          wthv_sec_zi,&                       ! Input
          w3)                                 ! Output

  ! perform clipping to prevent unrealistically large values from occuring
  call clipping_diag_third_shoc_moments(&
          nlevi,shcol,w_sec_zi,&    !Input
          w3)                       !Input/Output

  return

end subroutine diag_third_shoc_moments

subroutine compute_diag_third_shoc_moment(&
          shcol,nlev,nlevi, &                 ! Input
          w_sec,thl_sec, qw_sec, qwthl_sec,&  ! Input
          wthl_sec, tke, dz_zt, dz_zi,&       ! Input
          zt_grid,zi_grid, isotropy_zi,&      ! Input
          brunt_zi,w_sec_zi,thetal_zi,&       ! Input
          wthv_sec_zi,&                       ! Input
          w3)                                 ! Output

  implicit none
! INPUT VARIABLES
  ! number of SHOC columns
  integer, intent(in) :: shcol
  ! number of midpoint levels
  integer, intent(in) :: nlev
  ! number of interface levels
  integer, intent(in) :: nlevi
  ! second order vertical velocity [m2/s2]
  real(rtype), intent(in) :: w_sec(shcol,nlev)
  ! second order liquid wat. potential temperature [K^2]
  real(rtype), intent(in) :: thl_sec(shcol,nlevi)
  ! second order total water mixing ratio [kg2/kg2]
  real(rtype), intent(in) :: qw_sec(shcol,nlevi)
  ! covariance of temp and moisture [K kg/kg]
  real(rtype), intent(in) :: qwthl_sec(shcol,nlevi)
  ! vertical flux of heat [K m/s]
  real(rtype), intent(in) :: wthl_sec(shcol,nlevi)
  ! turbulent kinetic energy [m2/s2]
  real(rtype), intent(in) :: tke(shcol,nlev)
  ! thickness centered on thermodynamic grid [m]
  real(rtype), intent(in) :: dz_zt(shcol,nlev)
  ! thickness centered on interface grid [m]
  real(rtype), intent(in) :: dz_zi(shcol,nlevi)
  ! heights of thermodynamics points [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! heights of interface points [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)

  !Interpolated varaibles
  real(rtype), intent(in) :: isotropy_zi(shcol,nlevi)
  real(rtype), intent(in) :: brunt_zi(shcol,nlevi)
  real(rtype), intent(in) :: w_sec_zi(shcol,nlevi)
  real(rtype), intent(in) :: thetal_zi(shcol,nlevi)
  real(rtype), intent(in) :: wthv_sec_zi(shcol,nlevi)
  ! third moment of vertical velocity
  real(rtype), intent(out) :: w3(shcol,nlevi)
! LOCAL VARIABLES
  integer i, j, k, kb, kc
  real(rtype) :: omega0, omega1, omega2
  real(rtype) :: X0, Y0, X1, Y1, AA0, AA1
  real(rtype) :: zvar, x5var, iso, thedz, thedz2
  real(rtype) :: theterm, tsign
  real(rtype) :: isosqrd
  real(rtype) :: buoy_sgs2, bet2
  real(rtype) :: f0, f1, f2, f3, f4, f5

  ! set lower condition
  w3(:,nlevi) = 0._rtype

  do k=2,nlev

     kb=k+1
     kc=k-1
     do i=1,shcol

        !Compute inputs for computing f0 to f5 terms
        call fterms_input_for_diag_third_shoc_moment(&
	     dz_zi(i,k), dz_zt(i,k), dz_zt(i,kc), &               ! Input
             isotropy_zi(i,k), brunt_zi(i,k), thetal_zi(i,k), &   ! Input
             thedz, thedz2, iso, isosqrd, buoy_sgs2, bet2)        ! Output

        !Compute f0 to f5 terms
        call f0_to_f5_diag_third_shoc_moment(&
	     thedz, thedz2, bet2, iso, isosqrd, &                 ! Input
	     wthl_sec (i,k), wthl_sec(i,kc), wthl_sec(i,kb), &    ! Input
	     thl_sec(i,k), thl_sec(i,kc), thl_sec(i,kb), &        ! Input
             w_sec(i,k), w_sec(i,kc), w_sec_zi(i,k), &            ! Input
             tke(i,k), tke(i,kc), &                               ! Input
             f0, f1, f2, f3, f4, f5)                              ! Output

        !Compute the omega terms
        call omega_terms_diag_third_shoc_moment(&
	     buoy_sgs2, f3, f4, &       ! Input
	     omega0, omega1, omega2)    ! Output

        !Compute the X0, Y0, X1, Y1 terms
        call x_y_terms_diag_third_shoc_moment(&
	     buoy_sgs2, f0, f1, f2, &   ! Input
	     x0, y0, x1, y1)            ! Output

        !Compute the AA0, AA1 terms
        call aa_terms_diag_third_shoc_moment(&
	     omega0, omega1, omega2, &  ! Input
	     x0, x1, y0, y1, &          ! Input
	     aa0, aa1)                  ! Output

        !Finally, we have the third moment of w
        w3(i,k) = w3_diag_third_shoc_moment(aa0, aa1, x0, x1, f5)

     enddo  ! end i loop (column loop)
  enddo  ! end k loop (vertical loop)

  ! set upper condition
  w3(:,1) = 0._rtype

end subroutine compute_diag_third_shoc_moment

subroutine fterms_input_for_diag_third_shoc_moment(&
     dz_zi, dz_zt, dz_zt_kc, &                      ! Input
     isotropy_zi, brunt_zi, thetal_zi, &            ! Input
     thedz, thedz2, iso, isosqrd, buoy_sgs2, bet2)  ! Output

  !Compute inputs for computing f0 to f5 terms

  implicit none

  !intent-ins
  real(rtype), intent(in) :: dz_zi, dz_zt, dz_zt_kc
  real(rtype), intent(in) :: isotropy_zi, brunt_zi, thetal_zi

  !intent-outs
  real(rtype), intent(out) :: thedz, thedz2, iso, isosqrd
  real(rtype), intent(out) :: buoy_sgs2, bet2

  thedz  = 1._rtype/dz_zi
  thedz2 = 1._rtype/(dz_zt+dz_zt_kc)

  iso       = isotropy_zi
  isosqrd   = iso**2
  buoy_sgs2 = isosqrd*brunt_zi
  bet2      = ggr/thetal_zi

  return
end subroutine fterms_input_for_diag_third_shoc_moment

subroutine f0_to_f5_diag_third_shoc_moment(&
     thedz, thedz2, bet2, iso, isosqrd, &    ! Input
     wthl_sec, wthl_sec_kc, wthl_sec_kb, &   ! Input
     thl_sec, thl_sec_kc, thl_sec_kb, &      ! Input
     w_sec, w_sec_kc,w_sec_zi, &             ! Input
     tke, tke_kc, &                          ! Input
     f0, f1, f2, f3, f4, f5)                 ! Output

  !Compute f0 to f5 terms

  implicit none

  !intent-ins
  real(rtype), intent(in) :: thedz, thedz2, bet2, iso, isosqrd
  real(rtype), intent(in) :: wthl_sec, wthl_sec_kc, wthl_sec_kb
  real(rtype), intent(in) :: thl_sec, thl_sec_kc, thl_sec_kb
  real(rtype), intent(in) :: w_sec, w_sec_kc, w_sec_zi, tke, tke_kc

  !intent-out
  real(rtype), intent(out) :: f0, f1, f2, f3, f4, f5

  !local variables
  real(rtype) :: thl_sec_diff, wthl_sec_diff, wsec_diff, tke_diff

  !Some common factors
  thl_sec_diff  = thl_sec_kc  - thl_sec_kb
  wthl_sec_diff = wthl_sec_kc - wthl_sec_kb
  wsec_diff     = w_sec_kc    - w_sec
  tke_diff      = tke_kc      - tke

  f0 = thedz2 * bet2**3 * iso**4 * wthl_sec * &
       thl_sec_diff

  f1 = thedz2 * bet2**2 * iso**3 * (wthl_sec * &
       wthl_sec_diff + 0.5_rtype * &
       w_sec_zi*thl_sec_diff)

  f2 = thedz * bet2 * isosqrd * wthl_sec * &
       wsec_diff+ 2._rtype * thedz2 * bet2 * &
       isosqrd * w_sec_zi * wthl_sec_diff

  f3 = thedz2 * bet2 * isosqrd * w_sec_zi * &
       wthl_sec_diff + thedz * &
       bet2 * isosqrd * (wthl_sec * tke_diff)

  f4 = thedz * iso * w_sec_zi * (wsec_diff + &
       tke_diff)

  f5 = thedz * iso * w_sec_zi * wsec_diff

  return
end subroutine f0_to_f5_diag_third_shoc_moment

subroutine omega_terms_diag_third_shoc_moment(&
           buoy_sgs2, f3, f4, &    ! Input
	   omega0, omega1, omega2) ! Output

  implicit none

  !Compute the omega terms

  !initent-ins
  real(rtype), intent(in) :: buoy_sgs2, f3, f4

  !intent-out
  real(rtype), intent(out) :: omega0, omega1, omega2

  real(rtype), parameter :: a4=2.4_rtype/(3._rtype*c_diag_3rd_mom+5._rtype)
  real(rtype), parameter :: a5=0.6_rtype/(c_diag_3rd_mom*(3._rtype+5._rtype*c_diag_3rd_mom))

  omega0 = a4 / (1._rtype - a5 * buoy_sgs2)
  omega1 = omega0/(2._rtype * c_diag_3rd_mom)
  omega2 = omega1 * f3 + (5._rtype/4._rtype) * omega0 * f4

  return
end subroutine omega_terms_diag_third_shoc_moment

subroutine x_y_terms_diag_third_shoc_moment(&
           buoy_sgs2, f0, f1, f2,&  ! Input
	   x0, y0, x1, y1)          ! Output

  implicit none

  !Compute the X0, Y0, X1, Y1 terms

  !intent-ins
  real(rtype), intent(in) :: buoy_sgs2, f0, f1, f2

  !intent-outs
  real(rtype), intent(out) :: x0, y0, x1, y1

  real(rtype), parameter :: a0=(0.52_rtype*c_diag_3rd_mom**(-2))/(c_diag_3rd_mom-2._rtype)
  real(rtype), parameter :: a1=0.87_rtype/(c_diag_3rd_mom**2)
  real(rtype), parameter :: a2=0.5_rtype/c_diag_3rd_mom
  real(rtype), parameter :: a3=0.6_rtype/(c_diag_3rd_mom*(c_diag_3rd_mom-2._rtype))

  x0 = (a2 * buoy_sgs2 * (1._rtype - a3 * buoy_sgs2)) / &
       (1._rtype - (a1 + a3) * buoy_sgs2)
  y0 = (2._rtype * a2 * buoy_sgs2 * x0) / (1._rtype - a3 * buoy_sgs2)
  x1 = (a0 * f0 + a1 * f1 + a2 * (1._rtype - a3 * buoy_sgs2) * f2) / &
       (1._rtype - (a1 + a3) * buoy_sgs2)
  y1 = (2._rtype * a2 * (buoy_sgs2 * x1 + (a0/a1) * f0 + f1)) / &
       (1._rtype - a3* buoy_sgs2)

  return
end subroutine x_y_terms_diag_third_shoc_moment

subroutine aa_terms_diag_third_shoc_moment(&
           omega0, omega1, omega2, & ! Input
	   x0, x1, y0, y1, &         ! Input
	   aa0, aa1)                 ! Output

  implicit none

  !Compute the AA0, AA1 terms

  !intent-ins
  real(rtype), intent(in) :: omega0, omega1, omega2, x0, x1, y0, y1

  !intent-outs
  real(rtype), intent(out) :: aa0, aa1

  aa0 = omega0 * x0 + omega1 * y0
  aa1 = omega0 * x1 + omega1 * y1 + omega2

  return
end subroutine aa_terms_diag_third_shoc_moment

pure function w3_diag_third_shoc_moment(aa0, aa1, x0, x1, f5) result(w3)

  implicit none

  !Compute third moment of w

  !intent-ins
  real(rtype), intent(in) :: aa0, aa1, x0, x1, f5

  !return type
  real(rtype) :: w3

  w3 = (aa1-1.2_rtype*x1-1.5_rtype*f5)/(c_diag_3rd_mom-1.2_rtype*x0+aa0)

  return
end function w3_diag_third_shoc_moment

subroutine clipping_diag_third_shoc_moments(&
           nlevi,shcol,w_sec_zi,& ! Input
	   w3)                    ! Input/Output

  ! perform clipping to prevent unrealistically large values from occuring

  implicit none

  integer, intent(in) :: nlevi
  integer, intent(in) :: shcol

  real(rtype), intent(in) :: w_sec_zi(shcol,nlevi)
  real(rtype), intent(inout) :: w3(shcol,nlevi)

  real(rtype) :: tsign
  real(rtype) :: cond
  real(rtype) :: theterm

  integer k, i

  do k=1, nlevi
    do i=1, shcol

      tsign = 1._rtype
      theterm = w_sec_zi(i,k)
      cond = w3clip * sqrt(2._rtype * theterm**3)
      if (w3(i,k) .lt. 0) tsign = -1._rtype
      if (tsign * w3(i,k) .gt. cond) w3(i,k) = tsign * cond

    enddo !end i loop (column loop)
  enddo ! end k loop (vertical loop)

end subroutine clipping_diag_third_shoc_moments

!==============================================================
! Assumed PDF closure for the SHOC scheme

subroutine shoc_assumed_pdf(&
         shcol,nlev,nlevi, &                ! Input
         thetal,qw,w_field,thl_sec,qw_sec,& ! Input
         wthl_sec,w_sec, &                  ! Input
         wqw_sec,qwthl_sec,w3,pres, &       ! Input
         zt_grid,zi_grid,&                  ! Input
         shoc_cldfrac,shoc_ql,&             ! Output
         wqls,wthv_sec,shoc_ql2)            ! Output

  ! Purpose of this subroutine is calculate the
  !  double Gaussian PDF of SHOC, which is the centerpiece
  !  of the scheme.  The main outputs are the SGS cloud
  !  fraction and liquid water amount, in addition to the
  !  SGS buoyancy flux which is needed to close the SGS
  !  TKE equation.  This code follows the appendix of
  !  Larson et al. (2002) for Analytic Double Gaussian 1

  implicit none

! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  ! number of midpoint layers
  integer, intent(in) :: nlev
  ! number of interface layers
  integer, intent(in) :: nlevi

  ! liquid water potential temperature [K]
  real(rtype), intent(in) :: thetal(shcol,nlev)
  ! total water mixing ratio [kg/kg]
  real(rtype), intent(in) :: qw(shcol,nlev)
  ! thetal variance [K^2]
  real(rtype), intent(in) :: thl_sec(shcol,nlevi)
  ! qw variance [kg/kg^2]
  real(rtype), intent(in) :: qw_sec(shcol,nlevi)
  ! vertical flux of heat [K m/s]
  real(rtype), intent(in) :: wthl_sec(shcol,nlevi)
  ! vertical velocity variance [m2/s2]
  real(rtype), intent(in) :: w_sec(shcol,nlev)
  ! vertical flux of moisture [kg/kg m/s]
  real(rtype), intent(in) :: wqw_sec(shcol,nlevi)
  ! qw and thetal correlation [K kg/kg]
  real(rtype), intent(in) :: qwthl_sec(shcol,nlevi)
  ! third moment vertical velocity [m^3/s^3]
  real(rtype), intent(in) :: w3(shcol,nlevi)
  ! large scale vertical velocity [m/s]
  real(rtype), intent(in) :: w_field(shcol,nlev)
  ! pressure [Pa]
  real(rtype), intent(in) :: pres(shcol,nlev)
  ! heights on midpoint grid [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! heights on interface grid [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)

! OUTPUT VARIABLES
  ! SGS cloud fraction [-]
  real(rtype), intent(out) :: shoc_cldfrac(shcol,nlev)
  ! SGS liquid water mixing ratio [kg/kg]
  real(rtype), intent(out) :: shoc_ql(shcol,nlev)
  ! SGS buoyancy flux [K m/s]
  real(rtype), intent(out) :: wthv_sec(shcol,nlev)
  ! SGS liquid water flux [kg/kg m/s]
  real(rtype), intent(out) :: wqls(shcol,nlev)
  ! SGS liquid water mixing ratio variance [kg/kg]
  real(rtype), intent(out) :: shoc_ql2(shcol,nlev)

! LOCAL VARIABLES
  integer i,k
  real(rtype) skew_w,a
  real(rtype) w1_1,w1_2,w2_1,w2_2,w3var
  real(rtype) thl1_1,thl1_2,thl2_1,thl2_2
  real(rtype) qw1_1,qw1_2,qw2_1,qw2_2
  real(rtype) r_qwthl_1
  real(rtype) s1,s2,std_s1,std_s2,C1,C2
  real(rtype) ql1,ql2
  real(rtype) thl_first,qw_first,w_first
  real(rtype) Tl1_1,Tl1_2,pval
  real(rtype) thlsec,qwsec,qwthlsec,wqwsec,wthlsec
  real(rtype) qn1,qn2
  real(rtype) beta1, beta2, qs1, qs2
  real(rtype) sqrtw2, sqrtthl, sqrtqt
  real(rtype) epsterm
  real(rtype) sqrtqw2_1, sqrtqw2_2, sqrtthl2_1, sqrtthl2_2
  real(rtype) thl_tol, rt_tol, w_tol_sqd, w_thresh
  character(len=200) :: err_msg

  ! variables on thermo grid
  real(rtype) :: wthl_sec_zt(shcol,nlev)
  real(rtype) :: wqw_sec_zt(shcol,nlev)
  real(rtype) :: w3_zt(shcol,nlev)
  real(rtype) :: thl_sec_zt(shcol,nlev)
  real(rtype) :: qwthl_sec_zt(shcol,nlev)
  real(rtype) :: qw_sec_zt(shcol,nlev)

  ! define these so they don't have to be computed more than once
  real(rtype), parameter :: sqrt2 = sqrt(2._rtype)
  real(rtype), parameter :: sqrtpi = sqrt(2._rtype*3.14_rtype)

  epsterm=rgas/rv

  thl_tol=1.e-2_rtype
  rt_tol=1.e-4_rtype
  w_tol_sqd=(2.e-2_rtype)**2
  w_thresh=0.0_rtype

  ! Initialize cloud variables to zero
  shoc_cldfrac(:,:)=0._rtype
  shoc_ql(:,1)=0._rtype
  shoc_ql2(:,:) = 0._rtype

  ! Interpolate many variables from interface grid to themo grid
  call linear_interp(zi_grid,zt_grid,w3,w3_zt,nlevi,nlev,shcol,largeneg)
  call linear_interp(zi_grid,zt_grid,thl_sec,thl_sec_zt,nlevi,nlev,shcol,0._rtype)
  call linear_interp(zi_grid,zt_grid,wthl_sec,wthl_sec_zt,nlevi,nlev,shcol,largeneg)
  call linear_interp(zi_grid,zt_grid,qwthl_sec,qwthl_sec_zt,nlevi,nlev,shcol,largeneg)
  call linear_interp(zi_grid,zt_grid,wqw_sec,wqw_sec_zt,nlevi,nlev,shcol,largeneg)
  call linear_interp(zi_grid,zt_grid,qw_sec,qw_sec_zt,nlevi,nlev,shcol,0._rtype)

  do k=1,nlev
    do i=1,shcol

      pval = pres(i,k)

      ! Get all needed input moments for the PDF
      !  at this particular point
      thl_first = thetal(i,k)
      w_first = w_field(i,k)
      qw_first = qw(i,k)

      w3var = w3_zt(i,k)
      thlsec = thl_sec_zt(i,k)
      qwsec = qw_sec_zt(i,k)
      qwthlsec = qwthl_sec_zt(i,k)
      wqwsec = wqw_sec_zt(i,k)
      wthlsec = wthl_sec_zt(i,k)

      ! Compute square roots of some variables so we don't
      !  have to compute these again
      sqrtw2 = sqrt(w_sec(i,k))
      sqrtthl = max(thl_tol,sqrt(thlsec))
      sqrtqt = max(rt_tol,sqrt(qwsec))

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  FIND PARAMETERS FOR VERTICAL VELOCITY

      call shoc_assumed_pdf_vv_parameters(&
         w_first,w_sec(i,k),w3var,&    ! Input
         Skew_w,w1_1,w1_2,w2_1,w2_2,a) ! Output

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  FIND PARAMETERS FOR THETAL

      call shoc_assumed_pdf_thl_parameters(&
         wthlsec,sqrtw2,sqrtthl,thlsec,thl_first,& ! Input
         w1_1,w1_2,Skew_w,a,dothetal_skew,&        ! Input
         thl1_1,thl1_2,thl2_1,thl2_2,sqrtthl2_1,&  ! Output
         sqrtthl2_2)                               ! Output

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  FIND PARAMETERS FOR TOTAL WATER MIXING RATIO

      call shoc_assumed_pdf_qw_parameters(&
         wqwsec,sqrtw2,Skew_w,sqrtqt,& ! Input
         qwsec,w1_2,w1_1,qw_first,a,&  ! Input
         qw1_1,qw1_2,qw2_1,&           ! Output
         qw2_2,sqrtqw2_1,sqrtqw2_2)    ! Output

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  CONVERT FROM TILDA VARIABLES TO "REAL" VARIABLES

      call shoc_assumed_pdf_tilda_to_real(&
         w_first,sqrtw2,& ! Input
         w1_1)            ! Output
      call shoc_assumed_pdf_tilda_to_real(&
         w_first,sqrtw2,& ! Input
         w1_2)            ! Output

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  FIND WITHIN-PLUME CORRELATIONS

      call shoc_assumed_pdf_inplume_correlations(&
        sqrtqw2_1,sqrtthl2_1,a,sqrtqw2_2,sqrtthl2_2,&           ! Input
        qwthlsec,qw1_1,qw_first,thl1_1,thl_first,qw1_2,thl1_2,& ! Input
        r_qwthl_1)                                              ! Output

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  BEGIN TO COMPUTE CLOUD PROPERTY STATISTICS

      call shoc_assumed_pdf_compute_temperature(&
        thl1_1,basepres,pval,& ! Input
        Tl1_1)                 ! Output
      call shoc_assumed_pdf_compute_temperature(&
        thl1_2,basepres,pval,& ! Input
        Tl1_2)                 ! Output

      ! Check to ensure Tl1_1 and Tl1_2 are not negative. endrun otherwise
      if (Tl1_1 .le. 0._rtype) then
         write(err_msg,*)'ERROR: Tl1_1 is .le. 0 before shoc_assumed_pdf_compute_qs in shoc. Tl1_1 is:',Tl1_1
         call endscreamrun(err_msg)
      endif

      if (Tl1_2 .le. 0._rtype) then
         write(err_msg,*)'ERROR: Tl1_2 is .le. 0 before shoc_assumed_pdf_compute_qs in shoc. Tl1_2 is:',Tl1_2
         call endscreamrun(err_msg)
      endif

      ! Now compute qs
      call shoc_assumed_pdf_compute_qs(&
        Tl1_1,Tl1_2,pval,&   ! Input
        qs1,beta1,qs2,beta2) ! Output

      !!!!!  Now compute cloud stuff
      !!!!!!  compute s term
      call shoc_assumed_pdf_compute_s(&
        qw1_1,qs1,beta1,pval,thl2_1,&          ! Input
        qw2_1,sqrtthl2_1,sqrtqw2_1,r_qwthl_1,& ! Input
        s1,std_s1,qn1,C1)                      ! Output

      !!!!! now compute non-precipitating cloud condensate

      ! If two plumes exactly equal, then just set many of these
      ! variables to themselves to save on computation.
      if (qw1_1 .eq. qw1_2 .and. thl2_1 .eq. thl2_2 .and. qs1 .eq. qs2) then
        s2=s1
        std_s2=std_s1
        C2=C1
        qn2=qn1
      else
        call shoc_assumed_pdf_compute_s(&
        qw1_2,qs2,beta2,pval,thl2_2,&          ! Input
        qw2_2,sqrtthl2_2,sqrtqw2_2,r_qwthl_1,& ! Input
        s2,std_s2,qn2,C2)                      ! Output
      endif

      ql1=min(qn1,qw1_1)
      ql2=min(qn2,qw1_2)

      ! Finally, compute SGS cloud fraction
      shoc_cldfrac(i,k) = min(1._rtype,a*C1+(1._rtype-a)*C2)

      ! Compute SGS liquid water mixing ratio
      call shoc_assumed_pdf_compute_sgs_liquid(&
        a,ql1,ql2,&   ! Input
        shoc_ql(i,k)) ! Output

      ! Compute cloud liquid variance (CLUBB formulation, adjusted to SHOC parameters based)
      call shoc_assumed_pdf_compute_cloud_liquid_variance(&
        a,s1,ql1,C1,std_s1,s2,ql2,C2,std_s2,shoc_ql(i,k),& ! Input
        shoc_ql2(i,k))                                     ! Output

      ! Compute liquid water flux
      call shoc_assumed_pdf_compute_liquid_water_flux(&
        a,w1_1,w_first,ql1,w1_2,ql2,& ! Input
        wqls(i,k))                    ! Output

      ! Compute the SGS buoyancy flux
      call shoc_assumed_pdf_compute_buoyancy_flux(&
        wthlsec, epsterm, wqwsec, pval, wqls(i,k),& ! Input
        wthv_sec(i,k))                              ! Output

    enddo  ! end i loop here
  enddo ! end k loop here

  return

end subroutine shoc_assumed_pdf

subroutine shoc_assumed_pdf_vv_parameters(&
   w_first,w_sec,w3var,&          ! Input
   Skew_w,w1_1,w1_2,w2_1,w2_2,a)  ! Output
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  FIND PARAMETERS FOR VERTICAL VELOCITY

  implicit none

  ! intent-ins
  real(rtype), intent(in) :: w_first
  real(rtype), intent(in) :: w_sec
  real(rtype), intent(in) :: w3var

  ! intent-out
  real(rtype), intent(out) :: Skew_w
  real(rtype), intent(out) :: w1_1
  real(rtype), intent(out) :: w1_2
  real(rtype), intent(out) :: w2_1
  real(rtype), intent(out) :: w2_2
  real(rtype), intent(out) :: a

  ! local vars
  real(rtype) :: sqrtw2t

  ! parameters
  real(rtype), parameter :: w_tol_sqd=(2.e-2_rtype)**2

  Skew_w=w3var/w_sec**(3./2.)

  if (w_sec .le. w_tol_sqd) then
    Skew_w=0._rtype
    w1_1=w_first
    w1_2=w_first
    w2_1=0._rtype
    w2_2=0._rtype
    a=0.5_rtype
  else

    w2_1=0.4_rtype
    w2_2=0.4_rtype

    a=max(0.01_rtype,min(0.5_rtype*(1._rtype-Skew_w*sqrt(1._rtype/(4._rtype*(1._rtype-w2_1)**3+Skew_w**2))),0.99_rtype))

    sqrtw2t=sqrt(1._rtype-w2_1)

    w1_1=sqrt((1._rtype-a)/a)*sqrtw2t
    w1_2=-1._rtype*sqrt(a/(1._rtype-a))*sqrtw2t

    w2_1=w2_1*w_sec
    w2_2=w2_2*w_sec

  endif


end subroutine shoc_assumed_pdf_vv_parameters

subroutine shoc_assumed_pdf_thl_parameters(&
  wthlsec,sqrtw2,sqrtthl,thlsec,thl_first,& ! Input
  w1_1,w1_2,Skew_w,a,dothetal_skew,&        ! Input
  thl1_1,thl1_2,thl2_1,thl2_2,sqrtthl2_1,&  ! Output
  sqrtthl2_2)                               ! Output

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  FIND PARAMETERS FOR TOTAL WATER MIXING RATIO
  implicit none

  ! intent-ins
  real(rtype), intent(in) :: wthlsec
  real(rtype), intent(in) :: sqrtw2
  real(rtype), intent(in) :: sqrtthl
  real(rtype), intent(in) :: thlsec
  real(rtype), intent(in) :: thl_first
  real(rtype), intent(in) :: w1_1
  real(rtype), intent(in) :: w1_2
  real(rtype), intent(in) :: Skew_w
  real(rtype), intent(in) ::  a
  logical(btype), intent(in)  :: dothetal_skew

  ! intent-outs
  real(rtype), intent(out) :: thl1_1
  real(rtype), intent(out) :: thl1_2
  real(rtype), intent(out) :: thl2_1
  real(rtype), intent(out) :: thl2_2
  real(rtype), intent(out) :: sqrtthl2_1
  real(rtype), intent(out) :: sqrtthl2_2

  ! local vars
  real(rtype) :: corrtest1, tsign, Skew_thl

  ! parameters
  real(rtype), parameter :: thl_tol = 1.e-2_rtype
  real(rtype), parameter :: w_thresh = 0.0_rtype

  corrtest1=max(-1._rtype,min(1._rtype,wthlsec/(sqrtw2*sqrtthl)))

  if (thlsec .le. thl_tol**2 .or. abs(w1_2-w1_1) .le. w_thresh) then
    thl1_1=thl_first
    thl1_2=thl_first
    thl2_1=0._rtype
    thl2_2=0._rtype
    sqrtthl2_1=0._rtype
    sqrtthl2_2=0._rtype
  else

    thl1_1=(-1._rtype*corrtest1)/w1_2
    thl1_2=(-1._rtype*corrtest1)/w1_1

    if (dothetal_skew) then
      tsign=abs(thl1_2-thl1_1)

      if (tsign .gt. 0.4_rtype) then
        Skew_thl=1.2_rtype*Skew_w
      else if (tsign .le. 0.2_rtype) then
        Skew_thl=0.0_rtype
      else
        Skew_thl=((1.2_rtype*Skew_w)/0.2_rtype)*(tsign-0.2_rtype)
      endif
    else
      Skew_thl = 0.0_rtype
    endif

    thl2_1=min(100._rtype,max(0._rtype,(3._rtype*thl1_2*(1._rtype-a*thl1_1**2-(1._rtype-a)*thl1_2**2) &
            -(Skew_thl-a*thl1_1**3-(1._rtype-a)*thl1_2**3))/ &
            (3._rtype*a*(thl1_2-thl1_1))))*thlsec

    thl2_2=min(100._rtype,max(0._rtype,(-3._rtype*thl1_1*(1._rtype-a*thl1_1**2-(1._rtype-a)*thl1_2**2) &
      +(Skew_thl-a*thl1_1**3-(1._rtype-a)*thl1_2**3))/ &
      (3._rtype*(1._rtype-a)*(thl1_2-thl1_1))))*thlsec


    thl1_1=thl1_1*sqrtthl+thl_first
    thl1_2=thl1_2*sqrtthl+thl_first

    sqrtthl2_1=sqrt(thl2_1)
    sqrtthl2_2=sqrt(thl2_2)

  endif

end subroutine shoc_assumed_pdf_thl_parameters

subroutine shoc_assumed_pdf_qw_parameters(&
  wqwsec, sqrtw2, Skew_w, sqrtqt, qwsec, w1_2, w1_1, qw_first, a, &
  qw1_1, qw1_2, qw2_1, qw2_2, sqrtqw2_1, sqrtqw2_2)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  FIND PARAMETERS FOR TOTAL WATER MIXING RATIO

  implicit none

  !intent-in
  real(rtype), intent(in) :: wqwsec
  real(rtype), intent(in) :: sqrtw2
  real(rtype), intent(in) :: Skew_w
  real(rtype), intent(in) :: sqrtqt
  real(rtype), intent(in) :: qwsec
  real(rtype), intent(in) :: w1_2
  real(rtype), intent(in) :: w1_1
  real(rtype), intent(in) :: qw_first
  real(rtype), intent(in) :: a

  ! intent-out
  real(rtype), intent(out) :: qw1_1
  real(rtype), intent(out) :: qw1_2
  real(rtype), intent(out) :: qw2_1
  real(rtype), intent(out) :: qw2_2
  real(rtype), intent(out) :: sqrtqw2_1
  real(rtype), intent(out) :: sqrtqw2_2

  ! local vars
  real(rtype) :: corrtest2, tsign, Skew_qw

  ! Parameters
  real(rtype), parameter :: rt_tol=1.e-4_rtype
  real(rtype), parameter :: w_thresh=0.0_rtype

  corrtest2=max(-1.0_rtype,min(1.0_rtype,wqwsec/(sqrtw2*sqrtqt)))


  if (qwsec .le. rt_tol**2 .or. abs(w1_2-w1_1) .le. w_thresh) then
    qw1_1=qw_first
    qw1_2=qw_first
    qw2_1=0._rtype
    qw2_2=0._rtype
    sqrtqw2_1=0._rtype
    sqrtqw2_2=0._rtype
  else

    qw1_1=(-1._rtype*corrtest2)/w1_2
    qw1_2=(-1._rtype*corrtest2)/w1_1

    tsign=abs(qw1_2-qw1_1)

    if (tsign .gt. 0.4_rtype) then
      Skew_qw=1.2_rtype*Skew_w
    else if (tsign .le. 0.2_rtype) then
      Skew_qw=0._rtype
    else
      Skew_qw=((1.2_rtype*Skew_w)/0.2_rtype)*(tsign-0.2_rtype)
    endif
     qw2_1=min(100._rtype,max(0._rtype,(3._rtype*qw1_2*(1._rtype-a*qw1_1**2-(1._rtype-a)*qw1_2**2) &
      -(Skew_qw-a*qw1_1**3-(1._rtype-a)*qw1_2**3))/ &
      (3._rtype*a*(qw1_2-qw1_1))))*qwsec

    qw2_2=min(100._rtype,max(0._rtype,(-3._rtype*qw1_1*(1._rtype-a*qw1_1**2-(1._rtype-a)*qw1_2**2) &
      +(Skew_qw-a*qw1_1**3-(1._rtype-a)*qw1_2**3))/ &
      (3._rtype*(1._rtype-a)*(qw1_2-qw1_1))))*qwsec

    qw1_1=qw1_1*sqrtqt+qw_first
    qw1_2=qw1_2*sqrtqt+qw_first

    sqrtqw2_1=sqrt(qw2_1)
    sqrtqw2_2=sqrt(qw2_2)

  endif

end subroutine shoc_assumed_pdf_qw_parameters

subroutine shoc_assumed_pdf_tilda_to_real(&
  w_first,sqrtw2,& ! intent-in
  w1)              ! intent-out

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  CONVERT FROM TILDA VARIABLES TO "REAL" VARIABLES
  implicit none

  ! intent-ins
  real(rtype), intent(in) :: w_first
  real(rtype), intent(in) :: sqrtw2

  !intent-inouts
  real(rtype), intent(inout) :: w1

  w1 = w1 * sqrtw2 + w_first

end subroutine shoc_assumed_pdf_tilda_to_real

subroutine shoc_assumed_pdf_inplume_correlations(&
  sqrtqw2_1, sqrtthl2_1, a, sqrtqw2_2, sqrtthl2_2,&              ! Input
  qwthlsec, qw1_1, qw_first, thl1_1, thl_first, qw1_2, thl1_2,&  ! Input
  r_qwthl_1)                                                     ! Output

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  FIND WITHIN-PLUME CORRELATIONS
  implicit none

  ! intent out
  real(rtype), intent(in)  :: sqrtqw2_1
  real(rtype), intent(in)  :: sqrtthl2_1
  real(rtype), intent(in)  :: a
  real(rtype), intent(in)  :: sqrtqw2_2
  real(rtype), intent(in)  :: sqrtthl2_2
  real(rtype), intent(in)  :: qwthlsec
  real(rtype), intent(in)  :: qw1_1
  real(rtype), intent(in)  :: qw_first
  real(rtype), intent(in)  :: thl1_1
  real(rtype), intent(in)  :: thl_first
  real(rtype), intent(in)  :: qw1_2
  real(rtype), intent(in)  :: thl1_2

  ! intent out
  real(rtype), intent(out) :: r_qwthl_1

  real(rtype) :: testvar
  testvar=(a*sqrtqw2_1*sqrtthl2_1+(1._rtype-a)*sqrtqw2_2*sqrtthl2_2)

  if (testvar .eq. 0._rtype) then
    r_qwthl_1=0._rtype
  else
    r_qwthl_1=max(-1.0_rtype,min(1.0_rtype,(qwthlsec-a*(qw1_1-qw_first) &
      *(thl1_1-thl_first)-(1._rtype-a)*(qw1_2-qw_first) &
      *(thl1_2-thl_first))/testvar))
  endif

end subroutine shoc_assumed_pdf_inplume_correlations


subroutine shoc_assumed_pdf_compute_temperature(&
  thl1,basepres,pval,& ! Input
  Tl1)                 ! Output

  implicit none
  ! intent-in
  real(rtype), intent(in)  :: thl1
  real(rtype), intent(in)  ::basepres
  real(rtype), intent(in)  ::pval

  ! intent-out
  real(rtype), intent(out) :: Tl1

  TL1 = thl1/((basepres/pval)**(rgas/cp))

end subroutine shoc_assumed_pdf_compute_temperature

subroutine shoc_assumed_pdf_compute_qs(&
  Tl1_1,Tl1_2,pval,&   ! Input
  qs1,beta1,qs2,beta2) ! Ouput

  use wv_sat_scream, only: MurphyKoop_svp
  implicit none

  ! intent-in
  real(rtype), intent(in) :: Tl1_1
  real(rtype), intent(in) :: Tl1_2
  real(rtype), intent(in) :: pval

  ! intent-out
  real(rtype), intent(out) :: qs1
  real(rtype), intent(out) ::   beta1
  real(rtype), intent(out) ::   qs2
  real(rtype), intent(out) ::   beta2

  ! local vars
  integer, parameter :: liquid = 0   ! liquid flag for MurphyKoop
  real(rtype) :: esval1_1
  real(rtype) :: esval1_2
  real(rtype) :: lstarn1
  real(rtype) :: lstarn2

  esval1_1=0._rtype
  esval1_2=0._rtype

  esval1_1=MurphyKoop_svp(Tl1_1, liquid)
  lstarn1=lcond

  qs1=0.622_rtype*esval1_1/max(esval1_1,pval-esval1_1)
  beta1=(rgas/rv)*(lstarn1/(rgas*Tl1_1))*(lstarn1/(cp*Tl1_1))

  ! Are the two plumes equal?  If so then set qs and beta
  ! in each column to each other to save computation
  lstarn2=lcond
  if (Tl1_1 .eq. Tl1_2) then
    qs2=qs1
    beta2=beta1
  else
    esval1_2=MurphyKoop_svp(Tl1_2, liquid)
    qs2=0.622_rtype*esval1_2/max(esval1_2,pval-esval1_2)
    beta2=(rgas/rv)*(lstarn2/(rgas*Tl1_2))*(lstarn2/(cp*Tl1_2))
  endif

end subroutine shoc_assumed_pdf_compute_qs

subroutine shoc_assumed_pdf_compute_s(&
  qw1,qs1,beta,pval,thl2,&        ! Input
  qw2,sqrtthl2,sqrtqw2,r_qwthl,&  ! Input
  s,std_s,qn,C)                   ! Ouput

  !!!!!!  compute s term
  implicit none

  ! intent-in
  real(rtype), intent(in)  :: qw1
  real(rtype), intent(in)  :: qs1
  real(rtype), intent(in)  :: beta
  real(rtype), intent(in)  :: pval
  real(rtype), intent(in)  :: thl2
  real(rtype), intent(in)  :: qw2
  real(rtype), intent(in)  :: sqrtthl2
  real(rtype), intent(in)  :: sqrtqw2
  real(rtype), intent(in)  :: r_qwthl

  ! intent-out
  real(rtype), intent(out) :: s
  real(rtype), intent(out) :: std_s
  real(rtype), intent(out) :: qn
  real(rtype), intent(out) :: C
  
  ! local variables
  real(rtype) :: cthl, cqt

  ! Parameters
  real(rtype), parameter :: sqrt2 = sqrt(2._rtype)
  real(rtype), parameter :: sqrtpi = sqrt(2._rtype*3.14_rtype)

  s=qw1-qs1*((1._rtype+beta*qw1)/(1._rtype+beta*qs1))
  cthl=((1._rtype+beta*qw1)/(1._rtype+beta*qs1)**2)*(cp/lcond) &
    *beta*qs1*(pval/basepres)**(rgas/cp)

  cqt=1._rtype/(1._rtype+beta*qs1)
  std_s=sqrt(max(0._rtype,cthl**2*thl2+cqt**2*qw2-2._rtype*cthl &
    *sqrtthl2*cqt*sqrtqw2*r_qwthl))

  qn=0._rtype
  C=0._rtype

  if (std_s .ne. 0.0_rtype) then
    C=0.5_rtype*(1._rtype+erf(s/(sqrt2*std_s)))
    IF (C .ne. 0._rtype) qn=s*C+(std_s/sqrtpi)*exp(-0.5_rtype*(s/std_s)**2)
  else
    if (s .gt. 0._rtype) then
      C=1.0_rtype
      qn=s
    endif
  endif

end subroutine shoc_assumed_pdf_compute_s

subroutine shoc_assumed_pdf_compute_sgs_liquid(&
  a, ql1, ql2,& ! Input
  shoc_ql)      ! Output

  ! Compute SGS liquid water mixing ratio
  implicit none

  ! intent-in
  real(rtype), intent(in)  :: a
  real(rtype), intent(in)  :: ql1
  real(rtype), intent(in)  :: ql2

  ! intent-out
  real(rtype), intent(out) :: shoc_ql

  shoc_ql = max(0._rtype,a*ql1+(1._rtype-a)*ql2)

end subroutine shoc_assumed_pdf_compute_sgs_liquid

subroutine shoc_assumed_pdf_compute_cloud_liquid_variance(&
  a,s1,ql1,C1,std_s1,s2,ql2,C2,std_s2,shoc_ql,& ! Input
  shoc_ql2)                                     ! Output

  ! Compute cloud liquid variance (CLUBB formulation, adjusted to SHOC parameters based)
  implicit none

  ! intent-in
  real(rtype), intent(in)  :: a
  real(rtype), intent(in)  :: s1
  real(rtype), intent(in)  :: ql1
  real(rtype), intent(in)  :: C1
  real(rtype), intent(in)  :: std_s1
  real(rtype), intent(in)  :: s2
  real(rtype), intent(in)  :: ql2
  real(rtype), intent(in)  :: C2
  real(rtype), intent(in)  :: std_s2
  real(rtype), intent(in)  :: shoc_ql

  ! intent-out
  real(rtype), intent(out) :: shoc_ql2

  shoc_ql2 = a * ( s1*ql1 + C1*std_s1**2.0 )                  &
    + ( 1._rtype-a ) * ( s2*ql2 + C2*std_s2**2.0 ) - shoc_ql**2.0
  shoc_ql2 = max( 0._rtype, shoc_ql2 )

end subroutine shoc_assumed_pdf_compute_cloud_liquid_variance

subroutine shoc_assumed_pdf_compute_liquid_water_flux(&
  a,w1_1,w_first,ql1,w1_2,ql2,& ! Input
  wqls)                         ! Output
  ! Compute liquid water flux
  implicit none

  ! intent-in
  real(rtype), intent(in)  :: a
  real(rtype), intent(in)  :: w1_1
  real(rtype), intent(in)  :: w_first
  real(rtype), intent(in)  :: ql1
  real(rtype), intent(in)  :: w1_2
  real(rtype), intent(in)  :: ql2

  ! intent-out
  real(rtype), intent(out) :: wqls

  wqls =a*((w1_1-w_first)*ql1)+(1._rtype-a)*((w1_2-w_first)*ql2)

end subroutine shoc_assumed_pdf_compute_liquid_water_flux

subroutine shoc_assumed_pdf_compute_buoyancy_flux(&
  wthlsec,epsterm,wqwsec,pval,wqls,& ! Input
  wthv_sec)                          ! Output
  ! Compute the SGS buoyancy flux
  implicit none

  ! intent-in
  real(rtype), intent(in) :: wthlsec
  real(rtype), intent(in) :: epsterm
  real(rtype), intent(in) :: wqwsec
  real(rtype), intent(in) :: pval
  real(rtype), intent(in) :: wqls

  ! intent-out
  real(rtype), intent(out) :: wthv_sec

  wthv_sec=wthlsec+((1._rtype-epsterm)/epsterm)*basetemp*wqwsec &
  +((lcond/cp)*(basepres/pval)**(rgas/cp)-(1._rtype/epsterm)*basetemp)*wqls

end subroutine shoc_assumed_pdf_compute_buoyancy_flux

!==============================================================
! Advance turbulent kinetic energy equation

subroutine shoc_tke(&
         shcol,nlev,nlevi,dtime,&    ! Input
         wthv_sec,shoc_mix,&         ! Input
         dz_zi,dz_zt,pres,&          ! Input
         u_wind,v_wind,brunt,obklen,&! Input
         zt_grid,zi_grid,pblh,&      ! Input
         tke,tk,tkh, &               ! Input/Output
         isotropy)                   ! Output

  ! Purpose of this subroutine is to advance the SGS
  !  TKE equation due to shear production, buoyant
  !  production, and dissipation processes.

! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  ! number of levels on midpoint grid
  integer, intent(in) :: nlev
  ! number of levels on interface grid
  integer, intent(in) :: nlevi

  ! timestep [s]
  real(rtype), intent(in) :: dtime
  ! SGS buoyancy flux [K m/s]
  real(rtype), intent(in) :: wthv_sec(shcol,nlev)
  ! Mixing length [m]
  real(rtype), intent(in) :: shoc_mix(shcol,nlev)
  ! Meridional wind [m/s]
  real(rtype), intent(in) :: u_wind(shcol,nlev)
  ! Zonal wind [m/s]
  real(rtype), intent(in) :: v_wind(shcol,nlev)
  ! Obukov length
  real(rtype), intent(in) :: obklen(shcol)
  ! thickness on interface grid [m]
  real(rtype), intent(in) :: dz_zi(shcol,nlevi)
  ! thickness on thermodynamic grid [m]
  real(rtype), intent(in) :: dz_zt(shcol,nlev)
  ! pressure [Pa]
  real(rtype), intent(in) :: pres(shcol,nlev)
  ! Brunt Vaisalla frequncy [/s]
  real(rtype), intent(in) :: brunt(shcol,nlev)
  ! heights on midpoint grid [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! heights on interface grid [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)
  ! PBLH height
  real(rtype), intent(in) :: pblh(shcol)

! INPUT/OUTPUT VARIABLES
  ! turbulent kinetic energy [m2/s2]
  real(rtype), intent(inout) :: tke(shcol,nlev)
  ! eddy coefficient for momentum [m2/s]
  real(rtype), intent(inout) :: tk(shcol,nlev)
  ! eddy coefficient for heat [m2/s]
  real(rtype), intent(inout) :: tkh(shcol,nlev)

! OUTPUT VARIABLES
  ! Return to isotropic timescale [s]
  real(rtype), intent(out) :: isotropy(shcol,nlev)

! LOCAL VARIABLES
  real(rtype) :: sterm(shcol,nlevi), sterm_zt(shcol,nlev)
  ! Dissipation term
  real(rtype) :: a_diss(shcol,nlev)
  !column integrated stability
  real(rtype) :: brunt_int(shcol)

  ! Compute integrated column stability in lower troposphere
  call integ_column_stability(nlev, shcol, dz_zt, pres, brunt, brunt_int)

  ! Compute shear production term, which is on interface levels
  ! This follows the methods of Bretheron and Park (2010)

  call compute_shr_prod(nlevi, nlev, shcol, dz_zi, u_wind, v_wind, sterm)

  ! Interpolate shear term from interface to thermo grid
  call linear_interp(zi_grid,zt_grid,sterm,sterm_zt,nlevi,nlev,shcol,0._rtype)

  !advance sgs TKE
  call adv_sgs_tke(nlev, shcol, dtime, shoc_mix, wthv_sec, &
       sterm_zt, tk, tke, a_diss)

  !Compute isotropic time scale [s]
  call isotropic_ts(nlev, shcol, brunt_int, tke, a_diss, brunt, isotropy)

  !Compute eddy diffusivity for heat and momentum
  call eddy_diffusivities(nlev, shcol, obklen, pblh, zt_grid, &
       shoc_mix, sterm_zt, isotropy, tke, tkh, tk)

  return

end subroutine shoc_tke

!==============================================================
! Compute column integrated stability in lower troposphere

subroutine integ_column_stability(nlev, shcol, dz_zt, pres, brunt, brunt_int)

  implicit none
  !intent-ins
  integer,     intent(in) :: nlev, shcol
  ! thickness on thermodynamic grid [m]
  real(rtype), intent(in) :: dz_zt(shcol,nlev)
  ! pressure [Pa]
  real(rtype), intent(in) :: pres(shcol,nlev)
  ! Brunt Vaisalla frequncy [/s]
  real(rtype), intent(in) :: brunt(shcol,nlev)

  !intent-out
  !column integrated stability
  real(rtype), intent(out) :: brunt_int(shcol)

  !local variables
  integer :: i, k

  brunt_int(1:shcol) = 0._rtype
  do k = 1, nlev
     do i = 1, shcol
        if (pres(i,k) .gt. troppres) then
           brunt_int(i) = brunt_int(i) + dz_zt(i,k)*brunt(i,k)
        endif
     enddo
  enddo

  return

end subroutine integ_column_stability

!==============================================================
! Compute shear production term, which is on interface levels
! This follows the methods of Bretheron and Park (2010)

subroutine compute_shr_prod(nlevi, nlev, shcol, dz_zi, u_wind, v_wind, sterm)

  implicit none

  integer,     intent(in)  :: nlevi, nlev, shcol
  ! thickness on interface grid [m]
  real(rtype), intent(in)  :: dz_zi(shcol,nlevi)
  ! Meridional wind [m/s]
  real(rtype), intent(in)  :: u_wind(shcol,nlev)
  ! Zonal wind [m/s]
  real(rtype), intent(in)  :: v_wind(shcol,nlev)

  !intent-outs
  real(rtype), intent(out) :: sterm(shcol,nlevi)

  !local variables
  integer :: i, k, km1
  real(rtype) :: grid_dz, u_grad, v_grad

  !compute shear production term
  do k = 2, nlev
     km1 = k - 1
     do i = 1, shcol
        grid_dz = 1._rtype/dz_zi(i,k)

        ! calculate vertical gradient of u&v wind
        u_grad = grid_dz*(u_wind(i,km1)-u_wind(i,k))
        v_grad = grid_dz*(v_wind(i,km1)-v_wind(i,k))
        sterm(i,k) = u_grad**2+v_grad**2
     enddo
  enddo

  ! Set lower and upper boundary for shear production
  ! Note that the lower bound for shear production has already
  ! been taken into account for the TKE boundary condition,
  ! thus zero out here
  sterm(1:shcol,1)     = 0._rtype
  sterm(1:shcol,nlevi) = 0._rtype

  return

end subroutine compute_shr_prod

!==============================================================
! Advance SGS TKE

subroutine adv_sgs_tke(nlev, shcol, dtime, shoc_mix, wthv_sec, &
     sterm_zt, tk, tke, a_diss)

  implicit none

  !intent -ins
  integer, intent(in) :: nlev, shcol

  ! timestep [s]
  real(rtype), intent(in) :: dtime
  ! Mixing length [m]
  real(rtype), intent(in) :: shoc_mix(shcol,nlev)
  ! SGS buoyancy flux [K m/s]
  real(rtype), intent(in) :: wthv_sec(shcol,nlev)
  ! Interpolate shear production to thermo grid
  real(rtype), intent(in) :: sterm_zt(shcol,nlev)
  ! eddy coefficient for momentum [m2/s]
  real(rtype), intent(in) :: tk(shcol,nlev)

  ! intent-inout
  ! turbulent kinetic energy [m2/s2]
  real(rtype), intent(inout) :: tke(shcol,nlev)

  !intent-out
  ! Dissipation term
  real(rtype), intent(out) :: a_diss(shcol,nlev)

  !local variables
  integer :: i, k
  real(rtype) :: a_prod_bu, a_prod_sh

  real(rtype) :: Ck, Cs, Ce, Ce1, Ce2, Cee

  Cs=0.15_rtype
  Ck=0.1_rtype
  Ce=Ck**3/Cs**4

  Ce1=Ce/0.7_rtype*0.19_rtype
  Ce2=Ce/0.7_rtype*0.51_rtype
  Cee=Ce1+Ce2

  do k = 1, nlev
     do i = 1, shcol

        ! Compute buoyant production term
        a_prod_bu=(ggr/basetemp)*wthv_sec(i,k)

        tke(i,k)=max(0._rtype,tke(i,k))

        ! Shear production term, use diffusivity from
        !  previous timestep
        a_prod_sh=tk(i,k)*sterm_zt(i,k)

        ! Dissipation term
        a_diss(i,k)=Cee/shoc_mix(i,k)*tke(i,k)**1.5

        ! March equation forward one timestep
        tke(i,k)=max(mintke,tke(i,k)+dtime* &
	  (max(0._rtype,a_prod_sh+a_prod_bu)-a_diss(i,k)))

        tke(i,k)=min(tke(i,k),maxtke)
     enddo
  enddo

  return

end subroutine adv_sgs_tke

subroutine isotropic_ts(nlev, shcol, brunt_int, tke, a_diss, brunt, isotropy)
  !------------------------------------------------------------
  ! Compute the return to isotropic timescale as per
  ! Canuto et al. 2004.  This is used to define the
  ! eddy coefficients as well as to diagnose higher
  ! moments in SHOC
  !------------------------------------------------------------

  implicit none

  !intent-ins
  integer, intent(in) :: nlev, shcol

  !column integrated stability
  real(rtype), intent(in) :: brunt_int(shcol)
  ! turbulent kinetic energy [m2/s2]
  real(rtype), intent(in) :: tke(shcol,nlev)
  ! Dissipation term
  real(rtype), intent(in) :: a_diss(shcol,nlev)
  ! Brunt Vaisalla frequncy [/s]
  real(rtype), intent(in) :: brunt(shcol,nlev)

  ! intent-out
  ! Return to isotropic timescale [s]
  real(rtype), intent(out) :: isotropy(shcol,nlev)

  !local vars
  integer     :: i, k
  real(rtype) :: tscale, lambda, buoy_sgs_save

  !Parameters
  real(rtype), parameter :: lambda_low   = 0.001_rtype
  real(rtype), parameter :: lambda_high  = 0.04_rtype
  real(rtype), parameter :: lambda_slope = 0.65_rtype
  real(rtype), parameter :: brunt_low    = 0.02_rtype
  real(rtype), parameter :: maxiso       = 20000.0_rtype ! Return to isotropic timescale [s]

  do k = 1, nlev
     do i = 1, shcol

        ! define the time scale
        tscale=(2.0_rtype*tke(i,k))/a_diss(i,k)

        ! define a damping term "lambda" based on column stability
        lambda=lambda_low+((brunt_int(i)/ggr)-brunt_low)*lambda_slope
        lambda=max(lambda_low,min(lambda_high,lambda))

        buoy_sgs_save=brunt(i,k)
        if (buoy_sgs_save .le. 0._rtype) lambda=0._rtype

        ! Compute the return to isotropic timescale
        isotropy(i,k)=min(maxiso,tscale/(1._rtype+lambda*buoy_sgs_save*tscale**2))
     enddo
  enddo

  return

end subroutine isotropic_ts

subroutine eddy_diffusivities(nlev, shcol, obklen, pblh, zt_grid, &
     shoc_mix, sterm_zt, isotropy, tke, tkh, tk)

  !------------------------------------------------------------
  ! Compute eddy diffusivity for heat and momentum
  !------------------------------------------------------------

  implicit none

  !intent-ins
  integer, intent(in) :: nlev, shcol

  ! Monin-Okbukov length [m]
  real(rtype), intent(in) :: obklen(shcol)
  ! PBL height [m]
  real(rtype), intent(in) :: pblh(shcol)
  ! Heights on the mid-point grid [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! Mixing length [m]
  real(rtype), intent(in) :: shoc_mix(shcol,nlev)
  ! Interpolate shear production to thermo grid
  real(rtype), intent(in) :: sterm_zt(shcol,nlev)
  ! Return to isotropic timescale [s]
  real(rtype), intent(in) :: isotropy(shcol,nlev)
  ! turbulent kinetic energy [m2/s2]
  real(rtype), intent(in) :: tke(shcol,nlev)

  !intent-outs
  ! eddy coefficient for heat [m2/s]
  real(rtype), intent(out) :: tkh(shcol,nlev)
  ! eddy coefficient for momentum [m2/s]
  real(rtype), intent(out) :: tk(shcol,nlev)

  !local vars
  integer     :: i, k
  real(rtype) :: z_over_L, zt_grid_1d(shcol)
  real(rtype) :: Ckh_s, Ckm_s

  !parameters
  ! Critical value of dimensionless Monin-Obukhov length,
  !  for which diffusivities are no longer damped
  real(rtype), parameter :: zL_crit_val = 100.0_rtype
  ! Transition depth [m] above PBL top to allow
  ! stability diffusivities
  real(rtype), parameter :: pbl_trans = 200.0_rtype
  ! Turbulent coefficients
  real(rtype), parameter :: Ckh = 0.1_rtype
  real(rtype), parameter :: Ckm = 0.1_rtype
  ! Default eddy coefficients for stable PBL diffusivities
  real(rtype), parameter :: Ckh_s_def = 1.0_rtype
  real(rtype), parameter :: Ckm_s_def = 1.0_rtype
  ! Minimum allowable value for stability diffusivities
  real(rtype), parameter :: Ck_s_min = 0.1_rtype

  !store zt_grid at nlev in 1d array
  zt_grid_1d(1:shcol) = zt_grid(1:shcol,nlev)

  do k = 1, nlev
     do i = 1, shcol

        ! Dimensionless Okukhov length considering only
        !  the lowest model grid layer height to scale
        z_over_L = zt_grid_1d(i)/obklen(i)

        if (z_over_L .gt. 0._rtype .and. (zt_grid(i,k) .lt. pblh(i)+pbl_trans)) then
           ! If surface layer is stable, based on near surface
           !  dimensionless Monin-Obukov use modified coefficients of
           !  tkh and tk that are primarily based on shear production
           !  and SHOC length scale, to promote mixing within the PBL
           !  and to a height slighty above to ensure smooth transition.

	   ! Compute diffusivity coefficient as function of dimensionless
           !  Obukhov, given a critical value
           Ckh_s = max(Ck_s_min,min(Ckh_s_def,z_over_L/zL_crit_val))
           Ckm_s = max(Ck_s_min,min(Ckm_s_def,z_over_L/zL_crit_val))

	   ! Compute stable PBL diffusivities
           tkh(i,k) = Ckh_s*(shoc_mix(i,k)**2)*sqrt(sterm_zt(i,k))
           tk(i,k)  = Ckm_s*(shoc_mix(i,k)**2)*sqrt(sterm_zt(i,k))
        else
           ! Default definition of eddy diffusivity for heat and momentum
           tkh(i,k) = Ckh*isotropy(i,k)*tke(i,k)
           tk(i,k)  = Ckm*isotropy(i,k)*tke(i,k)
        endif

     enddo
  enddo

  return

end subroutine eddy_diffusivities

!==============================================================
! Check the TKE

subroutine check_tke(&
             shcol,nlev,& ! Input
             tke)         ! Input/Output

  implicit none
  ! Make sure TKE falls within reasonable bounds
  ! If not, then clip

! INPUT VARIABLES
  integer, intent(in) :: shcol
  integer, intent(in) :: nlev

! IN/OUT VARIABLES
  real(rtype), intent(inout) :: tke(shcol,nlev)

! LOCAL VARIABLES
  integer :: i, k

  do k=1,nlev
    do i=1,shcol
      tke(i,k)=max(mintke,tke(i,k))
    enddo
  enddo

end subroutine check_tke

!==============================================================
! Compute the turbulent length scale

subroutine shoc_length(&
         shcol,nlev,nlevi,tke,&        ! Input
         host_dx,host_dy,pblh,&        ! Input
         zt_grid,zi_grid,dz_zt,dz_zi,& ! Input
         thetal,wthv_sec,thv,&         ! Input
         brunt,shoc_mix)               ! Output

  ! Purpose of this subroutine is to compute the SHOC
  !  mixing length scale, which is used to compute the
  !  turbulent dissipation in the SGS TKE equation

  implicit none

! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  !  number of midpoint levels
  integer, intent(in) :: nlev
  ! number of interface levels
  integer, intent(in) :: nlevi

  ! host model grid size [m]
  real(rtype), intent(in) :: host_dx(shcol)
  ! host model grid size [m]
  real(rtype), intent(in) :: host_dy(shcol)
  ! Planetary boundary layer (PBL) height [m]
  real(rtype), intent(in) :: pblh(shcol)
  ! turbulent kinetic energy [m^2/s^2]
  real(rtype), intent(in) :: tke(shcol,nlev)
  ! heights on midpoint grid [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! heights on interface grid [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)
  ! dz on midpoint grid [m]
  real(rtype), intent(in) :: dz_zt(shcol,nlev)
  ! dz on interface grid [m]
  real(rtype), intent(in) :: dz_zi(shcol,nlevi)
  ! SGS buoyancy flux [K m/s]
  real(rtype), intent(in) :: wthv_sec(shcol,nlev)
  ! liquid water potential temperature [K]
  real(rtype), intent(in) :: thetal(shcol,nlev)
  ! virtual potential temperature [K]
  real(rtype), intent(in) :: thv(shcol,nlev)

  ! OUTPUT VARIABLES
  ! Brunt-Vaisala frequency [/s]
  real(rtype), intent(out) :: brunt(shcol,nlev)
  ! SHOC mixing length [m]
  real(rtype), intent(out) :: shoc_mix(shcol,nlev)

  ! LOCAL VARIABLES
  real(rtype) :: conv_vel(shcol), tscale(shcol)
  real(rtype) :: thv_zi(shcol,nlevi)
  real(rtype) :: l_inf(shcol)

  ! Interpolate virtual potential temperature onto interface grid
  call linear_interp(zt_grid,zi_grid,thv,thv_zi,nlev,nlevi,shcol,0._rtype)

  ! Define the brunt vaisalia frequency
  call compute_brunt_shoc_length(nlev,nlevi,shcol,dz_zt,thv,thv_zi,brunt)

  ! Find L_inf
  call compute_l_inf_shoc_length(nlev,shcol,zt_grid,dz_zt,tke,l_inf)

  ! determine the convective velocity scale of
  !   the planetary boundary layer
  call compute_conv_vel_shoc_length(nlev,shcol,pblh,zt_grid,dz_zt,thv,wthv_sec,conv_vel)

  ! computed quantity above is wstar3
  ! clip, to avoid negative values and take the cubed
  !   root to get the convective velocity scale

  ! Compute eddy turnover timescale.  If
  !  convective velocity scale is zero then
  !  set to a minimum threshold
  call compute_conv_time_shoc_length(shcol,pblh,conv_vel,tscale)

  ! compute mixing-length
  call compute_shoc_mix_shoc_length(nlev,shcol,tke,brunt,tscale,zt_grid,l_inf,shoc_mix)

  ! Do checks on the length scale.  Make sure it is not
  !  larger than the grid mesh of the host model.
  call check_length_scale_shoc_length(nlev,shcol,host_dx,host_dy,shoc_mix)

  return

end subroutine shoc_length

!==============================================================
! Subroutine to determine superdiagonal and subdiagonal coeffs of
! the tridiagonal diffusion matrix.

subroutine vd_shoc_decomp( &
         shcol,nlev,nlevi,&          ! Input
         kv_term,tmpi,rdp_zt,dtime,& ! Input
         flux, &                     ! Input
         ca,cc,denom,ze)             ! Output

  implicit none

! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  ! number of mid-point levels
  integer, intent(in) :: nlev
  ! number of levels on the interface
  integer, intent(in) :: nlevi

  ! SHOC timestep [s]
  real(rtype), intent(in) :: dtime
  ! diffusion coefficent [m2/s]
  real(rtype), intent(in) :: kv_term(shcol,nlevi)
  ! dt*(g*rho)**2/dp at interfaces
  real(rtype), intent(in) :: tmpi(shcol,nlevi)
  ! 1/dp
  real(rtype), intent(in) :: rdp_zt(shcol,nlev)
  ! surface flux [varies]
  real(rtype), intent(in) :: flux(shcol)

! OUTPUT VARIABLES
  ! superdiagonal
  real(rtype), intent(out) :: ca(shcol,nlev)
  ! subdiagonal
  real(rtype), intent(out) :: cc(shcol,nlev)
  ! 1./(1.+ca(k)+cc(k)-cc(k)*ze(k-1))
  real(rtype), intent(out) :: denom(shcol,nlev)
  ! Term in tri-diag. matrix system
  real(rtype), intent(out) :: ze(shcol,nlev)

! LOCAL VARIABLES
  integer :: i, k

  ! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the
  ! tridiagonal diffusion matrix. The diagonal elements  (cb=1+ca+cc) are
  ! a combination of ca and cc; they are not required by the solver.

  do k=nlev-1,1,-1
    do i=1,shcol
      ca(i,k) = kv_term(i,k+1) * tmpi(i,k+1) * rdp_zt(i,k)
      cc(i,k+1) = kv_term(i,k+1) * tmpi(i,k+1) * rdp_zt(i,k+1)
    enddo
  enddo

  ! The bottom element of the upper diagonal (ca) is zero (not used).
  ! The subdiagonal (cc) is not needed in the solver.

  ca(:,nlev) = 0._rtype

  ! Calculate e(k). This term is required in the solution of the
  ! tridiagonal matrix as defined by the implicit diffusion equation.

  do i=1,shcol
    denom(i,nlev) = 1._rtype/ &
      (1._rtype + cc(i,nlev) + flux(i)*dtime*ggr*rdp_zt(i,nlev))
    ze(i,nlev) = cc(i,nlev) * denom(i,nlev)
  enddo

  do k=nlev-1,2,-1
    do i=1,shcol
      denom(i,k) = 1._rtype/ &
        (1._rtype + ca(i,k) + cc(i,k) - &
    ca(i,k) * ze(i,k+1))
      ze(i,k) = cc(i,k) * denom(i,k)
    enddo
  enddo

  do i=1,shcol
    denom(i,1) = 1._rtype/ &
      (1._rtype + ca(i,1) - ca(i,1) * ze(i,2))
  enddo

  return

end subroutine vd_shoc_decomp

!==============================================================
! Subroutine to solve the implicit vertical diffsion equation
! with zero flux boundary conditions.  Actual surface fluxes
! should be applied explicitly.  Procedure for solution of the
! implicit equation follows Richtmeyer and Morton (1967, pp 198-200).
! The equation solved is
!
!     -ca(k)*q(k+1) + cb(k)*q(k) - cc(k)*q(k-1) = d(k),
!
! where d(k) is the input profile and q(k) is the output profile
!
! The solution has the form
!
!     q(k) = ze(k)*q(k-1) + zf(k)
!
!     ze(k) = cc(k) * dnom(k)
!
!     zf(k) = [d(k) + ca(k)*zf(k+1)] * dnom(k)
!
!     dnom(k) = 1/[cb(k) - ca(k)*ze(k+1)]
!             = 1/[1 + ca(k) + cc(k) - ca(k)*ze(k+1)]
!
! Note that the same routine is used for temperature, momentum and
! tracers, and that input variables are replaced.
! ---------------------------------------------------------------

subroutine vd_shoc_solve(&
         shcol,nlev,nlevi,&   ! Input
         ca,cc,denom,ze,&     ! Input
         var)                 ! Input/Output

  implicit none

! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  ! number of mid-point levels
  integer, intent(in) :: nlev
  ! number of levels on the interface
  integer, intent(in) :: nlevi
  ! superdiagonal
  real(rtype), intent(in) :: ca(shcol,nlev)
  ! subdiagonal
  real(rtype), intent(in) :: cc(shcol,nlev)
  ! 1./(1.+ca(k)+cc(k)-cc(k)*ze(k-1))
  real(rtype), intent(in) :: denom(shcol,nlev)
  ! Term in tri-diag. matrix system
  real(rtype), intent(in) :: ze(shcol,nlev)

! IN/OUT VARIABLES
  real(rtype), intent(inout) :: var(shcol,nlev)

! LOCAL VARIABLES
  integer :: i, k
  ! Term in tri-diag solution
  real(rtype) :: zf(shcol,nlev)


  ! Calculate zf(k). Terms zf(k) and ze(k) are required in solution of
  ! tridiagonal matrix defined by implicit diffusion equation.
  ! Note that only levels ntop through nbot need be solved for.

  do i=1,shcol
    zf(i,nlev) = var(i,nlev) * denom(i,nlev)
  enddo

  do k=nlev-1,1,-1
    do i=1,shcol
      zf(i,k) = (var(i,k) + ca(i,k) * zf(i,k+1)) * denom(i,k)
    enddo
  enddo

  ! Perform back substitution

  do i=1,shcol
    var(i,1) = zf(i,1)
  enddo

  do k=2,nlev
    do i=1,shcol
      var(i,k) = zf(i,k) + ze(i,k)*var(i,k-1)
    enddo
  enddo

  return

end subroutine vd_shoc_solve

!==============================================================
! Subroutine to compute integrals for SHOC conservation
!  with host model

subroutine shoc_energy_integrals(&
         shcol,nlev,host_dse,pdel,&     ! Input
         rtm,rcm,u_wind,v_wind,&        ! Input
         se_int,ke_int,wv_int,wl_int)   ! Output

  implicit none

! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  ! number of levels
  integer, intent(in) :: nlev
  ! host model temperature [K]
  real(rtype), intent(in) :: host_dse(shcol,nlev)
  ! pressure layer thickness [Pa]
  real(rtype), intent(in) :: pdel(shcol,nlev)
  ! zonal wind [m/s]
  real(rtype), intent(in) :: u_wind(shcol,nlev)
  ! meridional wind [m/s]
  real(rtype), intent(in) :: v_wind(shcol,nlev)
  ! total water mixing ratio [kg/kg]
  real(rtype), intent(in) :: rtm(shcol,nlev)
  ! cloud liquid mixing ratio [kg/kg]
  real(rtype), intent(in) :: rcm(shcol,nlev)

! OUTPUT VARIABLES
  ! integrated static energy
  real(rtype), intent(out) :: se_int(shcol)
  ! integrated kinetic energy
  real(rtype), intent(out) :: ke_int(shcol)
  ! integrated water vapor
  real(rtype), intent(out) :: wv_int(shcol)
  ! integrated liquid water
  real(rtype), intent(out) :: wl_int(shcol)

! LOCAL VARIABLES
  integer :: i, k
  real(rtype) :: rvm

  se_int(:) = 0._rtype
  ke_int(:) = 0._rtype
  wv_int(:) = 0._rtype
  wl_int(:) = 0._rtype
  do k=1,nlev
    do i=1,shcol
       rvm = rtm(i,k) - rcm(i,k) ! compute water vapor
       se_int(i) = se_int(i) + host_dse(i,k)*pdel(i,k)/ggr
       ke_int(i) = ke_int(i) + 0.5_rtype*(u_wind(i,k)**2+v_wind(i,k)**2)*pdel(i,k)/ggr
       wv_int(i) = wv_int(i) + rvm*pdel(i,k)/ggr
       wl_int(i) = wl_int(i) + rcm(i,k)*pdel(i,k)/ggr
    enddo
  enddo

  return

end subroutine shoc_energy_integrals

!==============================================================
! Subroutine to update SHOC output to host model temperature

subroutine update_host_dse(&
         shcol,nlev,thlm,&                 ! Input
         shoc_ql,exner,zt_grid,phis,&      ! Input
         host_dse)                         ! Output

#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: update_host_dse_f
#endif

  implicit none

  ! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  ! number of levels
  integer, intent(in) :: nlev
  ! liquid water potential temperature [K]
  real(rtype), intent(in) :: thlm(shcol,nlev)
  ! cloud liquid water mixing ratio [kg/kg]
  real(rtype), intent(in) :: shoc_ql(shcol,nlev)
  ! exner function [-]
  real(rtype), intent(in) :: exner(shcol,nlev)
  ! heights at mid point [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! surface geopotential height of host model [m]
  real(rtype), intent(in) :: phis(shcol)

  ! OUTPUT VARIABLES
  ! host model temperature [K]
  real(rtype), intent(out) :: host_dse(shcol,nlev)

  ! LOCAL VARIABLES
  ! Temperature [K]
  real(rtype) :: temp

  integer :: i, k

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
      call update_host_dse_f(shcol,nlev,thlm,shoc_ql,exner,zt_grid,phis, &  ! Input
           host_dse)                              ! Input/Output)
      return
   endif
#endif

  do k=1,nlev
    do i=1,shcol
      temp = (thlm(i,k)+(lcond/cp)*shoc_ql(i,k))/exner(i,k)
      host_dse(i,k) = cp*temp+ggr*zt_grid(i,k)+phis(i)
    enddo
  enddo

  return

end subroutine update_host_dse

!==============================================================
! Subroutine for SHOC energy fixer with host model temp

subroutine shoc_energy_fixer(&
         shcol,nlev,nlevi,dtime,nadv,&  ! Input
         zt_grid,zi_grid,&              ! Input
         se_b,ke_b,wv_b,wl_b,&          ! Input
         se_a,ke_a,wv_a,wl_a,&          ! Input
         wthl_sfc,wqw_sfc,pdel,&        ! Input
         rho_zt,tke,pint,&              ! Input
         host_dse)                      ! Input/Output

  implicit none

  ! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  ! number of levels
  integer, intent(in) :: nlev
  ! number of levels on interface grid
  integer, intent(in) :: nlevi
  ! SHOC timestep
  real(rtype), intent(in) :: dtime
  ! number of SHOC iterations
  integer, intent(in) :: nadv
  ! integrated static energy before
  real(rtype), intent(in) :: se_b(shcol)
  ! integrated kinetic energy before
  real(rtype), intent(in) :: ke_b(shcol)
  ! integrated water vapor before
  real(rtype), intent(in) :: wv_b(shcol)
  ! integrated liquid water before
  real(rtype), intent(in) :: wl_b(shcol)
  ! integrated static energy after
  real(rtype), intent(in) :: se_a(shcol)
  ! integrated kinetic energy after
  real(rtype), intent(in) :: ke_a(shcol)
  ! integrated water vapor after
  real(rtype), intent(in) :: wv_a(shcol)
  ! integrated liquid water after
  real(rtype), intent(in) :: wl_a(shcol)
  ! Surface sensible heat flux [K m/s]
  real(rtype), intent(in) :: wthl_sfc(shcol)
  ! Surface latent heat flux [kg/kg m/s]
  real(rtype), intent(in) :: wqw_sfc(shcol)
  ! pressure differenes [Pa]
  real(rtype), intent(in) :: pdel(shcol,nlev)
  ! heights on midpoint grid [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! heights on interface grid [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)
  ! pressure on interface grid [Pa]
  real(rtype), intent(in) :: pint(shcol,nlevi)
  ! density on midpoint grid [kg/m^3]
  real(rtype), intent(in) :: rho_zt(shcol,nlev)
  !turbulent kinetic energy [m^2/s^2]
  real(rtype), intent(in) :: tke(shcol,nlev)

  ! INPUT VARIABLES
  !host temperature [K]
  real(rtype), intent(inout) :: host_dse(shcol,nlev)

  ! LOCAL VARIABLES
  real(rtype) :: se_dis(shcol), te_a(shcol), te_b(shcol)
  integer :: shoctop(shcol)

  call shoc_energy_total_fixer(&
         shcol,nlev,nlevi,dtime,nadv,&  ! Input
         zt_grid,zi_grid,&              ! Input
         se_b,ke_b,wv_b,wl_b,&          ! Input
         se_a,ke_a,wv_a,wl_a,&          ! Input
         wthl_sfc,wqw_sfc,rho_zt,&      ! Input
         te_a, te_b)                    ! Output

  call shoc_energy_threshold_fixer(&
         shcol,nlev,nlevi,&             ! Input
         pint,tke,te_a,te_b,&           ! Input
         se_dis,shoctop)                ! Output

  call shoc_energy_dse_fixer(&
         shcol,nlev,&                  ! Input
         se_dis,shoctop,   &           ! Input
         host_dse)                     ! Input/Output

  return

end subroutine shoc_energy_fixer

!==============================================================
! Subroutine foe SHOC energy fixer with host model temp

subroutine shoc_energy_total_fixer(&
         shcol,nlev,nlevi,dtime,nadv,&  ! Input
         zt_grid,zi_grid,&              ! Input
         se_b,ke_b,wv_b,wl_b,&          ! Input
         se_a,ke_a,wv_a,wl_a,&          ! Input
         wthl_sfc,wqw_sfc,rho_zt,&      ! Input
         te_a, te_b)                    ! Output

  implicit none

  ! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  ! number of levels
  integer, intent(in) :: nlev
  ! number of levels on interface grid
  integer, intent(in) :: nlevi
  ! SHOC timestep
  real(rtype), intent(in) :: dtime
  ! number of SHOC iterations
  integer, intent(in) :: nadv
  ! integrated static energy before
  real(rtype), intent(in) :: se_b(shcol)
  ! integrated kinetic energy before
  real(rtype), intent(in) :: ke_b(shcol)
  ! integrated water vapor before
  real(rtype), intent(in) :: wv_b(shcol)
  ! integrated liquid water before
  real(rtype), intent(in) :: wl_b(shcol)
  ! integrated static energy after
  real(rtype), intent(in) :: se_a(shcol)
  ! integrated kinetic energy after
  real(rtype), intent(in) :: ke_a(shcol)
  ! integrated water vapor after
  real(rtype), intent(in) :: wv_a(shcol)
  ! integrated liquid water after
  real(rtype), intent(in) :: wl_a(shcol)
  ! Surface sensible heat flux [K m/s]
  real(rtype), intent(in) :: wthl_sfc(shcol)
  ! Surface latent heat flux [kg/kg m/s]
  real(rtype), intent(in) :: wqw_sfc(shcol)
  ! heights on midpoint grid [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! heights on interface grid [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)
  ! density on midpoint grid [kg/m^3]
  real(rtype), intent(in) :: rho_zt(shcol,nlev)

  ! OUTPUT VARIABLES
  real(rtype), intent(out) :: te_a(shcol)
  real(rtype), intent(out) :: te_b(shcol)

  ! LOCAL VARIABLES
  ! density on interface grid [kg/m^3]
  real(rtype) :: rho_zi(shcol,nlevi)
  ! sensible and latent heat fluxes [W/m^2]
  real(rtype) :: shf, lhf, hdtime
  integer :: i, k

  ! compute the host timestep
  hdtime = dtime * float(nadv)

  call linear_interp(zt_grid,zi_grid,rho_zt,rho_zi,nlev,nlevi,shcol,0._rtype)

  ! Based on these integrals, compute the total energy before and after SHOC
  ! call
  do i=1,shcol
    ! convert shf and lhf to W/m^2
    shf=wthl_sfc(i)*cp*rho_zi(i,nlevi)
    lhf=wqw_sfc(i)*rho_zi(i,nlevi)
    te_a(i) = se_a(i) + ke_a(i) + (lcond+lice)*wv_a(i)+lice*wl_a(i)
    te_b(i) = se_b(i) + ke_b(i) + (lcond+lice)*wv_b(i)+lice*wl_b(i)
    te_b(i) = te_b(i)+(shf+(lhf)*(lcond+lice))*hdtime
  enddo

  return

end subroutine shoc_energy_total_fixer



!==============================================================
! Subroutine foe SHOC energy fixer with host model temp

subroutine shoc_energy_threshold_fixer(&
         shcol,nlev,nlevi,&             ! Input
         pint,tke,te_a,te_b,&           ! Input
         se_dis,shoctop)                ! Output

  implicit none

  ! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  ! number of levels
  integer, intent(in) :: nlev
  ! number of levels on interface grid
  integer, intent(in) :: nlevi
  ! pressure on interface grid [Pa]
  real(rtype), intent(in) :: pint(shcol,nlevi)
  !turbulent kinetic energy [m^2/s^2]
  real(rtype), intent(in) :: tke(shcol,nlev)

  real(rtype), intent(in) :: te_a(shcol)
  real(rtype), intent(in) :: te_b(shcol)


  ! OUTPUT VARIABLES
  real(rtype), intent(out) :: se_dis(shcol)
  integer, intent(out) :: shoctop(shcol)

  ! LOCAL VARIABLES
  ! sensible and latent heat fluxes [W/m^2]
  integer :: i, k

  ! Limit the energy fixer to find highest layer where SHOC is active
  ! Find first level where wp2 is higher than lowest threshold
  do i=1,shcol
    shoctop(i) = 1
    do while (tke(i,shoctop(i)) .eq. mintke .and. shoctop(i) .lt. nlev-1)
      shoctop(i) = shoctop(i) + 1
    enddo

    ! Compute the disbalance of total energy, over depth where SHOC is active
    se_dis(i) = (te_a(i) - te_b(i))/(pint(i,nlevi)-pint(i,shoctop(i)))
  enddo

  return

end subroutine shoc_energy_threshold_fixer

!==============================================================
! Subroutine foe SHOC energy fixer with host model temp

subroutine shoc_energy_dse_fixer(&
         shcol,nlev,&                  ! Input
         se_dis,shoctop,   &           ! Input
         host_dse)                     ! Input/Output

  implicit none

  ! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  ! number of levels
  integer, intent(in) :: nlev

  ! INPUT VARIABLES
  real(rtype), intent(in) :: se_dis(shcol)
  integer, intent(in) :: shoctop(shcol)

  !host temperature [K]
  real(rtype), intent(inout) :: host_dse(shcol,nlev)

  ! LOCAL VARIABLES
  integer :: i, k

  do i=1,shcol
    do k=shoctop(i),nlev
      host_dse(i,k) = host_dse(i,k) - se_dis(i)*ggr
    enddo
  enddo

  return

end subroutine shoc_energy_dse_fixer

!==============================================================
! Linear interpolation to get values on various grids

subroutine shoc_diag_obklen(&
         shcol,uw_sfc,vw_sfc,&      ! Input
         wthl_sfc,wqw_sfc,thl_sfc,& ! Input
         cldliq_sfc,qv_sfc,&        ! Input
         ustar,kbfs,obklen)         ! Ouput

  implicit none

! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  ! Surface sensible heat flux [K m/s]
  real(rtype), intent(in) :: wthl_sfc(shcol)
  ! Surface latent heat flux [kg/kg m/s]
  real(rtype), intent(in) :: wqw_sfc(shcol)
  ! Surface momentum flux (u-direction) [m2/s2]
  real(rtype), intent(in) :: uw_sfc(shcol)
  ! Surface momentum flux (v-direction) [m2/s2]
  real(rtype), intent(in) :: vw_sfc(shcol)
  ! Surface potential temperature [K]
  real(rtype), intent(in) :: thl_sfc(shcol)
  ! Surface cloud liquid water [kg /kg]
  real(rtype), intent(in) :: cldliq_sfc(shcol)
  ! Surface water vapor
  real(rtype), intent(in) :: qv_sfc(shcol)

! OUTPUT VARIABLES
  ! Obukhov length [m]
  real(rtype), intent(out) :: obklen(shcol)
  ! surface friction velocity [m/s]
  real(rtype), intent(out) :: ustar(shcol)
  ! surface kinematic buoyancy flux [m^s/s^3]
  real(rtype), intent(out) :: kbfs(shcol)

! LOCAL VARIABLES
  integer :: i
  real(rtype) :: th_sfc  ! potential temperature at surface
  real(rtype) :: thv_sfc ! virtual potential temperature at surface

  do i=1,shcol
    th_sfc = thl_sfc(i) + (lcond/cp)*cldliq_sfc(i)
    thv_sfc = th_sfc*(1._rtype+eps*qv_sfc(i)-cldliq_sfc(i))
    ustar(i) = max(sqrt(uw_sfc(i)**2 + vw_sfc(i)**2),ustar_min)
    kbfs(i) = wthl_sfc(i)+eps*th_sfc*wqw_sfc(i)
    obklen(i) = -thv_sfc*ustar(i)**3/(ggr*vk*(kbfs(i)+sign(1.e-10_rtype,kbfs(i))))
  enddo

  return

end subroutine shoc_diag_obklen

  !
  !===============================================================================
subroutine pblintd(&
       shcol,nlev,nlevi,&             ! Input
       z,zi,thl,ql,&                  ! Input
       q,u,v,&                        ! Input
       ustar,obklen,kbfs,cldn,&       ! Input
       pblh)                          ! Output

    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Diagnose standard PBL variables
    !
    ! Method:
    ! Diagnosis of PBL depth.
    ! The PBL depth follows:
    !    Holtslag, A.A.M., and B.A. Boville, 1993:
    !    Local versus Nonlocal Boundary-Layer Diffusion in a Global Climate
    !    Model. J. Clim., vol. 6., p. 1825--1842.
    !
    ! Updated by Holtslag and Hack to exclude the surface layer from the
    ! definition of the boundary layer Richardson number. Ri is now defined
    ! across the outer layer of the pbl (between the top of the surface
    ! layer and the pbl top) instead of the full pbl (between the surface and
    ! the pbl top). For simiplicity, the surface layer is assumed to be the
    ! region below the first model level (otherwise the boundary layer depth
    ! determination would require iteration).
    !
    ! Modified for boundary layer height diagnosis: Bert Holtslag, june 1994
    ! (Use ricr = 0.3 in this formulation)
    !
    ! Author: B. Stevens (extracted from pbldiff, August 2000)
    !
    !-----------------------------------------------------------------------
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    integer, intent(in) :: shcol                     ! number of atmospheric columns
    integer, intent(in) :: nlev                      ! number of mid-point layers
    integer, intent(in) :: nlevi                     ! number of interface layers

    real(rtype), intent(in)  :: thl(shcol,nlev)         ! liquid water potential temp [K]
    real(rtype), intent(in)  :: ql(shcol,nlev)          ! cloud liquid mixing ratio [kg/kg]
    real(rtype), intent(in)  :: q(shcol,nlev)           ! water vapor [kg/kg]
    real(rtype), intent(in)  :: z(shcol,nlev)           ! height above surface [m]
    real(rtype), intent(in)  :: u(shcol,nlev)           ! windspeed x-direction [m/s]
    real(rtype), intent(in)  :: v(shcol,nlev)           ! windspeed y-direction [m/s]
    real(rtype), intent(in)  :: ustar(shcol)            ! surface friction velocity [m/s]
    real(rtype), intent(in)  :: obklen(shcol)           ! Obukhov length
    real(rtype), intent(in)  :: kbfs(shcol)             ! sfc kinematic buoyancy flux [m^2/s^3]
    real(rtype), intent(in)  :: zi(shcol,nlevi)         ! height above surface [m]
    real(rtype), intent(in)  :: cldn(shcol,nlev)        ! new cloud fraction
    !
    ! Output arguments
    !
    real(rtype), intent(out) :: pblh(shcol)             ! boundary-layer height [m]
    !
    !---------------------------Local workspace-----------------------------
    !

    real(rtype) :: rino(shcol,nlev)        ! bulk Richardson no. from level to ref lev
    real(rtype) :: thv(shcol,nlev)         ! virtual potential temperature
    real(rtype) :: tlv(shcol)              ! ref. level pot tmp + tmp excess

    logical  :: check(shcol)            ! True=>chk if Richardson no.>critcal

    !
    ! Compute Obukhov length virtual temperature flux and various arrays for use later:
    !

    ! Compute virtual potential temperature
    call pblintd_init_pot(&
               shcol,nlev,&             ! Input
               thl,ql,q,&               ! Input
               thv)                     ! Output

    call pblintd_init(&
           shcol,nlev,&             ! Input
           z,&                      ! Input
           check,rino,pblh)         ! Output

    !
    ! PBL height calculation
    !
    call pblintd_height(&
       shcol,nlev,&                   ! Input
       z,u,v,ustar,&                  ! Input
       thv,thv(:,nlev),&              ! Input
       pblh,rino,check)               ! Output
    !
    ! Estimate an effective surface temperature to account for surface
    ! fluctuations
    !
    call pblintd_surf_temp(&
       shcol,nlev,nlevi,&          ! Input
       z,ustar,obklen,kbfs,thv,&   ! Input
       tlv,&                       ! Output
       pblh,check,rino)            ! InOutput
    !
    ! Improve pblh estimate for unstable conditions using the convective
    ! temperature excess as reference temperature:
    !
    call pblintd_height(&
       shcol,nlev,&             ! Input
       z,u,v,ustar,&            ! Input
       thv,tlv,&                ! Input
       pblh,rino,check)         ! Output
    !
    ! Check PBL height
    !
    call pblintd_check_pblh(&
       shcol,nlev,nlevi,&             ! Input
       z,ustar,check,&                ! Input
       pblh)                          ! Output
    !
    ! PBL check over ocean
    !
    call pblintd_cldcheck(      &
                   shcol,nlev,nlevi, &                  ! Input
                   zi,cldn,          &                  ! Input
                   pblh)                                ! InOutput

    return
end subroutine pblintd

subroutine pblintd_init_pot(&
       shcol,nlev,&             ! Input
       thl,ql,q,&               ! Input
       thv)                     ! Output
#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: shoc_pblintd_init_pot_f
#endif
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    integer, intent(in) :: shcol                     ! number of atmospheric columns
    integer, intent(in) :: nlev                      ! number of mid-point layers

    real(rtype), intent(in)  :: thl(shcol,nlev)         ! liquid water potential temp [K]
    real(rtype), intent(in)  :: ql(shcol,nlev)          ! cloud liquid mixing ratio [kg/kg]
    real(rtype), intent(in)  :: q(shcol,nlev)           ! water vapor [kg/kg]
    real(rtype), intent(out) :: thv(shcol,nlev)         ! virtual potential temperature
    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i                       ! longitude index
    integer  :: k                       ! level index
    real(rtype) :: th

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
      call shoc_pblintd_init_pot_f(shcol,nlev,thl,ql,q,&               ! Input
                                   thv)                     ! Output
      return
   endif
#endif
    ! Compute virtual potential temperature
    do k=1,nlev
      do i=1,shcol
        th=thl(i,k)+(lcond/cp)*ql(i,k)
        thv(i,k)=th+(1._rtype+eps*q(i,k)-ql(i,k))
      enddo
    enddo

end subroutine pblintd_init_pot

subroutine pblintd_init(&
       shcol,nlev,&             ! Input
       z,&                      ! Input
       check,rino,pblh)         ! Output
    !------------------------------Arguments--------------------------------
    ! Input arguments
    !
    integer, intent(in) :: shcol                     ! number of atmospheric columns
    integer, intent(in) :: nlev                      ! number of mid-point layers

    real(rtype), intent(in)  :: z(shcol,nlev)           ! height above surface [m]
    !
    ! Output arguments
    !
    real(rtype), intent(out) :: pblh(shcol)             ! boundary-layer height [m]
    real(rtype), intent(out) :: rino(shcol,nlev)        ! bulk Richardson no. from level to ref lev
    logical, intent(out)     :: check(shcol)            ! True=>chk if Richardson no.>critcal

    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i                       ! longitude index
    do i=1,shcol
       check(i)     = .true.
       rino(i,nlev) = 0.0_rtype
       pblh(i)      = z(i,nlev)
    end do
end subroutine pblintd_init

subroutine pblintd_height(&
       shcol,nlev,&              ! Input
       z,u,v,ustar,&             ! Input
       thv,thv_ref,&             ! Input
       pblh,rino,check)          ! Output
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    integer, intent(in) :: shcol                     ! number of atmospheric columns
    integer, intent(in) :: nlev                      ! number of mid-point layers

    real(rtype), intent(in)  :: z(shcol,nlev)           ! height above surface [m]
    real(rtype), intent(in)  :: u(shcol,nlev)           ! windspeed x-direction [m/s]
    real(rtype), intent(in)  :: v(shcol,nlev)           ! windspeed y-direction [m/s]
    real(rtype), intent(in)  :: ustar(shcol)            ! surface friction velocity [m/s]
    real(rtype), intent(in)  :: thv(shcol,nlev)         ! virtual potential temperature
    real(rtype), intent(in)  :: thv_ref(shcol)          ! ref. level pot tmp

    !
    ! Output arguments
    !
    real(rtype), intent(out)   :: pblh(shcol)             ! boundary-layer height [m]
    real(rtype), intent(inout) :: rino(shcol,nlev)                     ! bulk Richardson no. from level to ref lev
    logical, intent(inout)     :: check(shcol)            ! True=>chk if Richardson no.>critcal

    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i                       ! longitude index
    integer  :: k                       ! level index
    real(rtype) :: vvk                     ! velocity magnitude squared
    real(rtype) :: th
    !
    ! PBL height calculation:  Scan upward until the Richardson number between
    ! the first level and the current level exceeds the "critical" value.
    !
    do k=nlev-1,nlev-npbl+1,-1
       do i=1,shcol
          if (check(i)) then
             vvk = (u(i,k) - u(i,nlev))**2 + (v(i,k) - v(i,nlev))**2 +fac*ustar(i)**2
             vvk = max(vvk,tiny)
             rino(i,k) = ggr*(thv(i,k) -thv_ref(i))*(z(i,k)-z(i,nlev))/(thv(i,nlev)*vvk)
             if (rino(i,k) >= ricr) then
                pblh(i) = z(i,k+1) + (ricr - rino(i,k+1))/(rino(i,k) -rino(i,k+1)) * &
                     (z(i,k) - z(i,k+1))
                check(i) = .false.
             end if
          end if
       end do
    end do
    return
end subroutine pblintd_height

subroutine pblintd_surf_temp(&
       shcol,nlev,nlevi,&          ! Input
       z,ustar,obklen,kbfs,thv,&   ! Input
       tlv,&                       ! Output
       pblh,check,rino)            ! InOutput
    !------------------------------Arguments--------------------------------
    ! Input arguments
    !
    integer, intent(in) :: shcol                     ! number of atmospheric columns
    integer, intent(in) :: nlev                      ! number of mid-point layers
    integer, intent(in) :: nlevi                     ! number of interface layers

    real(rtype), intent(in)  :: z(shcol,nlev)           ! height above surface [m]
    real(rtype), intent(in)  :: ustar(shcol)            ! surface friction velocity [m/s]
    real(rtype), intent(in)  :: obklen(shcol)           ! Obukhov length
    real(rtype), intent(in)  :: kbfs(shcol)             ! sfc kinematic buoyancy flux [m^2/s^3]
    real(rtype), intent(in) :: thv(shcol,nlev)          ! virtual potential temperature

    real(rtype), intent(out) :: tlv(shcol)              ! ref. level pot tmp + tmp excess
    logical, intent(inout)  :: check(shcol)             ! True=>chk if Richardson no.>critcal
    real(rtype), intent(inout) :: rino(shcol,nlev)      ! bulk Richardson no. from level to ref lev
    real(rtype), intent(inout) :: pblh(shcol)              ! boundary-layer height [m]
    !
    !---------------------------Local workspace-----------------------------
    !
    real(rtype) :: phiminv
    integer  :: i                       ! longitude index

    !===================
    ! const parameter for Diagnosis of PBL depth
    real(rtype), parameter :: onet  = 1._rtype/3._rtype  ! 1/3 power in wind gradient expression
    real(rtype), parameter :: fak   =  8.5_rtype      ! Constant in surface temperature excess
    real(rtype), parameter :: betam = 15.0_rtype      ! Constant in wind gradient expression
    real(rtype), parameter :: sffrac=  0.1_rtype      ! Surface layer fraction of boundary layer
    real(rtype), parameter :: binm  = betam*sffrac ! betam * sffrac

    !
    ! Estimate an effective surface temperature to account for surface
    ! fluctuations
    !
    do i=1,shcol
       if (check(i)) pblh(i) = z(i,nlevi-npbl)
       check(i)  = (kbfs(i) > 0._rtype)
       if (check(i)) then
          phiminv      = (1._rtype - binm*pblh(i)/obklen(i))**onet
          rino(i,nlev) = 0.0_rtype
          tlv(i)       = thv(i,nlev) + kbfs(i)*fak/( ustar(i)*phiminv )
       end if
    end do
    return
end subroutine pblintd_surf_temp

subroutine pblintd_check_pblh(&
       shcol,nlev,nlevi,&             ! Input
       z,ustar,check,&                ! Input
       pblh)                          ! Output
    !------------------------------Arguments--------------------------------
    ! Input arguments
    !
    integer, intent(in) :: shcol                     ! number of atmospheric columns
    integer, intent(in) :: nlev                      ! number of mid-point layers
    integer, intent(in) :: nlevi                     ! number of interface layers

    real(rtype), intent(in)  :: z(shcol,nlev)           ! height above surface [m]
    real(rtype), intent(in)  :: ustar(shcol)            ! surface friction velocity [m/s]
    logical, intent(in)      :: check(shcol)            ! True=>chk if Richardson no.>critcal
    !
    ! Output arguments
    !
    real(rtype), intent(out) :: pblh(shcol)             ! boundary-layer height [m]
    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i                       ! longitude index
    !
    ! PBL height must be greater than some minimum mechanical mixing depth
    ! Several investigators have proposed minimum mechanical mixing depth
    ! relationships as a function of the local friction velocity, u*.  We
    ! make use of a linear relationship of the form h = c u* where c=700.
    ! The scaling arguments that give rise to this relationship most often
    ! represent the coefficient c as some constant over the local coriolis
    ! parameter.  Here we make use of the experimental results of Koracin
    ! and Berkowicz (1988) [BLM, Vol 43] for which they recommend 0.07/f
    ! where f was evaluated at 39.5 N and 52 N.  Thus we use a typical mid
    ! latitude value for f so that c = 0.07/f = 700.  Also, do not allow
    ! PBL to exceed some maximum (npbl) number of allowable points
    !
    do i=1,shcol
       if (check(i)) pblh(i) = z(i,nlevi-npbl)
       pblh(i) = max(pblh(i),700.0_rtype*ustar(i))
    end do
    return
end subroutine pblintd_check_pblh


subroutine pblintd_cldcheck(      &
                   shcol,nlev,nlevi, &                  ! Input
                   zi,cldn,          &                  ! Input
                   pblh)                                ! InOutput
    !------------------------------Arguments--------------------------------
    ! Input arguments
    !
    integer, intent(in) :: shcol                     ! number of atmospheric columns
    integer, intent(in) :: nlev                      ! number of mid-point layers
    integer, intent(in) :: nlevi                     ! number of interface layers

    real(rtype), intent(in)  :: zi(shcol,nlevi)         ! height above surface [m]
    real(rtype), intent(in)  :: cldn(shcol,nlev)        ! new cloud fraction
    !
    ! In/Output arguments
    !
    real(rtype), intent(inout) :: pblh(shcol)             ! boundary-layer height [m]
    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i                       ! longitude index
    logical  :: cldcheck(shcol)      ! True=>if cloud in lowest layer
    !
    ! Final requirement on PBL heightis that it must be greater than the depth
    ! of the lowest model level if there is any cloud diagnosed in
    ! the lowest model level.  This is to deal with the inadequacies of the
    ! current "dry" formulation of the boundary layer, where this test is
    ! used to identify circumstances where there is marine stratus in the
    ! lowest level, and to provide a weak ventilation of the layer to avoid
    ! a pathology in the cloud scheme (locking in low-level stratiform cloud)
    ! If  any cloud is diagnosed in the lowest level, set pblh to 50 meters
    ! higher than top interface of lowest level
    !
    do i=1,shcol
       cldcheck(i) = .false.
       if (cldn(i,nlev).ge.0.0_rtype) cldcheck(i) = .true.
       if (cldcheck(i)) pblh(i) = max(pblh(i),zi(i,nlev) + 50._rtype)
    end do
    return
end subroutine pblintd_cldcheck

  !==============================================================
  ! Linear interpolation to get values on various grids

subroutine linear_interp(x1,x2,y1,y2,km1,km2,ncol,minthresh)
    implicit none

    integer, intent(in) :: km1, km2
    integer, intent(in) :: ncol
    real(rtype), intent(in) :: x1(ncol,km1), y1(ncol,km1)
    real(rtype), intent(in) :: x2(ncol,km2)
    real(rtype), intent(in) :: minthresh
    real(rtype), intent(out) :: y2(ncol,km2)

    integer :: k1, k2, i

#if 1
    !i = check_grid(x1,x2,km1,km2,ncol)
    if (km1 .eq. km2+1) then
       do k2 = 1,km2
          k1 = k2+1
          do i = 1,ncol
             y2(i,k2) = y1(i,k1-1) + (y1(i,k1)-y1(i,k1-1))*(x2(i,k2)-x1(i,k1-1))/(x1(i,k1)-x1(i,k1-1))
          end do
       end do
    elseif (km2 .eq. km1+1) then
       k2 = 1
       do i = 1,ncol
          y2(i,k2) = y1(i,1) + (y1(i,2)-y1(i,1))*(x2(i,k2)-x1(i,1))/(x1(i,2)-x1(i,1))
       end do
       do k2 = 2, km2-1
          k1 = k2
          do i = 1,ncol
             y2(i,k2) = y1(i,k1-1) + (y1(i,k1)-y1(i,k1-1))*(x2(i,k2)-x1(i,k1-1))/(x1(i,k1)-x1(i,k1-1))
          end do
       end do
       k2 = km2
       do i = 1,ncol
          y2(i,k2) = y1(i,km1) + (y1(i,km1)-y1(i,km1-1))*(x2(i,k2)-x1(i,km1))/(x1(i,km1)-x1(i,km1-1))
       end do
    else
       print *,km1,km2
    end if
    do k2 = 1,km2
       do i = 1,ncol
          if (y2(i,k2) .lt. minthresh) then
             y2(i,k2) = minthresh
          endif
       end do
    end do
#else
    do i=1,ncol
       do k2=1,km2
          if( x2(i,k2) <= x1(i,1) ) then
             y2(i,k2) = y1(i,1) + (y1(i,2)-y1(i,1))*(x2(i,k2)-x1(i,1))/(x1(i,2)-x1(i,1))
          elseif( x2(i,k2) >= x1(i,km1) ) then
             y2(i,k2) = y1(i,km1) + (y1(i,km1)-y1(i,km1-1))*(x2(i,k2)-x1(i,km1))/(x1(i,km1)-x1(i,km1-1))
          else
             do k1 = 2,km1
                if( (x2(i,k2)>=x1(i,k1-1)).and.(x2(i,k2)<x1(i,k1)) ) then
                   y2(i,k2) = y1(i,k1-1) + (y1(i,k1)-y1(i,k1-1))*(x2(i,k2)-x1(i,k1-1))/(x1(i,k1)-x1(i,k1-1))
                endif
             enddo ! end k1 loop
          endif

          if (y2(i,k2) .lt. minthresh) then
             y2(i,k2) = minthresh
          endif

       enddo ! end k2 loop
    enddo ! i loop
#endif

    return

end subroutine linear_interp

subroutine compute_brunt_shoc_length(nlev,nlevi,shcol,dz_zt,thv,thv_zi,brunt)

  !=========================================================
  !
  ! Computes the brunt_visala frequency

  implicit none
  integer, intent(in) :: nlev, nlevi, shcol
  ! Grid difference centereted on thermo grid [m]
  real(rtype), intent(in) :: dz_zt(shcol,nlev)
  ! virtual potential temperature [K]
  real(rtype), intent(in) :: thv(shcol,nlev)
  ! virtual potential temperature [K] at interface
  real(rtype), intent(in) :: thv_zi(shcol,nlevi)
  ! brunt vaisala frequency [s-1]
  real(rtype), intent(out) :: brunt(shcol, nlev)
  integer k, i

  do k=1,nlev
    do i=1,shcol
      brunt(i,k) = (ggr/thv(i,k)) * (thv_zi(i,k) - thv_zi(i,k+1))/dz_zt(i,k)
    enddo
  enddo

end subroutine compute_brunt_shoc_length

subroutine compute_l_inf_shoc_length(nlev,shcol,zt_grid,dz_zt,tke,l_inf)

  !=========================================================
  !

  implicit none
  integer, intent(in) :: nlev, shcol
  real(rtype), intent(in) :: zt_grid(shcol,nlev), dz_zt(shcol,nlev), tke(shcol,nlev)
  real(rtype), intent(inout) :: l_inf(shcol)
  real(rtype) :: tkes, numer(shcol), denom(shcol)
  integer k, i

  numer(:) = 0._rtype
  denom(:) = 0._rtype

  do k=1,nlev
    do i=1,shcol
        tkes=sqrt(tke(i,k))
        numer(i)=numer(i)+tkes*zt_grid(i,k)*dz_zt(i,k)
        denom(i)=denom(i)+tkes*dz_zt(i,k)
    enddo
  enddo

  do i=1,shcol
    l_inf(i)=0.1_rtype*(numer(i)/denom(i))
  enddo

end subroutine compute_l_inf_shoc_length

subroutine compute_conv_vel_shoc_length(nlev,shcol,pblh,zt_grid,dz_zt,thv,wthv_sec,conv_vel)

  !=========================================================
  ! determine the convective velocity scale of
  !   the planetary boundary layer

  implicit none
  integer, intent(in) :: nlev, shcol
! Planetary boundary layer (PBL) height [m]
  real(rtype), intent(in) :: pblh(shcol)
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  real(rtype), intent(in) :: dz_zt(shcol,nlev)
  real(rtype), intent(in) :: thv(shcol,nlev)
  real(rtype), intent(in) :: wthv_sec(shcol,nlev)
  real(rtype), intent(inout) :: conv_vel(shcol)
  integer k, i
  conv_vel(:) = 0._rtype

  do k=nlev,1,-1
    do i=1,shcol
      if (zt_grid(i,k) .lt. pblh(i)) then
        conv_vel(i) = conv_vel(i)+2.5_rtype*dz_zt(i,k)*(ggr/thv(i,k))*wthv_sec(i,k)
      endif
    enddo ! end i loop (column loop)
  enddo ! end k loop (vertical loop)

end subroutine compute_conv_vel_shoc_length

subroutine compute_conv_time_shoc_length(shcol,pblh,conv_vel,tscale)
  ! Compute eddy turnover timescale.  If
  !  convective velocity scale is zero then
  !  set to a minimum threshold

  implicit none
  integer, intent(in) :: shcol
! Planetary boundary layer (PBL) height [m]
  real(rtype), intent(in) :: pblh(shcol)
  ! Convective velocity scale
  real(rtype), intent(inout) :: conv_vel(shcol)
  ! Convective time scale
  real(rtype), intent(inout) ::  tscale(shcol)

  integer i

  do i=1,shcol
    conv_vel(i) = max(0._rtype,conv_vel(i))**(1._rtype/3._rtype)

    if (conv_vel(i) .gt. 0._rtype) then
      tscale(i)=pblh(i)/conv_vel(i)
    else
      tscale(i)=100._rtype
    endif
  enddo

end subroutine compute_conv_time_shoc_length

subroutine compute_shoc_mix_shoc_length(nlev,shcol,tke,brunt,tscale,zt_grid,l_inf,shoc_mix)

#ifdef SCREAM_CONFIG_IS_CMAKE
  use shoc_iso_f, only: compute_shoc_mix_shoc_length_f
#endif

  implicit none

  integer, intent(in) :: nlev, shcol
  ! turbulent kinetic energy [m^2/s^2]
  real(rtype), intent(in) :: tke(shcol,nlev)
  ! brunt vaisala frequency [s-1]
  real(rtype), intent(in) :: brunt(shcol,nlev)
  ! convective time scale
  real(rtype), intent(in) :: tscale(shcol)
  ! heights, for thermo grid [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  real(rtype), intent(in) :: l_inf(shcol)

  ! Turbulent length scale [m]
  real(rtype), intent(out) :: shoc_mix(shcol,nlev)

  !  LOCAL VARIABLES
  real(rtype) :: brunt2(shcol,nlev)
  integer k, i
  real(rtype) :: tkes

#ifdef SCREAM_CONFIG_IS_CMAKE
  if (use_cxx) then
    call compute_shoc_mix_shoc_length_f(nlev,shcol,tke,brunt,tscale,zt_grid,l_inf,& !Input
                                        shoc_mix) ! Ouptut
    return
  endif
#endif

  brunt2(:,:) = 0.0

  do k=1,nlev
    do i=1,shcol

      tkes = sqrt(tke(i,k))

      if(brunt(i,k) .ge. 0) brunt2(i,k) = brunt(i,k)

      shoc_mix(i,k)=min(maxlen,(2.8284_rtype*sqrt(1._rtype/((1._rtype/(tscale(i)*tkes*vk*zt_grid(i,k)))&
        +(1._rtype/(tscale(i)*tkes*l_inf(i)))+0.01_rtype*(brunt2(i,k)/tke(i,k)))))/length_fac)
    enddo ! end i loop (column loop)
  enddo ! end k loop (vertical loop)

end subroutine compute_shoc_mix_shoc_length

subroutine check_length_scale_shoc_length(nlev,shcol,host_dx,host_dy,shoc_mix)
  ! Do checks on the length scale.  Make sure it is not
  !  larger than the grid mesh of the host model.

  implicit none
  integer, intent(in) :: nlev, shcol
  real(rtype), intent(in) :: host_dx(shcol), host_dy(shcol)
  ! Turbulent length scale [m]
  real(rtype), intent(inout) :: shoc_mix(shcol, nlev)
  integer k, i

  do k=1,nlev
    do i=1,shcol
      shoc_mix(i,k)=min(maxlen,shoc_mix(i,k))
      shoc_mix(i,k)=max(minlen,shoc_mix(i,k))
      shoc_mix(i,k)=min(sqrt(host_dx(i)*host_dy(i)),shoc_mix(i,k))
    enddo
  enddo

end subroutine check_length_scale_shoc_length

end module

!==============================================================
! This is the end of the SHOC parameterization
! We hope you have enjoyed your time here
