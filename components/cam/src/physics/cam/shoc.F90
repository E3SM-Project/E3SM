module shoc

!--------------------------------------------------------------
! SHOC parameterization
!   SHOC = Simplified Higher Order Closure 
!   reference, Bogenschutz and Krueger 2013
!
! PDF-based parameterization for low clouds and turbulence   
! email: bogenschutz1@llnl.gov
!--------------------------------------------------------------

implicit none

! define r8 for double precision 
integer, parameter, public :: r8 = selected_real_kind(12)
real(r8), parameter, public :: largeneg = -99999999.99_r8

!=========================================================
! Physical constants used in SHOC 
!=========================================================

! These are set in initialization and should be set to 
!  to the values used in whatever host model SHOC is
!  implemented in
real(r8) :: ggr   ! gravity [m/s^2]
real(r8) :: rgas  ! dry air gas constant [J/kg.K]
real(r8) :: rv    ! water vapor gas constant [J/kg.K]
real(r8) :: cp    ! specific heat of dry air [J/kg.K]
real(r8) :: lcond ! latent heat of vaporization [J/kg]

!=========================================================
! Private module parameters
!=========================================================

! These are adjustable tunable parameters for SHOC
!  for certain variables and turbulent moments.
!  All are unitless

! temperature variance
real(r8), parameter :: thl2tune=1.0_r8 
! moisture variance
real(r8), parameter :: qw2tune=1.0_r8 
! temp moisture covariance
real(r8), parameter :: qwthl2tune=1.0_r8
! vertical velocity variance 
real(r8), parameter :: w2tune=1.0_r8 
! third moment of vertical velocity
real(r8), parameter :: w3clip=1.2_r8 

! =========
! Below are options to activate certain features in SHOC

! Allow temperature skewness to be independent of moisture 
!  variance 
logical, parameter :: dothetal_skew = .false.

! Use an implicit diffusion solver for SHOC
!  If running with long timesteps (dt > 20 s), then 
!  this should be set to true. If set to false then 
!  an explicit solver will be used 
logical, parameter :: do_implicit = .true.

! ========
! Below define some parameters for SHOC

! Reference temperature [K]
real(r8), parameter :: basetemp = 300._r8
! Reference pressure [Pa]
real(r8), parameter :: basepres = 100000._r8

! ========
! Set upper limits for certain SHOC quantities
! Note that these upper limits are quite high 
! and they are rarely reached in a stable simulation

! Return to isotropic timescale [s]
real(r8), parameter :: maxiso = 20000.0_r8
! Mixing length [m]
real(r8), parameter :: maxlen = 20000.0_r8
! Maximum TKE [m2/s2]
real(r8), parameter :: maxtke = 50.0_r8
! Minimum TKE [m2/s2]
real(r8), parameter :: mintke = 0.0004_r8

!==============================================================
! Begin SHOC parameterization code!
contains
!==============================================================

subroutine shoc_init( &
             gravit, rair, rh2o, cpair, &
	     latvap)

  implicit none
  
  ! Purpose:  Initialize constants for SHOC
  ! These should be set to the same values used in 
  ! whatever host model SHOC is implemented in
	     
  real(r8), intent(in)  :: gravit ! gravity
  real(r8), intent(in)  :: rair   ! dry air gas constant 
  real(r8), intent(in)  :: rh2o   ! water vapor gas constant
  real(r8), intent(in)  :: cpair  ! specific heat of dry air
  real(r8), intent(in)  :: latvap ! latent heat of vaporization
  
  ggr = gravit   ! [m/s2] 
  rgas = rair    ! [J/kg.K]
  rv = rh2o      ! [J/kg.K]
  cp = cpair     ! [J/kg.K]
  lcond = latvap ! [J/kg]
  
end subroutine shoc_init	     				   

!==============================================================
! Main driver for the SHOC scheme
! Host models should call the following routine to call SHOC 

subroutine shoc_main ( &
     shcol, nlev, nlevi, dtime, &         ! Input
     host_dx, host_dy,thv, cldliq, &      ! Input
     zt_grid,zi_grid,pres,pdel,&          ! Input
     wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, & ! Input
     wtracer_sfc,num_qtracers,w_field, &  ! Input     
     tke, thetal, qw, &                   ! Input/Output
     u_wind, v_wind,qtracers,&            ! Input/Output
     wthv_sec,tkh,tk,&                    ! Input/Output
     shoc_cldfrac,shoc_ql,&               ! Output
     shoc_mix, isotropy,&                 ! Output (diagnostic)
     w_sec, thl_sec, qw_sec, qwthl_sec,&  ! Output (diagnostic)
     wthl_sec, wqw_sec, wtke_sec,&        ! Output (diagnostic)
     uw_sec, vw_sec, w3,&                 ! Output (diagnostic)    
     wqls_sec, brunt)                     ! Output (diagnostic)

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
  
  ! SHOC timestep [s]
  real(r8), intent(in) :: dtime	
  ! grid spacing of host model in x direction [m]
  real(r8), intent(in) :: host_dx(shcol)
  ! grid spacing of host model in y direction [m] 
  real(r8), intent(in) :: host_dy(shcol) 
  ! heights, for thermo grid [m]
  real(r8), intent(in) :: zt_grid(shcol,nlev)	
  ! heights, for interface grid [m]
  real(r8), intent(in) :: zi_grid(shcol,nlevi) 
  ! pressure levels on thermo grid [Pa]  
  real(r8), intent(in) :: pres(shcol,nlev)
  ! Differences in pressure levels [Pa] 
  real(r8), intent(in) :: pdel(shcol,nlev)
  ! virtual potential temperature [K] 
  real(r8), intent(in) :: thv(shcol,nlev) 
  ! cloud liquid mixing ratio [kg/kg]
  real(r8), intent(in) :: cldliq(shcol,nlev) 
  ! large scale vertical velocity [m/s]
  real(r8), intent(in) :: w_field(shcol,nlev) 
  ! Surface sensible heat flux [K m/s]
  real(r8), intent(in) :: wthl_sfc(shcol) 
  ! Surface latent heat flux [kg/kg m/s]
  real(r8), intent(in) :: wqw_sfc(shcol)
  ! Surface momentum flux (u-direction) [m2/s2] 
  real(r8), intent(in) :: uw_sfc(shcol) 
  ! Surface momentum flux (v-direction) [m2/s2]
  real(r8), intent(in) :: vw_sfc(shcol) 
  ! Surface flux for tracers [varies]
  real(r8), intent(in) :: wtracer_sfc(shcol,num_qtracers) 

! INPUT/OUTPUT VARIABLES  
  ! turbulent kinetic energy [m2/s2]
  real(r8), intent(inout) :: tke(shcol,nlev)  
  ! liquid water potential temperature [K]
  real(r8), intent(inout) :: thetal(shcol,nlev) 
  ! total water mixing ratio [kg/kg]
  real(r8), intent(inout) :: qw(shcol,nlev) 
  ! u wind component [m/s]					    
  real(r8), intent(inout) :: u_wind(shcol,nlev)
  ! v wind component [m/s] 
  real(r8), intent(inout) :: v_wind(shcol,nlev) 
  ! buoyancy flux [K m/s]
  real(r8), intent(inout) :: wthv_sec(shcol,nlev) 
  ! tracers [varies]
  real(r8), intent(inout) :: qtracers(shcol,nlev,num_qtracers) 
  ! eddy coefficient for momentum [m2/s]
  real(r8), intent(inout) :: tk(shcol,nlev) 
  ! eddy coefficent for heat [m2/s]
  real(r8), intent(inout) :: tkh(shcol,nlev)

! OUTPUT VARIABLES
  ! Cloud fraction [-]
  real(r8), intent(out) :: shoc_cldfrac(shcol,nlev)
  ! cloud liquid mixing ratio [kg/kg] 
  real(r8), intent(out) :: shoc_ql(shcol,nlev) 
  
  ! also output variables, but part of the SHOC diagnostics
  !  to be output to history file by host model (if desired)
  
  ! Turbulent length scale [m]  
  real(r8) :: shoc_mix(shcol,nlev)
  ! vertical velocity variance [m2/s2] 
  real(r8) :: w_sec(shcol,nlev) 
  ! temperature variance [K^2]
  real(r8) :: thl_sec(shcol,nlevi) 
  ! moisture variance [kg2/kg2]
  real(r8) :: qw_sec(shcol,nlevi) 
  ! temp moisture covariance [K kg/kg]
  real(r8) :: qwthl_sec(shcol,nlevi)
  ! vertical heat flux [K m/s] 
  real(r8) :: wthl_sec(shcol,nlevi)
  ! vertical moisture flux [K m/s] 
  real(r8) :: wqw_sec(shcol,nlevi) 
  ! vertical tke flux [m3/s3]
  real(r8) :: wtke_sec(shcol,nlevi)
  ! vertical zonal momentum flux [m2/s2] 
  real(r8) :: uw_sec(shcol,nlevi) 
  ! vertical meridional momentum flux [m2/s2]
  real(r8) :: vw_sec(shcol,nlevi)
  ! third moment vertical velocity [m3/s3] 
  real(r8) :: w3(shcol,nlevi) 
  ! liquid water flux [kg/kg m/s]
  real(r8) :: wqls_sec(shcol,nlev)
  ! brunt vaisala frequency [s-1] 
  real(r8) :: brunt(shcol,nlev) 
  ! return to isotropic timescale [s]
  real(r8) :: isotropy(shcol,nlev) 

  !============================================================================
! LOCAL VARIABLES
  
  ! vertical flux of tracers [varies]
  real(r8) :: wtracer_sec(shcol,nlevi,num_qtracers) 
  ! air density on thermo grid [kg/m3]
  real(r8) :: rho_zt(shcol,nlev) 
 
  ! Grid difference centereted on thermo grid [m] 
  real(r8) :: dz_zt(shcol,nlev)
  ! Grid difference centereted on interface grid [m] 
  real(r8) :: dz_zi(shcol,nlev) 

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
         dz_zt,dz_zi,rho_zt)  		      ! Output

  ! Update the turbulent length scale	 
  call shoc_length(&
         shcol,nlev,nlevi, tke,&              ! Input
         host_dx, host_dy, cldliq,&           ! Input
         zt_grid,zi_grid,dz_zt,dz_zi,&        ! Input
	 thetal,wthv_sec,thv,&                ! Input
	 brunt,shoc_mix)  		      ! Output

  ! Advance the SGS TKE equation	 
  call shoc_tke(&
         shcol,nlev,nlevi,dtime,&             ! Input
         wthv_sec,shoc_mix,&                  ! Input
	 dz_zi,dz_zt,pres,&                   ! Input
	 u_wind,v_wind,brunt,&                ! Input
	 uw_sfc,vw_sfc,&                      ! Input
	 zt_grid,zi_grid,&                    ! Input
         tke,tk,tkh,&		              ! Input/Output
	 isotropy)                            ! Output
  
  ! If implicit diffusion solver is used, 
  !  update SHOC prognostic variables here
  if (do_implicit) then
    call update_prognostics_implicit(&        ! Input
           shcol,nlev,nlevi,num_qtracers,&    ! Input
	   dtime,dz_zt,dz_zi,rho_zt,&         ! Input
	   zt_grid,zi_grid,tk,tkh,&           ! Input
           uw_sfc,vw_sfc,wthl_sfc,wqw_sfc,&   ! Input
           thetal,qw,qtracers,tke,&           ! Input/Output
	   u_wind,v_wind)                     ! Input/Output
  endif	 

  ! Diagnose the second order moments, needed
  !  for the PDF closure
  call diag_second_shoc_moments(&
         shcol,nlev,nlevi,&                   ! Input
         num_qtracers,thetal,qw,&             ! Input
	 u_wind,v_wind,qtracers,tke,&         ! Input
	 isotropy,tkh,tk,&                    ! Input
	 dz_zi,zt_grid,zi_grid,&              ! Input
	 wthl_sfc,wqw_sfc,uw_sfc,vw_sfc,&     ! Input
	 wtracer_sfc,shoc_mix,&               ! Input
         w_sec,thl_sec,qw_sec,&               ! Output
	 wthl_sec,wqw_sec,&                   ! Output
	 qwthl_sec,uw_sec,vw_sec,wtke_sec,&   ! Output
	 wtracer_sec)

  ! Diagnose the third moment of vertical velocity, 
  !  needed for the PDF closure	 
  call diag_third_shoc_moments(&
         shcol,nlev,nlevi,&                   ! Input
         w_sec,thl_sec,qw_sec,qwthl_sec,&     ! Input
	 wthl_sec,isotropy,brunt,&            ! Input
	 thetal,tke,wthv_sec,shoc_mix,&       ! Input
	 dz_zt,dz_zi,&                        ! Input
	 zt_grid,zi_grid,&                    ! Input
	 w3)                                  ! Output
  
  ! Update thetal, qw, tracers, and wind components
  !   based on SGS mixing, if explicit scheme is used
  if (.not. do_implicit) then
    call update_prognostics(&
           shcol,nlev,nlevi,num_qtracers,&    ! Input
           dtime,dz_zt,wthl_sec,&             ! Input
           wqw_sec,wtke_sec,uw_sec,&          ! Input
	   vw_sec,wtracer_sec,&               ! Input
           rho_zt,zt_grid,zi_grid,&           ! Input
           thetal,qw,qtracers,tke,&           ! Input/Output	
	   u_wind,v_wind)                     ! Input/Output
  endif
	 
  ! Call the PDF to close on SGS cloud and turbulence
  call shoc_assumed_pdf(&
         shcol,nlev,nlevi,&                   ! Input
         thetal,qw,w_field,thl_sec,qw_sec,&   ! Input
	 wthl_sec,w_sec,&                     ! Input
	 wqw_sec,qwthl_sec,w3,pres,&          ! Input
	 zt_grid,zi_grid,&                    ! Input
	 shoc_cldfrac,shoc_ql,&               ! Output
         wqls_sec,wthv_sec)                   ! Output

  ! Check TKE to make sure values lie within acceptable 
  !  bounds after vertical advection, etc.
  call check_tke(shcol,nlev,tke)
	
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
  real(r8), intent(in) :: zt_grid(shcol,nlev) 
  ! interface grid heights [m] 
  real(r8), intent(in) :: zi_grid(shcol,nlevi) 
  ! pressure differences centered on mid-point grid [Pa]
  real(r8), intent(in) :: pdel(shcol,nlev) 

! OUTPUT VARIABLES
  ! thickness (dz) on the thermo grid [m] 
  real(r8), intent(out) :: dz_zt(shcol,nlev)
  ! thickness (dz) on the interface grid [m] 
  real(r8), intent(out) :: dz_zi(shcol,nlev)
  ! air density on the thermo grid [kg/m3] 
  real(r8), intent(out) :: rho_zt(shcol,nlev)  

  ! local variables
  integer :: i, k
 
  do k=1,nlev
    do i=1,shcol
      ! define thickness of the thermodynamic gridpoints
      dz_zt(i,k) = zi_grid(i,k) - zi_grid(i,k+1)
 
      ! define thickness of the interface grid points
      if (k .eq. nlev) then
        dz_zi(i,k) = zt_grid(i,k)
      else
        dz_zi(i,k) = zt_grid(i,k) - zt_grid(i,k+1)
      endif

      ! Define the air density on the thermo grid
      rho_zt(i,k) = (1._r8/ggr)*(pdel(i,k)/dz_zt(i,k))
      
    enddo ! end i loop (column loop)
  enddo ! end k loop (vertical loop)
  
  return
  
end subroutine shoc_grid

!==============================================================
! Update T, q, tracers, tke, u, and v based on SGS mixing
! using explicit diffusion solver.  Note that this routine 
! should only be called if using very small time steps 
! (< 20 s), otherwise the implicit diffusion solver should be used

subroutine update_prognostics( &
             shcol,nlev,nlevi,num_tracer,&    ! Input
             dtime,dz_zt,wthl_sec,&           ! Input
	     wqw_sec,wtke_sec,uw_sec,&        ! Input
	     vw_sec,wtracer_sec,&             ! Input
             rho_zt,zt_grid,zi_grid,&         ! Input
             thetal,qw,tracer,tke,&           ! Input/Output
	     u_wind,v_wind)                   ! Input/Output

! Purpose of this subroutine is to update T, q, u, v, tke, and 
!  tracers based on SGS mixing due to SHOC, using 
!  explicit diffusion solver
	     
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
  ! time step [s]
  real(r8), intent(in) :: dtime
  ! thickness of grid centered on thermo points [m]
  real(r8), intent(in) :: dz_zt(shcol,nlev)
  ! vertical flux of heat [K m/s]
  real(r8), intent(in) :: wthl_sec(shcol,nlevi)
  ! vertical flux of moisture [kg/kg m/s]
  real(r8), intent(in) :: wqw_sec(shcol,nlevi)
  ! vertical flux of tracers [varies]
  real(r8), intent(in) :: wtracer_sec(shcol,nlevi,num_tracer)
  ! vertical zonal momentum flux [m2/s2]
  real(r8), intent(in) :: uw_sec(shcol,nlevi)
  ! vertical meridional momentum flux [m2/s2]
  real(r8), intent(in) :: vw_sec(shcol,nlevi)
  ! vertical flux of TKE [m3/s3]
  real(r8), intent(in) :: wtke_sec(shcol,nlevi)
  ! air density [kg/m3]
  real(r8), intent(in) :: rho_zt(shcol,nlev)
  ! heights centered on thermo points [m]
  real(r8), intent(in) :: zt_grid(shcol,nlev)
  ! heights centered on interface points [m]
  real(r8), intent(in) :: zi_grid(shcol,nlevi)

! IN/OUT VARIABLES
  ! liquid water potential temperature [K]  
  real(r8), intent(inout) :: thetal(shcol,nlev)
  ! total water mixing ratio [kg/kg]
  real(r8), intent(inout) :: qw(shcol,nlev)   
  ! tracers [varies]
  real(r8), intent(inout) :: tracer(shcol,nlev,num_tracer)
  ! zonal wind [m/s]
  real(r8), intent(inout) :: u_wind(shcol,nlev)
  ! meridional wind [m/s]
  real(r8), intent(inout) :: v_wind(shcol,nlev)
  ! turbulent kinetic energy [m2/s2]
  real(r8), intent(inout) :: tke(shcol,nlev)
  
! LOCAL VARIABLES  
  integer :: kb, kt, k, i, p  
  real(r8) :: thedz, r1, r2, r3
  real(r8) :: rho_zi(shcol,nlevi)
 
  ! linearly interpolate air density from thermo to interface grid
  call linear_interp(zt_grid,zi_grid,rho_zt,rho_zi,nlev,nlevi,shcol,0._r8)

  do i=1,shcol
    do k=1,nlev 
      kb=k+1

      ! define air densities on various levels for mass weighted 
      !  diffusion for conservation of mass
      r1=rho_zi(i,k)
      r2=rho_zi(i,kb)
      r3=rho_zt(i,k)

      ! mass weighted 1/dz
      thedz=1._r8/(dz_zt(i,k)*r3)

      ! Update temperature via vertical diffusion
      thetal(i,k)=thetal(i,k)-dtime*(r1*wthl_sec(i,k)-r2*wthl_sec(i,kb))*thedz
     
      ! Update total water mixing ratio via vertical diffusion
      qw(i,k)=qw(i,k)-dtime*(r1*wqw_sec(i,k)-r2*wqw_sec(i,kb))*thedz

      ! Update turbulent kinetic energy via vertical diffusion
      tke(i,k)=tke(i,k)-dtime*(r1*wtke_sec(i,k)-r2*wtke_sec(i,kb))*thedz

      ! Update tracers via vertical diffusion
      do p=1,num_tracer
        tracer(i,k,p)=tracer(i,k,p)-dtime*(r1*wtracer_sec(i,k,p)-r2*wtracer_sec(i,kb,p))*thedz
      enddo

      ! Update the u and v wind components via vertical diffusion
      u_wind(i,k)=u_wind(i,k)-dtime*(r1*uw_sec(i,k)-r2*uw_sec(i,kb))*thedz
      v_wind(i,k)=v_wind(i,k)-dtime*(r1*vw_sec(i,k)-r2*vw_sec(i,kb))*thedz

    enddo ! end i loop (column loop)
  enddo ! end k loop (vertical loop)

  return
  
end subroutine update_prognostics

!==============================================================
! Update T, q, tracers, tke, u, and v based on implicit diffusion
! If running with time steps longer than ~ 20 s then to preserve
! numerical stability, an implicit diffusion solver will need to be
! used.  Here we use a backward Euler scheme.  This is the default
! diffusion solver for SHOC.  Switching to an explicit scheme is possible
! by setting do_implicit = .false.

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
  real(r8), intent(in) :: dtime
  ! Eddy coefficient for momentum [m2/s]
  real(r8), intent(in) :: tk(shcol,nlev)
  ! Eddy coefficient for heat [m2/s]
  real(r8), intent(in) :: tkh(shcol,nlev)
  ! Air density on thermo grid [kg/m3]
  real(r8), intent(in) :: rho_zt(shcol,nlev)
  ! height thickness centered on thermo grid [m]
  real(r8), intent(in) :: dz_zt(shcol,nlev)
  ! height thickness centered on interface grid [m]
  real(r8), intent(in) :: dz_zi(shcol,nlev)
  ! vertical zonal momentum flux at surface [m3/s3]
  real(r8), intent(in) :: uw_sfc(shcol)
  ! vertical meridional momentum flux at surface [m3/s3]
  real(r8), intent(in) :: vw_sfc(shcol)
  ! vertical heat flux at surface [K m/s]
  real(r8), intent(in) :: wthl_sfc(shcol)
  ! vertical moisture flux at surface [kg/kg m/s]
  real(r8), intent(in) :: wqw_sfc(shcol)
  ! heights of mid-point [m]
  real(r8), intent(in) :: zt_grid(shcol,nlev) 
  ! heights at interfaces [m]
  real(r8), intent(in) :: zi_grid(shcol,nlevi)  
	     
! IN/OUT VARIABLES
  ! liquid water potential temperature [K] 
  real(r8), intent(inout) :: thetal(shcol,nlev)
  ! total water mixing ratio [kg/kg]
  real(r8), intent(inout) :: qw(shcol,nlev) 
  ! tracers [varies]  
  real(r8), intent(inout) :: tracer(shcol,nlev,num_tracer)
  ! zonal wind [m/s]
  real(r8), intent(inout) :: u_wind(shcol,nlev)
  ! meridional wind [m/s]
  real(r8), intent(inout) :: v_wind(shcol,nlev)
  ! turbulent kinetic energy [m2/s2]
  real(r8), intent(inout) :: tke(shcol,nlev)
 
! LOCAL VARIABLES         
  integer :: i, k, p
  real(r8) :: rdp_zt(shcol,nlev)
  real(r8) :: tmpi(shcol,nlevi)
  real(r8) :: tkh_zi(shcol,nlevi)
  real(r8) :: tk_zi(shcol,nlevi)
  real(r8) :: rho_zi(shcol,nlevi)
 
  real(r8) :: flux_dummy(shcol)
  real(r8) :: ws(shcol)
  real(r8) :: tau(shcol), taux(shcol), tauy(shcol)
  real(r8) :: ksrf(shcol), ustar, wtke_flux(shcol)
  
  real(r8) :: ca(shcol,nlev) ! superdiagonal for solver
  real(r8) :: cc(shcol,nlev) ! subdiagonal for solver
  real(r8) :: denom(shcol,nlev) ! denominator in solver
  real(r8) :: ze(shcol,nlev)

  real(r8) :: wsmin, ksrfmin
  real(r8) :: timeres 
 
  wsmin = 1._r8      ! Minimum wind speed for ksrfturb computation [ m/s ]
  ksrfmin = 1.e-4_r8 ! Minimum surface drag coefficient  [ kg/s/m^2 ]
  timeres = 7200._r8 ! Relaxation time scale of residual stress ( >= dt ) [s]

  ! linearly interpolate tkh, tk, and air density onto the interface grids
  call linear_interp(zt_grid,zi_grid,tkh,tkh_zi,nlev,nlevi,shcol,0._r8)
  call linear_interp(zt_grid,zi_grid,tk,tk_zi,nlev,nlevi,shcol,0._r8)
  call linear_interp(zt_grid,zi_grid,rho_zt,rho_zi,nlev,nlevi,shcol,0._r8)
 
  ! Define the tmpi variable, which is really dt*(g*rho)**2/dp 
  !  at interfaces. Sub dp = g*rho*dz
  do k=1,nlev
    do i=1,shcol
       tmpi(i,k) = dtime * (ggr*rho_zi(i,k)) / dz_zi(i,k)
    enddo
  enddo

  ! compute 1/dp term, needed in diffusion solver
  do k=1,nlev
    do i=1,shcol 
      rdp_zt(i,k) = 1._r8/(ggr*rho_zt(i,k)*dz_zt(i,k))
    enddo
  enddo  

  ! define terms needed for the implicit surface stress
  do i=1,shcol
    taux(i) = rho_zi(i,nlevi)*uw_sfc(i) ! stress in N/m2
    tauy(i) = rho_zi(i,nlevi)*vw_sfc(i) ! stress in N/m2
    ! compute the wind speed
    ws(i) = max(sqrt(u_wind(i,nlev)**2._r8 + v_wind(i,nlev)**2._r8),wsmin)
    tau(i) = sqrt( taux(i)**2._r8 + tauy(i)**2._r8 )
    ksrf(i) = max(tau(i) / ws(i), ksrfmin)
    ustar=max(sqrt(sqrt(uw_sfc(i)**2 + vw_sfc(i)**2)),0.01_r8)
    wtke_flux(i) = ustar**3
  enddo

  ! Apply the surface fluxes explicitly for temperature and moisture
  thetal(:,nlev) = thetal(:,nlev) + dtime * (ggr * rho_zi(:,nlev) * rdp_zt(:,nlev)) * wthl_sfc(:)  
  qw(:,nlev) = qw(:,nlev) + dtime * (ggr * rho_zi(:,nlevi) * rdp_zt(:,nlev)) * wqw_sfc(:)
  tke(:,nlev) = tke(:,nlev) + dtime * (ggr * rho_zi(:,nlevi) * rdp_zt(:,nlev)) * wtke_flux(:)

  ! Call decomp for momentum variables
  call vd_shoc_decomp(shcol,nlev,nlevi,tk_zi,tmpi,rdp_zt,dtime,&
	 ksrf,ca,cc,denom,ze)

  ! march u_wind one step forward using implicit solver
  call vd_shoc_solve(shcol,nlev,nlevi,ca,cc,denom,ze,u_wind)

  ! march v_wind one step forward using implicit solver
  call vd_shoc_solve(shcol,nlev,nlevi,ca,cc,denom,ze,v_wind)

! Call decomp for thermo variables
  flux_dummy(:) = 0._r8 ! fluxes applied explicitly, so zero fluxes out
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

!==============================================================
! SHOC Diagnose the second order moments

subroutine diag_second_shoc_moments(&
             shcol,nlev,nlevi, &                    ! Input
             num_tracer,thetal,qw, &                ! Input
	     u_wind,v_wind,tracer,tke, &            ! Input
	     isotropy,tkh,tk,&                      ! Input
	     dz_zi,zt_grid,zi_grid,&                ! Input
	     wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, &   ! Input   
	     wtracer_sfc,shoc_mix, &                ! Input 
             w_sec, thl_sec, qw_sec,&               ! Output
	     wthl_sec,wqw_sec,&                     ! Output
	     qwthl_sec, uw_sec, vw_sec, wtke_sec, & ! Output
	     wtracer_sec)                           ! Output

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
  real(r8), intent(in) :: thetal(shcol,nlev) 
  ! total water mixing ratio [kg/kg]
  real(r8), intent(in) :: qw(shcol,nlev) 
  ! zonal wind component [m/s]
  real(r8), intent(in) :: u_wind(shcol,nlev) 
  ! meridional wind component [m/s]
  real(r8), intent(in) :: v_wind(shcol,nlev) 
  ! turbulent kinetic energy [m2/s2]
  real(r8), intent(in) :: tke(shcol,nlev) 
  ! return to isotropy timescale [s]
  real(r8), intent(in) :: isotropy(shcol,nlev) 
  ! eddy coefficient for heat [m2/s]
  real(r8), intent(in) :: tkh(shcol,nlev) 
  ! eddy coefficient for momentum [m2/s]
  real(r8), intent(in) :: tk(shcol,nlev) 
  ! tracers [varies]
  real(r8), intent(in) :: tracer(shcol,nlev,num_tracer) ! tracers
  ! heights of mid-point grid [m]
  real(r8), intent(in) :: zt_grid(shcol,nlev) 
  ! heights of interface grid [m]
  real(r8), intent(in) :: zi_grid(shcol,nlevi) 
  ! thickness centered on interface grid [m]
  real(r8), intent(in) :: dz_zi(shcol,nlev) 
  
  ! Surface sensible heat flux [K m/s]
  real(r8), intent(in) :: wthl_sfc(shcol)
  ! Surface latent heat flux [kg/kg m/s] 
  real(r8), intent(in) :: wqw_sfc(shcol) 
  ! Surface momentum flux (u-direction) [m3/s3]
  real(r8), intent(in) :: uw_sfc(shcol) 
  ! Surface momentum flux (v-direction) [m3/s3]
  real(r8), intent(in) :: vw_sfc(shcol) 
  ! Tracer flux [varies m/s]
  real(r8), intent(in) :: wtracer_sfc(shcol,num_tracer)
  ! Mixing length [m]
  real(r8), intent(in) :: shoc_mix(shcol,nlev)

! OUTPUT VARIABLES
  ! second order vertical velocity [m2/s2]
  real(r8), intent(out) :: w_sec(shcol,nlev)	
  ! second order liquid wat. potential temp. [K^2]
  real(r8), intent(out) :: thl_sec(shcol,nlevi) 
  ! second order total water mixing rat. [kg^2/kg^2] 
  real(r8), intent(out) :: qw_sec(shcol,nlevi)  
  ! covariance of temp and moisture [K kg/kg] 
  real(r8), intent(out) :: qwthl_sec(shcol,nlevi) 
  ! vertical flux of heat [K m/s]
  real(r8), intent(out) :: wthl_sec(shcol,nlevi)
  ! vertical flux of total water [kg/kg m/s] 
  real(r8), intent(out) :: wqw_sec(shcol,nlevi)
  ! vertical flux of zonal wind [m2/s2] 
  real(r8), intent(out) :: uw_sec(shcol,nlevi) 
  ! vertical flux of meridional wind [m2/s2]
  real(r8), intent(out) :: vw_sec(shcol,nlevi) 
  ! vertical flux of tke [m3/s3]
  real(r8), intent(out) :: wtke_sec(shcol,nlevi) 
  ! vertical flux of tracer [varies m/s]
  real(r8), intent(out) :: wtracer_sec(shcol,nlevi,num_tracer) 

! LOCAL VARIABLES  
  integer :: kb, kt, k, i, p
  real(r8) :: grid_dz2,grid_dz,grid_dzw
  real(r8) :: gr1,gr2,grw1
  real(r8) :: isotropy_zi(shcol,nlevi)
  real(r8) :: tkh_zi(shcol,nlevi)
  real(r8) :: tk_zi(shcol,nlevi)
  real(r8) :: shoc_mix_zi(shcol,nlevi)
  real(r8) :: sm ! Mixing coefficient
  real(r8) :: ustar2, wstar, uf
  
  ! Constants to parameterize surface variances
  real(r8), parameter :: a_const = 1.8_r8
  real(r8), parameter :: z_const = 1.0_r8 
  real(r8), parameter :: ufmin = 0.01_r8
  
  ! Interpolate some variables from the midpoint grid to the interface grid
  call linear_interp(zt_grid,zi_grid,isotropy,isotropy_zi,nlev,nlevi,shcol,0._r8)
  call linear_interp(zt_grid,zi_grid,tkh,tkh_zi,nlev,nlevi,shcol,0._r8)
  call linear_interp(zt_grid,zi_grid,tk,tk_zi,nlev,nlevi,shcol,0._r8)
  call linear_interp(zt_grid,zi_grid,shoc_mix,shoc_mix_zi,nlev,nlevi,shcol,0._r8)
 
  ! Vertical velocity variance is assumed to be propotional
  !  to the TKE
  w_sec = w2tune*(2._r8/3._r8)*tke
  
  ! apply the surface conditions to diagnose turbulent 
  !  moments at the surface
  do i=1,shcol
    
    ! Parameterize thermodyanmics variances via Andre et al. 1978
    ustar2 = sqrt(uw_sfc(i) * uw_sfc(i) + vw_sfc(i) * vw_sfc(i)) 
    if (wthl_sfc(i) > 0._r8) then
      wstar = (1._r8/basetemp * ggr * wthl_sfc(i) * z_const)**(1._r8/3._r8)
    else
      wstar = 0._r8
    endif
    uf = sqrt(ustar2 + 0.3_r8 * wstar * wstar)
    uf = max(ufmin,uf)
    
    ! Diagnose thermodynamics variances and covariances
    thl_sec(i,nlevi) = 0.4_r8 * a_const * (wthl_sfc(i)/uf)**2
    qw_sec(i,nlevi) = 0.4_r8 * a_const * (wqw_sfc(i)/uf)**2
    qwthl_sec(i,nlevi) = 0.2_r8 * a_const * (wthl_sfc(i)/uf) * &
                     (wqw_sfc(i)/uf)
   
    ! Vertical fluxes of heat and moisture, simply 
    !  use the surface fluxes given by host model
    wthl_sec(i,nlevi) = wthl_sfc(i)
    wqw_sec(i,nlevi) = wqw_sfc(i)
    uw_sec(i,nlevi) = uw_sfc(i)
    vw_sec(i,nlevi) = vw_sfc(i)
    wtke_sec(i,nlevi) = max(sqrt(ustar2),0.01_r8)**3
    do p=1,num_tracer
      wtracer_sec(i,nlevi,p) = wtracer_sfc(i,p)
    enddo

  enddo ! end i loop (column loop)
  
  ! Calculate the second moments, which are on the
  !  interface grid.
  do k=2,nlev
    do i=1,shcol
    
      kt=k-1 ! define upper grid point indicee
      grid_dz=1._r8/dz_zi(i,k) ! Define the vertical grid difference
      grid_dz2=(1._r8/dz_zi(i,k))**2 !squared
    
      sm=isotropy_zi(i,k)*tkh_zi(i,k) ! coefficient for variances
      
      ! Compute variance of thetal
      thl_sec(i,k)=thl2tune*sm*grid_dz2*(thetal(i,kt)-thetal(i,k))**2
      ! Compute variance of total water mixing ratio
      qw_sec(i,k)=qw2tune*sm*grid_dz2*(qw(i,kt)-qw(i,k))**2
      ! Compute correlation betwen total water and thetal
      qwthl_sec(i,k)=qwthl2tune*sm*grid_dz2*((thetal(i,kt)-thetal(i,k))* &
        (qw(i,kt)-qw(i,k)))
	
      ! diagnose vertical heat flux
      wthl_sec(i,k)=-1._r8*tkh_zi(i,k)*grid_dz*(thetal(i,kt)-thetal(i,k))
      ! diagnose vertical moisture flux
      wqw_sec(i,k)=-1._r8*tkh_zi(i,k)*grid_dz*(qw(i,kt)-qw(i,k))

      ! The below fluxes are used for diagnostic purposes
      ! diagnose vertical TKE flux      
      wtke_sec(i,k)=-1._r8*tkh_zi(i,k)*grid_dz*(tke(i,kt)-tke(i,k))
 
      ! diagnose vertical momentum transport
      uw_sec(i,k)=-1._r8*tk_zi(i,k)*grid_dz*(u_wind(i,kt)-u_wind(i,k))
      vw_sec(i,k)=-1._r8*tk_zi(i,k)*grid_dz*(v_wind(i,kt)-v_wind(i,k))

      ! diagnose vertical flux of tracers      
      do p=1,num_tracer
        wtracer_sec(i,k,p)=-1._r8*tkh_zi(i,k)*grid_dz*(tracer(i,kt,p)-tracer(i,k,p))
      enddo
    
    enddo ! end i loop (column loop)
  enddo  ! end k loop (vertical loop)

  ! apply the upper boundary condition
  do i=1,shcol
    wthl_sec(i,1) = 0._r8
    wqw_sec(i,1) = 0._r8
    uw_sec(i,1) = 0._r8
    vw_sec(i,1) = 0._r8
    wtracer_sec(i,1,:) = 0._r8
    wtke_sec(i,1) = 0._r8
    
    thl_sec(i,1) = 0._r8
    qw_sec(i,1) = 0._r8
    qwthl_sec(i,1) = 0._r8    
  enddo ! end i loop (column loop)

  return

end subroutine diag_second_shoc_moments

!==============================================================
! SHOC Diagnose the third order moment of vertical velocity

subroutine diag_third_shoc_moments(&
              shcol,nlev,nlevi, &                 ! Input 
	      w_sec, thl_sec, qw_sec, qwthl_sec,& ! Input
	      wthl_sec, isotropy, brunt,&         ! Input
	      thetal,tke,wthv_sec,shoc_mix,&      ! Input
	      dz_zt, dz_zi,&                      ! Input
	      zt_grid,zi_grid,&                   ! Input
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
  real(r8), intent(in) :: w_sec(shcol,nlev)
  ! second order liquid wat. potential temperature [K^2]	
  real(r8), intent(in) :: thl_sec(shcol,nlevi)  
  ! second order total water mixing ratio [kg2/kg2]
  real(r8), intent(in) :: qw_sec(shcol,nlevi) 
  ! covariance of temp and moisture [K kg/kg]
  real(r8), intent(in) :: qwthl_sec(shcol,nlevi)
  ! vertical flux of heat [K m/s]
  real(r8), intent(in) :: wthl_sec(shcol,nlevi)
  ! return to isotropy timescale [s] 
  real(r8), intent(in) :: isotropy(shcol,nlev)
  ! brunt vaisallia frequency [s]
  real(r8), intent(in) :: brunt(shcol,nlev) 
  ! liquid water potential temperature [K]
  real(r8), intent(in) :: thetal(shcol,nlev) 
  ! turbulent kinetic energy [m2/s2]
  real(r8), intent(in) :: tke(shcol,nlev)
  ! buoyancy flux [K m/s] 
  real(r8), intent(in) :: wthv_sec(shcol,nlev)  
  ! mixing length [m]
  real(r8), intent(in) :: shoc_mix(shcol,nlev)
  ! thickness centered on thermodynamic grid [m]
  real(r8), intent(in) :: dz_zt(shcol,nlev) 
  ! thickness centered on interface grid [m]
  real(r8), intent(in) :: dz_zi(shcol,nlev) 	
  ! heights of thermodynamics points [m]
  real(r8), intent(in) :: zt_grid(shcol,nlev)
  ! heights of interface points [m]
  real(r8), intent(in) :: zi_grid(shcol,nlevi) 
  
! OUTPUT VARIABLES
  ! third moment of vertical velocity
  real(r8), intent(out) :: w3(shcol,nlevi) 
  
! LOCAL VARIABLES
  integer i, j, k, kb, kc
  real(r8) :: omega0, omega1, omega2
  real(r8) :: X0, Y0, X1, Y1, AA0, AA1
  real(r8) :: zvar, x5var, iso, thedz, thedz2
  real(r8) :: theterm, cond, tsign
  real(r8) :: isosqrt, dthl2, dwthl, dtke, dw2, aw2  
  real(r8) :: a0, a1, a2, a3, a4, a5 
  real(r8) :: buoy_sgs2, c, grd, bet2, bet
  real(r8) :: f0, f1, f2, f3, f4, f5
  
  real(r8) :: w_sec_zi(shcol,nlevi)	! second order vertical velocity 
  real(r8) :: isotropy_zi(shcol,nlevi)
  real(r8) :: brunt_zi(shcol,nlevi)  
  real(r8) :: thetal_zi(shcol,nlevi)
  real(r8) :: wthv_sec_zi(shcol,nlevi)
  real(r8) :: shoc_mix_zi(shcol,nlevi)

  ! Interpolate variables onto the interface levels
  call linear_interp(zt_grid,zi_grid,isotropy,isotropy_zi,nlev,nlevi,shcol,0._r8)
  call linear_interp(zt_grid,zi_grid,brunt,brunt_zi,nlev,nlevi,shcol,largeneg)
  call linear_interp(zt_grid,zi_grid,w_sec,w_sec_zi,nlev,nlevi,shcol,(2._r8/3._r8)*mintke)
  call linear_interp(zt_grid,zi_grid,thetal,thetal_zi,nlev,nlevi,shcol,0._r8)
  call linear_interp(zt_grid,zi_grid,wthv_sec,wthv_sec_zi,nlev,nlevi,shcol,largeneg)
  call linear_interp(zt_grid,zi_grid,shoc_mix,shoc_mix_zi,nlev,nlevi,shcol,10._r8)
 
  c=7.0_r8
  a0=(0.52_r8*c**(-2))/(c-2._r8)
  a1=0.87_r8/(c**2)
  a2=0.5_r8/c
  a3=0.6_r8/(c*(c-2._r8))
  a4=2.4_r8/(3._r8*c+5._r8)
  a5=0.6_r8/(c*(3._r8+5._r8*c))    
  
  ! set lower condition
  w3(:,nlevi) = 0._r8

  do k=2,nlev  
    do i=1,shcol

     kb=k+1
     kc=k-1
    
     thedz=dz_zi(i,k)
     thedz2=dz_zt(i,k)+dz_zt(i,kc)
     thedz=1._r8/thedz
     thedz2=1._r8/thedz2    
      
      iso=isotropy_zi(i,k)
      isosqrt=iso**2
      buoy_sgs2=isosqrt*brunt_zi(i,k)
      bet2=ggr/thetal_zi(i,k)
      
      f0=thedz2 * bet2**3 * iso**4 * wthl_sec(i,k) * &
        (thl_sec(i,kc)-thl_sec(i,kb))
	
      f1=thedz2 * bet2**2 * iso**3 * (wthl_sec(i,k) * &
        (wthl_sec(i,kc)-wthl_sec(i,kb)) + 0.5_r8 * &
	w_sec_zi(i,k)*(thl_sec(i,kc)-thl_sec(i,kb))) ! bug here
	
      f2=thedz * bet2 * isosqrt * wthl_sec(i,k) * &
        (w_sec(i,kc)-w_sec(i,k))+ 2._r8 * thedz2 * bet2 * &   
	isosqrt * w_sec_zi(i,k) * (wthl_sec(i,kc) - wthl_sec(i,kb)) 
	
      f3=thedz2 * bet2 * isosqrt * w_sec_zi(i,k) * &
        (wthl_sec(i,kc) - wthl_sec(i,kb)) + thedz * &
	bet2 * isosqrt * (wthl_sec(i,k) * (tke(i,kc) - tke(i,k))) 
	
      f4=thedz * iso * w_sec_zi(i,k) * ((w_sec(i,kc) - w_sec(i,k) + &
        (tke(i,kc) - tke(i,k))))
	
      f5=thedz * iso * w_sec_zi(i,k) * (w_sec(i,kc) - w_sec(i,k))
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Compute the omega terms     
      
      omega0 = a4 / (1._r8 - a5 * buoy_sgs2)
      omega1 = omega0/(2._r8 * c)
      omega2 = omega1 * f3 + (5._r8/4._r8) * omega0 * f4
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Compute the X0, Y0, X1, Y1 terms  
      
      X0 = (a2 * buoy_sgs2 * (1._r8 - a3 * buoy_sgs2)) / &
        (1._r8 - (a1 + a3) * buoy_sgs2)
      Y0 = (2._r8 * a2 * buoy_sgs2 * X0) / (1._r8 - a3 * buoy_sgs2)
      X1 = (a0 * f0 + a1 * f1 + a2 * (1._r8 - a3 * buoy_sgs2) * f2) / &
        (1._r8 - (a1 + a3) * buoy_sgs2)
      Y1 = (2._r8 * a2 * (buoy_sgs2 * X1 + (a0/a1) * f0 + f1)) / &   ! bug here!
        (1._r8 - a3* buoy_sgs2)     

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
      ! Compute the A0, A1 terms	
      
      AA0 = omega0 * X0 + omega1 * Y0
      AA1 = omega0 * X1 + omega1 * Y1 + omega2      
      
      ! Finally, we have the third moment of w
      w3(i,k)=(AA1-1.2_r8*X1-1.5_r8*f5)/(c-1.2_r8*X0+AA0)      
      
    enddo  ! end i loop (column loop)
  enddo  ! end k loop (vertical loop)
  
  ! set upper condition
  w3(:,1) = 0._r8

  ! perform clipping to prevent unrealistically large values from occuring
  do k=1,nlevi
    do i=1,shcol
    
      tsign = 1._r8
      theterm = w_sec_zi(i,k)
      cond = w3clip * sqrt(2._r8 * theterm**3)
      if (w3(i,k) .lt. 0) tsign = -1._r8
      if (tsign * w3(i,k) .gt. cond) w3(i,k) = tsign * cond
    
    enddo ! end i loop (column loop)
  enddo ! end k loop (vertical loop)
 
  return
  
end subroutine diag_third_shoc_moments

!==============================================================
! Assumed PDF closure for the SHOC scheme

subroutine shoc_assumed_pdf(&
             shcol,nlev,nlevi, &                ! Input
             thetal,qw,w_field,thl_sec,qw_sec,& ! Input
	     wthl_sec,w_sec, &                  ! Input
	     wqw_sec,qwthl_sec,w3,pres, &       ! Input
	     zt_grid,zi_grid,&                  ! Input
	     shoc_cldfrac,shoc_ql,&             ! Output
             wqls,wthv_sec)                     ! Output

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
  real(r8), intent(in) :: thetal(shcol,nlev) 
  ! total water mixing ratio [kg/kg]
  real(r8), intent(in) :: qw(shcol,nlev) 
  ! thetal variance [K^2]
  real(r8), intent(in) :: thl_sec(shcol,nlevi) 
  ! qw variance [kg/kg^2]
  real(r8), intent(in) :: qw_sec(shcol,nlevi) 
  ! vertical flux of heat [K m/s]
  real(r8), intent(in) :: wthl_sec(shcol,nlevi) 
  ! vertical velocity variance [m2/s2]
  real(r8), intent(in) :: w_sec(shcol,nlev) 
  ! vertical flux of moisture [kg/kg m/s]
  real(r8), intent(in) :: wqw_sec(shcol,nlevi) 
  ! qw and thetal correlation [K kg/kg]
  real(r8), intent(in) :: qwthl_sec(shcol,nlevi) 
  ! third moment vertical velocity [m^3/s^3]
  real(r8), intent(in) :: w3(shcol,nlevi) 
  ! large scale vertical velocity [m/s]
  real(r8), intent(in) :: w_field(shcol,nlev)
  ! pressure [Pa] 
  real(r8), intent(in) :: pres(shcol,nlev) 
  ! heights on midpoint grid [m] 
  real(r8), intent(in) :: zt_grid(shcol,nlev) 
  ! heights on interface grid [m]
  real(r8), intent(in) :: zi_grid(shcol,nlevi) 
  
! OUTPUT VARIABLES
  ! SGS cloud fraction [-]
  real(r8), intent(out) :: shoc_cldfrac(shcol,nlev) 
  ! SGS liquid water mixing ratio [kg/kg]
  real(r8), intent(out) :: shoc_ql(shcol,nlev) 
  ! SGS buoyancy flux [K m/s]
  real(r8), intent(out) :: wthv_sec(shcol,nlev) 
  ! SGS liquid water flux [kg/kg m/s]
  real(r8), intent(out) :: wqls(shcol,nlev)

! LOCAL VARIABLES
  integer i,j,k,dothis,nmicro_fields
  real(r8) skew_w,skew_thl,skew_qw,a
  real(r8) w1_1,w1_2,w2_1,w2_2,w3var
  real(r8) thl1_1,thl1_2,thl2_1,thl2_2
  real(r8) qw1_1,qw1_2,qw2_1,qw2_2
  real(r8) r_qwthl_1,r_wqw_1,r_wthl_1
  real(r8) testvar,s1,s2,std_s1,std_s2,C1,C2
  real(r8) ql1,ql2,flux_ws1,flux_ws2
  real(r8) w_ql1,w_ql2,thl_first,qw_first,w_first
  real(r8) Tl1_1,Tl1_2,betatest,pval
  real(r8) w2thl,w2qw,w2ql,w2ql_1,w2ql_2
  real(r8) s,thec,thlsec,qwsec,qwthlsec,wqwsec,wthlsec
  real(r8) thestd,erf,exp,dum
  real(r8) cqt1,cthl1,cqt2,cthl2
  real(r8) qn1,qn2,qi1,qi2,omn1,omn2, omn
  real(r8) lstarn, test
  real(r8) m, bigm, basetemp2
  real(r8) beta1, beta2, qs1, qs2
  real(r8) esval1_1, esval2_1, esval1_2, esval2_2, om1, om2
  real(r8) lstarn1, lstarn2, sqrtw2, sqrtthl, sqrtqt
  real(r8) sqrtstd1, sqrtstd2, tsign, tvar, sqrtw2t
  real(r8) wqis, skip, epsterm
  real(r8) sqrtqw2_1, sqrtqw2_2, sqrtthl2_1, sqrtthl2_2
  real(r8) corrtest1, corrtest2, thl_tol, rt_tol, w_tol_sqd, w_thresh
  
  ! variables on thermo grid
  real(r8) :: wthl_sec_zt(shcol,nlev)
  real(r8) :: wqw_sec_zt(shcol,nlev)
  real(r8) :: w3_zt(shcol,nlev)
  real(r8) :: thl_sec_zt(shcol,nlev)
  real(r8) :: qwthl_sec_zt(shcol,nlev)
  real(r8) :: qw_sec_zt(shcol,nlev)
  real(r8) :: w_field_zt(shcol,nlev)

  ! define these so they don't have to be computed more than once
  real(r8), parameter :: sqrt2 = sqrt(2._r8)
  real(r8), parameter :: sqrtpi = sqrt(2._r8*3.14_r8)
  
  epsterm=rgas/rv 
  
  thl_tol=1.e-2_r8
  rt_tol=1.e-4_r8
  w_tol_sqd=(2.e-2_r8)**2
  w_thresh=0.0_r8
  
  ! Initialize cloud variables to zero  
  shoc_cldfrac(:,:)=0._r8
  shoc_ql(:,1)=0._r8 

  ! Interpolate many variables from interface grid to themo grid
  call linear_interp(zi_grid,zt_grid,w3,w3_zt,nlevi,nlev,shcol,largeneg)  
  call linear_interp(zi_grid,zt_grid,thl_sec,thl_sec_zt,nlevi,nlev,shcol,0._r8)
  call linear_interp(zi_grid,zt_grid,wthl_sec,wthl_sec_zt,nlevi,nlev,shcol,largeneg) !Alert
  call linear_interp(zi_grid,zt_grid,qwthl_sec,qwthl_sec_zt,nlevi,nlev,shcol,largeneg)
  call linear_interp(zi_grid,zt_grid,wqw_sec,wqw_sec_zt,nlevi,nlev,shcol,largeneg) !Alert
  call linear_interp(zi_grid,zt_grid,qw_sec,qw_sec_zt,nlevi,nlev,shcol,0._r8)  
  call linear_interp(zi_grid,zt_grid,w_field,w_field_zt,nlevi,nlev,shcol,largeneg)
  
  do k=1,nlev
    do i=1,shcol

      pval = pres(i,k) 

      ! Get all needed input moments for the PDF
      !  at this particular point
      thl_first = thetal(i,k)
      w_first = w_field_zt(i,k)
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
      
      Skew_w=w3var/w_sec(i,k)**(3./2.)

      if (w_sec(i,k) .le. w_tol_sqd) then
        Skew_w=0._r8
        w1_1=w_first
        w1_2=w_first
        w2_1=0._r8
        w2_2=0._r8
        a=0.5_r8
      else

        w2_1=0.4_r8
        w2_2=0.4_r8

        a=max(0.01_r8,min(0.5_r8*(1._r8-Skew_w*sqrt(1._r8/(4._r8*(1._r8-w2_1)**3+Skew_w**2))),0.99_r8))

        sqrtw2t=sqrt(1._r8-w2_1)

        w1_1=sqrt((1._r8-a)/a)*sqrtw2t
        w1_2=-1._r8*sqrt(a/(1._r8-a))*sqrtw2t

        w2_1=w2_1*w_sec(i,k)
        w2_2=w2_2*w_sec(i,k)


      endif  
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  FIND PARAMETERS FOR THETAL
    
      corrtest1=max(-1._r8,min(1._r8,wthlsec/(sqrtw2*sqrtthl)))

      if (thlsec .le. thl_tol**2 .or. abs(w1_2-w1_1) .le. w_thresh) then
        thl1_1=thl_first
        thl1_2=thl_first
        thl2_1=0._r8
        thl2_2=0._r8
        sqrtthl2_1=0._r8
        sqrtthl2_2=0._r8
      else

        thl1_1=(-1._r8*corrtest1)/w1_2
        thl1_2=(-1._r8*corrtest1)/w1_1
      
        if (dothetal_skew) then
          tsign=abs(thl1_2-thl1_1)
	
	  if (tsign .gt. 0.4_r8) then
	    Skew_thl=1.2_r8*Skew_w
	  else if (tsign .le. 0.2_r8) then
	    Skew_thl=0.0_r8
	  else
	    Skew_thl=((1.2_r8*Skew_w)/0.2_r8)*(tsign-0.2_r8)
	  endif 
        else
          Skew_thl = 0.0_r8
        endif
	
        thl2_1=min(100._r8,max(0._r8,(3._r8*thl1_2*(1._r8-a*thl1_1**2-(1._r8-a)*thl1_2**2) &     
	    -(Skew_thl-a*thl1_1**3-(1._r8-a)*thl1_2**3))/ &
	    (3._r8*a*(thl1_2-thl1_1))))*thlsec
      
        thl2_2=min(100._r8,max(0._r8,(-3._r8*thl1_1*(1._r8-a*thl1_1**2-(1._r8-a)*thl1_2**2) &
          +(Skew_thl-a*thl1_1**3-(1._r8-a)*thl1_2**3))/ &
          (3._r8*(1._r8-a)*(thl1_2-thl1_1))))*thlsec

        thl1_1=thl1_1*sqrtthl+thl_first
        thl1_2=thl1_2*sqrtthl+thl_first

        sqrtthl2_1=sqrt(thl2_1)
        sqrtthl2_2=sqrt(thl2_2)

      endif	                  
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  FIND PARAMETERS FOR TOTAL WATER MIXING RATIO

      corrtest2=max(-1.0_r8,min(1.0_r8,wqwsec/(sqrtw2*sqrtqt)))

      if (qwsec .le. rt_tol**2 .or. abs(w1_2-w1_1) .le. w_thresh) then
        qw1_1=qw_first
        qw1_2=qw_first
        qw2_1=0._r8
        qw2_2=0._r8
        sqrtqw2_1=0._r8
        sqrtqw2_2=0._r8
      else

        qw1_1=(-1._r8*corrtest2)/w1_2
        qw1_2=(-1._r8*corrtest2)/w1_1      

        tsign=abs(qw1_2-qw1_1)

        if (tsign .gt. 0.4_r8) then
          Skew_qw=1.2_r8*Skew_w
        else if (tsign .le. 0.2_r8) then
          Skew_qw=0._r8
        else
          Skew_qw=((1.2_r8*Skew_w)/0.2_r8)*(tsign-0.2_r8)
        endif

        qw2_1=min(100._r8,max(0._r8,(3._r8*qw1_2*(1._r8-a*qw1_1**2-(1._r8-a)*qw1_2**2) &
          -(Skew_qw-a*qw1_1**3-(1._r8-a)*qw1_2**3))/ &
          (3._r8*a*(qw1_2-qw1_1))))*qwsec

        qw2_2=min(100._r8,max(0._r8,(-3._r8*qw1_1*(1._r8-a*qw1_1**2-(1._r8-a)*qw1_2**2) &
          +(Skew_qw-a*qw1_1**3-(1._r8-a)*qw1_2**3))/ &
          (3._r8*(1._r8-a)*(qw1_2-qw1_1))))*qwsec

        qw1_1=qw1_1*sqrtqt+qw_first
        qw1_2=qw1_2*sqrtqt+qw_first

        sqrtqw2_1=sqrt(qw2_1)
        sqrtqw2_2=sqrt(qw2_2)

      endif
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  CONVERT FROM TILDA VARIABLES TO "REAL" VARIABLES

      w1_1=w1_1*sqrtw2+w_first
      w1_2=w1_2*sqrtw2+w_first
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  FIND WITHIN-PLUME CORRELATIONS 

      testvar=(a*sqrtqw2_1*sqrtthl2_1+(1._r8-a)*sqrtqw2_2*sqrtthl2_2)

      if (testvar .eq. 0._r8) then
        r_qwthl_1=0._r8
      else
        r_qwthl_1=max(-1.0_r8,min(1.0_r8,(qwthlsec-a*(qw1_1-qw_first) &
	  *(thl1_1-thl_first)-(1._r8-a)*(qw1_2-qw_first) &
	  *(thl1_2-thl_first))/testvar))
      endif    
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  BEGIN TO COMPUTE CLOUD PROPERTY STATISTICS

      Tl1_1=thl1_1/((basepres/pval)**(rgas/cp))
      Tl1_2=thl1_2/((basepres/pval)**(rgas/cp))

      ! Now compute qs

      esval1_1=0._r8
      esval1_2=0._r8
      om1=1._r8
      om2=1._r8 
    
      esval1_1=esatw_shoc(Tl1_1)*100._r8
      lstarn1=lcond 
	
      qs1=0.622_r8*esval1_1/max(esval1_1,pval-esval1_1)
      beta1=(rgas/rv)*(lstarn1/(rgas*Tl1_1))*(lstarn1/(cp*Tl1_1))
      
      ! Are the two plumes equal?  If so then set qs and beta
      ! in each column to each other to save computation
      lstarn2=lcond
      if (Tl1_1 .eq. Tl1_2) then
        qs2=qs1     
        beta2=beta1
      else
        esval1_2=esatw_shoc(Tl1_2)*100._r8  
        qs2=0.622_r8*esval1_2/max(esval1_2,pval-esval1_2)
        beta2=(rgas/rv)*(lstarn2/(rgas*Tl1_2))*(lstarn2/(cp*Tl1_2)) 
      endif      

      !!!!!  Now compute cloud stuff
      !!!!!!  compute s term    
      
      s1=qw1_1-qs1*((1._r8+beta1*qw1_1)/(1._r8+beta1*qs1))
      cthl1=((1._r8+beta1*qw1_1)/(1._r8+beta1*qs1)**2)*(cp/lcond) &
        *beta1*qs1*(pval/basepres)**(rgas/cp)

      cqt1=1._r8/(1._r8+beta1*qs1)
      std_s1=sqrt(max(0._r8,cthl1**2*thl2_1+cqt1**2*qw2_1-2._r8*cthl1 &
        *sqrtthl2_1*cqt1*sqrtqw2_1*r_qwthl_1))
	
      qn1=0._r8
      C1=0._r8
      
      if (std_s1 .ne. 0.0_r8) then
        C1=0.5_r8*(1._r8+erf(s1/(sqrt2*std_s1)))
	IF (C1 .ne. C1) C1 = 0._r8
        IF (C1 .ne. 0._r8) qn1=s1*C1+(std_s1/sqrtpi)*exp(-0.5_r8*(s1/std_s1)**2)
      else
        if (s1 .gt. 0._r8) then
          C1=1.0_r8
          qn1=s1
        endif
      endif
      
      !!!!! now compute non-precipitating cloud condensate 

      ! If two plumes exactly equal, then just set many of these 
      ! variables to themselves to save on computation.
      if (qw1_1 .eq. qw1_2 .and. thl2_1 .eq. thl2_2 .and. qs1 .eq. qs2) then
        s2=s1
        cthl2=cthl1
        cqt2=cqt1
        std_s2=std_s1
        C2=C1
        qn2=qn1
      else
        
        s2=qw1_2-qs2*((1._r8+beta2*qw1_2)/(1._r8+beta2*qs2))
        cthl2=((1._r8+beta2*qw1_2)/(1._r8+beta2*qs2)**2)*(cp/lcond) &
	  *beta2*qs2*(pval/basepres)**(rgas/cp)
        cqt2=1._r8/(1._r8+beta2*qs2)
        std_s2=sqrt(max(0._r8,cthl2**2*thl2_2+cqt2**2*qw2_2-2._r8*cthl2* &
	  sqrtthl2_2*cqt2*sqrtqw2_2*r_qwthl_1))

        qn2=0._r8
        C2=0._r8
	
        if (std_s2 .ne. 0._r8) then
          C2=0.5_r8*(1.+erf(s2/(sqrt2*std_s2)))
	  if (C2 .ne. C2) C2 = 0._r8
          if (C2 .ne. 0._r8) qn2=s2*C2+(std_s2/sqrtpi)*exp(-0.5_r8*(s2/std_s2)**2)
        else
          if (s2 .gt. 0._r8) then
            C2=1.0_r8
            qn2=s2
          endif
        endif

      endif
 
      ! Finally, compute SGS cloud fraction
      shoc_cldfrac(i,k) = min(1._r8,a*C1+(1._r8-a)*C2)
 
      ql1=min(qn1,qw1_1)
      ql2=min(qn2,qw1_2)
      
      ! Compute SGS liquid water mixing ratio
      shoc_ql(i,k) = max(0._r8,a*ql1+(1._r8-a)*ql2)
      
      ! Compute liquid water flux
      wqls(i,k)=a*((w1_1-w_first)*ql1)+(1._r8-a)*((w1_2-w_first)*ql2)
      ! Compute the SGS buoyancy flux
      wthv_sec(i,k)=wthlsec+((1._r8-epsterm)/epsterm)*basetemp*wqwsec &
        +((lcond/cp)*(basepres/pval)**(rgas/cp)-(1._r8/epsterm)*basetemp)*wqls(i,k)  
     	
    enddo  ! end i loop here
  enddo	  ! end k loop here
	                     
  return                              
     
end subroutine shoc_assumed_pdf

!==============================================================
! Advance turbulent kinetic energy equation

subroutine shoc_tke(&
             shcol,nlev,nlevi,dtime,&    ! Input
             wthv_sec,shoc_mix,&         ! Input
	     dz_zi,dz_zt,pres,&          ! Input  
	     u_wind,v_wind,brunt,&       ! Input
	     uw_sfc,vw_sfc,&             ! Input
	     zt_grid,zi_grid,&           ! Input
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
  real(r8), intent(in) :: dtime 
  ! SGS buoyancy flux [K m/s]
  real(r8), intent(in) :: wthv_sec(shcol,nlev) 
  ! Mixing length [m]
  real(r8), intent(in) :: shoc_mix(shcol,nlev)
  ! Meridional wind [m/s] 
  real(r8), intent(in) :: u_wind(shcol,nlev) 
  ! Zonal wind [m/s]
  real(r8), intent(in) :: v_wind(shcol,nlev) 
  ! Zonal momentum flux at sfc [m2/s2]
  real(r8), intent(in) :: uw_sfc(shcol) 
  ! Meridional momentum flux at sfc [m2/s2]
  real(r8), intent(in) :: vw_sfc(shcol) 
  ! thickness on interface grid [m]
  real(r8), intent(in) :: dz_zi(shcol,nlev)
  ! thickness on thermodynamic grid [m]
  real(r8), intent(in) :: dz_zt(shcol,nlev)
  ! pressure [Pa]
  real(r8), intent(in) :: pres(shcol,nlev)
  ! Brunt Vaisalla frequncy [/s] 
  real(r8), intent(in) :: brunt(shcol,nlev)  
  ! heights on midpoint grid [m]
  real(r8), intent(in) :: zt_grid(shcol,nlev) 
  ! heights on interface grid [m]
  real(r8), intent(in) :: zi_grid(shcol,nlevi) 
  
! INPUT/OUTPUT VARIABLES
  ! turbulent kinetic energy [m2/s2]
  real(r8), intent(inout) :: tke(shcol,nlev)
  ! eddy coefficient for momentum [m2/s]
  real(r8), intent(inout) :: tk(shcol,nlev)
  ! eddy coefficient for heat [m2/s]
  real(r8), intent(inout) :: tkh(shcol,nlev)
  
! OUTPUT VARIABLES
  ! Return to isotropic timescale [s]
  real(r8), intent(out) :: isotropy(shcol,nlev) 
  
! LOCAL VARIABLES	     
  real(r8) :: shear_prod(shcol,nlevi)
  real(r8) :: shear_prod_zt(shcol,nlev), tk_zi(shcol,nlevi)
  real(r8) :: grd,betdz,Ck,Ckh,Ckm,Ce,Ces,Ce1,Ce2,smix,Cee,Cs
  real(r8) :: buoy_sgs,ratio,a_prod_sh,a_prod_bu,a_diss
  real(r8) :: lstarn, lstarp, bbb, omn, omp, ustar
  real(r8) :: qsatt,dqsat,tk_in, u_grad, v_grad
  real(r8) :: tscale1,lambda,buoy_sgs_save,grid_dzw,grw1,grid_dz
  real(r8) :: lambda_low,lambda_high,lambda_slope, brunt_low
  real(r8) :: brunt_int(shcol)
  integer i,j,k,kc,kb	 
  
  lambda_low=0.001_r8
  lambda_high=0.04_r8
  lambda_slope=0.65_r8
  brunt_low=0.02
 
  ! Turbulent coefficients
  Cs=0.15_r8
  Ck=0.1_r8
  Ckh=0.1_r8
  Ckm=0.1_r8
  Ce=Ck**3/Cs**4 
  
  Ce1=Ce/0.7_r8*0.19_r8
  Ce2=Ce/0.7_r8*0.51_r8 

  Cee=Ce1+Ce2
  
  ! Compute integrated column stability in lower troposphere
  brunt_int(:)=0._r8
  do k=1,nlev
    do i=1,shcol
      if (pres(i,k) .gt. 80000._r8) then
        brunt_int(i) = brunt_int(i) + dz_zt(i,k)*brunt(i,k)
      endif
    enddo
  enddo

  ! Interpolate tk onto interface grid
  call linear_interp(zt_grid,zi_grid,tk,tk_zi,nlev,nlevi,shcol,0._r8)

  ! Compute shear production term, which is on interface levels
  ! This follows the methods of Bretheron and Park (2010)
  do k=1,nlev-1
    do i=1,shcol
    
      kb=k+1     
      grid_dz = 1._r8/dz_zi(i,k)
  
      tk_in=tk_zi(i,k)
      ! calculate vertical gradient of u&v wind
      u_grad=grid_dz*(u_wind(i,k)-u_wind(i,kb))
      v_grad=grid_dz*(v_wind(i,k)-v_wind(i,kb))  
      shear_prod(i,k)=tk_in*(u_grad**2+v_grad**2) 
    enddo
  enddo
  
  ! Set lower and upper boundary for shear production
  ! Note that the lower bound for shear production has already 
  !  been taken into account for the TKE boundary condition, 
  !  thus zero out here
  shear_prod(:,1) = 0._r8
  shear_prod(:,nlevi) = 0._r8
  
  ! Interpolate shear production from interface to thermo grid
  call linear_interp(zi_grid,zt_grid,shear_prod,shear_prod_zt,nlevi,nlev,shcol,largeneg)

  do k=1,nlev
    do i=1,shcol
    
      smix=shoc_mix(i,k)
      ! Compute buoyant production term
      a_prod_bu=(ggr/basetemp)*wthv_sec(i,k)
           
      tke(i,k)=max(0._r8,tke(i,k))
      
      ! Shear production term
      a_prod_sh=shear_prod_zt(i,k)
      
      ! Dissipation term
      a_diss=Cee/shoc_mix(i,k)*tke(i,k)**1.5
      
      ! March equation forward one timestep
      tke(i,k)=max(0._r8,tke(i,k)+dtime*(max(0._r8,a_prod_sh+a_prod_bu)-a_diss))  
      
      tke(i,k)=min(tke(i,k),maxtke)
   
      ! Now compute the return to isotropic timescale as per
      ! Canuto et al. 2004.  This is used to define the 
      ! eddy coefficients as well as to diagnose higher 
      ! moments in SHOC
      
      ! define the time scale
      tscale1=(2.0_r8*tke(i,k))/a_diss
      
      ! define a damping term "lambda" based on column stability
      lambda=lambda_low+((brunt_int(i)/ggr)-brunt_low)*lambda_slope
      lambda=max(lambda_low,min(lambda_high,lambda))
       
      buoy_sgs_save=brunt(i,k)
      if (buoy_sgs_save .le. 0._r8) lambda=0._r8
      
      ! Compute the return to isotropic timescale
      isotropy(i,k)=min(maxiso,tscale1/(1._r8+lambda*buoy_sgs_save*tscale1**2))

      ! Define the eddy coefficients for heat and momentum
      tkh(i,k)=Ckh*isotropy(i,k)*tke(i,k)  
      tk(i,k)=Ckm*isotropy(i,k)*tke(i,k)            

      tke(i,k) = max(mintke,tke(i,k))
 
    enddo ! end i loop
  enddo ! end k loop

  return
  
end subroutine shoc_tke

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
  real(r8), intent(inout) :: tke(shcol,nlev)

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
             host_dx,host_dy,cldin,&       ! Input
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
  real(r8), intent(in) :: host_dx(shcol)
  ! host model grid size [m] 
  real(r8), intent(in) :: host_dy(shcol) 
  ! turbulent kinetic energy [m^2/s^2]
  real(r8), intent(in) :: tke(shcol,nlev) 
  ! cloud liquid water mixing ratio [kg/kg]
  real(r8), intent(in) :: cldin(shcol,nlev) 
  ! heights on midpoint grid [m]
  real(r8), intent(in) :: zt_grid(shcol,nlev)
  ! heights on interface grid [m] 
  real(r8), intent(in) :: zi_grid(shcol,nlev)
  ! dz on midpoint grid [m] 
  real(r8), intent(in) :: dz_zt(shcol,nlev) 
  ! dz on interface grid [m]
  real(r8), intent(in) :: dz_zi(shcol,nlev) 
  ! SGS buoyancy flux [K m/s]
  real(r8), intent(in) :: wthv_sec(shcol,nlev)
  ! liquid water potential temperature [K] 
  real(r8), intent(in) :: thetal(shcol,nlev) 
  ! virtual potential temperature [K]
  real(r8), intent(in) :: thv(shcol,nlev)

! OUTPUT VARIABLES
  ! brunt vailsailla frequency [/s]
  real(r8), intent(out) :: brunt(shcol,nlev) 
  ! SHOC mixing length [m]
  real(r8), intent(out) :: shoc_mix(shcol,nlev) 

! LOCAL VARIABLES
  integer i, j, k, kk, kt
  integer kl, ku, kb, kc, dothis, kli, kui
  real(r8) :: deep_thresh, deep_thick, cloud_thick, lstarn, thresh
  real(r8) :: cldmix, vonk, thedel, depth
  real(r8) :: omn, betdz, bbb, term, qsatt, dqsat, bet
  real(r8) :: thv_up, thv_dn, thedz, tscale, thefac, thecoef, thegam, norm
  real(r8) :: stabterm, conv_var, tkes, mmax, cldthresh
  logical lf, indexr
  logical doclouddef
  real(r8) :: conv_vel(shcol,nlev)
  real(r8) :: thv_zi(shcol,nlevi)
  
  real(r8) :: numer(shcol)
  real(r8) :: denom(shcol)
  real(r8) :: cldarr(shcol) 
  real(r8) :: l_inf(shcol)
  real(r8) :: brunt2(shcol,nlev)
 
  doclouddef = .true.
 
  vonk=0.35_r8   ! Vonkarman constnt
  tscale=400._r8 ! time scale set based on similarity results
  brunt2(:,:) = 0._r8
  numer(:) = 0._r8
  denom(:) = 0._r8

  ! Interpolate virtual potential temperature onto interface grid
  call linear_interp(zt_grid,zi_grid,thv,thv_zi,nlev,nlevi,shcol,0._r8)

  ! Define the brunt vaisalia frequency 
  do k=1,nlev
    do i=1,shcol
      brunt(i,k) = (ggr/thv(i,k)) * (thv_zi(i,k) - thv_zi(i,k+1))/dz_zt(i,k) 
    enddo
  enddo 
 
  ! Find length scale outside of clouds
  do k=1,nlev
    do i=1,shcol
    
      if (cldin(i,k) .eq. 0 .or. .not. doclouddef) then
        tkes=sqrt(tke(i,k))
	numer(i)=numer(i)+tkes*zt_grid(i,k)*dz_zt(i,k)
	denom(i)=denom(i)+tkes*dz_zt(i,k)
      else
        cldarr(i)=1
      endif
    
    enddo
  enddo
  
  do i=1,shcol
    if (denom(i) .gt. 0._r8) then
      l_inf(i)=0.1_r8*(numer(i)/denom(i))
    else
      l_inf(i)=100._r8
    endif
  enddo
 
  do k=1,nlev
    do i=1,shcol
    
      tkes=sqrt(tke(i,k))
      
      if (brunt(i,k) .ge. 0) brunt2(i,k) = brunt(i,k) 

      shoc_mix(i,k)=min(maxlen,(2.8284_r8*sqrt(1._r8/((1._r8/(tscale*tkes*vonk*zt_grid(i,k))) &
        +(1._r8/(tscale*tkes*l_inf(i)))+0.01_r8*(brunt2(i,k)/tke(i,k)))))/0.3_r8)     
      
    enddo  ! end i loop (column loop)
  enddo ! end k loop (vertical loop)
  
  ! Now find length scale in clouds
  cldthresh=0._r8
  
  ! determine the convective velocity scale at
  !   the top of the cloud
  
  conv_vel(:,nlev)=0._r8
  do k=nlev-1,1,-1
    do i=1,shcol
      conv_vel(i,k) = conv_vel(i,k+1)+2.5_r8*dz_zt(i,k)*(ggr/thv(i,k))*wthv_sec(i,k)
    enddo ! end i loop (column loop)
  enddo ! end k loop (vertical loop)
 
  if (doclouddef) then
 
  do i=1,shcol
  
    if (cldarr(i) .eq. 1) then
      
      kl=0
      ku=0
      
      do k=nlev-2,2,-1
        
	! Look for cloud base in this column
	if (cldin(i,k) .gt. cldthresh .and. kl .eq. 0) then
          kl=k
        endif
      
        ! Look for cloud top in this column
	if (cldin(i,k) .gt. cldthresh .and. cldin(i,k-1) .le. cldthresh) then 
	  ku=k
	  conv_var=conv_vel(i,k)**(1._r8/3._r8)
	endif
	
	! Compute the mixing length for the layer just determined
	if (kl .gt. 0 .and. ku .gt. 0 .and. kl-ku .gt. 1) then
	
	  if (conv_var .gt. 0) then
	    
	    depth=(zt_grid(i,ku) - zt_grid(i,kl)) + dz_zt(i,kl)
	    mmax=maxlen
	    if (zt_grid(i,ku) .gt. maxlen) mmax=maxlen
	    
            shoc_mix(i,ku:kl)=min(mmax,sqrt(1._r8/(((conv_var)/ &
              (depth*sqrt(tke(i,ku:kl))))**2+0.01_r8* &
              (brunt2(i,ku:kl)/tke(i,ku:kl))))/0.3_r8)	
	      
	  endif
	    
	  kl=0
	  ku=0
	
	endif 
	
      enddo ! end k loop
	
    endif ! end cldarr conditional  
	
  enddo ! end i loop

  endif 
  
  ! Do checks on the length scale.  Make sure it is not
  !  larger than the grid mesh of the host model.
  do k=1,nlev
    do i=1,shcol
     
      shoc_mix(i,k)=min(maxlen,shoc_mix(i,k))
      shoc_mix(i,k)=max(20._r8,shoc_mix(i,k))
      shoc_mix(i,k)=min(sqrt(host_dx(i)*host_dy(i)),shoc_mix(i,k))

    enddo
  enddo 
  
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
  real(r8), intent(in) :: dtime
  ! diffusion coefficent [m2/s]
  real(r8), intent(in) :: kv_term(shcol,nlevi)
  ! dt*(g*rho)**2/dp at interfaces
  real(r8), intent(in) :: tmpi(shcol,nlevi)
  ! 1/dp 
  real(r8), intent(in) :: rdp_zt(shcol,nlev)
  ! surface flux [varies]
  real(r8), intent(in) :: flux(shcol)
  
! OUTPUT VARIABLES
  ! superdiagonal
  real(r8), intent(out) :: ca(shcol,nlev)
  ! subdiagonal
  real(r8), intent(out) :: cc(shcol,nlev)
  ! 1./(1.+ca(k)+cc(k)-cc(k)*ze(k-1))
  real(r8), intent(out) :: denom(shcol,nlev)
  ! Term in tri-diag. matrix system
  real(r8), intent(out) :: ze(shcol,nlev)
 
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
  
  ca(:,nlev) = 0._r8
  
  ! Calculate e(k). This term is required in the solution of the
  ! tridiagonal matrix as defined by the implicit diffusion equation.   
  
  do i=1,shcol
    denom(i,nlev) = 1._r8/ &
      (1._r8 + cc(i,nlev) + flux(i)*dtime*ggr*rdp_zt(i,nlev))
    ze(i,nlev) = cc(i,nlev) * denom(i,nlev)
  enddo
  
  do k=nlev-1,2,-1
    do i=1,shcol
      denom(i,k) = 1._r8/ &
        (1._r8 + ca(i,k) + cc(i,k) - &
	ca(i,k) * ze(i,k+1))
      ze(i,k) = cc(i,k) * denom(i,k)
    enddo 
  enddo
  
  do i=1,shcol
    denom(i,1) = 1._r8/ &
      (1._r8 + ca(i,1) - ca(i,1) * ze(i,2))
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
  real(r8), intent(out) :: ca(shcol,nlev)
  ! subdiagonal
  real(r8), intent(out) :: cc(shcol,nlev)
  ! 1./(1.+ca(k)+cc(k)-cc(k)*ze(k-1))
  real(r8), intent(out) :: denom(shcol,nlev)
  ! Term in tri-diag. matrix system
  real(r8), intent(out) :: ze(shcol,nlev)
  
! IN/OUT VARIABLES
  real(r8), intent(inout) :: var(shcol,nlev)
  
! LOCAL VARIABLES
  integer :: i, k  
  ! Term in tri-diag solution
  real(r8) :: zf(shcol,nlev)
 
 
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
  ! Linear interpolation to get values on various grids

  subroutine linear_interp(x1,x2,y1,y2,km1,km2,ncol,minthresh)
    implicit none

    integer, intent(in) :: km1, km2
    integer, intent(in) :: ncol
    real(r8), intent(in) :: x1(ncol,km1), y1(ncol,km1)
    real(r8), intent(in) :: x2(ncol,km2)
    real(r8), intent(in) :: minthresh 
    real(r8), intent(out) :: y2(ncol,km2)

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

!==============================================================
! 
! Saturation vapor pressure and mixing ratio. 
! Based on Flatau et.al, (JAM, 1992:1507) - valid for T > -80C
! sat. vapor over ice below -80C - used Murphy and Koop (2005)
! For water below -80C simply assumed esw/esi = 2.
! des/dT below -80C computed as a finite difference of es

real(r8) function esatw_shoc(t)
implicit none
real(r8) t	! temperature (K)
real(r8) a0,a1,a2,a3,a4,a5,a6,a7,a8 
data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
	6.105851_r8, 0.4440316_r8, 0.1430341e-1_r8, &
        0.2641412e-3_r8, 0.2995057e-5_r8, 0.2031998e-7_r8, &
        0.6936113e-10_r8, 0.2564861e-13_r8,-0.3704404e-15_r8/
!     	6.11239921, 0.443987641, 0.142986287e-1, &
!       0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
!       0.640689451e-10, -0.952447341e-13,-0.976195544e-15/
real(r8) dt
 dt = t-273.16_r8
if(dt.gt.-80._r8) then
 esatw_shoc = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))) 
else
 esatw_shoc = 2._r8*0.01_r8*exp(9.550426_r8 - 5723.265_r8/t + 3.53068_r8*Log(t) - 0.00728332_r8*t)
end if
end     

end module shoc

!==============================================================
! This is the end of the SHOC parameterization
! We hope you have enjoyed your time here 
