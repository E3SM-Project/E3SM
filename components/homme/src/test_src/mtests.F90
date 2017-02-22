MODULE mtests

  ! M tests 1,2,3 as per Klemp, Skamarock, Park 2015, source modified from
  ! dcmip tests 2.

  !  20 - Impact of orography on a steady-state at rest
  !  21 and 22 - Non-Hydrostatic Mountain Waves Over A Schaer-Type Mountain without and with vertical wind shear
  !
  !  Given a point specified by: 
  !  	lon	longitude (radians) 
  ! 	lat	latitude (radians) 
  ! 	p/z	pressure/height
  !  the functions will return:
  !	u	zonal wind (m s^-1)
  !	v	meridional wind (m s^-1)
  !	w	vertical velocity (m s^-1)
  !	t	temperature (K)
  !	phis	surface geopotential (m^2 s^-2)
  !	ps	surface pressure (Pa)
  !	rho	density (kj m^-3)
  !	q	specific humidity (kg/kg)
  !	qi	tracers (kg/kg)
  !     p       pressure if height based (Pa)

  ! Use physical constants consistent with HOMME
  !use physical_constants, only: a=>rearth_std, Rd => Rgas, g, cp, pi=>dd_pi, p0

  IMPLICIT NONE

!-----------------------------------------------------------------------
!     Physical Parameters
!-----------------------------------------------------------------------

	real(8), parameter ::	a	= 6371220.0d0,	&	! Earth's Radius (m)
				Rd 	= 287.0d0,	&	! Ideal gas const dry air (J kg^-1 K^1)
				g	= 9.80616d0,	&	! Gravity (m s^2)
				cp	= 1004.5d0,	&	! Specific heat capacity (J kg^-1 K^1)
				pi	= 4.d0*atan(1.d0)       ! pi

! parameters of mountain and u_eq that dcmip2012_2 sets in the namelist, 
!but let's not plan on changing these for m tests.
        real(8),parameter :: mtest_ueq     = 20.d0          ! wind speed at equator (m/s)
        real(8),parameter :: mtest_h0      = 250.0d0        ! mountain height (m)
        real(8),parameter :: mtest_d       = 5000.0d0       ! mountain half width   (m)
        real(8),parameter :: mtest_xi      = 4000.0d0       ! mountain wavelength   (m)

!-----------------------------------------------------------------------
!     Additional constants
!-----------------------------------------------------------------------

	real(8), parameter :: p0 = 100000.d0 ! reference pressure (Pa)


CONTAINS

!=====================================================================================
! Tests 2-1 and 2-2:  Non-hydrostatic Mountain Waves over a Schaer-type Mountain
!=====================================================================================

SUBROUTINE mtest_state (lon,lat,p,z,hyam,hybm,u,v,w,t,phis,ps,rho,testid)

IMPLICIT NONE
!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------

  real(8), intent(in)   :: lon          ! Longitude (radians)
  real(8), intent(in)   :: lat          ! Latitude (radians)
  real(8), intent(inout):: z            ! Height (m)
  real(8), intent(in)   :: hyam         ! A coefficient for hybrid-eta coordinate, at model level midpoint
  real(8), intent(in)   :: hybm         ! B coefficient for hybrid-eta coordinate, at model level midpoint
  real(8), intent(inout):: p            ! Pressure  (Pa)
  integer, intent(in)   :: testid       ! 
  real(8), intent(out)  :: u            ! Zonal wind (m s^-1)
  real(8), intent(out)  :: v            ! Meridional wind (m s^-1)
  real(8), intent(out)  :: w            ! Vertical Velocity (m s^-1)
  real(8), intent(out)  :: T            ! Temperature (K)
  real(8), intent(out)  :: phis         ! Surface Geopotential (m^2 s^-2)
  real(8), intent(out)  :: ps           ! Surface Pressure (Pa)
  real(8), intent(out)  :: rho          ! density (kg m^-3)
!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
  real(8), parameter ::   &
  X       = 500.d0,     &   ! Reduced Earth reduction factor
  Om      = 0.d0,       &   ! Rotation Rate of Earth
  as      = a/X,        &   ! New Radius of small Earth
  Teq     = 300.d0,     &   ! Temperature at Equator
  Peq     = 100000.d0,  &   ! Reference PS at Equator


!Z TOP is defined in 3 different places!
  ztop    = 30000.d0,   &   ! Model Top
  lambdac = pi/4.d0,    &   ! Lon of Schar Mountain Center
  phic    = 0.d0,       &   ! Lat of Schar Mountain Center
  cs      = 0.00025d0       ! Wind Shear (shear=1)
                            
  real(8) :: height           ! Model level heights
  real(8) :: sin_tmp, cos_tmp ! Calculation of great circle distance
  real(8) :: r                ! Great circle distance
  real(8) :: zs               ! Surface height
  real(8) :: c                ! Shear

  real(8) :: ueq              ! Reference Velocity
  real(8) :: h0               ! Height of Mountain
  real(8) :: d                ! Mountain Half-Width
  real(8) :: xi               ! Mountain Wavelength

  ueq = mtest_ueq
  h0  = mtest_h0
  d   = mtest_d
  xi  = mtest_xi

  !-----------------------------------------------------------------------
  !    PHIS (surface geopotential)
  !-----------------------------------------------------------------------

  sin_tmp = sin(lat) * sin(phic)
  cos_tmp = cos(lat) * cos(phic)
 
  ! great circle distance with 'a/X'  
  r    = as * ACOS (sin_tmp + cos_tmp*cos(lon-lambdac))
!if mountain
!if ridge, insert ridge

  zs   = h0 * exp(-(r**2)/(d**2))*(cos(pi*r/xi)**2)
  !zs   = h0 * exp(-(r**2)/(d**2))

  phis = g*zs

  !-----------------------------------------------------------------------
  !    SHEAR FLOW OR CONSTANT FLOW
  !-----------------------------------------------------------------------

  if (testid .eq. 3) then
    c = cs
  else
    c = 0.d0
  endif

  !-----------------------------------------------------------------------
  !    TEMPERATURE 
  !-----------------------------------------------------------------------
!modified
  T = Teq *exp( - (c*ueq*ueq/g)*(sin(lat)**2) )

  !-----------------------------------------------------------------------
  !    PS (surface pressure)
  !-----------------------------------------------------------------------
!modified, t-> Teq
  ps = peq*exp( -(ueq*ueq/(2.d0*Rd*Teq))*(sin(lat)**2) - phis/(Rd*Teq)    )

  !-----------------------------------------------------------------------
  !    HEIGHT AND PRESSURE 
  !-----------------------------------------------------------------------

    ! compute the pressure based on the surface pressure and hybrid coefficients
  p = hyam*p0 + hybm*ps
!modified
  height = (Rd*Teq/(g))*log(peq/p) - (ueq*ueq/(2.d0*g))*(sin(lat)**2)
  z = height

  !-----------------------------------------------------------------------
  !    THE VELOCITIES
  !-----------------------------------------------------------------------

  ! Zonal Velocity
!modified
  u = ueq * cos(lat) * sqrt( 1.0d0 + 2.d0*c*height )

  ! Meridional Velocity
  v = 0.d0

  ! Vertical Velocity = Vertical Pressure Velocity = 0
  w = 0.d0

  !-----------------------------------------------------------------------
  !    RHO (density)
  !-----------------------------------------------------------------------

  rho = p/(Rd*T)


END SUBROUTINE mtest_state

END MODULE

