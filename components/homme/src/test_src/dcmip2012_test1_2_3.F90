MODULE dcmip2012_test1_2_3

  !=======================================================================
  !
  !  Functions for setting up initial conditions for the dynamical core tests:
  !
  !  11 - Deformational Advection Test
  !  12 - Hadley Cell Advection Test
  !  13 - Orography Advection Test
  !  20 - Impact of orography on a steady-state at rest
  !  21 and 22 - Non-Hydrostatic Mountain Waves Over A Schaer-Type Mountain without and with vertical wind shear
  !  31 - Non-Hydrostatic Gravity Waves
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
  !
  !
  !  Authors: James Kent, Paul Ullrich, Christiane Jablonowski 
  !		(University of Michigan, dcmip@ucar.edu)
  !          version 4
  !          July/8/2012
  !
  ! Change log: (v3, June/8/2012, v4 July/8/2012, v5 July/20/2012)
  !
  ! v2: bug fixes in the tracer initialization for height-based models
  ! v3: test 3-1: the density is now initialized with the unperturbed background temperature (not the perturbed temperature)
  ! v3: constants converted to double precision
  ! v4: modified tracers in test 1-1, now with cutoff altitudes. Outside of the vertical domain all tracers are set to 0
  ! v4: modified cos-term in vertical velocity (now cos(2 pi t/tau)) in test 1-1, now completing one full up and down cycle
  ! v4: added subroutine test1_advection_orography for test 1-3
  ! v4: added subroutine test2_steady_state_mountain for test 2-0
  ! v4: modified parameter list for tests 2-1 and 2-2 (routine test2_schaer_mountain): addition of parameters hybrid_eta, hyam, hybm
  !     if the logical flag hybrid_eta is true then the pressure in pressure-based model with hybrid sigma-p (eta) coordinates is 
  !     computed internally. In that case the hybrid coefficients hyam and hybm need to be supplied via the parameter list, 
  !     otherwise they are not used.
  ! v5: Change in test 11 - change to u and w, cutoff altitudes (introduced in v4) are removed again
  ! v5: Change in test 12 - velocities multiplies by rho0/rho, different w0 and vertical location of the initial tracer
  !
  !
  !=======================================================================

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

!-----------------------------------------------------------------------
!     Additional constants
!-----------------------------------------------------------------------

	real(8), parameter ::	p0	= 100000.d0		! reference pressure (Pa)


CONTAINS

!==========================================================================================
! TEST CASE 11 - PURE ADVECTION - 3D DEFORMATIONAL FLOW
!==========================================================================================

! The 3D deformational flow test is based on the deformational flow test of Nair and Lauritzen (JCP 2010), 
! with a prescribed vertical wind velocity which makes the test truly 3D. An unscaled planet 
! (with scale parameter X = 1) is selected.

SUBROUTINE test1_advection_deformation (time,lon,lat,p,z,zcoords,u,v,w,t,phis,ps,rho,q1,q2,q3,q4)

IMPLICIT NONE
!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------

  real(8), intent(in)     :: time       ! simulation time (s)
	real(8), intent(in)     :: lon        ! Longitude (radians)
  real(8), intent(in)     :: lat        ! Latitude (radians)
  real(8), intent(in)     :: z          ! Height (m)
	real(8), intent(inout)  :: p          ! Pressure  (Pa)
	integer, intent(in)     :: zcoords    ! 0 or 1 see below
	real(8), intent(out)    :: u          ! Zonal wind (m s^-1)
	real(8), intent(out)    :: v          ! Meridional wind (m s^-1)
	real(8), intent(out)    :: w          ! Vertical Velocity (m s^-1)
	real(8), intent(out)    :: T          ! Temperature (K)
	real(8), intent(out)    :: phis       ! Surface Geopotential (m^2 s^-2)
	real(8), intent(out)    :: ps         ! Surface Pressure (Pa)
	real(8), intent(out)    :: rho        ! density (kg m^-3)
  real(8), intent(out)    :: q1         ! Tracer q1 (kg/kg)
  real(8), intent(out)    :: q2         ! Tracer q2 (kg/kg)
  real(8), intent(out)    :: q3         ! Tracer q3 (kg/kg)
  real(8), intent(out)    :: q4         ! Tracer q4 (kg/kg)

	! if zcoords = 1, then we use z and output p
	! if zcoords = 0, then we use p 

!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
	real(8), parameter ::           &
    tau     = 12.d0 * 86400.d0,   &	! period of motion 12 days
    u0      = (2.d0*pi*a)/tau,    &	! 2 pi a / 12 days
    k0      = (10.d0*a)/tau,      &	! velocity magnitude
    omega0	= (23000.d0*pi)/tau,	&	! velocity magnitude
    T0      = 300.d0,             &	! temperature
    H       = Rd * T0 / g,        &	! scale height
    RR      = 1.d0/2.d0,          &	! horizontal half width divided by 'a'
    ZZ      = 1000.d0,            &	! vertical half width
    z0      = 5000.d0,            &	! center point in z
    lambda0 = 5.d0*pi/6.d0,       &	! center point in longitudes
    lambda1 = 7.d0*pi/6.d0,       &	! center point in longitudes
    phi0    = 0.d0,               &	! center point in latitudes
    phi1    = 0.d0
                            
  real(8) :: height                                                     ! The height of the model levels
  real(8) :: ptop                                                       ! model top in p
  real(8) :: sin_tmp, cos_tmp, sin_tmp2, cos_tmp2                       ! Calculate great circle distances
  real(8) :: d1, d2, r, r2, d3, d4                                      ! For tracer calculations
  real(8) :: s, bs                                                      ! Shape function, and parameter
  real(8) :: lonp                                                       ! Translational longitude, depends on time
  real(8) :: ud                                                         ! Divergent part of u

  !---------------------------------------------------------------------
  !    HEIGHT AND PRESSURE
  !---------------------------------------------------------------------
	
	! height and pressure are aligned (p = p0 exp(-z/H))
	if (zcoords .eq. 1) then
		height = z
		p = p0 * exp(-z/H)

	else
		height = H * log(p0/p)

	endif

	! model top in p
	ptop = p0*exp(-12000.d0/H)

  !---------------------------------------------------------------------
  !    THE VELOCITIES ARE TIME DEPENDENT AND THEREFORE MUST BE UPDATED
  !    IN THE DYNAMICAL CORE
  !---------------------------------------------------------------------

	! translational longitude
	lonp = lon - 2.d0*pi*time/tau

	! shape function
	bs = 0.2
	s = 1.0 + exp( (ptop-p0)/(bs*ptop) ) - exp( (p-p0)/(bs*ptop)) - exp( (ptop-p)/(bs*ptop))

	! zonal velocity
	ud = (omega0*a)/(bs*ptop) * cos(lonp) * (cos(lat)**2.0) * cos(2.0*pi*time/tau) * &
		( - exp( (p-p0)/(bs*ptop)) + exp( (ptop-p)/(bs*ptop))  )

	u = k0*sin(lonp)*sin(lonp)*sin(2.d0*lat)*cos(pi*time/tau) + u0*cos(lat) + ud

	! meridional velocity
	v = k0*sin(2.d0*lonp)*cos(lat)*cos(pi*time/tau)

	! vertical velocity - can be changed to vertical pressure velocity by
	! omega = -(g*p)/(Rd*T0)*w

  w = -((Rd*T0)/(g*p))*omega0*sin(lonp)*cos(lat)*cos(2.0*pi*time/tau)*s

  !-----------------------------------------------------------------------
  !    TEMPERATURE IS CONSTANT 300 K
  !-----------------------------------------------------------------------
	t = T0

  !-----------------------------------------------------------------------
  !    PHIS (surface geopotential) 
  !-----------------------------------------------------------------------
	phis = 0.d0

  !-----------------------------------------------------------------------
  !    PS (surface pressure)
  !-----------------------------------------------------------------------
	ps = p0

  !-----------------------------------------------------------------------
  !    RHO (density)
  !-----------------------------------------------------------------------
	rho = p/(Rd*t)

  !-----------------------------------------------------------------------
  !     initialize Q, set to zero 
  !-----------------------------------------------------------------------
  !	q = 0.d0

  !-----------------------------------------------------------------------
  !     initialize tracers
  !-----------------------------------------------------------------------
	! tracer 1 - Cosine Bells

	sin_tmp = sin(lat) * sin(phi0)
	cos_tmp = cos(lat) * cos(phi0)
	sin_tmp2 = sin(lat) * sin(phi1)
	cos_tmp2 = cos(lat) * cos(phi1)

  ! great circle distance without 'a'
	r  = ACOS (sin_tmp + cos_tmp*cos(lon-lambda0)) 
	r2  = ACOS (sin_tmp2 + cos_tmp2*cos(lon-lambda1)) 
	d1 = min( 1.d0, (r/RR)**2 + ((height-z0)/ZZ)**2 )
	d2 = min( 1.d0, (r2/RR)**2 + ((height-z0)/ZZ)**2 )
	
	q1 = 0.5d0 * (1.d0 + cos(pi*d1)) + 0.5d0 * (1.d0 + cos(pi*d2))

	! tracer 2 - correlated cosine bells
	q2 = 0.9d0 - 0.8d0*q1**2

	! tracer 3 - slotted ellipse

  ! make the ellipse
	if (d1 .le. RR) then
		q3 = 1.d0
	elseif (d2 .le. RR) then
		q3 = 1.d0
	else
		q3 = 0.1d0
	endif

  ! put in the slot
	if (height .gt. z0 .and. abs(lat) .lt. 0.125d0) then
		q3 = 0.1d0
	endif

	! tracer 4: q4 is chosen so that, in combination with the other three tracer
  !           fields with weight (3/10), the sum is equal to one
	q4 = 1.d0 - 0.3d0*(q1+q2+q3)

  !************
  ! change in version 4: added cutoff altitudes, tracers are set to zero outside this region
  ! prevents tracers from being trapped near the bottom and top of the domain
  !************

  ! use a cutoff altitude for
	! tracer 2 3 and 4
	! Set them to zero outside `buffer zone'

  !************
  ! change in version 5: change from v4 reversed, no cutoff altitudes due to continuity equation being satisfied
  !       commented out below
  !************

	!if (height .gt. (z0+1.25*ZZ) .or. height .lt. (z0-1.25*ZZ)) then
	!	q2 = 0.0
	!	q3 = 0.0
	!	q4 = 0.0
	!endif

END SUBROUTINE test1_advection_deformation


!==========================================================================================
! TEST CASE 12 - PURE ADVECTION - 3D HADLEY-LIKE FLOW
!==========================================================================================

SUBROUTINE test1_advection_hadley (time, lon,lat,p,z,zcoords,u,v,w,t,phis,ps,rho,q,q1)

IMPLICIT NONE
!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------

  real(8),  intent(in)    :: time     ! time (s)
	real(8),  intent(in)    :: lon      ! Longitude (radians)
	real(8),  intent(in)    :: lat      ! Latitude (radians)
	real(8),  intent(inout) :: z        ! Height (m)
	real(8),  intent(inout) :: p        ! Pressure  (Pa)
	integer,  intent(in)    :: zcoords 	! 0 or 1 see below
	real(8), intent(out)    :: u        ! Zonal wind (m s^-1)
	real(8), intent(out)    :: v        ! Meridional wind (m s^-1)
	real(8), intent(out)    :: w        ! Vertical Velocity (m s^-1)
	real(8), intent(out)    :: T        ! Temperature (K)
	real(8), intent(out)    :: phis     ! Surface Geopotential (m^2 s^-2)
	real(8), intent(out)    :: ps       ! Surface Pressure (Pa)
  real(8), intent(out)    :: rho      ! density (kg m^-3)
  real(8), intent(out)    :: q        ! Specific Humidity (kg/kg)
	real(8), intent(out)    :: q1       ! Tracer q1 (kg/kg)

	! if zcoords = 1, then we use z and output p
	! if zcoords = 0, then we use p

!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 

	real(8), parameter :: &
    tau     = 1.d0 * 86400.d0,	&	! period of motion 1 day (in s)
    u0      = 40.d0,            &	! Zonal velocity magnitude (m/s)
    w0      = 0.15d0,           &	! Vertical velocity magnitude (m/s)
    T0      = 300.d0,           &	! temperature (K)
    H       = Rd * T0 / g,      &	! scale height
    K       = 5.d0,             &	! number of Hadley-like cells
    z1      = 2000.d0,          &	! position of lower tracer bound (m)
    z2      = 5000.d0,          &	! position of upper tracer bound (m)
    z0      = 0.5d0*(z1+z2),    &	! midpoint (m)
    ztop    = 12000.d0            ! model top (m)
                            
  real(8) :: rho0                 ! reference density at z=0 m
  real(8) :: height               ! Model level heights

!-----------------------------------------------------------------------
!    HEIGHT AND PRESSURE
!-----------------------------------------------------------------------

	! Height and pressure are aligned (p = p0 exp(-z/H))

	if (zcoords .eq. 1) then
		height = z
		p = p0 * exp(-z/H)
!print *,'using z to get p'

	else
		height = H * log(p0/p)
    z = height
!print *,'using p to get z'

	endif

!-----------------------------------------------------------------------
!    TEMPERATURE IS CONSTANT 300 K
!-----------------------------------------------------------------------

	t = T0

!-----------------------------------------------------------------------
!    PHIS (surface geopotential)
!-----------------------------------------------------------------------
      
	phis = 0.d0

!-----------------------------------------------------------------------
!    PS (surface pressure)
!-----------------------------------------------------------------------

	ps = p0

!-----------------------------------------------------------------------
!    RHO (density)
!-----------------------------------------------------------------------

	rho = p/(Rd*T)
	rho0 = p0/(Rd*T)

!-----------------------------------------------------------------------
!    THE VELOCITIES ARE TIME DEPENDENT AND THEREFORE MUST BE UPDATED
!    IN THE DYNAMICAL CORE
!-----------------------------------------------------------------------

	! Zonal Velocity
	u = u0*cos(lat)

	! Meridional Velocity
	v = -(rho0/rho) * (a*w0*pi)/(K*ztop) *cos(lat)*sin(K*lat)*cos(pi*height/ztop)*cos(pi*time/tau)

	! Vertical Velocity - can be changed to vertical pressure velocity by omega = -g*rho*w
	w = (rho0/rho) *(w0/K)*(-2.d0*sin(K*lat)*sin(lat) + K*cos(lat)*cos(K*lat)) &
		*sin(pi*height/ztop)*cos(pi*time/tau)

!-----------------------------------------------------------------------
!     initialize Q, set to zero 
!-----------------------------------------------------------------------

	q = 0.d0

!-----------------------------------------------------------------------
!     initialize tracers
!-----------------------------------------------------------------------

	! Tracer 1 - Layer

	if (height .lt. z2 .and. height .gt. z1) then
	
		q1 = 0.5d0 * (1.d0 + cos( 2.d0*pi*(height-z0)/(z2-z1) ) )

	else

		q1 = 0.d0

	endif

END SUBROUTINE test1_advection_hadley



!==========================================================================================
! TEST CASE 13 - HORIZONTAL ADVECTION OF THIN CLOUD-LIKE TRACERS IN THE PRESENCE OF OROGRAPHY
!==========================================================================================

SUBROUTINE test1_advection_orography (lon,lat,p,z,zcoords,cfv,hybrid_eta,hyam,hybm,gc,u,v,w,t,phis,ps,rho,q1,q2,q3,q4)

IMPLICIT NONE

	real(8), intent(in)  :: lon             ! Longitude (radians)
  real(8), intent(in)  :: lat             ! Latitude (radians)
  real(8), intent(in)  :: hyam            ! A coefficient for hybrid-eta coordinate, at model level midpoint
  real(8), intent(in)  :: hybm            ! B coefficient for hybrid-eta coordinate, at model level midpoint
  real(8), intent(in)  :: gc              ! bar{z} for Gal-Chen coordinate

	logical, intent(in)  :: hybrid_eta      ! flag to indicate whether the hybrid sigma-p (eta) coordinate is used
                                          ! if set to .true., then the pressure will be computed via the
                                          !    hybrid coefficients hyam and hybm, they need to be initialized
                                          ! if set to .false. (for pressure-based models): the pressure is already pre-computed
                                          !    and is an input value for this routine
                                          ! for height-based models: pressure will always be computed based on the height and
                                          !    hybrid_eta is not used
                                          ! Note that we only use hyam and hybm for the hybrid-eta coordiantes, and we only use
                                          ! gc for the Gal-Chen coordinates. If not required then they become dummy variables

	real(8), intent(inout)  :: p            ! Pressure  (Pa)
  real(8), intent(inout)  :: z               ! Height (m)

	integer, intent(in)     :: zcoords      ! 0 or 1 see below
	integer, intent(in)     :: cfv          ! 0, 1 or 2 see below
	real(8), intent(out)    :: u            ! Zonal wind (m s^-1)
  real(8), intent(out)    :: v            ! Meridional wind (m s^-1)
  real(8), intent(out)    :: w            ! Vertical Velocity (m s^-1)
  real(8), intent(out)    :: t            ! Temperature (K)
  real(8), intent(out)    :: phis         ! Surface Geopotential (m^2 s^-2)
  real(8), intent(out)    :: ps           ! Surface Pressure (Pa)
  real(8), intent(out)    :: rho          ! density (kg m^-3)
  real(8), intent(out)    :: q1           ! Tracer q1 (kg/kg)
  real(8), intent(out)    :: q2           ! Tracer q2 (kg/kg)
  real(8), intent(out)    :: q3           ! Tracer q3 (kg/kg)
  real(8), intent(out)    :: q4           ! Tracer q4 (kg/kg)

	! if zcoords = 1, then we use z and output p
	! if zcoords = 0, then we use p

	! if cfv = 0 we assume that our horizontal velocities are not coordinate following
	! if cfv = 1 then our velocities follow hybrid eta coordinates and we need to specify w
	! if cfv = 2 then our velocities follow Gal-Chen coordinates and we need to specify w

	! In hybrid-eta coords: p = hyam p0 + hybm ps
	! In Gal-Chen coords: z = zs + (gc/ztop)*(ztop - zs)

  ! if other orography-following coordinates are used, the w wind needs to be newly derived for them

!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
	real(8), parameter :: 	tau     = 12.d0 * 86400.d0,	&	! period of motion 12 days (s)
    u0      = 2.d0*pi*a/tau,    &	! Velocity Magnitude (m/s)
    T0      = 300.d0,           &	! temperature (K)
    H       = Rd * T0 / g,      &	! scale height (m)
    alpha   = pi/6.d0,          &	! rotation angle (radians), 30 degrees
    lambdam = 3.d0*pi/2.d0,     &	! mountain longitude center point (radians)
    !lambdam = 2.d0*pi/2.d0,     &	! mountain longitude center point (radians)
    phim    = 0.d0,             &	! mountain latitude center point (radians)
    h0      = 2000.d0,          &	! peak height of the mountain range (m)
    !h0      = 4000.d0,          &	! peak height of the mountain range (m)
    Rm      = 3.d0*pi/4.d0,     &	! mountain radius (radians)
    zetam   = pi/16.d0,         &	! mountain oscillation half-width (radians)
    !zetam   = pi,               &	! mountain oscillation half-width (radians)
    lambdap = pi/2.d0,          &	! cloud-like tracer longitude center point (radians)
    phip    = 0.d0,             &	! cloud-like tracer latitude center point (radians)
    Rp      = pi/4.d0,          &	! cloud-like tracer radius (radians)
    zp1     = 3050.d0,          &	! midpoint of first (lowermost) tracer (m)
    zp2     = 5050.d0,          &	! midpoint of second tracer (m)
    zp3     = 8200.d0,          &	! midpoint of third (topmost) tracer (m)
    dzp1    = 1000.d0,          &	! thickness of first (lowermost) tracer (m)
    dzp2    = 1000.d0,          &	! thickness of second tracer (m)
    dzp3    = 400.d0,           &	! thickness of third (topmost) tracer (m)
    ztop    = 12000.d0            ! model top (m)
                            
    real(8) :: height							! Model level heights (m)
    real(8) :: r                  ! Great circle distance (radians)
    real(8) :: rz                 ! height differences
    real(8) :: zs                 ! Surface elevation (m)

!-----------------------------------------------------------------------
!    PHIS (surface geopotential)
!-----------------------------------------------------------------------
      
	r = acos( sin(phim)*sin(lat) + cos(phim)*cos(lat)*cos(lon - lambdam) )

	if (r .lt. Rm) then

		zs = (h0/2.d0)*(1.d0+cos(pi*r/Rm))*cos(pi*r/zetam)**2.d0

	else

		zs = 0.d0

	endif

	phis = g*zs

!-----------------------------------------------------------------------
!    PS (surface pressure)
!-----------------------------------------------------------------------

	ps = p0 * exp(-zs/H)

!-----------------------------------------------------------------------
!    HEIGHT AND PRESSURE
!-----------------------------------------------------------------------

	! Height and pressure are aligned (p = p0 exp(-z/H))

	if (zcoords .eq. 1) then

		height = z
		p = p0 * exp(-z/H)

	else

    if (hybrid_eta) p = hyam*p0 + hybm*ps   ! compute the pressure based on the surface pressure and hybrid coefficients
		height = H * log(p0/p)
    z = height

	endif

!-----------------------------------------------------------------------
!    THE VELOCITIES ARE TIME INDEPENDENT 
!-----------------------------------------------------------------------

	! Zonal Velocity

	u = u0*(cos(lat)*cos(alpha)+sin(lat)*cos(lon)*sin(alpha))

	! Meridional Velocity

	v = -u0*(sin(lon)*sin(alpha))

!-----------------------------------------------------------------------
!    TEMPERATURE IS CONSTANT 300 K
!-----------------------------------------------------------------------

	t = T0

!-----------------------------------------------------------------------
!    RHO (density)
!-----------------------------------------------------------------------

	rho = p/(Rd*t)

!-----------------------------------------------------------------------
!     initialize Q, set to zero 
!-----------------------------------------------------------------------

	!q = 0.d0

!-----------------------------------------------------------------------
!    VERTICAL VELOCITY IS TIME INDEPENDENT 
!-----------------------------------------------------------------------

	! Vertical Velocity - can be changed to vertical pressure velocity by 
	! omega = -(g*p)/(Rd*T0)*w

	! NOTE that if orography-following coordinates are used then the vertical 
	! velocity needs to be translated into the new coordinate system due to
	! the variation of the height along coordinate surfaces
	! See section 1.3 and the appendix of the test case document

	if (cfv .eq. 0) then

		! if the horizontal velocities do not follow the vertical coordinate

		w = 0.d0

	elseif (cfv .eq. 1) then

		! if the horizontal velocities follow hybrid eta coordinates then
		! the perceived vertical velocity is

		call test1_advection_orograph_hybrid_eta_velocity(w)

	elseif (cfv .eq. 2) then

		! if the horizontal velocities follow Gal Chen coordinates then
		! the perceived vertical velocity is	

		call test1_advection_orograph_Gal_Chen_velocity(w)	

!        else
!               compute your own vertical velocity if other orography-following 
!               vertical coordinate is used
!
	endif



!-----------------------------------------------------------------------
!     initialize tracers
!-----------------------------------------------------------------------

	! Tracer 1 - Cloud Layer
      
	r = acos( sin(phip)*sin(lat) + cos(phip)*cos(lat)*cos(lon - lambdap) )

	rz = abs(height - zp1)

	if (rz .lt. 0.5d0*dzp1 .and. r .lt. Rp) then

		q1 = 0.25d0*(1.d0+cos(2.d0*pi*rz/dzp1))*(1.d0+cos(pi*r/Rp))

	else

		q1 = 0.d0

	endif

	rz = abs(height - zp2)

	if (rz .lt. 0.5d0*dzp2 .and. r .lt. Rp) then

		q2 = 0.25d0*(1.d0+cos(2.d0*pi*rz/dzp2))*(1.d0+cos(pi*r/Rp))

	else

		q2 = 0.d0

	endif

	rz = abs(height - zp3)

	if (rz .lt. 0.5d0*dzp3 .and. r .lt. Rp) then

		q3 = 1.d0

	else

		q3 = 0.d0

	endif

	q4 = q1 + q2 + q3


	CONTAINS

	!-----------------------------------------------------------------------
	!    SUBROUTINE TO CALCULATE THE PERCEIVED VERTICAL VELOCITY 
	!    		UNDER HYBRID-ETA COORDINATES
	!-----------------------------------------------------------------------

	SUBROUTINE test1_advection_orograph_hybrid_eta_velocity(w)
	IMPLICIT NONE
	real(8), intent(out) ::	w

	real(8) :: 	press,  &		! hyam *p0 + hybm *ps
			r,              &		! Great Circle Distance
			dzsdx,          &		! Part of surface height derivative
			dzsdlambda,     & 	! Derivative of zs w.r.t lambda
			dzsdphi,        &   ! Derivative of zs w.r.t phi
			dzdlambda,      &   ! Derivative of z w.r.t lambda
			dzdphi,         &   ! Derivative of z w.r.t phi
			dpsdlambda,     &   ! Derivative of ps w.r.t lambda
			dpsdphi             ! Derivative of ps w.r.t phi

		! Calculate pressure and great circle distance to mountain center

		press = hyam*p0 + hybm*ps

		r = acos( sin(phim)*sin(lat) + cos(phim)*cos(lat)*cos(lon - lambdam) )

		! Derivatives of surface height

		if (r .lt. Rm) then
			dzsdx = -h0*pi/(2.d0*Rm)*sin(pi*r/Rm)*cos(pi*r/zetam)**2 - &
				(h0*pi/zetam)*(1.d0+cos(pi*r/Rm))*cos(pi*r/zetam)*sin(pi*r/zetam)
		else
			dzsdx = 0.d0
		endif

		! Prevent division by zero

		if (1.d0-cos(r)**2 .gt. 0.d0) then
			dzsdlambda = dzsdx * (cos(phim)*cos(lat)*sin(lon-lambdam)) &
					/sqrt(1.d0-cos(r)**2)
			dzsdphi    = dzsdx * (-sin(phim)*cos(lat) + cos(phim)*sin(lat)*cos(lon-lambdam)) & 
					/sqrt(1.d0-cos(r)**2)
		else
			dzsdlambda = 0.d0
			dzsdphi    = 0.d0
		endif

		! Derivatives of surface pressure

		dpsdlambda = -(g*p0/(Rd*T0))*exp(-g*zs/(Rd*T0))*dzsdlambda
		dpsdphi    = -(g*p0/(Rd*T0))*exp(-g*zs/(Rd*T0))*dzsdphi

		! Derivatives of coordinate height

		dzdlambda = -(Rd*T0/(g*press))*hybm*dpsdlambda
		dzdphi    = -(Rd*T0/(g*press))*hybm*dpsdphi

		! Prevent division by zero

		if (abs(lat) .lt. pi/2.d0) then
			w = - (u/(a*cos(lat)))*dzdlambda - (v/a)*dzdphi
		else
			w = 0.d0
		endif

	END SUBROUTINE test1_advection_orograph_hybrid_eta_velocity



	!-----------------------------------------------------------------------
	!    SUBROUTINE TO CALCULATE THE PERCEIVED VERTICAL VELOCITY 
	!    		UNDER GAL-CHEN COORDINATES
	!-----------------------------------------------------------------------

	SUBROUTINE test1_advection_orograph_Gal_Chen_velocity(w)
	IMPLICIT NONE
	real(8), intent(out) ::	w

	real(8) :: 	r, &			! Great Circle Distance
			dzsdx, &		! Part of surface height derivative
			dzsdlambda, & 		! Derivative of zs w.r.t lambda
			dzsdphi, & 		! Derivative of zs w.r.t phi
			dzdlambda, & 		! Derivative of z w.r.t lambda
			dzdphi	 		! Derivative of z w.r.t phi

		! Calculate great circle distance to mountain center

		r = acos( sin(phim)*sin(lat) + cos(phim)*cos(lat)*cos(lon - lambdam) )

		! Derivatives of surface height

		if (r .lt. Rm) then
			dzsdx = -h0*pi/(2.d0*Rm)*sin(pi*r/Rm)*cos(pi*r/zetam)**2 - &
				(h0*pi/zetam)*(1.d0+cos(pi*r/Rm))*cos(pi*r/zetam)*sin(pi*r/zetam)
		else
			dzsdx = 0.d0
		endif

		! Prevent division by zero

		if (1.d0-cos(r)**2 .gt. 0.d0) then
			dzsdlambda = dzsdx * (cos(phim)*cos(lat)*sin(lon-lambdam)) &
					/sqrt(1.d0-cos(r)**2)
			dzsdphi    = dzsdx * (-sin(phim)*cos(lat) + cos(phim)*sin(lat)*cos(lon-lambdam)) & 
					/sqrt(1.d0-cos(r)**2)
		else
			dzsdlambda = 0.d0
			dzsdphi    = 0.d0
		endif

		! Derivatives of coordinate height

		dzdlambda = (1.d0-gc/ztop)*dzsdlambda
		dzdphi    = (1.d0-gc/ztop)*dzsdphi

		! Prevent division by zero

		if (abs(lat) .lt. pi/2.d0) then
			w = - (u/(a*cos(lat)))*dzdlambda - (v/a)*dzdphi
		else
			w = 0.d0
		endif

	END SUBROUTINE test1_advection_orograph_Gal_Chen_velocity

END SUBROUTINE test1_advection_orography

!==========================================================================================
! TEST CASE 2X - IMPACT OF OROGRAPHY ON A NON-ROTATING PLANET
!==========================================================================================
! The tests in section 2-x examine the impact of 3D Schaer-like circular mountain profiles on an 
! atmosphere at rest (2-0), and on flow fields with wind shear (2-1) and without vertical wind shear (2-2).  
! A non-rotating planet is used for all configurations. Test 2-0 is conducted on an unscaled regular-size 
! planet and primarily examines the accuracy of the pressure gradient calculation in a steady-state 
! hydrostatically-balanced atmosphere at rest. This test is especially appealing for models with 
! orography-following vertical coordinates. It increases the complexity of test 1-3, that investigated 
! the impact of the same Schaer-type orographic profile on the accuracy of purely-horizontal passive 
! tracer advection.
!
! Tests 2-1 and 2-2 increase the complexity even further since non-zero flow fields are now prescribed 
! with and without vertical wind shear. In order to trigger non-hydrostatic responses the two tests are 
! conducted on a reduced-size planet with reduction factor $X=500$ which makes the horizontal and 
! vertical grid spacing comparable. This test clearly discriminates between non-hydrostatic and hydrostatic 
! models since the expected response is in the non-hydrostatic regime. Therefore, the flow response is 
! captured differently by hydrostatic models.

!=========================================================================
! Test 2-0:  Steady-State Atmosphere at Rest in the Presence of Orography
!=========================================================================
SUBROUTINE test2_steady_state_mountain (lon,lat,p,z,zcoords,hybrid_eta,hyam,hybm,u,v,w,t,phis,ps,rho,q)

IMPLICIT NONE

	real(8),  intent(in)    :: lon        ! Longitude (radians)
  real(8),  intent(in)    :: lat        ! Latitude (radians)
  real(8),  intent(in)    :: z          ! Height (m)
  real(8),  intent(in)    :: hyam       ! A coefficient for hybrid-eta coordinate, at model level midpoint
  real(8),  intent(in)    :: hybm       ! B coefficient for hybrid-eta coordinate, at model level midpoint
	logical,  intent(in)    :: hybrid_eta ! if true,p from hybrid coeff. if false, p is precomputed
  real(8),  intent(inout) :: p          ! Pressure  (Pa)
	integer,  intent(in)    :: zcoords    ! 0 or 1 see below
  real(8),  intent(out)   :: u          ! Zonal wind (m s^-1)
  real(8),  intent(out)   :: v          ! Meridional wind (m s^-1)
  real(8),  intent(out)   :: w          ! Vertical Velocity (m s^-1)
  real(8),  intent(out)   :: t          ! Temperature (K)
  real(8),  intent(out)   :: phis       ! Surface Geopotential (m^2 s^-2)
  real(8),  intent(out)   :: ps         ! Surface Pressure (Pa)
  real(8),  intent(out)   :: rho        ! density (kg m^-3)
  real(8),  intent(out)   :: q          ! Specific Humidity (kg/kg)

	! if zcoords = 1, then we use z and output p
	! if zcoords = 0, then we compute or use p
	! In hybrid-eta coords: p = hyam p0 + hybm ps

  !-----------------------------------------------------------------------
  !     test case parameters
  !----------------------------------------------------------------------- 
	real(8), parameter ::&
      T0      = 300.d0,       &	! temperature (K)
      gamma   = 0.0065d0,     & ! temperature lapse rate (K/m)
      lambdam = 3.d0*pi/2.d0,	&	! mountain longitude center point (radians)
      phim    = 0.d0,         &	! mountain latitude center point (radians)
      h0      = 2000.d0,      &	! peak height of the mountain range (m)
      Rm      = 3.d0*pi/4.d0,	&	! mountain radius (radians)
      zetam   = pi/16.d0,     &	! mountain oscillation half-width (radians)
      !zetam   = Rm,           &	! mountain oscillation half-width (radians)
      ztop    = 12000.d0        ! model top (m)
                            
  real(8) :: height         ! Model level heights (m)
  real(8) :: r							! Great circle distance (radians)
  real(8) :: zs							! Surface elevation (m)
  real(8) :: exponent       ! exponent: g/(Rd * gamma)
  real(8) :: exponent_rev   ! reversed exponent

  exponent     = g/(Rd*gamma)
  exponent_rev = 1.d0/exponent

  !-----------------------------------------------------------------------
  !    PHIS (surface geopotential)
  !-----------------------------------------------------------------------

	r   = acos( sin(phim)*sin(lat) + cos(phim)*cos(lat)*cos(lon - lambdam) )

  zs  = 0.d0
  if (r < Rm) zs=(h0/2.d0)*(1.d0+cos(pi*r/Rm))*cos(pi*r/zetam)**2.d0    ! mountain height

	phis = g*zs

  !-----------------------------------------------------------------------
  !    PS (surface pressure)
  !-----------------------------------------------------------------------

	ps = p0 * (1.d0 - gamma/T0*zs)**exponent

  !-----------------------------------------------------------------------
  !    HEIGHT AND PRESSURE
  !-----------------------------------------------------------------------

	! Height and pressure are aligned (p = p0 * (1.d0 - gamma/T0*z)**exponent)

	if (zcoords .eq. 1) then
		height = z
		p = p0 * (1.d0 - gamma/T0*z)**exponent

	else

    if (hybrid_eta) p = hyam*p0 + hybm*ps              ! compute the pressure based on the surface pressure and hybrid coefficients
		height = T0/gamma * (1.d0 - (p/p0)**exponent_rev)  ! compute the height at this pressure

	endif

!-----------------------------------------------------------------------
!    THE VELOCITIES ARE ZERO (STATE AT REST)
!-----------------------------------------------------------------------

	u = 0.d0    ! Zonal Velocity
	v = 0.d0    ! Meridional Velocity
  w = 0.d0    ! Vertical Velocity

!-----------------------------------------------------------------------
!    TEMPERATURE WITH CONSTANT LAPSE RATE
!-----------------------------------------------------------------------

	t = T0 - gamma*height

!-----------------------------------------------------------------------
!    RHO (density)
!-----------------------------------------------------------------------

	rho = p/(Rd*t)

!-----------------------------------------------------------------------
!     initialize Q, set to zero 
!-----------------------------------------------------------------------

	q = 0.d0

END SUBROUTINE test2_steady_state_mountain



!=====================================================================================
! Tests 2-1 and 2-2:  Non-hydrostatic Mountain Waves over a Schaer-type Mountain
!=====================================================================================

SUBROUTINE test2_schaer_mountain (lon,lat,p,z,zcoords,hybrid_eta,hyam,hybm,shear,u,v,w,t,phis,ps,rho,q)

 ! use control_mod, only: dcmip2_ueq, dcmip2_h0, dcmip2_d, dcmip2_xi

IMPLICIT NONE
!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------

	real(8), intent(in)   :: lon          ! Longitude (radians)
  real(8), intent(in)   :: lat          ! Latitude (radians)
  real(8), intent(inout):: z            ! Height (m)
  real(8), intent(in)   :: hyam         ! A coefficient for hybrid-eta coordinate, at model level midpoint
  real(8), intent(in)   :: hybm         ! B coefficient for hybrid-eta coordinate, at model level midpoint
	logical, intent(in)   :: hybrid_eta   ! flag to indicate whether the hybrid sigma-p (eta) coordinate is used
  real(8), intent(inout):: p            ! Pressure  (Pa)
	integer, intent(in)   :: zcoords      ! 0 or 1 see below
  integer, intent(in)   :: shear        ! 0 or 1 see below
	real(8), intent(out)  :: u            ! Zonal wind (m s^-1)
  real(8), intent(out)  :: v            ! Meridional wind (m s^-1)
  real(8), intent(out)  :: w            ! Vertical Velocity (m s^-1)
  real(8), intent(out)  :: T            ! Temperature (K)
  real(8), intent(out)  :: phis         ! Surface Geopotential (m^2 s^-2)
  real(8), intent(out)  :: ps           ! Surface Pressure (Pa)
  real(8), intent(out)  :: rho          ! density (kg m^-3)
  real(8), intent(out)  :: q            ! Specific Humidity (kg/kg)

	! if zcoords  = 1, then we use z and output p
	! if zcoords  = 0, then we either compute or use p
	! if shear    = 1, then we use shear flow
	! if shear    = 0, then we use constant u

!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
	real(8), parameter ::   &
    X       = 500.d0,     &   ! Reduced Earth reduction factor
    Om      = 0.d0,       &   ! Rotation Rate of Earth
    as      = a/X,        &   ! New Radius of small Earth
    Teq     = 300.d0,     &   ! Temperature at Equator
    Peq     = 100000.d0,	&   ! Reference PS at Equator
    ztop    = 30000.d0,		&   ! Model Top
    lambdac = pi/4.d0, 		&   ! Lon of Schar Mountain Center
    phic    = 0.d0,       & 	! Lat of Schar Mountain Center
    cs      = 0.00025d0 		 	! Wind Shear (shear=1)
                            
  real(8) :: height						! Model level heights
  real(8) :: sin_tmp, cos_tmp ! Calculation of great circle distance
  real(8) :: r                ! Great circle distance
  real(8) :: zs               ! Surface height
  real(8) :: c                ! Shear

  real(8) :: ueq              ! Reference Velocity
  real(8) :: h0               ! Height of Mountain
  real(8) :: d                ! Mountain Half-Width
  real(8) :: xi               ! Mountain Wavelength

stop
  ! Set orography parameters from namelist variables
  !ueq = dcmip2_ueq
  !h0  = dcmip2_h0
  !d   = dcmip2_d
  !xi  = dcmip2_xi

  !-----------------------------------------------------------------------
  !    PHIS (surface geopotential)
  !-----------------------------------------------------------------------

	sin_tmp = sin(lat) * sin(phic)
	cos_tmp = cos(lat) * cos(phic)
	
	! great circle distance with 'a/X'  

	r    = as * ACOS (sin_tmp + cos_tmp*cos(lon-lambdac))
	zs   = h0 * exp(-(r**2)/(d**2))*(cos(pi*r/xi)**2)
  !zs   = h0 * exp(-(r**2)/(d**2))

	phis = g*zs

  !-----------------------------------------------------------------------
  !    SHEAR FLOW OR CONSTANT FLOW
  !-----------------------------------------------------------------------

	if (shear .eq. 1) then

		c = cs

	else

		c = 0.d0

	endif

  !-----------------------------------------------------------------------
  !    TEMPERATURE 
  !-----------------------------------------------------------------------

	T = Teq *(1.d0 - (c*ueq*ueq/(g))*(sin(lat)**2) )

  !-----------------------------------------------------------------------
  !    PS (surface pressure)
  !-----------------------------------------------------------------------

	ps = peq*exp( -(ueq*ueq/(2.d0*Rd*Teq))*(sin(lat)**2) - phis/(Rd*t)    )

  !-----------------------------------------------------------------------
  !    HEIGHT AND PRESSURE 
  !-----------------------------------------------------------------------

	if (zcoords .eq. 1) then

		height = z
		p = peq*exp( -(ueq*ueq/(2.d0*Rd*Teq))*(sin(lat)**2) - g*height/(Rd*T)    )

	else

    ! compute the pressure based on the surface pressure and hybrid coefficients
    if (hybrid_eta) p = hyam*p0 + hybm*ps
		height = (Rd*T/(g))*log(peq/p) - (T*ueq*ueq/(2.d0*Teq*g))*(sin(lat)**2)
    z      = height
	endif

  !-----------------------------------------------------------------------
  !    THE VELOCITIES
  !-----------------------------------------------------------------------

	! Zonal Velocity

	u = ueq * cos(lat) * sqrt( (2.d0*Teq/T)*c*height + T/(Teq) )

	! Meridional Velocity

	v = 0.d0

	! Vertical Velocity = Vertical Pressure Velocity = 0

	w = 0.d0

  !-----------------------------------------------------------------------
  !    RHO (density)
  !-----------------------------------------------------------------------

	rho = p/(Rd*T)

  !-----------------------------------------------------------------------
  !     initialize Q, set to zero 
  !-----------------------------------------------------------------------

	q = 0.d0

END SUBROUTINE test2_schaer_mountain

!==========================================================================================
! TEST CASE 3 - GRAVITY WAVES
!==========================================================================================

! The non-hydrostatic gravity wave test examines the response of models to short time-scale wavemotion
! triggered by a localized perturbation. The formulation presented in this document is new,
! but is based on previous approaches by Skamarock et al. (JAS 1994), Tomita and Satoh (FDR 2004), and
! Jablonowski et al. (NCAR Tech Report 2008) 

!==========
! Test 3-1
!==========
SUBROUTINE test3_gravity_wave (lon,lat,p,z,zcoords,hybrid_eta,hyam,hybm,u,v,w,t,t_mean,phis,ps,rho,rho_mean,q)

IMPLICIT NONE
!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------

	real(8), intent(in)     :: lon        ! Longitude (radians)
  real(8), intent(in)     :: lat        ! Latitude (radians)
  real(8), intent(inout)  :: z          ! Height (m)
  real(8), intent(in)     :: hyam       ! A coefficient for hybrid-eta coordinate, at model level midpoint
  real(8), intent(in)     :: hybm       ! B coefficient for hybrid-eta coordinate, at model level midpoint
  logical, intent(in)     :: hybrid_eta ! flag to indicate whether the hybrid sigma-p (eta) coordinate is used
	real(8), intent(inout)  :: p          ! Pressure  (Pa)
	integer, intent(in)     :: zcoords    ! 0 or 1 see below
	real(8), intent(out)    :: u          ! Zonal wind (m s^-1)
  real(8), intent(out)    :: v          ! Meridional wind (m s^-1)
  real(8), intent(out)    :: w          ! Vertical Velocity (m s^-1)
  real(8), intent(out)    :: t          ! Temperature (K)
  real(8), intent(out)    :: t_mean     ! Temperature (K)
  real(8), intent(out)    :: phis       ! Surface Geopotential (m^2 s^-2)
  real(8), intent(out)    :: ps         ! Surface Pressure (Pa)
  real(8), intent(out)    :: rho        ! density (kg m^-3)
  real(8), intent(out)    :: rho_mean   ! density (kg m^-3)
  real(8), intent(out)    :: q          ! Specific Humidity (kg/kg)

	! if zcoords = 1, then we use z and output z
	! if zcoords = 0, then we use p

!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
	real(8), parameter :: &
    X       = 125.d0,       &	! Reduced Earth reduction factor
    Om      = 0.d0,         &	! Rotation Rate of Earth
    as      = a/X,          &	! New Radius of small Earth
    u0      = 20.d0,        &	! Reference Velocity
    Teq     = 300.d0,       &	! Temperature at Equator
    Peq     = 100000.d0,		&	! Reference PS at Equator
    ztop    = 10000.d0,     &	! Model Top
    lambdac = 2.d0*pi/3.d0, & ! Lon of Pert Center
    d       = 5000.d0,      & ! Width for Pert
    phic    = 0.d0,         & ! Lat of Pert Center
    delta_theta = 1.d0,     & ! Max Amplitude of Pert
    Lz      = 20000.d0, 		& ! Vertical Wavelength of Pert
    N       = 0.01d0,       & ! Brunt-Vaisala frequency
    N2      = N*N,          &	! Brunt-Vaisala frequency Squared
    bigG    = (g*g)/(N2*cp)   ! Constant
                            
  real(8) :: height           ! Model level height
  real(8) :: sin_tmp, cos_tmp ! Calculation of great circle distance
  real(8) :: r, s							! Shape of perturbation
  real(8) :: Ts 							! Surface temperature
  real(8) :: t_pert           ! temperature perturbation
  real(8) :: theta_pert       ! Pot-temp perturbation

!-----------------------------------------------------------------------
!    THE VELOCITIES
!-----------------------------------------------------------------------

	! Zonal Velocity 

	u = u0 * cos(lat)

	! Meridional Velocity

	v = 0.d0

	! Vertical Velocity = Vertical Pressure Velocity = 0

	w = 0.d0

!-----------------------------------------------------------------------
!    PHIS (surface geopotential)
!-----------------------------------------------------------------------

	phis = 0.d0

!-----------------------------------------------------------------------
!    SURFACE TEMPERATURE 
!-----------------------------------------------------------------------

	Ts = bigG + (Teq-bigG)*exp( -(u0*N2/(4.d0*g*g))*(u0+2.d0*Om*as)*(cos(2.d0*lat)-1.d0)    )

!-----------------------------------------------------------------------
!    PS (surface pressure)
!-----------------------------------------------------------------------

	ps = peq*exp((u0/(4.0*bigG*Rd))*(u0+2.0*Om*as)*(cos(2.0*lat)-1.0))*(Ts/Teq)**(cp/Rd)

!-----------------------------------------------------------------------
!    HEIGHT AND PRESSURE AND MEAN TEMPERATURE
!-----------------------------------------------------------------------

	if (zcoords .eq. 1) then

		height = z
		p = ps*( (bigG/TS)*exp(-N2*height/g)+1.d0 - (bigG/TS)  )**(cp/Rd)

	else

    if (hybrid_eta) p = hyam*p0 + hybm*ps
		height = (-g/N2)*log( (TS/bigG)*( (p/ps)**(Rd/cp) - 1.d0  ) + 1.d0 )
    z      = height

	endif

	t_mean = bigG*(1.d0 - exp(N2*height/g))+ TS*exp(N2*height/g)

!-----------------------------------------------------------------------
!    rho (density), unperturbed using the background temperature t_mean
!-----------------------------------------------------------------------

!***********
!       change in version 3: density is now initialized with unperturbed background temperature,
!                            temperature perturbation is added afterwards
!***********
	rho_mean = p/(Rd*t_mean)

!-----------------------------------------------------------------------
!    POTENTIAL TEMPERATURE PERTURBATION, 
!    here: converted to temperature and added to the temperature field
!    models with a prognostic potential temperature field can utilize
!    the potential temperature perturbation theta_pert directly and add it
!    to the background theta field (not included here)
!-----------------------------------------------------------------------

	sin_tmp = sin(lat) * sin(phic)
	cos_tmp = cos(lat) * cos(phic)

  ! great circle distance with 'a/X'

	r  = as * ACOS (sin_tmp + cos_tmp*cos(lon-lambdac)) 

	s = (d**2)/(d**2 + r**2)

	theta_pert = delta_theta*s*sin(2.d0*pi*height/Lz)

	t_pert = theta_pert*(p/p0)**(Rd/cp)

	t = t_mean + t_pert

  rho = p/(Rd*t)

!-----------------------------------------------------------------------
!     initialize Q, set to zero 
!-----------------------------------------------------------------------

	q = 0.d0

END SUBROUTINE test3_gravity_wave


END MODULE

