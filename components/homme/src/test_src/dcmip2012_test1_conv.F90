module dcmip2012_test1_conv
  implicit none

contains

SUBROUTINE test1_conv_advection_deformation (time,lon,lat,p,z,zcoords,u,v,w,t,phis,ps,rho,q1,q2,q3,q4)
!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------

  ! Use physical constants consistent with HOMME
  use physical_constants, only: a=>rearth0, Rd => Rgas, g, cp, pi=>dd_pi, p0

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
    omega0	= (2*23000.d0*pi)/tau,	&	! velocity magnitude
    T0      = 300.d0,             &	! temperature
    H       = Rd * T0 / g,        &	! scale height
    RR      = 1.d0/2.d0,          &	! horizontal half width divided by 'a'
    ZZ      = 1000.d0,            &	! vertical half width
    z0      = 5000.d0,            &	! center point in z
    lambda0 = 5.d0*pi/6.d0,       &	! center point in longitudes
    lambda1 = 7.d0*pi/6.d0,       &	! center point in longitudes
    phi0    = 0.d0,               &	! center point in latitudes
    phi1    = 0.d0, &
    ztop    = 12000.d0
                            
  real(8) :: height                                                     ! The height of the model levels
  real(8) :: ptop                                                       ! model top in p
  real(8) :: sin_tmp, cos_tmp, sin_tmp2, cos_tmp2                       ! Calculate great circle distances
  real(8) :: d1, d2, r, r2, d3, d4                                      ! For tracer calculations
  real(8) :: s, bs, s_p                                                 ! Shape function, and parameter
  real(8) :: lonp                                                       ! Translational longitude, depends on time
  real(8) :: ud                                                         ! Divergent part of u
  real(8) :: x,y,zeta,tmp

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
	ptop = p0*exp(-ztop/H)

  !---------------------------------------------------------------------
  !    THE VELOCITIES ARE TIME DEPENDENT AND THEREFORE MUST BE UPDATED
  !    IN THE DYNAMICAL CORE
  !---------------------------------------------------------------------

	! shape function
	bs = 1.0d0
	s = 1.0 + exp((ptop-p0)/(bs*ptop)) - exp((p-p0)/(bs*ptop)) - exp((ptop-p)/(bs*ptop))
  s_p = (-exp((p-p0)/(bs*ptop)) + exp((ptop-p)/(bs*ptop)))/(bs*ptop)

	! translational longitude
	lonp = lon - 2.d0*pi*time/tau

	! zonal velocity
	ud = (omega0*a) * cos(lonp) * (cos(lat)**2.0) * cos(pi*time/tau) * s_p
		
	u = k0*sin(lonp)*sin(lonp)*sin(2.d0*lat)*cos(pi*time/tau) + u0*cos(lat) + ud

	! meridional velocity
	v = k0*sin(2.d0*lonp)*cos(lat)*cos(pi*time/tau)

	! vertical velocity - can be changed to vertical pressure velocity by
	! omega = -(g*p)/(Rd*T0)*w

  w = -((Rd*T0)/(g*p))*omega0*sin(lonp)*cos(lat)*cos(pi*time/tau)*s

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
	! tracer 1 - a C^inf tracer field for order of accuracy analysis

  x = cos(lat)*cos(lon)
  y = cos(lat)*sin(lon)
  zeta = sin(lat)
  q1 = 0.3d0*(1.1 + sin(0.25d0*pi*x)*sin(0.3d0*pi*y)*sin(0.25d0*pi*zeta)*sin(pi*(p-ptop)/(p0-ptop)))

	! tracer 2 - correlated with 1
	q2 = 0.9d0 - 0.8d0*q1**2

	! tracer 3 - slotted ellipse

	sin_tmp = sin(lat) * sin(phi0)
	cos_tmp = cos(lat) * cos(phi0)
	sin_tmp2 = sin(lat) * sin(phi1)
	cos_tmp2 = cos(lat) * cos(phi1)

  ! great circle distance without 'a'
	r  = ACOS (sin_tmp + cos_tmp*cos(lon-lambda0)) 
	r2 = ACOS (sin_tmp2 + cos_tmp2*cos(lon-lambda1)) 
	d1 = min( 1.d0, (r/RR)**2 + ((height-z0)/ZZ)**2 )
	d2 = min( 1.d0, (r2/RR)**2 + ((height-z0)/ZZ)**2 )

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

END SUBROUTINE test1_conv_advection_deformation

end module dcmip2012_test1_conv
