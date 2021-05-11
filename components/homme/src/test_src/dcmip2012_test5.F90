MODULE dcmip2012_test5

  !=======================================================================
  !
  !  Function for setting up idealized tropical cyclone initial conditions
  !
  !  Given a point specified by: 
  !  	longitude (radians) 
  ! 	latitude (radians) 
  ! 	pressure/height
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
  !  Initial data are currently identical to:
  !
  !                 Reed, K. A., and C. Jablonowski, 2011: An analytic
  !                 vortex initialization technique for idealized tropical
  !                 cyclone studies in AGCMs. Mon. Wea. Rev., 139, 689-710. 
  !
  !  Author: Kevin A. Reed (University of Michigan, kareed@umich.edu)
  !          version 1
  !          5/21/2012
  !
  !=======================================================================

  IMPLICIT NONE

!-----------------------------------------------------------------------
!     Physical Parameters (or use a module?)
!-----------------------------------------------------------------------

	real(8), parameter ::	a	= 6371220.d0,	&	! Earth's Radius (m)
				Rd 	= 287.0d0,	&	! Ideal gas const dry air (J kg^-1 K^1)
				g	= 9.80616d0,	&	! Gravity (m s^2)
                                omega   = 7.292115d-5,  &       ! angular velocity 1/s
                                pi      = 4.d0*atan(1.d0)       ! pi

!-----------------------------------------------------------------------
!     Additional constants
!-----------------------------------------------------------------------
	real(8), parameter ::	convert = 180.d0/pi             ! conversion factor: radians to degrees

!-----------------------------------------------------------------------
!     Tropical Cyclone Test Case Tuning Parameters
!-----------------------------------------------------------------------
      real(8), parameter :: rp         = 282000.d0,  & ! Radius for calculation of PS
                            dp         = 1115.d0,    & ! Delta P for calculation of PS
                            zp         = 7000.d0,    & ! Height for calculation of P
                            q0         = 0.021d0,    & ! q at surface from Jordan
                            gamma      = 0.007d0,    & ! lapse rate
                            Ts0        = 302.15d0,   & ! Surface temperature (SST)
                            p00        = 101500.d0,  & ! global mean surface pressure
                            p0         = 100000.d0,  & ! p for model level calculation
                            cen_lat    = 10.d0,      & ! Center latitude of initial vortex
                            cen_lon    = 180.d0,     & ! Center longitufe of initial vortex
                            zq1        = 3000.d0,    & ! Height 1 for q calculation
                            zq2        = 8000.d0,    & ! Height 2 for q calculation
                            exppr      = 1.5d0,      & ! Exponent for r dependence of p
                            exppz      = 2.d0,       & ! Exponent for z dependence of p
                            ztrop      = 15000.d0,   & ! Tropopause Height
                            qtrop      = 1.d-11,     & ! Tropopause specific humidity
                            rfpi       = 1000000.d0, & ! Radius within which to use fixed-point iter.
                            constTv    = 0.608d0,    & ! Constant for Virtual Temp Conversion
                            deltaz     = 2.d-13,     & ! Small number to ensure convergence in FPI
                            epsilon    = 1.d-25,     & ! Small number to aviod dividing by zero in wind calc
                            exponent = Rd*gamma/g,   & ! exponent
                            T0    = Ts0*(1.d0+constTv*q0),             & ! Surface temp
                            Ttrop = T0 - gamma*ztrop,                  & ! Tropopause temp
                            ptrop = p00*(Ttrop/T0)**(1.d0/exponent)      ! Tropopause pressure

CONTAINS 


!==========================================================================================
! TEST CASE 5 - Tropical Cyclone 
!==========================================================================================

SUBROUTINE test5_tropical_cyclone (lon,lat,p,z,zcoords,u,v,w,t,phis,ps,rho,q)

IMPLICIT NONE

!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------

	real(8), intent(in)  :: lon, &		! Longitude (radians)
				lat, &		! Latitude (radians)
				z		! Height (m)

	real(8), intent(inout) :: p		! Pressure  (Pa)

	integer,  intent(in) :: zcoords 	! 0 or 1 see below
	real(8), intent(out) :: u, & 		! Zonal wind (m s^-1)
				v, &		! Meridional wind (m s^-1)
				w, &		! Vertical velocity (m s^-1)
				t, & 		! Temperature (K)
				phis, & 	! Surface Geopotential (m^2 s^-2)
				ps, & 		! Surface Pressure (Pa)
				rho, & 		! density (kg m^-3)
				q  		! Specific Humidity (kg/kg)

	! if zcoords = 1, then we use z and output p
	! if zcoords = 0, then we use p

!-----------------------------------------------------------------------
!    Additional parameters
!-----------------------------------------------------------------------
    real(8)  :: d1, d2, d, vfac, ufac, height, zhere, gr, f, zn

    integer  n

!-----------------------------------------------------------------------
!    Define Great circle distance (gr) and Coriolis parameter (f)
!-----------------------------------------------------------------------

    f  = 2.d0*omega*sin(cen_lat/convert)           ! Coriolis parameter
    gr = a*acos(sin(cen_lat/convert)*sin(lat) + &  ! Great circle distance
         (cos(cen_lat/convert)*cos(lat)*cos(lon-cen_lon/convert)))

!-----------------------------------------------------------------------
!    initialize PS (surface pressure of moist air)
!-----------------------------------------------------------------------
    ps = p00-dp*exp(-(gr/rp)**exppr) 

!-----------------------------------------------------------------------
!    initialize height field if provided pressure or pressure if provided z
!-----------------------------------------------------------------------

    if (zcoords .eq. 1) then

       height = z
 
       if (height > ztrop) then
          p = ptrop*exp(-(g*(height-ztrop))/(Rd*Ttrop))
       else
          p = (p00-dp*exp(-(gr/rp)**exppr)*exp(-(height/zp)**exppz)) &
              * ((T0-gamma*height)/T0)**(1/exponent)
       end if
 
    else

       height = (T0/gamma)*(1.d0-(p/ps)**exponent)

   ! If inside a certain distance of the center of the storm
   ! perform a Fixed-point iteration to calculate the height
   ! more accurately

       if (gr < rfpi ) then
          zhere = height 
          n = 1
          20 continue
          n = n+1
          zn = zhere - fpiF(p,gr,zhere)/fpidFdz(gr,zhere)
          if (n.gt.20) then
              PRINT *,'FPI did not converge after 20 interations in q & T!!!'
          else if ( abs(zn-zhere)/abs(zn) > deltaz) then
              zhere = zn
              goto 20
          end if
          height = zn
       end if

    end if

!-----------------------------------------------------------------------
!    initialize U and V (wind components)
!-----------------------------------------------------------------------

    d1 = sin(cen_lat/convert)*cos(lat) - &
         cos(cen_lat/convert)*sin(lat)*cos(lon-cen_lon/convert)
    d2 = cos(cen_lat/convert)*sin(lon-cen_lon/convert)
    d  = max(epsilon, sqrt(d1**2.d0 + d2**2.d0))
    ufac = d1/d
    vfac = d2/d
    if (height > ztrop) then
        u = 0.d0
        v = 0.d0
    else
        v = vfac*(-f*gr/2.d0+sqrt((f*gr/2.d0)**(2.d0) &
            - exppr*(gr/rp)**exppr*Rd*(T0-gamma*height) &
            /(exppz*height*Rd*(T0-gamma*height)/(g*zp**exppz) &
            +(1.d0-p00/dp*exp((gr/rp)**exppr)*exp((height/zp)**exppz)))))
        u = ufac*(-f*gr/2.d0+sqrt((f*gr/2.d0)**(2.d0) &
            - exppr*(gr/rp)**exppr*Rd*(T0-gamma*height) &
            /(exppz*height*Rd*(T0-gamma*height)/(g*zp**exppz) &
            +(1.d0-p00/dp*exp((gr/rp)**exppr)*exp((height/zp)**exppz)))))
    end if

!-----------------------------------------------------------------------
!    set the vertical velocity to zero (only required for non-hydrostatic models)
!-----------------------------------------------------------------------

    w = 0.d0

!-----------------------------------------------------------------------
!    tracer q (specific humidity)
!-----------------------------------------------------------------------

    if (height > ztrop) then
        q = qtrop
    else
        q = q0*exp(-height/zq1)*exp(-(height/zq2)**exppz)
    end if

!-----------------------------------------------------------------------
!    initialize T (temperature)
!-----------------------------------------------------------------------

    if (height > ztrop) then
        t = Ttrop
    else
        t = (T0-gamma*height)/(1.d0+constTv*q)/(1.d0+exppz*Rd*(T0-gamma*height)*height &
            /(g*zp**exppz*(1.d0-p00/dp*exp((gr/rp)**exppr)*exp((height/zp)**exppz))))
    end if

!-----------------------------------------------------------------------
!    initialize PHIS (surface geopotential)
!-----------------------------------------------------------------------

    phis = 0.d0  ! constant

!-----------------------------------------------------------------------
!    initialize RHO (density of moist air)
!-----------------------------------------------------------------------

    rho = p/(Rd*t*(1.d0+constTv*q))

  END SUBROUTINE test5_tropical_cyclone


!********************************************************************
! Function for fixed point iterations
!********************************************************************
  REAL*8 FUNCTION fpiF(phere,gr,zhere)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: phere, gr, zhere

      fpiF = phere-(p00-dp*exp(-(gr/rp)**exppr)*exp(-(zhere/zp)**exppz)) &
             *((T0-gamma*zhere)/T0)**(g/(Rd*gamma))

  END FUNCTION fpiF

!********************************************************************
! Function for fixed point iterations 
!********************************************************************
  REAL*8 FUNCTION fpidFdz(gr,zhere) 
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: gr, zhere

      fpidFdz =-exppz*zhere*dp*exp(-(gr/rp)**exppr)*exp(-(zhere/zp)**exppz)/(zp*zp)*((T0-gamma*zhere)/T0)**(g/(Rd*gamma)) &
               +g/(Rd*T0)*(p00-dp*exp(-(gr/rp)**exppr)*exp(-(zhere/zp)**exppz))*((T0-gamma*zhere)/T0)**(g/(Rd*gamma)-1.d0)

  END FUNCTION fpidFdz


END MODULE 
