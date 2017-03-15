MODULE tropical_cyclone

!=======================================================================
!
!  Date:  July 29, 2015
!
!  Function for setting up idealized tropical cyclone initial conditions
!
!  SUBROUTINE tropical_cyclone_sample(
!    lon,lat,p,z,zcoords,u,v,t,thetav,phis,ps,rho,q)
!
!  Given a point specified by: 
!      lon    longitude (radians) 
!      lat    latitude (radians) 
!      p/z    pressure (Pa) / height (m)
!  zcoords    1 if z is specified, 0 if p is specified
!
!  the functions will return:
!        p    pressure if z is specified (Pa)
!        z    geopotential height if p is specified (m)
!        u    zonal wind (m s^-1)
!        v    meridional wind (m s^-1)
!        t    temperature (K)
!   thetav    virtual potential temperature (K)
!     phis    surface geopotential (m^2 s^-2)
!       ps    surface pressure (Pa)
!      rho    density (kj m^-3)
!        q    specific humidity (kg/kg)
!
!  Initial data are currently identical to:
!
!       Reed, K. A., and C. Jablonowski, 2011: An analytic
!       vortex initialization technique for idealized tropical
!       cyclone studies in AGCMs. Mon. Wea. Rev., 139, 689-710. 
!
!  Author: Kevin A. Reed
!          Stony Brook University
!          Email: kevin.a.reed@stonybrook.edu
!
!=======================================================================

  IMPLICIT NONE

!=======================================================================
!    Physical constants
!=======================================================================

  REAL(8), PARAMETER ::               &
       a     = 6371220.0d0,           & ! Reference Earth's Radius (m)
       Rd    = 287.0d0,               & ! Ideal gas const dry air (J kg^-1 K^1)
       g     = 9.80616d0,             & ! Gravity (m s^2)
       cp    = 1004.5d0,              & ! Specific heat capacity (J kg^-1 K^1)
       Lvap  = 2.5d6,                 & ! Latent heat of vaporization of water
       Rvap  = 461.5d0,               & ! Ideal gas constnat for water vapor
       Mvap  = 0.608d0,               & ! Ratio of molar mass of dry air/water
       pi    = 3.14159265358979d0,    & ! pi
       p0    = 100000.0d0,            & ! surface pressure (Pa)
       kappa = 2.d0/7.d0,             & ! Ratio of Rd to cp
       omega = 7.29212d-5,            & ! Reference rotation rate of the Earth (s^-1)
       deg2rad  = pi/180.d0             ! Conversion factor of degrees to radians

!=======================================================================
!    Test case parameters
!=======================================================================
  REAL(8), PARAMETER ::         &
       rp         = 282000.d0,  & ! Radius for calculation of PS
       dp         = 1115.d0,    & ! Delta P for calculation of PS
       zp         = 7000.d0,    & ! Height for calculation of P
       q0         = 0.021d0,    & ! q at surface from Jordan
       gamma      = 0.007d0,    & ! lapse rate
       Ts0        = 302.15d0,   & ! Surface temperature (SST)
       p00        = 101500.d0,  & ! global mean surface pressure
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

!=======================================================================
!    Evaluate the tropical cyclone initial conditions
!=======================================================================
  SUBROUTINE tropical_cyclone_test(lon,lat,p,z,zcoords,u,v,t,thetav,phis,ps,rho,q) &
    BIND(c, name = "tropical_cyclone_test")

    IMPLICIT NONE

    !------------------------------------------------
    !   Input / output parameters
    !------------------------------------------------

    REAL(8), INTENT(IN) ::     &
              lon,             &     ! Longitude (radians)
              lat                    ! Latitude (radians)

    REAL(8), INTENT(INOUT) ::  &
              p,               &     ! Pressure (Pa)
              z                      ! Height (m)

    INTEGER, INTENT(IN) :: zcoords     ! 1 if z coordinates are specified
                                     ! 0 if p coordinates are specified

    REAL(8), INTENT(OUT) ::    &
              u,               &     ! Zonal wind (m s^-1)
              v,               &     ! Meridional wind (m s^-1)
              t,               &     ! Temperature (K)
              thetav,          &     ! Virtual potential temperature (K)
              phis,            &     ! Surface Geopotential (m^2 s^-2)
              ps,              &     ! Surface Pressure (Pa)
              rho,             &     ! Density (kg m^-3)
              q                      ! Specific Humidity (kg/kg)

    !------------------------------------------------
    !   Local variables
    !------------------------------------------------
    real(8)  :: d1, d2, d, vfac, ufac, height, zhere, gr, f, zn

    integer  n

    !------------------------------------------------
    !   Define Great circle distance (gr) and
    !   Coriolis parameter (f)
    !------------------------------------------------
    f  = 2.d0*omega*sin(cen_lat*deg2rad)           ! Coriolis parameter
    gr = a*acos(sin(cen_lat*deg2rad)*sin(lat) + &  ! Great circle radius
         (cos(cen_lat*deg2rad)*cos(lat)*cos(lon-cen_lon*deg2rad)))

    !------------------------------------------------
    !   Initialize PS (surface pressure)
    !------------------------------------------------
    ps = p00-dp*exp(-(gr/rp)**exppr) 

    !------------------------------------------------
    !   Initialize altitude (z) if pressure provided
    !   or pressure if altitude (z) is provided
    !------------------------------------------------
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

    !------------------------------------------------
    !   Initialize U and V (wind components)
    !------------------------------------------------
    d1 = sin(cen_lat*deg2rad)*cos(lat) - &
         cos(cen_lat*deg2rad)*sin(lat)*cos(lon-cen_lon*deg2rad)
    d2 = cos(cen_lat*deg2rad)*sin(lon-cen_lon*deg2rad)
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

    !------------------------------------------------
    !   Initialize water vapor mixing ratio (q)
    !------------------------------------------------
    if (height > ztrop) then
        q = qtrop
    else
        q = q0*exp(-height/zq1)*exp(-(height/zq2)**exppz)
    end if

    !------------------------------------------------
    !   Initialize temperature (T)
    !------------------------------------------------
    if (height > ztrop) then
        t = Ttrop
    else
        t = (T0-gamma*height)/(1.d0+constTv*q)/(1.d0+exppz*Rd*(T0-gamma*height)*height &
            /(g*zp**exppz*(1.d0-p00/dp*exp((gr/rp)**exppr)*exp((height/zp)**exppz))))
    end if

    !-----------------------------------------------------
    !   Initialize virtual potential temperature (thetav)
    !-----------------------------------------------------
    thetav = t * (1.d0+constTv*q) * (p0/p)**(Rd/cp)

    !-----------------------------------------------------
    !   Initialize surface geopotential (PHIS)
    !-----------------------------------------------------
    phis = 0.d0  ! constant

    !-----------------------------------------------------
    !   Initialize density (rho)
    !-----------------------------------------------------
    rho = p/(Rd*t*(1.d0+constTv*q))

  END SUBROUTINE tropical_cyclone_test

!-----------------------------------------------------------------------
!    First function for fixed point iterations
!-----------------------------------------------------------------------
  REAL(8) FUNCTION fpiF(phere, gr, zhere)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: phere, gr, zhere

      fpiF = phere-(p00-dp*exp(-(gr/rp)**exppr)*exp(-(zhere/zp)**exppz)) &
             *((T0-gamma*zhere)/T0)**(g/(Rd*gamma))

  END FUNCTION fpiF

!-----------------------------------------------------------------------
!    Second function for fixed point iterations
!-----------------------------------------------------------------------
  REAL(8) FUNCTION fpidFdz(gr, zhere) 
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: gr, zhere

      fpidFdz =-exppz*zhere*dp*exp(-(gr/rp)**exppr)*exp(-(zhere/zp)**exppz)/(zp*zp)*((T0-gamma*zhere)/T0)**(g/(Rd*gamma)) &
               +g/(Rd*T0)*(p00-dp*exp(-(gr/rp)**exppr)*exp(-(zhere/zp)**exppz))*((T0-gamma*zhere)/T0)**(g/(Rd*gamma)-1.d0)

  END FUNCTION fpidFdz

END MODULE tropical_cyclone 
