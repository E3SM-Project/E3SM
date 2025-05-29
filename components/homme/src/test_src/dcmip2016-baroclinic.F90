MODULE baroclinic_wave

!=======================================================================
!
!  Date:  July 29, 2015
!
!  Functions for setting up idealized initial conditions for the
!  Ullrich, Melvin, Staniforth and Jablonowski baroclinic instability.
!
!  SUBROUTINE baroclinic_wave_sample(
!    deep,moist,pertt,X,lon,lat,p,z,zcoords,u,v,w,t,phis,ps,rho,q)
!
!  Options:
!     deep    deep atmosphere (1 = yes or 0 = no)
!    moist    include moisture (1 = yes or 0 = no)
!    pertt    type of perturbation (0 = exponential, 1 = stream function)
!        X    Earth scaling factor
!
!  Given a point specified by: 
!      lon    longitude (radians) 
!      lat    latitude (radians) 
!      p/z    pressure (Pa) / height (m)
!  zcoords    1 if z is specified, 0 if p is specified
!
!  the functions will return:
!        p    pressure if z is specified and zcoords = 1 (Pa)
!        u    zonal wind (m s^-1)
!        v    meridional wind (m s^-1)
!        t    temperature (K)
!   thetav    virtual potential temperature (K)
!     phis    surface geopotential (m^2 s^-2)
!       ps    surface pressure (Pa)
!      rho    density (kj m^-3)
!        q    water vapor mixing ratio (kg/kg)
!
!
!  Author: Paul Ullrich
!          University of California, Davis
!          Email: paullrich@ucdavis.edu
!
!=======================================================================

use physical_constants, only: g0=>g,kappa0=>kappa,Rgas,Cp0=>Cp,Rwater_vapor,rearth0,omega0, dd_pi

  IMPLICIT NONE


!=======================================================================
!    Physical constants
!=======================================================================

 ! REAL(8), PARAMETER ::               &
 !      a     = 6371220.0d0,           & ! Reference Earth's Radius (m)
 !      Rd    = 287.0d0,               & ! Ideal gas const dry air (J kg^-1 K^1)
 !      g     = 9.80616d0,             & ! Gravity (m s^2)
 !      cp    = 1004.5d0,              & ! Specific heat capacity (J kg^-1 K^1)
 !      Lvap  = 2.5d6,                 & ! Latent heat of vaporization of water
 !      Rvap  = 461.5d0,               & ! Ideal gas constnat for water vapor
 !      Mvap  = 0.608d0,               & ! Ratio of molar mass of dry air/water
 !      pi    = 3.14159265358979d0,    & ! pi
 !      p0    = 100000.0d0,            & ! surface pressure (Pa)
 !      kappa = 2.d0/7.d0,             & ! Ratio of Rd to cp
 !      omega = 7.29212d-5,            & ! Reference rotation rate of the Earth (s^-1)
 !      deg2rad  = pi/180.d0             ! Conversion factor of degrees to radians

 REAL(8), PARAMETER ::                &
       a     = rearth0,               & ! Reference Earth's Radius (m)
       Rd    = Rgas,                  & ! Ideal gas const dry air (J kg^-1 K^1)
       g     = g0,                    & ! Gravity (m s^2)
       cp    = Cp0,                   & ! Specific heat capacity (J kg^-1 K^1)
       Lvap  = 2.5d6,                 & ! Latent heat of vaporization of water
       Rvap  = Rwater_vapor,          & ! Ideal gas constnat for water vapor
       Mvap  = Rvap/Rd-1.0d0,         & ! Ratio of molar mass of dry air/water
       pi    = dd_pi,                 & ! pi
       p0    = 100000.0d0,            & ! surface pressure (Pa)
       kappa = kappa0,                & ! Ratio of Rd to cp
       omega = omega0,                & ! Reference rotation rate of the Earth (s^-1)
       deg2rad  = pi/180.d0             ! Conversion factor of degrees to radians

!=======================================================================
!    Test case parameters
!=======================================================================
  REAL(8), PARAMETER ::               &
       T0E        = 310.d0     ,      & ! temperature at equatorial surface (K)
       T0P        = 240.d0     ,      & ! temperature at polar surface (K)
       B          = 2.d0       ,      & ! jet half-width parameter
       K          = 3.d0       ,      & ! jet width parameter
       lapse      = 0.005d0             ! lapse rate parameter

  REAL(8), PARAMETER ::               &
       pertu0     = 0.5d0      ,      & ! SF Perturbation wind velocity (m/s)
       pertr      = 1.d0/6.d0  ,      & ! SF Perturbation radius (Earth radii)
       pertup     = 1.0d0      ,      & ! Exp. perturbation wind velocity (m/s)
       pertexpr   = 0.1d0      ,      & ! Exp. perturbation radius (Earth radii)
       pertlon    = pi/9.d0    ,      & ! Perturbation longitude
       pertlat    = 2.d0*pi/9.d0,     & ! Perturbation latitude
       pertz      = 15000.d0   ,      & ! Perturbation height cap
       dxepsilon  = 1.d-5               ! Small value for numerical derivatives
 
  REAL(8), PARAMETER ::               &
       moistqlat  = 2.d0*pi/9.d0,     & ! Humidity latitudinal width
       moistqp    = 34000.d0,         & ! Humidity vertical pressure width
       moisttr    = 0.1d0,            & ! Vertical cut-off pressure for humidity
       moistqs    = 1.d-12,           & ! Humidity above cut-off
       moistq0    = 0.018d0,          & ! Maximum specific humidity
       moistqr    = 0.9d0,            & ! Maximum saturation ratio
       moisteps   = 0.622d0,          & ! Ratio of gas constants
       moistT0    = 273.16d0,         & ! Reference temperature (K)
       moistE0Ast = 610.78d0            ! Saturation vapor pressure at T0 (Pa) 

  real(8), parameter, private ::  &
       mountain_height = 2000,                &  ! Peak height of mountain (m)
       ! HOMME shift mountains by pi*4/5 to put them over CONUS
       mountain_lon_0 = (7.d0/9.d0 + 4d0/5d0) * pi,&  ! Longitudinal center of mountain 1(rad) 
       mountain_lon_1 = (0.4d0     + 4d0/5d0) * pi,&  ! Longitudinal center of mountain 2 (rad)
       mountain_lon_width =  7.0d0 * deg2rad, &  ! Distance between 10th percentile of mountain in longitude
       mountain_lat_0 = pi / 4.0d0,           &  ! Position of center of mountain (rad) 
       mountain_lat_width = 40.0d0 * deg2rad,   &  ! Distance between 10th percentile of mountain in latitude
       mountain_lat_scale = mountain_lat_width / (2.0d0*(-log(0.1d0))**(1.0d0/6.0d0) ),  & 
       mountain_lon_scale = mountain_lon_width / (2.0d0*(-log(0.1d0))**(1.0d0/2.0d0) )
 

CONTAINS

!=======================================================================
!    Generate the baroclinic instability initial conditions
!=======================================================================
  SUBROUTINE baroclinic_wave_test(deep,moist,pertt,X,lon,lat,p,z,zcoords,u,v,t,thetav,phis,ps,rho,q) &
    BIND(c, name = "baroclinic_wave_test")
 
    IMPLICIT NONE

!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------
    INTEGER, INTENT(IN)  :: &
                deep,       & ! Deep (1) or Shallow (0) test case
                moist,      & ! Moist (1) or Dry (0) test case
                pertt         ! Perturbation type

    REAL(8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                X             ! Earth scaling parameter

    REAL(8), INTENT(INOUT) :: &
                p,            & ! Pressure (Pa)
                z               ! Altitude (m)

    INTEGER, INTENT(IN) :: zcoords     ! 1 if z coordinates are specified
                                       ! 0 if p coordinates are specified

    REAL(8), INTENT(OUT) :: &
                u,          & ! Zonal wind (m s^-1)
                v,          & ! Meridional wind (m s^-1)
                t,          & ! Temperature (K)
                thetav,     & ! Virtual potential temperature (K)
                phis,       & ! Surface Geopotential (m^2 s^-2)
                ps,         & ! Surface Pressure (Pa)
                rho,        & ! density (kg m^-3)
                q             ! water vapor mixing ratio (kg/kg)

    !------------------------------------------------
    !   Local variables
    !------------------------------------------------
    REAL(8) :: aref, omegaref
    REAL(8) :: T0, constH, constC, scaledZ, inttau2, rratio
    REAL(8) :: inttermU, bigU, rcoslat, omegarcoslat
    REAL(8) :: eta, qratio, qnum, qden

    !------------------------------------------------
    !   Pressure and temperature
    !------------------------------------------------
    if (zcoords .eq. 1) then
      CALL evaluate_pressure_temperature(deep, X, lon, lat, z, p, t)
    else
      CALL evaluate_z_temperature(deep, X, lon, lat, p, z, t)
    end if

    !------------------------------------------------
    !   Compute test case constants
    !------------------------------------------------
    aref = a / X
    omegaref = omega * X

    T0 = 0.5d0 * (T0E + T0P)

    constH = Rd * T0 / g

    constC = 0.5d0 * (K + 2.d0) * (T0E - T0P) / (T0E * T0P)

    scaledZ = z / (B * constH)

    inttau2 = constC * z * exp(- scaledZ**2)

    ! radius ratio
    if (deep .eq. 0) then
      rratio = 1.d0
    else
      rratio = (z + aref) / aref;
    end if

    !-----------------------------------------------------
    !   Initialize surface pressure
    !-----------------------------------------------------
    ps = p0

    !-----------------------------------------------------
    !   Initialize velocity field
    !-----------------------------------------------------
    inttermU = (rratio * cos(lat))**(K - 1.d0) - (rratio * cos(lat))**(K + 1.d0)
    bigU = g / aref * K * inttau2 * inttermU * t
    if (deep .eq. 0) then
      rcoslat = aref * cos(lat)
    else
      rcoslat = (z + aref) * cos(lat)
    end if

    omegarcoslat = omegaref * rcoslat
    
    u = - omegarcoslat + sqrt(omegarcoslat**2 + rcoslat * bigU)
    v = 0.d0

    !-----------------------------------------------------
    !   Add perturbation to the velocity field
    !-----------------------------------------------------

    ! Exponential type
    if (pertt .eq. 0) then
      u = u + evaluate_exponential(lon, lat, z)

    ! Stream function type
    elseif (pertt .eq. 1) then
      u = u - 1.d0 / (2.d0 * dxepsilon) *                       &
          ( evaluate_streamfunction(lon, lat + dxepsilon, z)    &
          - evaluate_streamfunction(lon, lat - dxepsilon, z))

      v = v + 1.d0 / (2.d0 * dxepsilon * cos(lat)) *            &
          ( evaluate_streamfunction(lon + dxepsilon, lat, z)    &
          - evaluate_streamfunction(lon - dxepsilon, lat, z))
    end if

    !-----------------------------------------------------
    !   Initialize surface geopotential
    !-----------------------------------------------------
    phis = 0.d0

    !-----------------------------------------------------
    !   Initialize density
    !-----------------------------------------------------
    rho = p / (Rd * t)

    !-----------------------------------------------------
    !   Initialize specific humidity
    !-----------------------------------------------------
    if (moist .eq. 1) then
      eta = p/p0

      if (eta .gt. moisttr) then
        q = moistq0 * exp(- (lat/moistqlat)**4)          &
                    * exp(- ((eta-1.d0)*p0/moistqp)**2)
      else
        q = moistqs
      end if

      ! Convert virtual temperature to temperature
      t = t / (1.d0 + Mvap * q)

    else
      q = 0.d0
    end if

    !-----------------------------------------------------
    !   Initialize virtual potential temperature
    !-----------------------------------------------------
    thetav = t * (1.d0 + Mvap * q) * (p0 / p)**(Rd / cp)

  END SUBROUTINE baroclinic_wave_test
!=======================================================================
!    Generate the baroclinic instability initial conditions
!=======================================================================
  SUBROUTINE baroclinic_topo_test(deep,moist,pertt,X,lon,lat,p,z,zcoords,u,v,t,thetav,phis,ps,rho,q) &
    BIND(c, name = "baroclinic_topo_test")
    use parallel_mod,         only: abortmp  
    IMPLICIT NONE

!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------
    INTEGER, INTENT(IN)  :: &
                deep,       & ! Deep (1) or Shallow (0) test case
                moist,      & ! Moist (1) or Dry (0) test case
                pertt         ! Perturbation type

    REAL(8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                X             ! Earth scaling parameter

    REAL(8), INTENT(INOUT) :: &
                p,            & ! Pressure (Pa)
                z               ! Altitude (m)

    INTEGER, INTENT(IN) :: zcoords     ! 1 if z coordinates are specified
                                       ! 0 if p coordinates are specified

    REAL(8), INTENT(OUT) :: &
                u,          & ! Zonal wind (m s^-1)
                v,          & ! Meridional wind (m s^-1)
                t,          & ! Temperature (K)
                thetav,     & ! Virtual potential temperature (K)
                phis,       & ! Surface Geopotential (m^2 s^-2)
                ps,         & ! Surface Pressure (Pa)
                rho,        & ! density (kg m^-3)
                q             ! water vapor mixing ratio (kg/kg)

    !------------------------------------------------
    !   Local variables
    !------------------------------------------------
    REAL(8) :: aref, omegaref
    REAL(8) :: T0, constH, constC, scaledZ, inttau2, rratio
    REAL(8) :: inttermU, bigU, rcoslat, omegarcoslat
    REAL(8) :: eta, qratio, qnum, qden

    !------------------------------------------------
    !   Pressure and temperature
    !------------------------------------------------
    if (zcoords .eq. 1) then
        call abortmp('ERROR: baroclinic topo case is not implemented for z coordinates')
    else
       call evaluate_phis_ps(deep,X,lon,lat,ps,phis)
       CALL evaluate_z_temperature(deep, X, lon, lat, p, z, t)
    end if

    !------------------------------------------------
    !   Compute test case constants
    !------------------------------------------------
    aref = a / X
    omegaref = omega * X

    T0 = 0.5d0 * (T0E + T0P)

    constH = Rd * T0 / g

    constC = 0.5d0 * (K + 2.d0) * (T0E - T0P) / (T0E * T0P)

    scaledZ = z / (B * constH)

    inttau2 = constC * z * exp(- scaledZ**2)

    ! radius ratio
    if (deep .eq. 0) then
      rratio = 1.d0
    else
      rratio = (z + aref) / aref;
    end if

    !-----------------------------------------------------
    !   Initialize surface pressure
    !-----------------------------------------------------
    !ps = p0

    !-----------------------------------------------------
    !   Initialize velocity field
    !-----------------------------------------------------
    inttermU = (rratio * cos(lat))**(K - 1.d0) - (rratio * cos(lat))**(K + 1.d0)
    bigU = g / aref * K * inttau2 * inttermU * t
    if (deep .eq. 0) then
      rcoslat = aref * cos(lat)
    else
      rcoslat = (z + aref) * cos(lat)
    end if

    omegarcoslat = omegaref * rcoslat
    
    u = - omegarcoslat + sqrt(omegarcoslat**2 + rcoslat * bigU)
    v = 0.d0

    !-----------------------------------------------------
    !   Add perturbation to the velocity field
    !-----------------------------------------------------

    !! Exponential type
    !if (pertt .eq. 0) then
    !  u = u + evaluate_exponential(lon, lat, z)

    ! Stream function type
    !elseif (pertt .eq. 1) then
    !  u = u - 1.d0 / (2.d0 * dxepsilon) *                       &
    !      ( evaluate_streamfunction(lon, lat + dxepsilon, z)    &
    !      - evaluate_streamfunction(lon, lat - dxepsilon, z))

    !  v = v + 1.d0 / (2.d0 * dxepsilon * cos(lat)) *            &
    !      ( evaluate_streamfunction(lon + dxepsilon, lat, z)    &
    !      - evaluate_streamfunction(lon - dxepsilon, lat, z))
    !end if

    !-----------------------------------------------------
    !   Initialize surface geopotential
    !-----------------------------------------------------
    !phis = 0.d0

    !-----------------------------------------------------
    !   Initialize density
    !-----------------------------------------------------
    rho = p / (Rd * t)

    !-----------------------------------------------------
    !   Initialize specific humidity
    !-----------------------------------------------------
    if (moist .eq. 1) then
      eta = p/p0

      if (eta .gt. moisttr) then
        q = moistq0 * exp(- (lat/moistqlat)**4)          &
                    * exp(- ((eta-1.d0)*p0/moistqp)**2)
      else
        q = moistqs
      end if

      ! Convert virtual temperature to temperature
      t = t / (1.d0 + Mvap * q)

    else
      q = 0.d0
    end if

    !-----------------------------------------------------
    !   Initialize virtual potential temperature
    !-----------------------------------------------------
    thetav = t * (1.d0 + Mvap * q) * (p0 / p)**(Rd / cp)

  END SUBROUTINE baroclinic_topo_test



!-----------------------------------------------------------------------
!    Calculate pointwise pressure and temperature
!-----------------------------------------------------------------------
  SUBROUTINE evaluate_pressure_temperature(deep, X, lon, lat, z, p, t)

    INTEGER, INTENT(IN)  :: deep ! Deep (1) or Shallow (0) test case

    REAL(8), INTENT(IN)  :: &
                X,          & ! Earth scaling ratio
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                z             ! Altitude (m)

    REAL(8), INTENT(OUT) :: &
                p,          & ! Pressure (Pa)
                t             ! Temperature (K)

    REAL(8) :: aref, omegaref
    REAL(8) :: T0, constA, constB, constC, constH, scaledZ
    REAL(8) :: tau1, tau2, inttau1, inttau2
    REAL(8) :: rratio, inttermT

    !--------------------------------------------
    ! Constants
    !--------------------------------------------
    aref = a / X
    omegaref = omega * X

    T0 = 0.5d0 * (T0E + T0P)
    constA = 1.d0 / lapse
    constB = (T0 - T0P) / (T0 * T0P)
    constC = 0.5d0 * (K + 2.d0) * (T0E - T0P) / (T0E * T0P)
    constH = Rd * T0 / g

    scaledZ = z / (B * constH)

    !--------------------------------------------
    !    tau values
    !--------------------------------------------
    tau1 = constA * lapse / T0 * exp(lapse * z / T0) &
         + constB * (1.d0 - 2.d0 * scaledZ**2) * exp(- scaledZ**2)
    tau2 = constC * (1.d0 - 2.d0 * scaledZ**2) * exp(- scaledZ**2)

    inttau1 = constA * (exp(lapse * z / T0) - 1.d0) &
            + constB * z * exp(- scaledZ**2)
    inttau2 = constC * z * exp(- scaledZ**2)

    !--------------------------------------------
    !    radius ratio
    !--------------------------------------------
    if (deep .eq. 0) then
      rratio = 1.d0
    else
      rratio = (z + aref) / aref;
    end if

    !--------------------------------------------
    !    interior term on temperature expression
    !--------------------------------------------
    inttermT = (rratio * cos(lat))**K &
             - K / (K + 2.d0) * (rratio * cos(lat))**(K + 2.d0)

    !--------------------------------------------
    !    temperature
    !--------------------------------------------
    t = 1.d0 / (rratio**2 * (tau1 - tau2 * inttermT))

    !--------------------------------------------
    !    hydrostatic pressure
    !--------------------------------------------
    p = p0 * exp(- g / Rd * (inttau1 - inttau2 * inttermT))

  END SUBROUTINE evaluate_pressure_temperature

  SUBROUTINE evaluate_phis_ps(deep,X,lon,lat,ps,phis)
    integer, intent(in) :: deep ! ignored

    real(8), intent(in) :: &
                 X,        & ! earth reduction factor
                 lon,      & ! longitude in radians
                 lat         ! latitude in radians
    real(8), intent(inout), optional :: &
                 ps,       & ! surface pressure
                 phis          ! geometric surface height
    
    ! local variables

    real(8)            :: ztmp, ttmp,d0,d1
    if (present(ps) .or. present(phis)) then
       d0=mod(lon - mountain_lon_0,2*pi)
       d0=min(d0,2*pi-d0)
       d1=mod(lon - mountain_lon_1,2*pi)
       d1=min(d1,2*pi-d0)
       
       ztmp = exp(-1.0d0 * (((lat - mountain_lat_0)/ mountain_lat_scale)**6 +&
            (d0/ mountain_lon_scale)**2))
       ztmp = ztmp + exp(-1.0d0 * (((lat - mountain_lat_0)/ mountain_lat_scale)**6 +&
            (d1/mountain_lon_scale)**2))
        ztmp = ztmp *  mountain_height/X
    end if
    if (present(phis)) then
        phis = g * ztmp
    end if
    if (present(ps)) then
        call evaluate_pressure_temperature(deep, X, lon, lat, ztmp, ps, ttmp)
    end if
  END SUBROUTINE evaluate_phis_ps


!-----------------------------------------------------------------------
!    Calculate pointwise z and temperature given pressure
!-----------------------------------------------------------------------
  SUBROUTINE evaluate_z_temperature(deep, X, lon, lat, p, z, t)
    
    INTEGER, INTENT(IN)  :: deep ! Deep (1) or Shallow (0) test case

    REAL(8), INTENT(IN)  :: &
                X,          & ! Earth scaling ratio
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                p             ! Pressure (Pa)

    REAL(8), INTENT(OUT) :: &
                z,          & ! Altitude (m)
                t             ! Temperature (K)

    INTEGER :: ix

    REAL(8) :: z0, z1, z2
    REAL(8) :: p0, p1, p2

    z0 = 0.d0
    z1 = 10000.d0

    CALL evaluate_pressure_temperature(deep, X, lon, lat, z0, p0, t)
    CALL evaluate_pressure_temperature(deep, X, lon, lat, z1, p1, t)

    DO ix = 1, 1000
      z2 = z1 - (p1 - p) * (z1 - z0) / (p1 - p0)

      CALL evaluate_pressure_temperature(deep, X, lon, lat, z2, p2, t)

      IF (ABS((p2 - p)/p) .lt. 1.0d-13) THEN
        EXIT
      END IF

      z0 = z1
      p0 = p1

      z1 = z2
      p1 = p2
    END DO

    z = z2

    CALL evaluate_pressure_temperature(deep, X, lon, lat, z, p0, t)

  END SUBROUTINE evaluate_z_temperature

!-----------------------------------------------------------------------
!    Exponential perturbation function
!-----------------------------------------------------------------------
  REAL(8) FUNCTION evaluate_exponential(lon, lat, z)

    REAL(8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                z             ! Altitude (meters)

    REAL(8) :: greatcircler, perttaper

    ! Great circle distance
    greatcircler = 1.d0 / pertexpr &
      * acos(sin(pertlat) * sin(lat) + cos(pertlat) * cos(lat) * cos(lon - pertlon))

    ! Vertical tapering of stream function
    if (z < pertz) then
      perttaper = 1.d0 - 3.d0 * z**2 / pertz**2 + 2.d0 * z**3 / pertz**3
    else
      perttaper = 0.d0
    end if

    ! Zonal velocity perturbation
    if (greatcircler < 1.d0) then
      evaluate_exponential = pertup * perttaper * exp(- greatcircler**2)
    else
      evaluate_exponential = 0.d0
    end if

  END FUNCTION evaluate_exponential

!-----------------------------------------------------------------------
!    Stream function perturbation function
!-----------------------------------------------------------------------
  REAL(8) FUNCTION evaluate_streamfunction(lon, lat, z)

    REAL(8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                z             ! Altitude (meters)

    REAL(8) :: greatcircler, perttaper, cospert

    ! Great circle distance
    greatcircler = 1.d0 / pertr &
      * acos(sin(pertlat) * sin(lat) + cos(pertlat) * cos(lat) * cos(lon - pertlon))

    ! Vertical tapering of stream function
    if (z < pertz) then
      perttaper = 1.d0 - 3.d0 * z**2 / pertz**2 + 2.d0 * z**3 / pertz**3
    else
      perttaper = 0.d0
    end if

    ! Horizontal tapering of stream function
    if (greatcircler .lt. 1.d0) then
      cospert = cos(0.5d0 * pi * greatcircler)
    else
      cospert = 0.d0
    end if

    evaluate_streamfunction = &
        (- pertu0 * pertr * perttaper * cospert**4)

  END FUNCTION evaluate_streamfunction

END MODULE baroclinic_wave
