MODULE dcmip2012_test4

  !=======================================================================
  !
  !  Functions for setting up idealized and perturbed initial conditions
  !  for the Jablonowski-Williamson baroclinic instability with invariant
  !  tracers (Jablonowski and Williamson, QJ 2006).
  !
  !  Options:
  !   moist    1 for the moist baroclinic instability 
  !            0 for the dry baroclinic instability test
  !       X    scale factor of the Earth
  !
  !  Given a point specified by: 
  !     lon    longitude (radians) 
  !     lat    latitude (radians) 
  !     p/z    pressure (Pa)/height (m)
  ! zcoords    1 if z is specified, 0 if p is specified
  !
  !  the functions will return:
  !       p    pressure if z is specified and zcoords = 1 (Pa)
  !       u    zonal wind (m s^-1)
  !       v    meridional wind (m s^-1)
  !       w    vertical velocity (m s^-1)
  !       t    temperature (K)
  !    phis    surface geopotential (m^2 s^-2)
  !      ps    surface pressure (Pa)
  !     rho    density (kg m^-3)
  !       q    specific humidity (kg/kg)
  !      q1    potential temperature tracer (K)
  !      q2    absolute value of Ertel's potential vorticity (EPV) tracer,
  !            in PVU units for X=1: PVU = 10^-6 K m^2 kg^-1 s^-1,
  !
  !  Authors: Paul Ullrich, Christiane Jablonowski, James Kent
  !           (University of Michigan, dcmip@ucar.edu)
  !           version 3
  !           corrected on July/20/2012
  !
  !           Change Log:
  !           version 2 (v2), June/8/2012
  !             (1) correction of a typo (removal of 1/a) in the dthetadphi equation
  !             (2) initialization of q2 with the absolute value of EPV
  !           version 3 (v3), July/20/2012
  !            - newly added if construct prevents division by zero in relative vorticity (zeta) 
  !              calculation in case the perturbation center point, its antipode, or the north or 
  !              south poles are part of the computational grid
  !=======================================================================

  IMPLICIT NONE

!=======================================================================
! use physical constants
!=======================================================================

  !REAL(8), PARAMETER ::                               &
  !     a     = 6371220.0d0,                           & ! Reference Earth's Radius (m)
  !     Rd    = 287.0d0,                               & ! Ideal gas const dry air (J kg^-1 K^1)
  !     g     = 9.80616d0,                             & ! Gravity (m s^2)
  !     cp    = 1004.5d0,                              & ! Specific heat capacity (J kg^-1 K^1)
  !     pi    = 4.d0*atan(1.d0),                       & ! pi
  !     kappa = 2.d0/7.d0,                             & ! Ratio of Rd to cp
  !     omega = 7.29212d-5                              ! Reference rotation rate of the Earth (s^-1)

  ! use physical constants consistent with HOMME
  REAL(8), PARAMETER ::                               &
       a     = 6.376D6 ,                              & ! Reference Earth's Radius (m)
       Rd    = 287.04D0,                              & ! Ideal gas const dry air (J kg^-1 K^1)
       g     = 9.80616D0,                             & ! Gravity (m s^2)
       cp    = 1005.0D0,                              & ! Specific heat capacity (J kg^-1 K^1)
       pi    = 4.d0*atan(1.d0),                       & ! pi
       kappa = Rd/Cp,                                 & ! Ratio of Rd to cp
       omega = 7.29212d-5                               ! Reference rotation rate of the Earth (s^-1)


!=======================================================================
! additional constants
!=======================================================================

  REAL(8), PARAMETER ::                               &
       p0       = 100000.0d0,                         & ! surface pressure (Pa)
       deg2rad  = pi/180.d0                             ! Conversion factor of degrees to radians

!-----------------------------------------------------------------------
! steady-state and baroclinic wave tuning parameter
!----------------------------------------------------------------------- 
  REAL(8), PARAMETER ::                            &
       eta_tropo  = 0.2d0     ,                    & ! tropopause level
       u0         = 35.d0     ,                    & ! 35 m/s
       T0         = 288.d0    ,                    & ! horizontal mean T at surface
       eta0       = 0.252d0   ,                    & ! center of jets (hybrid)
       !
       radius                 = 10.d0,             & ! reciprocal radius of the perturbation without 'a'
       perturbation_amplitude =  1.d0,             & ! amplitude of u perturbation 1 m/s
       perturbation_longitude = 20.d0,             & ! longitudinal position, 20E
       perturbation_latitude  = 40.d0,             & ! latitudinal position, 40N
       eta_sfc                = 1.d0,              & ! hybrid value at surface
       delta_T                = 480000.d0,         & ! in K, for T mean calculation
       gamma                  = 0.005d0,           & ! lapse rate
       !
       perturbation_latitude_tracer = 55.d0,       &        
       a_omega                = a*omega,           &
       exponent               = Rd*gamma/g,        &
       q0                     = 0.021d0,           & ! maximum specific humidity in kg/kg
       lat_hw                 = 2.d0*pi/9.d0,      & ! halfwidth for q: 40 degrees in radians
       p_hw                   = 34000.d0             ! halfwidth for q: 34000 Pa   

CONTAINS

  SUBROUTINE test4_baroclinic_wave (moist,X,lon,lat,p,z,zcoords,u,v,w,t,phis,ps,rho,q,q1,q2)
 
    IMPLICIT NONE

!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------
    INTEGER, INTENT(IN)  :: moist      ! Moist (1) or non-moist (0) test case

    REAL(8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                X             ! Scale factor, not used in this version since unscaled EPV is selected

    REAL(8), INTENT(INOUT) :: p,z        ! Pressure  (Pa), and height(m)

    INTEGER, INTENT(IN) :: zcoords     ! 0 or 1 see below

    REAL(8), INTENT(OUT) :: &
                u,          & ! Zonal wind (m s^-1)
                v,          & ! Meridional wind (m s^-1)
                w,          & ! Vertical Velocity (m s^-1)
                t,          & ! Temperature (K)
                phis,       & ! Surface Geopotential (m^2 s^-2)
                ps,         & ! Surface Pressure (Pa)
                rho,        & ! density (kg m^-3)
                q,          & ! Specific Humidity (kg/kg)
                q1,         & ! Tracer q1 - Potential temperature (kg/kg)
                q2            ! Tracer q2 - Ertel's potential vorticity (kg/kg)

    REAL(8) :: eta

!-----------------------------------------------------------------------
!    initialize PS (surface pressure)
!-----------------------------------------------------------------------
     ps = p0  ! constant

!-----------------------------------------------------------------------
!    calculate eta
!-----------------------------------------------------------------------
    if (zcoords .eq. 0) then
      eta = p / ps
      z = 0 
      z = geopotential(lon,lat,eta)/g
    else
      eta = eta_from_z(lon,lat,z)
      p = eta * ps
    end if

!-----------------------------------------------------------------------
!    initialize wind field with perturbation
!-----------------------------------------------------------------------
     u = u_wind(lon,lat,eta,.TRUE.)
     v = v_wind(lon,lat,eta,.TRUE.)
     w = 0.d0

!-----------------------------------------------------------------------
!    initialize temperature field (virtual temperature for moist conditions)
!-----------------------------------------------------------------------
    t = temperature(lon,lat,eta)

!-----------------------------------------------------------------------
!    initialize surface geopotential
!-----------------------------------------------------------------------
    phis = surface_geopotential(lon,lat)  

!-----------------------------------------------------------------------
!    initialize density from ideal gas law, t is the virtual temperature in the moist case
!-----------------------------------------------------------------------
    rho = p/(t*Rd)

!-----------------------------------------------------------------------
!      tracer q (specific humidity), moist or dry, and dry temperature
!-----------------------------------------------------------------------
    if (moist == 1) then
       q = q0 * exp(-(lat/lat_hw)**4) * exp(-((eta-1.d0)*p0/p_hw)**2)  ! moist
       t = t/(1.d0+0.608d0*q)                                          ! convert virtual temperature to temperature
    else
       q = 0.d0        ! dry
    endif

!-----------------------------------------------------------------------
!     tracer q1, potential temperature
!-----------------------------------------------------------------------
    q1 = theta(lon,lat,eta)

!-----------------------------------------------------------------------
!     tracer q2, absolute value of EPV
!-----------------------------------------------------------------------
    q2 = abs(epv(lon,lat,eta)) * X ! the value of |EPV| scales with X
!
!   The absolute value of Ertel's potential vorticity
!   is selected to avoid the negative EPV values in the
!   Southern Hemisphere. Such negative values might interfere with positive-definite 
!   constraints in the tracer advection algorithm.

  end subroutine test4_baroclinic_wave


!********************************************************************
! Temperature = mean + perturbation temperature
!********************************************************************
  REAL(8) FUNCTION temperature(lon,lat,eta)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: lon, lat, eta

    temperature  = t_mean(eta) + t_deviation(lon,lat,eta)

  END FUNCTION temperature



!***********************************************************************************************
! Horizontally averaged temperature
!***********************************************************************************************
  REAL(8) FUNCTION t_mean(eta)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: eta
    
    if (eta .ge. eta_tropo) then
      t_mean = T0*eta**exponent       ! mean temperature at each level
    else
      t_mean = T0*eta**exponent + delta_T*(eta_tropo-eta)**5
    endif
  END FUNCTION t_mean



!***********************************************************************************************
! Temperature deviation from the horizontal mean 
!***********************************************************************************************
  REAL(8) FUNCTION t_deviation(lon,lat,eta)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: eta, lon, lat
    REAL(8)             :: factor, phi_vertical

    factor       = eta*pi*u0/Rd             
    phi_vertical = (eta - eta0) * 0.5d0*pi

    t_deviation = factor * 1.5d0 * SIN(phi_vertical) * (cos(phi_vertical))**0.5d0 *                &
                  ((-2.d0*(SIN(lat))**6 * ((COS(lat))**2 + 1.d0/3.d0) + 10.d0/63.d0)*              &
                  u0 * (COS(phi_vertical))**1.5d0  +                                               &
                  (8.d0/5.d0*(COS(lat))**3 * ((SIN(lat))**2 + 2.d0/3.d0) - pi/4.d0)*a_omega*0.5d0 )
  END FUNCTION t_deviation



!**************************************************************************
! Surface geopotential
!**************************************************************************  
  REAL(8) FUNCTION surface_geopotential(lon,lat)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: lon, lat
    REAL(8)             :: cos_tmp

    cos_tmp    = u0 * (cos((eta_sfc-eta0)*pi*0.5d0))**1.5d0

    surface_geopotential = ((-2.d0*(SIN(lat))**6 * ((COS(lat))**2 + 1.d0/3.d0) + 10.d0/63.d0)*COS_tmp   &
                 + (8.d0/5.d0*(COS(lat))**3 * ((SIN(lat))**2 + 2.d0/3.d0) - pi/4.d0)*a_omega)*COS_tmp

  END FUNCTION surface_geopotential

!**************************************************************************
! 3D geopotential
!**************************************************************************  
  REAL(8) FUNCTION geopotential(lon,lat,eta)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: lon, lat, eta
    REAL(8)             :: cos_tmp

    cos_tmp    = u0 * (cos((eta-eta0)*pi*0.5d0))**1.5d0

    geopotential = ((-2.d0*(SIN(lat))**6 * ((COS(lat))**2 + 1.d0/3.d0) + 10.d0/63.d0)*COS_tmp   &
                 + (8.d0/5.d0*(COS(lat))**3 * ((SIN(lat))**2 + 2.d0/3.d0) - pi/4.d0)*a_omega)*COS_tmp

    geopotential = geopotential + horiz_mean_geopotential(eta)

  END FUNCTION geopotential

!**************************************************************************
! mean geopotential
!**************************************************************************  
  REAL(8) FUNCTION horiz_mean_geopotential(eta)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: eta
    REAL(8)             :: delta_phi

    horiz_mean_geopotential = T0 * g / gamma * (1.0d0 - eta**(Rd * gamma / g))

    if (eta < eta_tropo) then
      delta_phi = Rd * delta_T * (                         &
        (log(eta/eta_tropo) + 137.d0/60.d0) * eta_tropo**5 &
        - 5.d0 * eta_tropo**4 * eta                        &
        + 5.d0 * eta_tropo**3 * eta**2                     &
        - 10.d0/3.d0 * eta_tropo**2 * eta**3               &
        + 5.d0/4.d0 * eta_tropo * eta**4                   &
        - 1.d0/5.d0 * eta**5)
 
      horiz_mean_geopotential = horiz_mean_geopotential - delta_phi
    end if

  END FUNCTION horiz_mean_geopotential

!********************************************************************
! u wind component
!********************************************************************  
  REAL(8) FUNCTION u_wind(lon,lat,eta,lperturb)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: lon,lat,eta
    LOGICAL, INTENT(IN)  :: lperturb        ! if .true. perturbation is added
    REAL(8) :: phi_vertical, sin_tmp, cos_tmp, r, u_perturb
    REAL(8) :: perturb_lon, perturb_lat


!---------------
!   unperturbed u wind
!---------------
    phi_vertical = (eta - eta0) *0.5d0*pi
    u_wind = (COS(phi_vertical))**1.5d0 * 4.d0 * u0 * (sin(lat))**2 * (cos(lat))**2

!---------------
! if needed (lperturb == .TRUE.) add a Gaussian hill perturbation
!---------------
    IF (lperturb) THEN

       !  center of perturbation, convert from degrees to radians
       perturb_lon = perturbation_longitude*deg2rad
       perturb_lat = perturbation_latitude*deg2rad

       sin_tmp = SIN(perturb_lat)*SIN(lat)
       cos_tmp = COS(perturb_lat)*COS(lat)
                  
       r = ACOS( sin_tmp + cos_tmp*COS(lon-perturb_lon) )         ! great circle distance without radius 'a'
       u_perturb = perturbation_amplitude*EXP(- (r*radius)**2 )   ! perturbation_amplitude determines strength
       u_wind    = u_perturb + u_wind                             ! perturbation + unperturbed u wind
    ENDIF

  END FUNCTION u_wind




!********************************************************************
! v wind component
!********************************************************************  
  REAL(8) FUNCTION v_wind(lon,lat,eta,lperturb)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: lon,lat,eta
    LOGICAL, INTENT(IN)  :: lperturb

    v_wind = 0.0d0

  END FUNCTION v_wind

!********************************************************************
! Calculate eta from height
!********************************************************************
  REAL(8) FUNCTION eta_from_z(lon,lat,z)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: lon, lat, z
    REAL(8) :: eta_new, f, df
    INTEGER :: iter

    REAL(8), PARAMETER ::    &
        initial_eta = 1.0d-7, &
        convergence = 1.0d-13

    eta_from_z = initial_eta

    do iter=0, 25
      f  = - g * z + geopotential(lon,lat,eta_from_z)
      df = - Rd / eta_from_z * temperature(lon,lat,eta_from_z)

      eta_new = eta_from_z - f / df

      if (ABS(eta_from_z - eta_new) < convergence) then
        exit
      end if

      eta_from_z = eta_new
    enddo

  END FUNCTION eta_from_z

!********************************************************************
! Tracers
!********************************************************************  
!-----------------------------------------------------------------------
! Potential temperature
!-----------------------------------------------------------------------
  REAL(8) FUNCTION theta(lon,lat,eta)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: eta, lon, lat
    REAL(8) :: eta_nu, cos_tmp, Y

    eta_nu     = (eta-eta0)*pi*0.5d0
    cos_tmp    = u0 * (cos(eta_nu))**1.5d0

    Y = (-2.d0*(SIN(lat))**6 * ((COS(lat))**2 + 1.d0/3.d0) + 10.d0/63.d0)*2.d0*cos_tmp   &
        + (8.d0/5.d0*(COS(lat))**3 * ((SIN(lat))**2 + 2.d0/3.d0) - pi/4.d0)*a_omega

    theta = t_mean(eta) * eta**(-kappa) + 3.d0/4.d0*pi*u0/Rd * eta**(1.d0-kappa) &
          * SIN(eta_nu) * SQRT(COS(eta_nu)) * Y

  END FUNCTION theta

!-----------------------------------------------------------------------
! Ertel's potential vorticity
!-----------------------------------------------------------------------
  REAL(8) FUNCTION epv(lon,lat,eta)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: eta, lon, lat
    REAL(8) :: perturb_lon, perturb_lat
    REAL(8) :: eta_nu, cos_tmp, Y, r, zeta, K, DK, F
    REAL(8) :: dudeta, dmeanthetadeta, dthetadeta
    REAL(8) :: dthetadphi1, dthetadphi2, dthetadphi
    REAL(8) :: epsilon

    !  center of perturbation, convert from degrees to radians
    perturb_lon = perturbation_longitude*deg2rad
    perturb_lat = perturbation_latitude*deg2rad

    ! initialize epsilon, used to judge the distance of a point to prescribed points
    epsilon = 1.d-6

    eta_nu     = (eta-eta0)*pi*0.5d0
    cos_tmp    = u0 * (cos(eta_nu))**1.5d0

    Y = (-2.d0*(SIN(lat))**6 * ((COS(lat))**2 + 1.d0/3.d0) + 10.d0/63.d0)*2.d0*COS_tmp   &
        + (8.d0/5.d0*(COS(lat))**3 * ((SIN(lat))**2 + 2.d0/3.d0) - pi/4.d0)*a_omega

    ! great circle distance without radius 'a'
    K  = SIN(perturb_lat)*SIN(lat) + COS(perturb_lat)*COS(lat)*COS(lon-perturb_lon)
    DK = SIN(perturb_lat)*COS(lat) - COS(perturb_lat)*SIN(lat)*COS(lon-perturb_lon)
    r  = ACOS(K)         

    ! relative vorticity

    ! version v3 (7/20/2012):
    ! if construct prevents DIVISION BY ZERO in zeta calculation

    if (abs(lat-perturb_lat).le.epsilon .and. abs(lon-perturb_lon).le.epsilon) then  ! grid point is at center position
      zeta = -4.d0/a*cos_tmp*SIN(lat)*COS(lat)*(2.d0-5.d0*(SIN(lat))**2) + perturbation_amplitude/a*TAN(lat)

    else if ( (abs(lat+perturb_lat).le.epsilon .and. abs(lon-(perturb_lon+pi)).le.epsilon) & ! antipode
          .or. abs(lat-pi*0.5d0).le.epsilon                                                & ! north pole
          .or. abs(lat+pi*0.5d0).le.epsilon) then                                            ! south pole
      zeta = -4.d0/a*cos_tmp*SIN(lat)*COS(lat)*(2.d0-5.d0*(SIN(lat))**2)

    else                                                                                     ! all other positions
      zeta = -4.d0/a*cos_tmp*SIN(lat)*COS(lat)*(2.d0-5.d0*(SIN(lat))**2) &
           + perturbation_amplitude/a*EXP(- (r*radius)**2 )              &
           * (TAN(lat) - 2.d0*radius**2*ACOS(K)*DK/SQRT(1.d0-K**2))
    endif

    ! derivative of u with respect to eta
    dudeta = -u0*(SIN(2.d0*lat))**2*3.d0*pi/4.d0*SQRT(COS(eta_nu))*SIN(eta_nu)

    ! derivative of mean theta with respect to eta
    dmeanthetadeta = T0*(Rd*gamma/g - kappa)*eta**(Rd*gamma/g-kappa-1.d0)

    if (eta < eta_tropo) then
      dmeanthetadeta = dmeanthetadeta                          &
        - delta_T * (                                          &
          5.d0*(eta_tropo - eta)**4*eta**(-kappa)              &
          + kappa * (eta_tropo - eta)**5 * eta**(-kappa-1.d0))
    end if

    ! derivative of theta with respect to eta
    dthetadeta = dmeanthetadeta                                                             &
      + 3.d0/4.d0*pi*u0/Rd*(1.d0-kappa)*eta**(-kappa)*SIN(eta_nu)*SQRT(COS(eta_nu))*Y       &
      + 3.d0/8.d0*pi*pi/Rd*eta**(1.d0-kappa)*cos_tmp*Y                                      &
      - 3.d0/16.d0*pi*pi*u0/Rd*eta**(1.d0-kappa)*(SIN(eta_nu))**2*(COS(eta_nu))**(-0.5d0)*Y &
      - 9.d0/8.d0*pi*pi*u0*u0/Rd*eta**(1.d0-kappa)*(SIN(eta_nu))**2*COS(eta_nu)             &
        * (-2.d0*(SIN(lat))**6*((COS(lat))**2+1.d0/3.d0) + 10.d0/63.d0)

    ! derivative of theta with respect to phi
    dthetadphi1 = 2.d0*cos_tmp                                    &
      * (- 12.d0*COS(lat)*(SIN(lat))**5*((COS(lat))**2+1.d0/3.d0) &
          + 4.d0*COS(lat)*(SIN(lat))**7)

    dthetadphi2 = a_omega                                             &
      * (-24.d0/5.d0*SIN(lat)*(COS(lat))**2*((SIN(lat))**2+2.d0/3.d0) &
        + 16.d0/5.d0*(COS(lat))**4*SIN(lat))

!   corrected on June/6/2012
    dthetadphi = 3.d0/4.d0*pi*u0/Rd*eta**(1.d0-kappa)*SIN(eta_nu)*SQRT(COS(eta_nu)) &
      * (dthetadphi1 + dthetadphi2)

    ! Coriolis parameter
    f = 2.d0 * omega * SIN(lat)

    ! EPV, unscaled since real-Earth 'a' and 'Omega' were used, EPV will be scaled in calling routine
    epv = g/p0*(-dudeta*dthetadphi/a - (zeta+f) * dthetadeta)

  END FUNCTION epv

END MODULE
