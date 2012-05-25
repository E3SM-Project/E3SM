#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!
!  This file contains the initial condititions for the ASP 2008 summer school test cases
!
MODULE jw 
  IMPLICIT NONE
  !
  !  Functions for setting up initial conditions for the Jablonowski-Williamson test case.
  !
  !  Given longitude (radians), latitude (radians), eta (pressure) and rotation_angle (degrees)
  !  the functions will return temperature, surface geopotential, zonal and meridional wind
  !  components, respectively.
  !
  !  lperturb=.FALSE. result in initial conditions for the steady-state test case.
  !  lperturb=.TRUE.  result in initial conditions for the baroclinic wave test case.
  !
  !     T   : FUNCTION temperature         (lon,lat,eta,rotation_angle)
  !     PHIS: FUNCTION surface_geopotential(lon,lat,rotation_angle)
  !     U   : FUNCTION u_wind              (lon,lat,eta,lperturb,rotation_angle)
  !     V   : FUNCTION v_wind              (lon,lat,eta,lperturb,rotation_angle)
  !
  !  The non-rotated (rotation_angle=0) version of the test cases is described in: 
  !
  !                 Jablonowski, C., and D. L. Williamson, 2006: A baroclinic instability 
  !                 test case for atmospheric model dynamical cores. 
  !                 Quart. J. Roy. Meteor. Soc., 132, 2943-2975.
  !
  !                 Jablonowski, C., and D. L. Williamson, 2006: A Baroclinic Wave Test Case 
  !                 for Dynamical Cores of General Circulation Models: Model Intercomparisons, 
  !                 NCAR Technical Note, NCAR/TN-469+STR, 89 pp. 
  !  
  !  The rotated version simply rotates the initial conditions so that the spherical coordinate
  !  poles do not conicide with the earth's rotation axis. Thereby the Coriolis parameter is
  !  a function of latitude and longitude:
  !
  !      f = 2*Omega*(-cos(lon)*cos(lat)*sin(rotation_angle)+sin(lat)*cos(rotation_angle))
  !
  !  where Omega = 7.292 x 10E-5/s and rotation_angle is the angle between the flow direction
  !  and equator.
  !
  !  Author: Peter Hjort Lauritzen (NCAR, pel@ucar.edu)
  !
  INTEGER,PARAMETER :: r8 = SELECTED_REAL_KIND(12) ! 8 byte real

  REAL(r8), PARAMETER ::       &
       eta_strato = 0.2d0     ,&! tropopause level
       u0         = 35.d0     ,&! 35 m/s
       T0         = 288.d0    ,&! global mean T at surface
       p0         = 100000.d0 ,&! global mean surface pressure
       eta0       = 0.252d0   ,&! center of jets (hybrid)
       !
       radius                 = 10.d0,             & ! radius of the perturbation
       perturbation_amplitude =  1.d0,             & ! amplitude of u perturbation 1 m/s
       perturbation_longitude = 20.d0,             & ! longitudinal position, 20E
       perturbation_latitude  = 40.d0,             & ! latitudinal position, 40N
       !
       !
       !
       Rd         = 287.d0    ,&! gas constant J/(K kg)
       g          = 9.80616d0 ,&! gravitational acceleration (m/s^2)
       a          = 6371229.d0,&! Earth's radius in m
       omega      = 7.29212d-5,&! angular velocity 1/s
       gamma      = 0.005d0   ,&! lapse rate
       pi         = 3.14159265358979323846_R8,&  ! pi
       deg2rad    = pi/180.d0, &
       perturbation_latitude_tracer = 55.d0

CONTAINS
  !
  !********************************************************************
  !
  ! Temperature (equation (6) in Jablonowski and Williamson, 2006)
  !
  !********************************************************************
  !
  REAL(r8) FUNCTION temperature(lon,lat,eta,rotation_angle)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: eta, lon, lat, rotation_angle
    REAL(r8)             :: rot_lon, rot_lat

    IF (ABS(rotation_angle)<1.0E-8) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)
    ENDIF

    temperature  = t_mean(eta) + t_deviation(rot_lon,rot_lat,eta)
  END FUNCTION temperature
  !
  ! Horizontally averaged temperature (equation (4) and (5) in Jablonowski and Williamson (2006))
  !
  REAL(r8) FUNCTION t_mean(eta)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: eta
    REAL(r8)             :: exponent, delta_T
    
    exponent = Rd*gamma/g
    delta_T  = 480000.d0  ! in K, for T mean calculation
    IF (eta.gt.(eta_strato)) THEN
       t_mean = T0*eta**exponent                                ! mean temperature at each level (troposphere)
    ELSE
       t_mean = T0*eta**exponent + delta_T*(eta_strato-eta)**5  ! mean temperature at each level (stratosphere)
    ENDIF
  END FUNCTION t_mean
  !
  ! Temperature deviation from the horizontal mean 
  ! (equation (6) minus horizontally averaged temperature)
  !
  REAL(r8) FUNCTION t_deviation(lon,lat,eta)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: eta, lon, lat
    REAL(r8)             :: factor, phi_vertical, rot_lon, rot_lat, a_omega

    factor       = eta*pi*u0/Rd             
    phi_vertical = (eta - eta0) * 0.5d0*pi
    a_omega      = a*omega

    rot_lon = lon
    rot_lat = lat

    t_deviation = factor * 1.5d0 * SIN(phi_vertical) * (cos(phi_vertical))**0.5d0 *                        &
                  ((-2.d0*(SIN(rot_lat))**6 * ((COS(rot_lat))**2 + 1.d0/3.d0) + 10.d0/63.d0)*              &
                  u0 * (COS(phi_vertical))**1.5d0  +                                                       &
                  (8.d0/5.d0*(COS(rot_lat))**3 * ((SIN(rot_lat))**2 + 2.d0/3.d0) - pi/4.d0)*a_omega*0.5d0 )
  END FUNCTION t_deviation
  !
  !**************************************************************************
  !
  ! Surface geopotential (equaiton (7) in Jablonowski and Williamson, 2006)
  !
  !**************************************************************************
  !  
  REAL(r8) FUNCTION surface_geopotential(lon,lat,rotation_angle)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: lon, lat, rotation_angle
    REAL(r8)             :: cos_tmp, eta_sfc, rot_lon, rot_lat, a_omega

    IF (ABS(rotation_angle)<1.0E-8) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)
    ENDIF

    eta_sfc    = 1.d0
    cos_tmp    = u0 * (cos((eta_sfc-eta0)*pi*0.5d0))**1.5d0
    a_omega    = a*omega

    surface_geopotential = ((-2.d0*(SIN(rot_lat))**6 * ((COS(rot_lat))**2 + 1.d0/3.d0) + 10.d0/63.d0)*COS_tmp   &
                 + (8.d0/5.d0*(COS(rot_lat))**3 * ((SIN(rot_lat))**2 + 2.d0/3.d0) - pi/4.d0)*a_omega)*COS_tmp
  END FUNCTION surface_geopotential
  !
  !********************************************************************
  !
  ! wind components (equation 2 in Jablonowski and Williamson, 2006)
  !
  !********************************************************************
  !  
  REAL(r8) FUNCTION u_wind(lon,lat,eta,lperturb,rotation_angle)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: lon,lat,eta,rotation_angle
    LOGICAL, INTENT(IN)  :: lperturb
    REAL(r8) :: cos_lat, u_lat, phi_vertical, rot_lon, rot_lat, sin_tmp, cos_tmp, r, u_perturb, v_lat
    REAL(r8) :: perturb_lon, perturb_lat, v_tmp

    perturb_lon = perturbation_longitude*deg2rad
    perturb_lat = perturbation_latitude*deg2rad

    IF (ABS(rotation_angle)<1.0E-8) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)
    ENDIF

    phi_vertical = (eta - eta0) *0.5d0*pi
    u_lat = (COS(phi_vertical))**1.5d0 * 4.d0 * u0 * (sin(rot_lat))**2 * (cos(rot_lat))**2
    u_wind = u_lat

    IF (lperturb) THEN

       sin_tmp = SIN(perturb_lat)*SIN(rot_lat)
       cos_tmp = COS(perturb_lat)*COS(rot_lat)
                  
       r = ACOS( sin_tmp + cos_tmp*COS(rot_lon-perturb_lon) )    ! great circle distance
       u_perturb = perturbation_amplitude*EXP(- (r*radius)**2 )
       u_lat     = u_perturb + u_lat
    ENDIF
    IF (ABS(rotation_angle)<1.0E-8) THEN
       u_wind = u_lat
    ELSE
       v_lat = 0.0d0
       !
       ! rotate wind components
       !
       CALL turnwi(u_lat,v_lat, u_wind,v_tmp,lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,-1)
       IF (ABS(u_wind)<1.0E-10) u_wind=0.0d0
    ENDIF
  END FUNCTION u_wind

  REAL(r8) FUNCTION v_wind(lon,lat,eta,lperturb,rotation_angle)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: lon,lat,eta,rotation_angle
    LOGICAL, INTENT(IN)  :: lperturb
    REAL(r8) :: cos_lat, u_lat, phi_vertical, rot_lon, rot_lat, sin_tmp, cos_tmp, r, u_perturb, v_lat
    REAL(r8) :: perturb_lon, perturb_lat, u_tmp

    perturb_lon = perturbation_longitude*deg2rad
    perturb_lat = perturbation_latitude*deg2rad

    IF (ABS(rotation_angle)<1.0E-8) THEN
       v_wind = 0.0d0
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)


       phi_vertical = (eta - eta0) *0.5d0*pi
       u_lat = (COS(phi_vertical))**1.5d0 * 4.d0 * u0 * (sin(rot_lat))**2 * (cos(rot_lat))**2
 
       IF (lperturb) THEN
          
          sin_tmp = SIN(perturb_lat)*SIN(rot_lat)
          cos_tmp = COS(perturb_lat)*COS(rot_lat)
          
          r = ACOS( sin_tmp + cos_tmp*COS(rot_lon-perturb_lon) )    ! great circle distance
          u_perturb = perturbation_amplitude*EXP(- (r*radius)**2 )
          u_lat     = u_perturb + u_lat
       ENDIF

       v_lat = 0.0d0
       !
       ! pole point velocities are not well-defined
       !
       IF (ABS(pi*0.5d0-lat)<1.0E-8.OR.ABS(pi*0.5d0+lat)<1.0E-8) THEN
          v_wind = 0.0d0
       ELSE
          !
          ! rotate wind components
          !
          CALL turnwi(u_lat,v_lat, u_tmp,v_wind,lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,-1)
       ENDIF
    ENDIF
  END FUNCTION v_wind




!-----------------------------------------------------------------------
! Tracer q1 and q2
!-----------------------------------------------------------------------
  REAL(r8) FUNCTION tracer_q1_q2(lon,lat,eta,rotation_angle, eta_c)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: eta, lon, lat, rotation_angle, eta_c
    REAL(r8) :: rot_lon, rot_lat, sin_tmp, cos_tmp, r
    REAL(r8) :: rot_perturb_lon, rot_perturb_lat, tmp
    
    rot_perturb_lon = perturbation_longitude*deg2rad
    rot_perturb_lat = perturbation_latitude_tracer *deg2rad

    IF (ABS(rotation_angle)<1.0E-8) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)
    ENDIF
    sin_tmp = SIN(rot_perturb_lat)*SIN(rot_lat)
    cos_tmp = COS(rot_perturb_lat)*COS(rot_lat)
    r = ACOS( sin_tmp + cos_tmp*COS(rot_lon-rot_perturb_lon) )    ! great circle distance

    tmp = EXP(- ((r*radius)**2 + ((eta-eta_c)/0.1_r8)**2))
    IF (ABS(tmp)<1.0E-8) tmp = 0.0
    tracer_q1_q2 = tmp
  END FUNCTION tracer_q1_q2
  
!-----------------------------------------------------------------------
! Tracer q3
!-----------------------------------------------------------------------
  REAL(r8) FUNCTION tracer_q3(lon,lat,eta,rotation_angle)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: eta, lon, lat, rotation_angle
    REAL(r8) :: rot_lon, rot_lat

    IF (ABS(rotation_angle)<1.0E-8) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)
    ENDIF
    tracer_q3 = 0.5_r8 * ( tanh( 3._r8*abs(rot_lat)-pi ) + 1._r8)

  END FUNCTION tracer_q3

!-----------------------------------------------------------------------
! Tracer q, absolute value of the relative vorticity of the unperturbed initial state
!           multiplied by 10^5
!-----------------------------------------------------------------------
  REAL(r8) FUNCTION tracer_q(lon,lat,eta,rotation_angle)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: eta, lon, lat, rotation_angle
    REAL(r8) :: rot_lon, rot_lat

    IF (ABS(rotation_angle)<1.0E-8) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)
    ENDIF
    tracer_q = abs(-4._r8 * u0/a * (cos((eta-eta0)*pi*0.5_r8))**1.5_r8 * sin(rot_lat) * &
               cos(rot_lat) * (2._r8-5._r8*(sin(rot_lat))**2)) * 1.e5_r8
    if (tracer_q < 1.e-9_r8) tracer_q = 0._r8  !  otherwise error in netcdf file

  END FUNCTION tracer_q






  !
  !******************************************************************************
  !
  ! Subroutines for rotation
  !
  !******************************************************************************
  !

  SUBROUTINE regrot(pxreg,pyreg,pxrot,pyrot,pxcen,pycen,kcall)
    IMPLICIT NONE
!
!----------------------------------------------------------------------
!
!*    conversion between regular and rotated spherical coordinates.
!*
!*    pxreg     longitudes of the regular coordinates
!*    pyreg     latitudes of the regular coordinates
!*    pxrot     longitudes of the rotated coordinates
!*    pyrot     latitudes of the rotated coordinates
!*              all coordinates given in degrees n (negative for s)
!*              and degrees e (negative values for w)
!*    pxcen     regular longitude of the south pole of the rotated grid
!*    pycen     regular latitude of the south pole of the rotated grid
!*
!*    kcall=-1: find regular as functions of rotated coordinates.
!*    kcall= 1: find rotated as functions of regular coordinates.
!
!-----------------------------------------------------------------------
!
      integer kxdim,kydim,kx,ky,kcall
      real(r8) :: pxreg,pyreg,&
                  pxrot,pyrot,&
                  pxcen,pycen
!
!-----------------------------------------------------------------------
!
      real(r8) zsycen,zcycen,zxmxc,zsxmxc,zcxmxc,zsyreg,zcyreg, &
               zsyrot,zcyrot,zcxrot,zsxrot,zpi,zpih
      integer jy,jx

      zpih = pi*0.5d0
!
      !----------------------------------------------------------------------
!
      zsycen = SIN((pycen+zpih))
      zcycen = COS((pycen+zpih))
!
      IF (kcall.eq.1) then
!
         zxmxc  = pxreg - pxcen
         zsxmxc = SIN(zxmxc)
         zcxmxc = COS(zxmxc)
         zsyreg = SIN(pyreg)
         zcyreg = COS(pyreg)
         zsyrot = zcycen*zsyreg - zsycen*zcyreg*zcxmxc
         zsyrot = max(zsyrot,-1.0_r8)
         zsyrot = min(zsyrot,+1.0_r8)
         !
         pyrot = ASIN(zsyrot)
         !
         zcyrot = COS(pyrot)
         zcxrot = (zcycen*zcyreg*zcxmxc +zsycen*zsyreg)/zcyrot
         zcxrot = max(zcxrot,-1.0_r8)
         zcxrot = min(zcxrot,+1.0_r8)
         zsxrot = zcyreg*zsxmxc/zcyrot
         !
         pxrot = ACOS(zcxrot)
         !
         IF (zsxrot<0.0) pxrot = -pxrot
               !
      ELSEIF (kcall.eq.-1) then
         !
         zsxrot = SIN(pxrot)
         zcxrot = COS(pxrot)
         zsyrot = SIN(pyrot)
         zcyrot = COS(pyrot)
         zsyreg = zcycen*zsyrot + zsycen*zcyrot*zcxrot
         zsyreg = max(zsyreg,-1.0_r8)
         zsyreg = min(zsyreg,+1.0_r8)
         !
         pyreg = ASIN(zsyreg)
         !
         zcyreg = COS(pyreg)
         zcxmxc = (zcycen*zcyrot*zcxrot -&
              zsycen*zsyrot)/zcyreg
         zcxmxc = max(zcxmxc,-1.0_r8)
         zcxmxc = min(zcxmxc,+1.0_r8)
         zsxmxc = zcyrot*zsxrot/zcyreg
         zxmxc  = ACOS(zcxmxc)
         IF (zsxmxc<0.0) zxmxc = -zxmxc
         !
         pxreg = zxmxc + pxcen
         !
      ELSE
         WRITE(6,'(1x,''invalid kcall in regrot'')')
         STOP
      ENDIF
    END SUBROUTINE regrot

    SUBROUTINE turnwi(puarg,pvarg,pures,pvres,         &
                      pxreg,pyreg,pxrot,pyrot,   &
                      pxcen,pycen,kcall)
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!*    turn horizontal velocity components between regular and
!*    rotated spherical coordinates.
!
!*    puarg : input u components
!*    pvarg : input v components
!*    pures : output u components
!*    pvres : output v components
!*    pa    : transformation coefficients
!*    pb    :    -"-
!*    pc    :    -"-
!*    pd    :    -"-
!*    pxreg : regular longitudes
!*    pyreg : regular latitudes
!*    pxrot : rotated longitudes
!*    pyrot : rotated latitudes
!*    kxdim              : dimension in the x (longitude) direction
!*    kydim              : dimension in the y (latitude) direction
!*    kx                 : number of gridpoints in the x direction
!*    ky                 : number of gridpoints in the y direction
!*    pxcen              : regular longitude of the south pole of the
!*                         transformed grid
!*    pycen              : regular latitude of the south pole of the
!*                         transformed grid
!*
!*    kcall < 0          : find wind components in regular coordinates
!*                         from wind components in rotated coordinates
!*    kcall > 0          : find wind components in rotated coordinates
!*                         from wind components in regular coordinates
!*    note that all coordinates are given in degrees n and degrees e.
!*       (negative values for s and w)
!
!-----------------------------------------------------------------------

      integer kxdim,kydim,kx,ky,kcall
      real(r8) puarg,pvarg,    &
               pures,pvres,    &
               pa,   pb,       &
               pc,   pd,       &
               pxreg,pyreg,    &
               pxrot,pyrot
      real(r8) pxcen,pycen
!
!-----------------------------------------------------------------------
!
      integer jy,jx
      real(r8) zpih,zsyc,zcyc,zsxreg,zcxreg,zsyreg,zcyreg,zxmxc,&
               zsxmxc,zcxmxc,zsxrot,zcxrot,zsyrot,zcyrot
!
!-----------------------------------------------------------------------
!
      IF (kcall.eq.1) then
         zpih = pi*0.5d0
         zsyc = SIN(pycen+zpih)
         zcyc = COS(pycen+zpih)
         !
         zsxreg = SIN(pxreg)
         zcxreg = COS(pxreg)
         zsyreg = SIN(pyreg)
         zcyreg = COS(pyreg)
         !
         zxmxc  = pxreg - pxcen
         zsxmxc = SIN(zxmxc)
         zcxmxc = COS(zxmxc)
         !
         zsxrot = SIN(pxrot)
         zcxrot = COS(pxrot)
         zsyrot = SIN(pyrot)
         zcyrot = COS(pyrot)
         !
         pa = zcyc*zsxmxc*zsxrot + zcxmxc*zcxrot
         pb = zcyc*zcxmxc*zsyreg*zsxrot - zsyc*zcyreg*zsxrot - &
              zsxmxc*zsyreg*zcxrot
         pc = zsyc*zsxmxc/zcyrot
         pd = (zsyc*zcxmxc*zsyreg + zcyc*zcyreg)/zcyrot
         !
         pures = pa*puarg + pb*pvarg
         pvres = pc*puarg + pd*pvarg
      ELSEIF (kcall.eq.-1) then
         zpih = pi*0.5d0
         zsyc = SIN(pycen+zpih)
         zcyc = COS(pycen+zpih)
         !
         zsxreg = SIN(pxreg)
         zcxreg = COS(pxreg)
         zsyreg = SIN(pyreg)
         zcyreg = COS(pyreg)
         !
         zxmxc  = pxreg - pxcen
         zsxmxc = SIN(zxmxc)
         zcxmxc = COS(zxmxc)
         !
         zsxrot = SIN(pxrot)
         zcxrot = COS(pxrot)
         zsyrot = SIN(pyrot)
         zcyrot = COS(pyrot)
         !
         pa = zcxmxc*zcxrot + zcyc*zsxmxc*zsxrot
         pb = zcyc*zsxmxc*zcxrot*zsyrot + zsyc*zsxmxc*zcyrot -&
              zcxmxc*zsxrot*zsyrot
         pc =-zsyc*zsxrot/zcyreg
         pd = (zcyc*zcyrot - zsyc*zcxrot*zsyrot)/zcyreg
         !
         pures = pa*puarg + pb*pvarg
         pvres = pc*puarg + pd*pvarg
      ELSE
         write(6,'(1x,''invalid kcall in turnwi'')')
         STOP
      ENDIF
    END SUBROUTINE turnwi  

END MODULE jw








MODULE testcases_3_4_5_6

  !=======================================================================
  !
  !  Functions for setting up initial conditions for the dynamical core test cases 3-6
  !  The input parameters depend on the test case (see comments for each test case below).
  !
  !  Given longitude (radians), latitude (radians), eta (pressure) and rotation_angle (degrees)
  !  the functions will return temperature, surface geopotential, zonal and meridional wind
  !  components, respectively.
  !
  !  Author: Christiane Jablonowski (University of Michigan, cjablono@umich.edu)
  !
  !  May/5/2008
  !
  !=======================================================================
  
implicit none
!=======================================================================
!  physical constants
!=======================================================================
  integer, parameter :: r8 = SELECTED_REAL_KIND(12) ! 8 byte real

  real(r8), parameter ::           &
       Rd         = 287.04_r8,     &             ! gas constant J/(K kg)
       cp         = 1004.64_r8,    &             ! specific heat at constant pressure J/(K kg)
       kappa      = Rd/cp,         &             ! kappa = 2/7
       g          = 9.80616_r8,    &             ! gravitational acceleration (m/s^2)
       a          = 6371229._r8,   &             ! Earth's radius in m
       pi         = 3.14159265358979323846_r8,&  ! pi
       omega      = 2._r8*pi/86164._r8, &        ! Earth's angular velocity 1/s
       pih        = pi*0.5_r8,     &             ! pi/2
       deg2rad    = pi/180._r8, &
       tau     = 4*86400,       &                ! period: 4 days expressed in s (used to be 3 days)
       omega_0 = pi*40000._r8 / tau              ! 0.4848 Pa/s

      real(r8),parameter :: p0      = 100000._r8  ! reference pressure 

CONTAINS
  SUBROUTINE advection_vertical(time,eta,eta_top,eta_dot)
  real(r8), intent(in)  :: time,eta,eta_top
  real(r8), intent(out)  :: eta_dot

  ! local
  real(r8) :: s

  s = 2*sqrt(  sin( pi* (eta-eta_top)/(1-eta_top) )  )
  if (s>1) s=1
     
  eta_dot = omega_0*cos(2*pi*time/tau)*sin( s*pi/2)/p0


end subroutine



  SUBROUTINE advection (tracer_variant, lon, lat, height, rotation_angle,  &
                        u_wind, v_wind, temperature, surface_geopotential, &
                        surface_pressure, q5, q6)
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      character*2, intent(in) :: tracer_variant                          ! identifies test variant 'yy', here tracers
                                                                         ! e.g. 0 : no tracer, set to zero
                                                                         !      5 : tracer q5 only
                                                                         !     56 : both tracers q5 and q6
      real(r8), intent(in)  :: lon,                                     & ! longitude in radians
                               lat,                                     & ! latitude in radians
                               height,                                  & ! height of the level in m
                               rotation_angle                             ! alpha in degrees
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      real(r8), intent(out) :: u_wind,                                  & ! zonal wind in m/s
                               v_wind,                                  & ! meridional wind in m/s
                               temperature,                             & ! temperature in K
                               surface_geopotential,                    & ! surface geopotential in m^2/s^2
                               surface_pressure,                        & ! surface pressure in Pa
                               q5,                                      & ! tracer q5
                               q6                                         ! tracer q6
!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
 real(r8),parameter  ::      u0      = (2._r8*pi*a)/(12._r8*86400._r8), & ! circumference / 12 days
                             T0      = 300._r8,                         & ! constant temperature
                             H       = Rd * T0 / g,                     & ! scale height
                             RR      = 1/3._r8,                         & ! horizontal half width divided by 'a'
                             ZZ      = 1000._r8,                        & ! vertical half width
                             z0      = 4500._r8,                        & ! center point in z
                             lambda0 = 1.5_r8*pi,                       & ! center point in longitudes
                             phi0    = 0._r8,                           & ! center point in latitudes
                             slot    = 1._r8/8._r8                        ! half width of the slot in radians
!----------------------------------------------------------------------- 
!     local variables
!-----------------------------------------------------------------------                             
      real(r8) :: alpha
      real(r8) :: sin_tmp, cos_tmp
      real(r8) :: d1, d2, s, r
      
      alpha = rotation_angle*deg2rad
!-----------------------------------------------------------------------
!    initialize the wind components
!-----------------------------------------------------------------------
      u_wind = u0 * (cos(lat)*cos(alpha) + sin(lat)*cos(lon)*sin(alpha))
      v_wind = -u0 *  sin(lon) * sin(alpha)
!-----------------------------------------------------------------------
!     initialization of the vertical velocity: 
!     must be implemented in the dynamical core
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     initialize T (temperature)
!-----------------------------------------------------------------------
      temperature = T0
!-----------------------------------------------------------------------
!     initialize surface geopotential
!-----------------------------------------------------------------------
      surface_geopotential = 0._r8
!-----------------------------------------------------------------------
!     initialize PS (surface pressure)
!-----------------------------------------------------------------------
      surface_pressure = p0
!-----------------------------------------------------------------------
!     Tracer variables
!-----------------------------------------------------------------------
      q5 = 0._r8   ! default
      q6 = 0._r8   ! default
!-----------------------------------------------------------------------
!     tracer q5
!-----------------------------------------------------------------------
      if (tracer_variant(1:1) == '5' .or. tracer_variant(2:2) == '5') then
        sin_tmp = sin(lat) * sin(phi0)
        cos_tmp = cos(lat) * cos(phi0)
        r  = ACOS (sin_tmp + cos_tmp*cos(lon-lambda0))       ! great circle distance without 'a'
        d1 = min( 1._r8, (r/RR)**2 + ((height-z0)/ZZ)**2 )
        q5 = 0.5_r8 * (1._r8 + cos(pi*d1))
      endif
!-----------------------------------------------------------------------
!     tracer q6
!-----------------------------------------------------------------------
      if (tracer_variant(1:1) == '6' .or. tracer_variant(2:2) == '6') then
        sin_tmp = sin(lat) * sin(phi0)
        cos_tmp = cos(lat) * cos(phi0)
        r  = ACOS (sin_tmp + cos_tmp*cos(lon-lambda0))       ! great circle distance without 'a'
        d2 = (r/RR)**2 + ((height-z0)/ZZ)**2
        if (d2 <= 1._r8) then
          q6 = 1._r8
        else
          q6 = 0._r8
        endif
        if ((height > z0) .and. ((phi0-slot) < lat .and. lat < (phi0+slot)) ) q6 = 0._r8   ! slotted ellipse               
      endif
  end subroutine advection
 
!==========================================================================================
! Rossby_Haurwitz wave, wavenumber 4
!==========================================================================================
  SUBROUTINE Rossby_Haurwitz (lon, lat, pressure,                                &
                              u_wind, v_wind, temperature, surface_geopotential, &
                              surface_pressure)
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      real(r8), intent(in)  :: lon,                                     & ! longitude in radians
                               lat,                                     & ! latitude in radians
                               pressure                                   ! pressure at full model level in Pa
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      real(r8), intent(out) :: u_wind,                                  & ! zonal wind in m/s
                               v_wind,                                  & ! meridional wind in m/s
                               temperature,                             & ! temperature in K
                               surface_geopotential,                    & ! surface geopotential in m^2/s^2
                               surface_pressure                           ! surface pressure in Pa

!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
      real(r8),parameter :: u0      = 50._r8,                          &   ! reference wind
                            T0      = 288._r8,                         &   ! reference temperature
                            n       = 4._r8,                           &   ! wavenumber
                            MM      = u0/(n*a),                        &   ! parameter M and p_ref=95500.
                            KK      = u0/(n*a),                        &   ! parameter K
                            gamma   = 0.0065_r8,                       &   ! lapse rate in K/m
                            p_ref   = 95500._r8                            ! reference pressure
                            
!----------------------------------------------------------------------- 
!     local
!----------------------------------------------------------------------- 
      real(r8) :: tmp1, tmp2, tmp3
      real(r8) :: sin_lat, cos_lat, sin_slat, cos_slat
      real(r8) :: exponent_1, exponent_2
      real(r8) :: AA, BB, CC
      real(r8) :: phis_perturb
      
!-----------------------------------------------------------------------
!     initialize the wind components
!-----------------------------------------------------------------------
      cos_lat = cos(lat)
      sin_lat = sin(lat)
      tmp1 = a * MM * cos_lat
      tmp2 = a * KK * cos_lat**(n-1._r8)*(n*sin_lat**2 - cos_lat**2)
      tmp3 = -a * KK * n * cos_lat**(n-1._r8) * sin_lat
      u_wind = tmp1 + tmp2 * cos(n*lon)
      v_wind = tmp3 * sin(n*lon)
!-----------------------------------------------------------------------
!     initialize surface geopotential
!-----------------------------------------------------------------------
      surface_geopotential = 0._r8     
!-----------------------------------------------------------------------
!     initialize surface pressure and temperature
!-----------------------------------------------------------------------
      tmp1       = gamma/(g*T0)
      tmp2       = a*a
      exponent_1 = g/(gamma*Rd)
      exponent_2 = (gamma*Rd)/g
      
      cos_lat = cos(lat)
      AA = tmp2 * (0.5_r8 * MM*(2._r8*omega+MM) * cos_lat**2 + 0.25_r8 * KK**2 * cos_lat**(2._r8*n) * &
                  ( (n+1._r8)*cos_lat**2 + (2._r8*n*n - n - 2._r8)) - 0.5_r8*n*n*KK**2 * cos_lat**(2._r8*(n-1)))
      BB = tmp2 * (2._r8*(omega+MM)*KK/((n+1._r8)*(n+2._r8)) * cos_lat**n * &
                   ( (n*n + 2._r8*n +2._r8) - (n+1._r8)**2 * cos_lat**2 ))
      CC = tmp2 * (0.25_r8 * KK**2 * cos_lat**(2._r8*n) * ( (n+1._r8)*cos_lat**2 - (n+2._r8)))
      phis_perturb = AA + BB * cos(n*lon) + CC * cos(2._r8*n*lon)
      surface_pressure = p_ref * (1._r8 + tmp1*phis_perturb)**exponent_1   ! surface pressure
      temperature      = T0 * (pressure/p_ref)**exponent_2                 ! temperature
      
  end subroutine Rossby_Haurwitz

!==========================================================================================
! Mountain induced Rossby wave
!==========================================================================================
  SUBROUTINE mountain_Rossby (lon, lat, pressure,                                &
                              u_wind, v_wind, temperature, surface_geopotential, &
                              surface_pressure)
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      real(r8), intent(in)  :: lon,                                     & ! longitude in radians
                               lat,                                     & ! latitude in radians
                               pressure                                   ! pressure at full model level in Pa

!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      real(r8), intent(out) :: u_wind,                                  & ! zonal wind in m/s
                               v_wind,                                  & ! meridional wind in m/s
                               temperature,                             & ! temperature in K
                               surface_geopotential,                    & ! surface geopotential in m^2/s^2
                               surface_pressure                           ! surface pressure in Pa
!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
      real(r8),parameter :: u0      = 20._r8,                          &   ! 20 m/s
                            T0      = 288._r8,                         &   ! temperature
                            N2      = g*g/(cp*T0),                     &   ! squared Brunt Vaisala frequency N^2
                            h0      = 2000._r8,                        &   ! amplitude of the mountain, 2km
                            d       = 1500.e3_r8,                      &   ! half width 1500 km
                            lambda0 = 0.5_r8*pi,                       &   ! center point in longitudes
                            phi0    = pi/6._r8,                        &   ! center point in latitudes
                            p_sp    = 93000._r8                            ! pressure at the South Pole in Pa
!-----------------------------------------------------------------------
!   local variables
!-----------------------------------------------------------------------
      real(r8) :: sin_tmp, cos_tmp
      real(r8) :: tmp1, tmp2, tmp3
      real(r8) :: r

!-----------------------------------------------------------------------
!    initialize the wind components
!-----------------------------------------------------------------------
      u_wind = u0 * cos(lat)
      v_wind = 0._r8
!-----------------------------------------------------------------------
!     initialize T (temperature)
!-----------------------------------------------------------------------
      temperature = T0
!-----------------------------------------------------------------------
!     initialize surface geopotential and surface pressure
!-----------------------------------------------------------------------
      tmp1 = (a * N2 * u0)/(2._r8 * g*g * kappa) * (u0/a + 2._r8 * omega)
      tmp2 = N2 / (g*g * kappa)
      sin_tmp = sin(lat) * sin(phi0)
      cos_tmp = cos(lat) * cos(phi0)
      tmp3 = tmp1*((sin(lat))**2 - 1._r8)
      r = a * ACOS (sin_tmp + cos_tmp*cos(lon-lambda0))   ! great circle distance with 'a'
      surface_geopotential = g*h0 * exp(-(r/d)**2)        ! Gaussian profile of the mountain
      surface_pressure     = p_sp * exp( -tmp3 - tmp2*surface_geopotential)

  end subroutine mountain_Rossby
  
!==========================================================================================
! gravity waves
!==========================================================================================
  SUBROUTINE gravity_wave (choice, lon, lat, height,                          &
                           u_wind, v_wind, temperature, surface_geopotential, &
                           surface_pressure)
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      integer, intent(in)   :: choice                                     ! identifies test variant 'x'
                                                                          ! e.g. 0 : no Coriolis, N=0.01 1/s, u0=0  m/s
                                                                          !      1 : no Coriolis, N=0.01 1/s, u0=40 m/s
      real(r8), intent(in)  :: lon,                                     & ! longitude in radians
                               lat,                                     & ! latitude in radians
                               height                                     ! height of the level in m
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      real(r8), intent(out) :: u_wind,                                  & ! zonal wind in m/s
                               v_wind,                                  & ! meridional wind in m/s
                               temperature,                             & ! temperature in K
                               surface_geopotential,                    & ! surface geopotential in m^2/s^2
                               surface_pressure

!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
    real(r8),parameter ::      T0      = 300._r8,                         & ! reference temperature
                            RR      = a/3._r8,                         & ! half width   
                            Lz      = 20.e3_r8,                        & ! vertical wave length, 20 km 
                            delta_theta = 10._r8                         ! potential temperature perturbation amplitude                   

!----------------------------------------------------------------------- 
!     local variables
!----------------------------------------------------------------------- 
      real(r8) :: sin_tmp, cos_tmp
      real(r8) :: tmp1, tmp2, tmp3
      real(r8) :: theta                                                    ! potential temperature
      real(r8) :: theta_mean
      real(r8) :: pres                                                     ! pressure
      real(r8) :: r                                                        ! great circle distance

!-----------------------------------------------------------------------
!     more test case parameters
!----------------------------------------------------------------------- 
      real(r8) :: N2,                                                  &   ! squared Brunt-Vaisala frequency
                  S,                                                   &   ! parameter
                  u0,                                                  &   ! background wind speed
                  lambda0,                                             &   ! center point in longitudes 
                  phi0,                                                &   ! center point in latitudes
                  gw_omega,                                            &   ! rotation
                  p_eq
                  
!-----------------------------------------------------------------------
!    initialize parameters
!-----------------------------------------------------------------------
      lambda0 = pi                                                         ! center point in longitudes
      p_eq    = p0                                                         ! surface pressure at the equator

      select case (choice)
      case (0) 
        N2 = 1.e-4_r8                                                      ! squared Brunt Vaisala frequency N^2
        u0 = 0._r8                                                         ! background wind speed
        phi0 = 0._r8                                                       ! center point in latitudes (0 deg)
        gw_omega = 0._r8                                                   ! no rotation
      case (1) 
        N2 = (g*g)/(cp*T0)                                                 ! squared Brunt Vaisala frequency N^2
        u0 = 0._r8                                                         ! background wind speed
        phi0 = 0._r8                                                       ! center point in latitudes (0 deg)
        gw_omega = 0._r8                                                   ! no rotation
      case (2) 
        N2 = (g*g)/(cp*T0)                                                 ! squared Brunt Vaisala frequency N^2
        u0 = 40._r8                                                        ! background wind speed
        phi0 = 0._r8                                                       ! center point in latitudes (0 deg)
        gw_omega = 0._r8                                                   ! no rotation
      case (3) 
        N2 = (g*g)/(cp*T0)                                                 ! squared  Brunt Vaisala frequency N^2
        u0 = 0._r8                                                         ! background wind speed
        phi0 = pi/4._r8                                                    ! center point in latitudes (45 deg)
        gw_omega = omega                                                   ! Earth's rotation
      end select
      S = g*g/(cp*N2)                                                      ! parameter

!-----------------------------------------------------------------------
!     initialize the wind components
!-----------------------------------------------------------------------
      u_wind = u0 * cos(lat)
      v_wind = 0._r8
!-----------------------------------------------------------------------
!     initialize surface geopotential
!-----------------------------------------------------------------------
      surface_geopotential = 0._r8
!-----------------------------------------------------------------------
!     initialize surface pressure
!-----------------------------------------------------------------------
      tmp1 = (a * N2 * u0)/(2._r8 * g*g * kappa) * (u0/a + 2._r8 * gw_omega)
      surface_pressure = p_eq * exp( - tmp1*((sin(lat))**2) )
!-----------------------------------------------------------------------
!     initialize temperature
!-----------------------------------------------------------------------                 
      pres    = p0 * ( (1._r8 - S/T0) + S/T0 * exp(- (N2*height)/g) )**(cp/Rd)
      sin_tmp = sin(lat) * sin(phi0)
      cos_tmp = cos(lat) * cos(phi0)
      theta_mean = T0 /( T0/S * ((pres/p0)**kappa - 1._r8) + 1._r8 )
      r = a * ACOS (sin_tmp + cos_tmp*cos(lon-lambda0))     ! great circle distance with radius
      if (r < RR) then
        tmp1 = 0.5_r8 * (1._r8 + cos(pi*r/RR))
      else
        tmp1 = 0._r8
      endif
      theta = theta_mean + delta_theta * tmp1 * sin(2._r8*pi*height/Lz)
      temperature = theta * (pres/p0)**kappa

  end subroutine gravity_wave


END MODULE testcases_3_4_5_6













module asp_tests
!
!  This module contains the initial condititions for the baroclinic
!  instability probelms in Jablonowski and Williamson, QJR (2006) 132 
!
use element_mod, only : element_t, timelevels
use cslam_control_volume_mod, only : cslam_struct
use hybrid_mod, only : hybrid_t
use hybvcoord_mod, only : hvcoord_t 
use kinds, only : real_kind

use physical_constants, only : p0, g
use dimensions_mod, only : nlev,np, qsize,nc,ntrac
use control_mod, only : test_case, u_perturb
use cube_mod, only : rotate_grid
use jw, only : u_wind, v_wind, temperature, surface_geopotential, tracer_q1_q2, tracer_q3, perturbation_longitude, perturbation_latitude, deg2rad ! _EXTERNAL
use parallel_mod, only : abortmp

implicit none
private

public :: asp_baroclinic
public :: asp_tracer, asp_advection_vertical, asp_gravity_wave, asp_rossby, asp_mountain



contains


subroutine asp_baroclinic(elem,hybrid,hvcoord,nets,nete, cslam)
!=======================================================================================================!
!  
! new version of JW Baroclinic test case which also allows for rotation
! Code from Jablonowski and Lauritzen
! 
!=======================================================================================================!
    use prim_si_mod, only : preq_hydrostatic

    type(element_t), intent(inout) :: elem(:)
    type(cslam_struct), optional, intent(inout) :: cslam(:)
    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t), intent(in) :: hybrid
    integer :: nets,nete

!   local
    real (kind=real_kind)  :: lat,lon,eta(nlev)
    integer :: i,j,k,ie,idex,t
    logical :: perturbation 
    real(kind=real_kind) :: eta_c


!og, added slotted ellipse from asp advection test as tracer 5
    real (kind=real_kind)  :: height, q6, d2, cos_tmp, sin_tmp, perturb_lat, perturb_lon, r
    real (kind=real_kind)  :: p(np,np,nlev),dp(np,np,nlev)
    real (kind=real_kind)  :: ph(np,np,nlev+1)
    real (kind=real_kind)  :: phi(np,np,nlev)
    real (kind=real_kind),parameter  ::    RR      = 1.0d0/3.0d0,                         & ! horizontal half width divided by 'a'
                             ZZ      = 1000.0d0,                        & ! vertical half width
                             z0      = 4500.0d0,                        & ! center point in z
                             slot    = 1.0d0/8.0d0                        ! half width of the slot in radians




    perturbation = (u_perturb/=0)

    if (perturbation) then
       if (u_perturb /= 1) then
          call abortmp("Error: baroclinc test only supports u_perturb = 0 or 1")
       endif
    endif

    !print *,'JW initial condition  u_perturb, rotation = ',u_perturb,rotate_grid,perturbation
    do ie=nets,nete
    do j=1,np
       do i=1,np
          lon = elem(ie)%spherep(i,j)%lon
          lat = elem(ie)%spherep(i,j)%lat
          do k=1,nlev
             elem(ie)%state%v(i,j,1,k,:) = u_wind(lon,lat,hvcoord%etam(k),perturbation,rotate_grid)
             elem(ie)%state%v(i,j,2,k,:) = v_wind(lon,lat,hvcoord%etam(k),perturbation,rotate_grid)

             elem(ie)%state%T(i,j,k,:)  = temperature(lon,lat,hvcoord%etam(k),rotate_grid)
          enddo
          elem(ie)%state%phis(i,j) = surface_geopotential(lon,lat,rotate_grid)
          elem(ie)%state%ps_v(i,j,:) = p0
          elem(ie)%state%lnps(i,j,:) = LOG(p0)
       enddo
    enddo
    enddo


! initialize passive tracers.  just for testing.  
! first tracer:  temperature at T=0

if (qsize>=1) then
   idex=1
   eta_c = 0.6
   do ie=nets,nete
      do j=1,np
      do i=1,np
         lon = elem(ie)%spherep(i,j)%lon
         lat = elem(ie)%spherep(i,j)%lat
         do k=1,nlev
            elem(ie)%state%Q(i,j,k,idex,:) = tracer_q1_q2(lon,lat,hvcoord%etam(k),rotate_grid, eta_c)
         enddo
      enddo
      enddo
   enddo
endif

if (qsize>=2) then
   idex=2
   eta_c = 1
   do ie=nets,nete
      do j=1,np
      do i=1,np
         lon = elem(ie)%spherep(i,j)%lon
         lat = elem(ie)%spherep(i,j)%lat
         do k=1,nlev
            elem(ie)%state%Q(i,j,k,idex,:) = tracer_q1_q2(lon,lat,hvcoord%etam(k),rotate_grid, eta_c)
         enddo
      enddo
      enddo
   enddo
endif

if (qsize>=3) then
   idex=3
   do ie=nets,nete
      do j=1,np
      do i=1,np
         lon = elem(ie)%spherep(i,j)%lon
         lat = elem(ie)%spherep(i,j)%lat
         do k=1,nlev
            elem(ie)%state%Q(i,j,k,idex,:) = tracer_q3(lon,lat,hvcoord%etam(k),rotate_grid)
         enddo
      enddo
      enddo
   enddo
endif

if (qsize>=4) then
   idex=4
   do ie=nets,nete
      elem(ie)%state%Q(:,:,:,idex,:) = 1
   enddo
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CSLAM tracers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (present(cslam)) then
if (ntrac>=1) then
   idex=1
   do ie=nets,nete
      cslam(ie)%c(:,:,:,idex,:) = 1
   enddo
endif
if (ntrac>=2) then
   idex=2
   eta_c = 1
   do ie=nets,nete
      do j=1,nc
      do i=1,nc
         lon = cslam(ie)%centersphere(i,j)%lon
         lat = cslam(ie)%centersphere(i,j)%lat
         do k=1,nlev
            cslam(ie)%c(i,j,k,idex,:) = tracer_q1_q2(lon,lat,hvcoord%etam(k),rotate_grid, eta_c)
         enddo
      enddo
      enddo
   enddo
endif
if (ntrac>=3) then
   idex=3
   do ie=nets,nete
      do j=1,nc
      do i=1,nc
         lon = cslam(ie)%centersphere(i,j)%lon
         lat = cslam(ie)%centersphere(i,j)%lat
         do k=1,nlev
            cslam(ie)%c(i,j,k,idex,:) = tracer_q3(lon,lat,hvcoord%etam(k),rotate_grid)
         enddo
      enddo
      enddo
   enddo
endif
if (ntrac>=4) then
   idex=4
   eta_c = 0.6
   do ie=nets,nete
      do j=1,nc
      do i=1,nc
         lon = cslam(ie)%centersphere(i,j)%lon
         lat = cslam(ie)%centersphere(i,j)%lat
         do k=1,nlev
            cslam(ie)%c(i,j,k,idex,:) = tracer_q1_q2(lon,lat,hvcoord%etam(k),rotate_grid, eta_c)
         enddo
      enddo
      enddo
   enddo
endif

if (ntrac>=5) then
   do idex=5,ntrac
   do ie=nets,nete
      do k=1,nlev
      do t=1,timelevels
         cslam(ie)%c(1:nc,1:nc,k,idex,t) = 1.0D0
      enddo
      enddo
   enddo
   enddo
endif
endif




if (qsize>=5) then
   idex=5

   perturb_lon = perturbation_longitude*deg2rad
   perturb_lat = perturbation_latitude*deg2rad

   ! now compute PHI, needed to init tracers:	
   do ie=nets,nete
      ! compute height
      do k=1,nlev+1
         do j=1,np
            do i=1,np
               ph(i,j,k)   = hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*elem(ie)%state%ps_v(i,j,1)
            end do
         end do
      end do
      do k=1,nlev
         do j=1,np
            do i=1,np
               p(i,j,k)   = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(i,j,1)
               dp(i,j,k)  = ph(i,j,k+1) - ph(i,j,k)
            end do
         end do
      end do
      call preq_hydrostatic(phi,elem(ie)%state%phis,elem(ie)%state%T(:,:,:,1),p,dp)

      ! init slotted ellipse 
      do j=1,np
      do i=1,np
         lon = elem(ie)%spherep(i,j)%lon
         lat = elem(ie)%spherep(i,j)%lat
         do k=1,nlev

	    height=phi(i,j,k)/g
            sin_tmp = sin(lat) * sin(perturb_lat)
            cos_tmp = cos(lat) * cos(perturb_lat)
            r  = ACOS (sin_tmp + cos_tmp*cos(lon-perturb_lon))       ! great circle distance without 'a'
            d2 = (r/RR)**2 + ((height-z0)/ZZ)**2
            if (d2 <= 1.0d0) then
               q6 = 1.0d0
            else
               q6 = 0.0d0
            endif
            if ((height > z0) .and. ((perturb_lat-slot) < lat .and. lat < (perturb_lat+slot)) ) q6 = 0.0d0   ! slotted ellipse     

            elem(ie)%state%Q(i,j,k,idex,:) = q6

         enddo
      enddo
      enddo
   enddo

endif

end subroutine




!==========================================================================================
! pure 3D advection, time-dependent
!==========================================================================================

subroutine asp_advection_vertical(time,hvcoord,ps_v,eta_dot_dpdn)
    use dimensions_mod, only : np,nlevp
    use testcases_3_4_5_6, only : advection_vertical ! _EXTERNAL
    real (kind=real_kind)  :: time
    real (kind=real_kind)  :: eta_dot_dpdn(np,np,nlevp)
    real (kind=real_kind)  :: ps_v(np,np)
    type (hvcoord_t)                  :: hvcoord


    ! local	
    real (kind=real_kind)  :: dpdn(np,np,nlevp)
    real (kind=real_kind)  :: eta_dot,p1,p2
    integer i,j,k

    ! compute dp/dn at interfaces	
    ! note that for this problem: ps = ps0 = 1000 mb.  
    ! thus  p = A p0 + B ps = (A + B) p0 = eta p0
    ! dp/dn = p0
    eta_dot_dpdn(:,:,1)=0
    eta_dot_dpdn(:,:,nlevp)=0
    do k=2,nlevp-1
       call advection_vertical(time,hvcoord%etai(k),hvcoord%etai(1),eta_dot)
       ! we need eta_dot_dpdn
       eta_dot_dpdn(:,:,k) = eta_dot*p0   ! dpdn(i,j,k)
    enddo
     
end subroutine


subroutine asp_tracer(elem,hybrid,hvcoord,nets,nete)
    use dimensions_mod, only : np,nlev
    use cube_mod, only : rotate_grid
    use prim_si_mod, only : preq_hydrostatic
    use testcases_3_4_5_6, only : advection,g ! _EXTERNAL

    type(element_t), intent(inout) :: elem(:)
    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t), intent(in) :: hybrid
    integer :: nets,nete

!   local
    real (kind=real_kind)  :: lat,lon,eta(nlev),height,q5,q6
    real (kind=real_kind)  :: u_wind,v_wind,temperature
    real (kind=real_kind)  :: surface_geopotential, surface_pressure
    real (kind=real_kind)  :: p(np,np,nlev),dp(np,np,nlev)
    real (kind=real_kind)  :: ph(np,np,nlev+1)
    real (kind=real_kind)  :: phi(np,np,nlev)
    integer :: i,j,k,ie,idex


    ! initialize phis, ps
    height=0   ! not needed for ps, phis
    do ie=nets,nete
    do j=1,np
       do i=1,np
          lon = elem(ie)%spherep(i,j)%lon
          lat = elem(ie)%spherep(i,j)%lat
          do k=1,nlev
             call advection("  ", lon, lat, height, rotate_grid,  &
                 u_wind, v_wind, temperature, surface_geopotential, &
                  surface_pressure, q5, q6)

             elem(ie)%state%phis(i,j) = surface_geopotential
             elem(ie)%state%ps_v(i,j,:) = surface_pressure
             elem(ie)%state%lnps(i,j,:) = LOG(surface_pressure)
             elem(ie)%state%T(i,j,k,:)  = temperature
             elem(ie)%state%v(i,j,1,k,:) = u_wind
             elem(ie)%state%v(i,j,2,k,:) = v_wind
          enddo
       enddo
    enddo
    enddo

    ! now compute PHI, needed to init tracers:	
    do ie=nets,nete
       ! compute height
       do k=1,nlev+1
          do j=1,np
             do i=1,np
                ph(i,j,k)   = hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*elem(ie)%state%ps_v(i,j,1)
             end do
          end do
       end do
       do k=1,nlev
          do j=1,np
             do i=1,np
                p(i,j,k)   = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(i,j,1)
                dp(i,j,k)  = ph(i,j,k+1) - ph(i,j,k)
             end do
          end do
       end do
       call preq_hydrostatic(phi,elem(ie)%state%phis,elem(ie)%state%T(:,:,:,1),p,dp)

       ! init tracers
       do j=1,np
       do i=1,np
          lon = elem(ie)%spherep(i,j)%lon
          lat = elem(ie)%spherep(i,j)%lat
          do k=1,nlev
             call advection("56", lon, lat, phi(i,j,k)/g, rotate_grid,  &
                 u_wind, v_wind, temperature, surface_geopotential, &
                  surface_pressure, q5, q6)


             if (qsize>=2) then
                idex=1
                elem(ie)%state%Q(i,j,k,idex,:) = q5
                idex=2
                elem(ie)%state%Q(i,j,k,idex,:) = q6
             endif
          enddo
       enddo
       enddo
    enddo


end subroutine




subroutine asp_rossby(elem,hybrid,hvcoord,nets,nete)
    use dimensions_mod, only : np,nlev
    use prim_si_mod, only : preq_hydrostatic
    use testcases_3_4_5_6, only : rossby_haurwitz ! _EXTERNAL

    type(element_t), intent(inout) :: elem(:)
    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t), intent(in) :: hybrid
    integer :: nets,nete

!   local
    real (kind=real_kind)  :: lat,lon,eta(nlev),height,q5,q6
    real (kind=real_kind)  :: u_wind,v_wind,temperature
    real (kind=real_kind)  :: surface_geopotential, surface_pressure
    real (kind=real_kind)  :: p(np,np,nlev),dp(np,np,nlev)
    real (kind=real_kind)  :: ph(np,np,nlev+1)
    real (kind=real_kind)  :: phi(np,np,nlev)
    integer :: i,j,k,ie,idex


    ! initialize surface pressure 
    p=0	
    do ie=nets,nete
    do j=1,np
       do i=1,np
          lon = elem(ie)%spherep(i,j)%lon
          lat = elem(ie)%spherep(i,j)%lat
          call Rossby_Haurwitz(lon, lat, p(i,j,1), &
                 u_wind, v_wind, temperature, surface_geopotential, &
                  surface_pressure)
          elem(ie)%state%phis(i,j) = surface_geopotential
          elem(ie)%state%ps_v(i,j,:) = surface_pressure
          elem(ie)%state%lnps(i,j,:) = LOG(surface_pressure)
       enddo
    enddo
    enddo

    ! now we can compute the pressure
    do ie=nets,nete
       do k=1,nlev
       do j=1,np
       do i=1,np
          p(i,j,k)   = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(i,j,1)
       end do
       end do
       end do

       do j=1,np
       do i=1,np
          lon = elem(ie)%spherep(i,j)%lon
          lat = elem(ie)%spherep(i,j)%lat
          do k=1,nlev
             call Rossby_Haurwitz(lon, lat, p(i,j,k), &
                 u_wind, v_wind, temperature, surface_geopotential, &
                  surface_pressure)
             elem(ie)%state%v(i,j,1,k,:) = u_wind
             elem(ie)%state%v(i,j,2,k,:) = v_wind
             elem(ie)%state%T(i,j,k,:)  = temperature
          enddo
       enddo
    enddo
    enddo


end subroutine



subroutine asp_mountain(elem,hybrid,hvcoord,nets,nete)
    use dimensions_mod, only : np,nlev
    use prim_si_mod, only : preq_hydrostatic
    use testcases_3_4_5_6, only : mountain_rossby ! _EXTERNAL

    type(element_t), intent(inout) :: elem(:)
    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t), intent(in) :: hybrid
    integer :: nets,nete

!   local
    real (kind=real_kind)  :: lat,lon,eta(nlev),height,q5,q6
    real (kind=real_kind)  :: u_wind,v_wind,temperature
    real (kind=real_kind)  :: surface_geopotential, surface_pressure, p
    integer :: i,j,k,ie,idex


    ! initialize surface pressure 
    p=0	 ! not used by initial condition
    do ie=nets,nete
       do j=1,np
       do i=1,np
          lon = elem(ie)%spherep(i,j)%lon
          lat = elem(ie)%spherep(i,j)%lat
          do k=1,nlev
             call mountain_Rossby(lon, lat, p, &
                 u_wind, v_wind, temperature, surface_geopotential, &
                  surface_pressure)
             if (k==1) then
                elem(ie)%state%phis(i,j) = surface_geopotential
                elem(ie)%state%ps_v(i,j,:) = surface_pressure
                elem(ie)%state%lnps(i,j,:) = LOG(surface_pressure)
             endif
             elem(ie)%state%v(i,j,1,k,:) = u_wind
             elem(ie)%state%v(i,j,2,k,:) = v_wind
             elem(ie)%state%T(i,j,k,:)  = temperature
          enddo
       enddo
       enddo

    enddo


end subroutine



subroutine asp_gravity_wave(elem,hybrid,hvcoord,nets,nete,choice)
    use dimensions_mod, only : np,nlev
    use cube_mod, only : rotate_grid
    use prim_si_mod, only : preq_hydrostatic
    use testcases_3_4_5_6, only : gravity_wave ! _EXTERNAL

    type(element_t), intent(inout) :: elem(:)
    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t), intent(in) :: hybrid
    integer :: nets,nete, choice

!   local
    real (kind=real_kind)  :: lat,lon,eta(nlev),height,dz,q5,q6
    real (kind=real_kind)  :: u_wind,v_wind,temperature
    real (kind=real_kind)  :: surface_geopotential, surface_pressure
    real (kind=real_kind)  :: p(np,np,nlev),dp(np,np,nlev)
    real (kind=real_kind)  :: ph(np,np,nlev+1)
    real (kind=real_kind)  :: phi(np,np,nlev)
    integer :: i,j,k,ie,idex


    ! initialize surface pressure, geopotential
    do ie=nets,nete
    do j=1,np
       do i=1,np
          lon = elem(ie)%spherep(i,j)%lon
          lat = elem(ie)%spherep(i,j)%lat

          do k=1,nlev
             dz = 500
             height = (nlev-k)*dz + dz/2
             call gravity_wave (choice, lon, lat, height,                          &
                           u_wind, v_wind, temperature, surface_geopotential, &
                           surface_pressure)
             if (k.eq.1) then
                elem(ie)%state%phis(i,j) = surface_geopotential
                elem(ie)%state%ps_v(i,j,:) = surface_pressure
                elem(ie)%state%lnps(i,j,:) = LOG(surface_pressure)
             endif
             elem(ie)%state%v(i,j,1,k,:) = u_wind
             elem(ie)%state%v(i,j,2,k,:) = v_wind
             elem(ie)%state%T(i,j,k,:)  = temperature
          enddo

       enddo
    enddo
    enddo


end subroutine




end module 






