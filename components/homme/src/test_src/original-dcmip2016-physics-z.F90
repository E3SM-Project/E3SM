!-----------------------------------------------------------------------
!
!  Version:  1.1
!
!  Date:  June 1st, 2016
!
!  Change log:
!
!  SUBROUTINE DCMIP2016_PHYSICS(test, u, v, p, qv, qc, qr, rho,
!                           dt, z, zi, lat, nz, precl, pbl_type, prec_type)
!
!  Input variables:
!     test      (IN) DCMIP2016 test id (1,2,3)
!     u      (INOUT) zonal velocity on model levels (m/s)
!     v      (INOUT) meridional velocity on model levels (m/s)
!     p      (INOUT) total pressure on model levels (Pa)
!     qv     (INOUT) water vapor mixing ratio on model levels (kg/kg)
!     qc     (INOUT) cloud water mixing ratio on model levels (kg/kg)
!     qr     (INOUT) rain water mixing ratio on model levels (kg/kg)
!     rho       (IN) dry air density on model levels (kg/m^3)
!     dt        (IN) time step (s)
!     z         (IN) heights of model levels in the grid column (m)
!     zi        (IN) heights of model interfaces levels in the grid column (m)
!     lat       (IN) latitude of column
!     nz        (IN) number of levels in the column
!     precl     (IN) large-scale precip rate (m/s)
!     pbl_type  (IN) type of planetary boundary layer to use (0,1)
!                    0 = Default Reed-Jablonowski boundary layer
!                    1 = Modified Bryan boundary layer
!     prec_type (IN) type of precipitation/microphysics to use (0,1)
!                    0 = Default Kessler physics routines
!                    1 = Reed-Jablonowski microphysics
!
!  Authors: Paul Ullrich
!           University of California, Davis
!           Email: paullrich@ucdavis.edu
!
!           Kevin Reed
!           Stony Brook University
!           Email: kevin.a.reed@stonybrook.edu
!
!           Kessler is based on a code by Joseph Klemp
!           (National Center for Atmospheric Research)
!
!  Reference:
!
!    Klemp, J. B., W. C. Skamarock, W. C., and S.-H. Park, 2015:
!    Idealized Global Nonhydrostatic Atmospheric Test Cases on a Reduced
!    Radius Sphere. Journal of Advances in Modeling Earth Systems. 
!    doi:10.1002/2015MS000435
!
!=======================================================================

SUBROUTINE DCMIP2016_PHYSICS(test, u, v, p, qv, qc, qr, rho, &
                             dt, z, zi, lat, nz, precl, pbl_type, prec_type)

  IMPLICIT NONE

  !------------------------------------------------
  !   Arguments
  !------------------------------------------------

  INTEGER, INTENT(IN) :: &
            test         ! DCMIP2016 test index

  REAL(8), DIMENSION(nz), INTENT(INOUT) :: &
            u       ,  & ! Zonal velocity on model levels (m/s)
            v       ,  & ! Meridional velocity on model levels (m/s)
            p       ,  & ! Pressure on model levels (Pa)
            qv      ,  & ! Water vapor mixing ratio (kg/kg)
            qc      ,  & ! Cloud water mixing ratio (kg/kg)
            qr           ! Rain water mixing ratio (kg/kg)

  REAL(8), DIMENSION(nz), INTENT(IN) :: &
            rho          ! Dry air density on model levels (kg/m^3)

  REAL(8), INTENT(IN) :: & 
            dt           ! Time step (s)

  REAL(8), DIMENSION(nz), INTENT(IN) :: &
            z       ,  & ! Heights of model levels (m)
            zi           ! Heights of model interfaces (m)

  REAL(8), INTENT(IN) :: &
            lat          ! Latitude of column (radians)

  INTEGER, INTENT(IN) :: &
            nz           ! Number of levels in grid column

  REAL(8), INTENT(INOUT) :: &
            precl        ! Large-scale precip beneath the grid column (mm)

  INTEGER, INTENT(IN) :: &
            pbl_type,    & ! Type of planetary boundary layer to use
            prec_type      ! Type of precipitation to use

  !------------------------------------------------
  ! Physical Constants - MAY BE MODEL DEPENDENT
  !------------------------------------------------
  REAL(8), PARAMETER ::      &
    one      = 1.d0,         & ! One
    gravit  = 9.80616d0,     & ! Gravity (m/s^2)
    rair    = 287.d0,        & ! Gas constant for dry air (J/kg/K)
    cpair   = 1004.5d0,      & ! Specific heat of dry air (J/kg/K)
    latvap  = 2.5d6,         & ! Latent heat of vaporization (J/kg)
    rh2o    = 461.5d0,       & ! Gas constant for water vapor (J/kg/K)
    epsilo  = rair/rh2o,     & ! Ratio of gas constants for dry air to vapor
    zvir    = (rh2o/rair) - 1.d0, & ! Constant for virtual temp. calc. (~0.608)
    a       = 6371220.0,     & ! Reference Earth's Radius (m)
    omega   = 7.29212d-5,    & ! Reference rotation rate of the Earth (s^-1)
    pi      = 4.d0*atan(1.d0)  ! Pi

  !------------------------------------------------
  ! Local Constants for Simple Physics
  !------------------------------------------------
  REAL(8), PARAMETER ::      &
    C        = 0.0011d0,     & ! From Simth and Vogl 2008
    SST_TC   = 302.15d0,     & ! Constant Value for SST
    T0       = 273.16d0,     & ! Control temp for calculation of qsat
    e0       = 610.78d0,     & ! Saturation vapor pressure at T0
    rhow     = 1000.0d0,     & ! Density of Liquid Water 
    Cd0      = 0.0007d0,     & ! Constant for Cd calc. Simth and Vogl 2008
    Cd1      = 0.000065d0,   & ! Constant for Cd calc. Simth and Vogl 2008
    Cm       = 0.002d0,      & ! Constant for Cd calc. Simth and Vogl 2008
    v20      = 20.0d0,       & ! Thresh. wind spd. for Cd Smith and Vogl 2008
    p0       = 100000.0d0,   & ! Constant for potential temp calculation
    pbltop   = 85000.d0,     & ! Top of boundary layer in p
    zpbltop  = 1000.d0,      & ! Top of boundary layer in z
    pblconst = 10000.d0,     & ! Constant for the decay of diffusivity
    T00      = 288.0d0,      & ! Horizontal mean T at surf. for moist baro test
    u0       = 35.0d0,       & ! Zonal wind constant for moist baro test
    latw     = 2.0d0*pi/9.0d0, & ! Halfwidth for  for baro test
    eta0     = 0.252d0,      & ! Center of jets (hybrid) for baro test
    etav     = (1.d0-eta0)*0.5d0*pi, & ! Auxiliary variable for baro test
    q0       = 0.021d0,      & ! Maximum specific humidity for baro test
    kappa    = 0.4d0           ! von Karman constant

  !------------------------------------------------
  ! Variables used in calculation
  !------------------------------------------------
  INTEGER :: k

  REAL(8) ::                  &
    za,                       & ! Altitude of lowest model level (m)
    Tsurf,                    & ! Sea surface temperature (K)
    ps,                       & ! Surface pressure (Pa)
    pres,                     & ! Pressure on model level (Pa)
    presi,                    & ! Pressure on model interface (Pa)
    rhomi,                    & ! Moist pressure on model interface (kg/m^3)
    thetav,                   & ! Virtual potential temperature (K)
    dz,                       & ! Thickness of model level (m)
    deltaqsv,                 & ! Change in specific humidity
    wind,                     & ! Wind speed in the lowest model level (m/s)
    Cd,                       & ! Drag coefficient for momentum
    qsat,                     & ! Saturation specific humidity
    qsats                       ! Saturation specific humidity at sea surface

  ! Matrix coefficients for PBL scheme
  REAL(8) ::                  &
    CA,                       & ! Matrix coefficients for PBL scheme
    CC,                       & ! Matrix coefficients for PBL scheme
    CAm,                      & ! Matrix coefficients for PBL scheme
    CAE,                      & ! Matrix coefficients for PBL scheme
    CBm,                      & ! Matrix coefficients for PBL scheme
    CBE,                      & ! Matrix coefficients for PBL scheme
    CCm,                      & ! Matrix coefficients for PBL scheme
    CCE                         ! Matrix coefficients for PBL scheme

  ! Variables defined on model levels
  REAL(8), DIMENSION(nz) ::   &
    rhom,                     & ! Moist density on model levels
    qsv,                      & ! Specific humidity (kg/kg)
    t,                        & ! Temperature (K)
    theta,                    & ! Potential temperature on model levels
    exner,                    & ! Exner function (p/p0)**(R/cp)
    CEm,                      & ! Matrix coefficients for PBL scheme
    CEE,                      & ! Matrix coefficients for PBL scheme
    CFu,                      & ! Matrix coefficients for PBL scheme
    CFv,                      & ! Matrix coefficients for PBL scheme
    CFt,                      & ! Matrix coefficients for PBL scheme
    CFq                         ! Matrix coefficients for PBL scheme

  ! Variables defined at interfaces
  REAL(8), DIMENSION(nz+1) :: &
    Km,                       & ! Eddy diffusivity for boundary layer
    Ke                          ! Eddy diffusivity for boundary layer

  !------------------------------------------------
  ! Store altitude of lowest model level
  !------------------------------------------------
  za = z(1)

  !------------------------------------------------
  ! Calculate sea surface temperature
  !------------------------------------------------

  ! Moist baroclinic wave test
  if (test .eq. 1) then 
    Tsurf = &
      (T00 + pi*u0/rair * 1.5d0 * sin(etav) * (cos(etav))**0.5d0 * &
      ((-2.d0*(sin(lat))**6 * ((cos(lat))**2 + 1.d0/3.d0) &
        + 10.d0/63.d0) * &
      u0 * (cos(etav))**1.5d0  + &
      (8.d0/5.d0*(cos(lat))**3 * ((sin(lat))**2 + 2.d0/3.d0) &
        - pi/4.d0)*a*omega*0.5d0 ))/ &
      (1.d0+zvir*q0*exp(-(lat/latw)**4))

  ! Tropical cyclone test
  elseif (test .eq. 2) then
    Tsurf = SST_TC

  ! Supercell test
  elseif (test .eq. 3) then
    Tsurf = 0.d0

  ! Invalid test
  else
    write(*,*) 'Invalid test specified in DCMIP2016_PHYSICS', test
    stop
  endif

  !------------------------------------------------
  ! Initialize precipitation rate to zero
  !------------------------------------------------
  precl = 0.d0

  !--------------------------------------------------------------
  ! Calculate moist density, specific humidity and temperature
  !--------------------------------------------------------------
  do k = 1, nz
    rhom(k) = rho(k) * (1.0 + qv(k))
    qsv(k) = qv(k) * rho(k) / rhom(k)
    t(k) = p(k) / (rhom(k) * rair * (one + zvir * qv(k)))
  enddo

  !------------------------------------------------
  ! Large-scale precipitation (Reed-Jablonowski)
  !------------------------------------------------
  if (prec_type .eq. 1) then
    precl = 0.d0

    do k=1, nz
      dz = zi(k+1) - zi(k)
      qsat = epsilo * e0 / p(k) * exp(-(latvap/rh2o) * ((one/t(k))-(one/T0)))
      if (qsv(k) > qsat) then
        deltaqsv = (qsv(k) - qsat) &
          / (one + (latvap/cpair) * epsilo * latvap * qsat / (rair*t(k)**2))
        t(k) = t(k) + latvap / cpair * deltaqsv
        qsv(k) = qsv(k) - deltaqsv
        precl = precl + deltaqsv * rhom(k) * dz / (dt * rhow)

        qv(k) = qsv(k) / (1.0 - qsv(k))
        rhom(k) = rho(k) / (1.0 - qsv(k))
        p(k) = rhom(k) * rair * t(k) * (one + zvir * qv(k))
      endif
    enddo

  !------------------------------------------------
  ! Large-scale precipitation (Kessler)
  !------------------------------------------------
  elseif (prec_type .eq. 0) then
    do k=1,nz
      exner(k) = (p(k) / p0)**(rair/cpair)
      theta(k) = t(k) / exner(k)
    enddo

    CALL KESSLER(   &
      theta,        &
      qv,           &
      qc,           &
      qr,           &
      rho,          &
      exner,        &
      dt,           &
      z,            &
      nz,           &
      precl)

    ! Convert qv to qsv and theta to pressure and temperature
    do k = 1,nz
      qsv(k) = qv(k) / (one + qv(k))
      rhom(k) = rho(k) / (one - qsv(k))

      thetav = theta(k) * (one + zvir * qv(k))
      p(k) = p0 * (rhom(k) * rair * thetav / p0)**(cpair/(cpair-rair))
      t(k) = p(k) / (rhom(k) * rair * (one + zvir * qv(k)))
    enddo

  else
    write(*,*) 'Invalid prec_type specified in DCMIP2016_PHYSICS', prec_type
    stop
  endif

  !----------------------------------------------------
  ! Do not apply surface fluxes or PBL for supercell
  !----------------------------------------------------
  if (test .eq. 3) then
    return
  endif

  !------------------------------------------------
  ! Turbulent mixing coefficients
  !------------------------------------------------
  wind = sqrt(u(1)**2 + v(1)**2)

  if (wind .lt. v20) then
    Cd = Cd0 + Cd1 * wind
  else
    Cd = Cm
  endif

  ! Reed-Jablonowski Boundary layer
  if (pbl_type .eq. 0) then
    do k = 2, nz
      presi = 0.5d0 * (p(k-1) + p(k))
      if (presi .ge. pbltop) then
        Km(k) = Cd * wind * za
        Ke(k) = C  * wind * za
      else
        Km(k) = Cd * wind * za * exp(-(pbltop-presi)**2/pblconst**2)
        Ke(k) = C  * wind * za * exp(-(pbltop-presi)**2/pblconst**2)
      endif
    enddo

  ! Bryan Planetary Boundary Layer
  elseif (pbl_type .eq. 1) then
    do k = 1, nz+1
      if (zi(k) .le. zpbltop) then
        Km(k) = kappa * sqrt(Cd) * wind * zi(k) &
              * (one - zi(k)/zpbltop) * (one - zi(k)/zpbltop)
        Ke(k) = kappa * sqrt(C) * wind * zi(k) &
              * (one - zi(k)/zpbltop) * (one - zi(k)/zpbltop)
      else
        Km(k) = 0.d0
        Ke(k) = 0.d0
      endif
    enddo

  ! Invalid PBL
  else
    write(*,*) 'Invalid pbl_type specified in DCMIP2016_PHYSICS', pbl_type
    stop
  endif

  !------------------------------------------------
  ! Surface fluxes
  !------------------------------------------------

  ! Hydrostatic surface pressure
  ps = 0.d0
  do k = 1, nz
    ps = ps + gravit * rhom(k) * (zi(k+1) - zi(k))
  enddo
  qsats = epsilo * e0 / ps * exp(-latvap / rh2o * ((one/Tsurf)-(one/T0)))

  u(1) = u(1) / (one + dt * Cd * wind / za)
  v(1) = v(1) / (one + dt * Cd * wind / za)
  qsv(1) = (qsv(1) + C * wind * qsats * dt / za) / (one + C * wind * dt / za)
  t(1) = (t(1) + C * wind * Tsurf * dt / za) / (one + C * wind * dt / za)

  qv(1) = qsv(1) / (1.0 - qsv(1))
  rhom(1) = rho(1) / (1.0 - qsv(1))
  p(1) = t(1) * rhom(1) * rair * (one + zvir * qv(1))

  !------------------------------------------------
  ! Boundary layer
  !------------------------------------------------
  do k = 1, nz
    theta(k) = t(k) * (p0 / p(k))**(rair/cpair)
  enddo

  do k = 1, nz
    CA = 0.d0
    CC = 0.d0

    if (k .ne. 1) then
      rhomi = 0.5d0 * (rhom(k-1) + rhom(k))
      CA = rhomi / rhom(k) * dt / (z(k) - z(k-1)) / (zi(k+1) - zi(k))
    endif
    if (k .ne. nz) then
      rhomi = 0.5d0 * (rhom(k) + rhom(k+1))
      CC = rhomi / rhom(k) * dt / (z(k+1) - z(k)) / (zi(k+1) - zi(k))
    endif

    CAm = CA * Km(k)
    CCm = CC * Km(k+1)

    CAE = CA * KE(k)
    CCE = CC * KE(k+1)

    CBm = one + CAm + CCm
    CBE = one + CAE + CCE

    if (k .eq. 1) then
      CEm(1) = CCm / CBm
      CEE(1) = CCE / CBE

      CFu(1) = u(1) / CBm
      CFv(1) = v(1) / CBm
      CFt(1) = theta(1) / CBE
      CFq(1) = qsv(1) / CBE

    else
      CEm(k) = CCm / (CBm - CAm * CEm(k-1))
      CEE(k) = CCE / (CBE - CAE * CEE(k-1))

      CFu(k) = (u(k) + CAm * CFu(k-1)) / (CBm - CAm * CEm(k-1))
      CFv(k) = (v(k) + CAm * CFv(k-1)) / (CBm - CAm * CEm(k-1))
      CFt(k) = (theta(k) + CAE * CFt(k-1)) / (CBE - CAE * CEE(k-1))
      CFq(k) = (qsv(k) + CAE * CFq(k-1)) / (CBE - CAE * CEE(k-1))
    endif
  enddo

  u(nz) = CFu(nz)
  v(nz) = CFv(nz)
  theta(nz) = CFt(nz)
  qsv(nz) = CFq(nz)

  do k=nz-1,1,-1
    u(k) = CEm(k) * u(k+1) + CFu(k)
    v(k) = CEm(k) * v(k+1) + CFv(k)
    theta(k) = CEE(k) * theta(k+1) + CFt(k)
    qsv(k) = CEE(k) * qsv(k+1) + CFq(k)
  enddo

  ! Convert theta to pressure
  do k = 1,nz
    qv(k) = qsv(k) / (one - qsv(k))
    rhom(k) = rho(k) / (one - qsv(k))

    thetav = theta(k) * (one + zvir * qv(k))
    p(k) = p0 * (rhom(k) * rair * thetav / p0)**(cpair/(cpair-rair))
  enddo

  return

END SUBROUTINE DCMIP2016_PHYSICS 

