  subroutine zeff(freq,D,N,nsizes,k2,tt,ice,xr,z_eff,z_ray,kr,qe,qs,rho_e)
  use math_lib
  use optics_lib
  implicit none
  
! Purpose:
!   Simulates radar return of a volume given DSD of spheres
!   Part of QuickBeam v1.03 by John Haynes
!   http://reef.atmos.colostate.edu/haynes/radarsim
!
! Inputs:
!   [freq]      radar frequency (GHz)
!   [D]         discrete drop sizes (um)
!   [N]         discrete concentrations (cm^-3 um^-1)
!   [nsizes]    number of discrete drop sizes
!   [k2]        |K|^2, -1=use frequency dependent default 
!   [tt]        hydrometeor temperature (C)
!   [ice]       indicates volume consists of ice
!   [xr]        perform Rayleigh calculations?
!   [qe]        if using a mie table, these contain ext/sca ...
!   [qs]        ... efficiencies; otherwise set to -1
!   [rho_e]     medium effective density (kg m^-3) (-1 = pure)
!
! Outputs:
!   [z_eff]     unattenuated effective reflectivity factor (mm^6/m^3)
!   [z_ray]     reflectivity factor, Rayleigh only (mm^6/m^3)
!   [kr]        attenuation coefficient (db km^-1)
!
! Created:
!   11/28/05  John Haynes (haynes@atmos.colostate.edu)

! ----- INPUTS -----  
  integer, intent(in) :: ice, xr
  integer, intent(in) :: nsizes
  real*8, intent(in) :: freq,D(nsizes),N(nsizes),tt,qe(nsizes), &
    qs(nsizes), rho_e(nsizes)
  real*8, intent(inout) :: k2
  
! ----- OUTPUTS -----
  real*8, intent(out) :: z_eff,z_ray,kr
    
! ----- INTERNAL -----
  integer :: &
  correct_for_rho        ! correct for density flag
  real*8, dimension(nsizes) :: &
  D0, &                  ! D in (m)
  N0, &                  ! N in m^-3 m^-1
  sizep, &               ! size parameter
  qext, &           ! extinction efficiency
  qbsca, &               ! backscatter efficiency
  rho_ice, &             ! bulk density ice (kg m^-3)
  f                 ! ice fraction
  real*8 :: &
  wl, &                  ! wavelength (m)
  cr                            ! kr(dB/km) = cr * kr(1/km)
  complex*16 :: &
  m                 ! complex index of refraction of bulk form
  complex*16, dimension(nsizes) :: &
  m0                ! complex index of refraction
  
  integer*4 :: i,one
  real*8 :: pi
  real*8 :: eta_sum, eta_mie, const, z0_eff, z0_ray, k_sum, &
            n_r, n_i, dqv(1), dqsc, dg, dph(1)
  integer*4 :: err
  complex*16 :: Xs1(1), Xs2(1)

  one=1
  pi = acos(-1.0)
  rho_ice(:) = 917
  z0_ray = 0.0

! // conversions
  D0 = d*1E-6            ! m
  N0 = n*1E12            ! 1/(m^3 m)
  wl = 2.99792458/(freq*10)   ! m
  
! // dielectric constant |k^2| defaults
  if (k2 < 0) then
    k2 = 0.933
    if (abs(94.-freq) < 3.) k2=0.75
    if (abs(35.-freq) < 3.) k2=0.88
    if (abs(13.8-freq) < 3.) k2=0.925
  endif  
  
  if (qe(1) < -9) then

!   // get the refractive index of the bulk hydrometeors
    if (ice == 0) then
      call m_wat(freq,tt,n_r,n_i)
    else
      call m_ice(freq,tt,n_r,n_i)
    endif
    m = cmplx(n_r,-n_i)
    m0(:) = m
    
    correct_for_rho = 0
    if ((ice == 1) .and. (minval(rho_e) >= 0)) correct_for_rho = 1
    
!   :: correct refractive index for ice density if needed
    if (correct_for_rho == 1) then
      f = rho_e/rho_ice
      m0 = ((2+m0**2+2*f*(m0**2-1))/(2+m0**2+f*(1-m0**2)))**(0.5)
    endif       
    
!   :: Mie calculations
    sizep = (pi*D0)/wl
    dqv(1) = 0.
    do i=1,nsizes
      call mieint(sizep(i), m0(i), one, dqv, qext(i), dqsc, qbsca(i), &
        dg, xs1, xs2, dph, err)
    end do
    
  else
!   // Mie table used
    
    qext = qe
    qbsca = qs
    
  endif
  
! // eta_mie = 0.25*sum[qbsca*pi*D^2*N(D)*deltaD]
!                   <--------- eta_sum --------->
! // z0_eff = (wl^4/!pi^5)*(1./k2)*eta_mie
  eta_sum = 0.
  if (size(D0) == 1) then
    eta_sum = qbsca(1)*(n(1)*1E6)*D0(1)**2
  else
    call avint(qbsca*N0*D0**2,D0,nsizes,D0(1),D0(size(D0,1)),eta_sum)
  endif
 
  eta_mie = eta_sum*0.25*pi
  const = (wl**4/pi**5)*(1./k2)
  z0_eff = const*eta_mie

! // kr = 0.25*cr*sum[qext*pi*D^2*N(D)*deltaD]
!                 <---------- k_sum --------->  
  k_sum = 0.
  if (size(D0) == 1) then
    k_sum = qext(1)*(n(1)*1E6)*D0(1)**2
  else
    call avint(qext*N0*D0**2,D0,nsizes,D0(1),D0(size(D0,1)),k_sum)
  endif
  cr = 10./log(10.)
  kr = k_sum*0.25*pi*(1000.*cr)
     
! // z_ray = sum[D^6*N(D)*deltaD]
  if (xr == 1) then
    z0_ray = 0.
    if (size(D0) == 1) then
      z0_ray = (n(1)*1E6)*D0(1)**6
    else
      call avint(N0*D0**6,D0,nsizes,D0(1),D0(size(D0)),z0_ray)
    endif
  endif
  
! // convert to mm^6/m^3
  z_eff = z0_eff*1E18 !  10.*alog10(z0_eff*1E18)
  z_ray = z0_ray*1E18 !  10.*alog10(z0_ray*1E18)
  
  end subroutine zeff
