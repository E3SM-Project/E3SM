  subroutine dsd(Q,Re,D,N,nsizes,dtype,rho_a,tc, &
             dmin,dmax,apm,bpm,rho_c,p1,p2,p3,fc,scaled)
  use array_lib
  use math_lib 
  implicit none

! Purpose:
!   Create a discrete drop size distribution
!   Part of QuickBeam v1.03 by John Haynes
!   http://reef.atmos.colostate.edu/haynes/radarsim
!
! Inputs:
!   [Q]        hydrometeor mixing ratio (g/kg)
!   [Re]       Optional Effective Radius (microns).  0 = use default.
!   [D]        discrete drop sizes (um)
!   [nsizes]   number of elements of [D]
!   [dtype]    distribution type
!   [rho_a]    ambient air density (kg m^-3)
!   [tc]       temperature (C)
!   [dmin]     minimum size cutoff (um)
!   [dmax]     maximum size cutoff (um)
!   [rho_c]    alternate constant density (kg m^-3)
!   [p1],[p2],[p3]  distribution parameters
!
! Input/Output:
!   [fc]       scaling factor for the distribution
!   [scaled]   has this hydrometeor type been scaled?
!   [apm]      a parameter for mass (kg m^[-bpm])
!   [bmp]      b params for mass
!
! Outputs:
!   [N]        discrete concentrations (cm^-3 um^-1)
!              or, for monodisperse, a constant (1/cm^3)
!
! Requires:
!   function infind
!
! Created:
!   11/28/05  John Haynes (haynes@atmos.colostate.edu)
! Modified:
!   01/31/06  Port from IDL to Fortran 90
!   07/07/06  Rewritten for variable DSD's
!   10/02/06  Rewritten using scaling factors (Roger Marchand and JMH)
 
! ----- INPUTS -----  
  
  integer, intent(in) :: nsizes
  integer, intent(in) :: dtype
  real*8, intent(in) :: Q,D(nsizes),rho_a,tc,dmin,dmax, &
    rho_c,p1,p2,p3
    
! ----- INPUT/OUTPUT -----

  real*8, intent(inout) :: fc(nsizes),apm,bpm,Re
  logical, intent(inout) :: scaled  
    
! ----- OUTPUTS -----

  real*8, intent(out) :: N(nsizes)
  
! ----- INTERNAL -----
  
  real*8 :: &
  N0,D0,vu,np,dm,ld, &             ! gamma, exponential variables
  dmin_mm,dmax_mm,ahp,bhp, &       ! power law variables
  rg,log_sigma_g, &           ! lognormal variables
  rho_e                       ! particle density (kg m^-3)
  
  real*8 :: tmp1, tmp2
  real*8 :: pi

  integer k,lidx,uidx

  pi = acos(-1.0)
  
! // if density is constant, store equivalent values for apm and bpm
  if ((rho_c > 0) .and. (apm < 0)) then
    apm = (pi/6)*rho_c
    bpm = 3.
  endif
  
  select case(dtype)
  
! ---------------------------------------------------------!
! // modified gamma                                        !
! ---------------------------------------------------------!
! :: N0 = total number concentration (m^-3)
! :: np = fixed number concentration (kg^-1)
! :: D0 = characteristic diameter (um)
! :: dm = mean diameter (um)
! :: vu = distribution width parameter

  case(1)  
    if (abs(p1+1) < 1E-8) then

!     // D0, vu are given  
      vu = p3 
      
      if(Re.le.0) then 
          dm = p2
     D0 = gamma(vu)/gamma(vu+1)*dm
      else
     D0 = 2.0*Re*gamma(vu+2)/gamma(vu+3)
      endif
     
      if (scaled .eqv. .false.) then
      
        fc = ( &
             ((D*1E-6)**(vu-1)*exp(-1*D/D0)) / &
             (apm*((D0*1E-6)**(vu+bpm))*gamma(vu+bpm)) &
          ) * 1E-12
     scaled = .true.

      endif       

      N = fc*rho_a*(Q*1E-3)
    
    elseif (abs(p2+1) < 1E-8) then

!     // N0, vu are given    
      np = p1
      vu = p3 
      tmp1 = (Q*1E-3)**(1./bpm)
      
      if (scaled .eqv. .false.) then

        fc = (D*1E-6 / (gamma(vu)/(apm*np*gamma(vu+bpm)))** &
             (1./bpm))**vu
          
        scaled = .true.

      endif

      N = ( &
          (rho_a*np*fc*(D*1E-6)**(-1.))/(gamma(vu)*tmp1**vu) * &
          exp(-1.*fc**(1./vu)/tmp1) &
       ) * 1E-12

    else

!     // vu isn't given
      print *, 'Error: Must specify a value for vu'
      stop
    
    endif
    
! ---------------------------------------------------------!
! // exponential                                           !
! ---------------------------------------------------------!
! :: N0 = intercept parameter (m^-4)
! :: ld = slope parameter (um)

  case(2)
    if (abs(p1+1) > 1E-8) then

!     // N0 has been specified, determine ld
      N0 = p1

      if(Re>0) then

     ! if Re is set and No is set than the distribution is fully defined.
     ! so we assume Re and No have already been chosen consistant with  
     ! the water content, Q.

     ! print *,'using Re pass ...'

     ld = 1.5/Re   ! units 1/um

     N = ( &
               N0*exp(-1*ld*D) &
        ) * 1E-12
    
      else

          tmp1 = 1./(1.+bpm)
      
          if (scaled .eqv. .false.) then
          fc = ((apm*gamma(1.+bpm)*N0)**tmp1)*(D*1E-6)
          scaled = .true.

          endif
     
          N = ( &
          N0*exp(-1.*fc*(1./(rho_a*Q*1E-3))**tmp1) &
     ) * 1E-12

      endif    

    elseif (abs(p2+1) > 1E-8) then

!     // ld has been specified, determine N0
      ld = p2

      if (scaled .eqv. .false.) then

        fc = (ld*1E6)**(1.+bpm)/(apm*gamma(1+bpm))* &
             exp(-1.*(ld*1E6)*(D*1E-6))*1E-12
        scaled = .true.

      endif

      N = fc*rho_a*(Q*1E-3)

    else

!     // ld will be determined from temperature, then N0 follows
      ld = 1220*10.**(-0.0245*tc)*1E-6
      N0 = ((ld*1E6)**(1+bpm)*Q*1E-3*rho_a)/(apm*gamma(1+bpm))
    
      N = ( &
          N0*exp(-1*ld*D) &
          ) * 1E-12
    
    endif
  
! ---------------------------------------------------------!
! // power law                                             !
! ---------------------------------------------------------!
! :: ahp = Ar parameter (m^-4 mm^-bhp)
! :: bhp = br parameter
! :: dmin_mm = lower bound (mm)
! :: dmax_mm = upper bound (mm)

  case(3)

!   :: br parameter
    if (abs(p1+2) < 1E-8) then
!     :: if p1=-2, bhp is parameterized according to Ryan (2000),
!     :: applicatable to cirrus clouds
      if (tc < -30) then
        bhp = -1.75+0.09*((tc+273)-243.16)
      elseif ((tc >= -30) .and. (tc < -9)) then
        bhp = -3.25-0.06*((tc+273)-265.66)
      else
        bhp = -2.15
      endif
    elseif (abs(p1+3) < 1E-8) then      
!     :: if p1=-3, bhp is parameterized according to Ryan (2000),
!     :: applicable to frontal clouds
      if (tc < -35) then
        bhp = -1.75+0.09*((tc+273)-243.16)
      elseif ((tc >= -35) .and. (tc < -17.5)) then
        bhp = -2.65+0.09*((tc+273)-255.66)
      elseif ((tc >= -17.5) .and. (tc < -9)) then
        bhp = -3.25-0.06*((tc+273)-265.66)
      else
        bhp = -2.15
      endif    
    else
!     :: otherwise the specified value is used
      bhp = p1
    endif

!   :: Ar parameter
    dmin_mm = dmin*1E-3
    dmax_mm = dmax*1E-3

!   :: commented lines are original method with constant density
      ! rc = 500.        ! (kg/m^3)
      ! tmp1 = 6*rho_a*(bhp+4)
      ! tmp2 = pi*rc*(dmax_mm**(bhp+4))*(1-(dmin_mm/dmax_mm)**(bhp+4))
      ! ahp = (Q*1E-3)*1E12*tmp1/tmp2

!   :: new method is more consistent with the rest of the distributions
!   :: and allows density to vary with particle size
      tmp1 = rho_a*(Q*1E-3)*(bhp+bpm+1)
      tmp2 = apm*(dmax_mm**bhp*dmax**(bpm+1)-dmin_mm**bhp*dmin**(bpm+1))
      ahp = tmp1/tmp2 * 1E24
      ! ahp = tmp1/tmp2 
 
      lidx = infind(D,dmin)
      uidx = infind(D,dmax)    
      do k=lidx,uidx
 
     N(k) = ( &
        ahp*(D(k)*1E-3)**bhp &
     ) * 1E-12    

      enddo

     ! print *,'test=',ahp,bhp,ahp/(rho_a*Q),D(100),N(100),bpm,dmin_mm,dmax_mm

! ---------------------------------------------------------!
! // monodisperse                                          !
! ---------------------------------------------------------!
! :: D0 = particle diameter (um)

  case(4)
  
    if (scaled .eqv. .false.) then
    
      D0 = p1
      rho_e = (6/pi)*apm*(D0*1E-6)**(bpm-3)
      fc(1) = (6./(pi*D0**3*rho_e))*1E12
      scaled = .true.
      
    endif
    
    N(1) = fc(1)*rho_a*(Q*1E-3)
    
! ---------------------------------------------------------!
! // lognormal                                             !
! ---------------------------------------------------------!
! :: N0 = total number concentration (m^-3)
! :: np = fixed number concentration (kg^-1)
! :: rg = mean radius (um)
! :: log_sigma_g = ln(geometric standard deviation)

  case(5)
    if (abs(p1+1) < 1E-8) then

!     // rg, log_sigma_g are given
      log_sigma_g = p3
      tmp2 = (bpm*log_sigma_g)**2.
      if(Re.le.0) then 
          rg = p2
      else
     rg =Re*exp(-2.5*(log_sigma_g**2))
      endif
 
      if (scaled .eqv. .false.) then
            
        fc = 0.5 * ( &
          (1./((2.*rg*1E-6)**(bpm)*apm*(2.*pi)**(0.5) * &
          log_sigma_g*D*0.5*1E-6)) * &
          exp(-0.5*((log(0.5*D/rg)/log_sigma_g)**2.+tmp2)) &
          ) * 1E-12
     scaled = .true.
          
      endif
             
      N = fc*rho_a*(Q*1E-3)
      
    elseif (abs(p2+1) < 1E-8) then

!     // N0, log_sigma_g are given    
      Np = p1
      log_sigma_g = p3
      N0 = np*rho_a
      tmp1 = (rho_a*(Q*1E-3))/(2.**bpm*apm*N0)
      tmp2 = exp(0.5*bpm**2.*(log_sigma_g))**2.      
      rg = ((tmp1/tmp2)**(1/bpm))*1E6
      
      N = 0.5*( &
        N0 / ((2.*pi)**(0.5)*log_sigma_g*D*0.5*1E-6) * &
     exp((-0.5*(log(0.5*D/rg)/log_sigma_g)**2.)) &
     ) * 1E-12      
      
    else

!     // vu isn't given
      print *, 'Error: Must specify a value for sigma_g'
      stop
    
    endif
    
  end select
  
  end subroutine dsd
