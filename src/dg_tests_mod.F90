#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


module dg_tests_mod
!=======================================================================================================! 
  ! ------------------------
  use kinds, only : real_kind
  ! ------------------------
  use physical_constants, only : omega, g, rearth, dd_pi
  ! ------------------------
  use dimensions_mod, only : nlev, np
  ! ------------------------
  use derivative_mod, only : derivative_t, derivative_stag_t, gradient_wk, divergence, vorticity
  ! ------------------------
  use edge_mod, only : edgebuffer_t, edgevpack, edgevunpack
  ! ------------------------
  use bndry_mod, only : bndry_exchangeV
  ! ------------------------
  use coordinate_systems_mod, only : spherical_polar_t
  ! ------------------------
  use quadrature_mod, only : gauss, gaussian_int, quadrature_t
  ! ------------------------
  use global_norms_mod, only : global_integral, linf_snorm, l1_snorm, l2_snorm
  ! ------------------------
  use time_mod, only : time_at, secpday, timelevel_t
  ! ------------------------
  use element_mod, only : element_t
  ! ------------------------
  use hybrid_mod, only : hybrid_t
  ! ------------------------
  use common_io_mod, only :output_prefix     ! Added to support output_prefix 
  ! ------------------------
!=======================================================================================================! 
 implicit none
 private 
!=======================================================================================================!
 integer, public :: alphatype 
 real (kind=real_kind), public :: alpha_dg
!=======================================================================================================!  
! Constant Parameter for SW-Test Case 1 and 2:								!
!=======================================================================================================!  
  real (kind=real_kind), private, parameter  :: h0     = 1000.0D0     ! height of cosine bell in meters
  real (kind=real_kind), private, parameter  :: gh0_sw1= g*h0         ! geopotential height of cosine bell
  real (kind=real_kind), private, parameter  :: rr     = 1.0/3.0      ! radius of cosine bell
  real (kind=real_kind), public, parameter   :: latc_dg= 0.0D0        ! init latitude of cosine bell center
  real (kind=real_kind), public, parameter   :: lonc_dg= 1.5D0*DD_PI  ! init longitude of cosine bell center
  real (kind=real_kind), private, parameter  :: u0     = 2.0D0*DD_PI*rearth/(12.0*86400.0) ! adv. wind velocity
  real (kind=real_kind), private, parameter  :: gh0_sw2= 2.94D4      ! m^2/s^2 (p. 218 of Williamson)
!=======================================================================================================!  
! Constant Parameter for SW-Test Case 5:								!
!=======================================================================================================! 
  real (kind=real_kind), private, parameter  :: h0_sw5  = 5960.0D0  ! height field
  real (kind=real_kind), private, parameter  :: gh0_sw5 = g*h0_sw5  ! test case 5 geopotential
  real (kind=real_kind), private, parameter  :: u0_sw5  = 20.0D0    ! wind velocity
  real (kind=real_kind), private, parameter  :: h_mtn   = 2000.0D0     ! height of mountain
  real (kind=real_kind), private, parameter  :: R_mtn   = DD_PI/9.0D0  ! radius of mountain
  real (kind=real_kind), private, parameter  :: lat_mtn = DD_PI/6.0D0  ! latitude of mountain
  real (kind=real_kind), private, parameter  :: lon_mtn = 1.5D0*DD_PI  ! longitude of mountain
!=======================================================================================================!
! Constant Parameter for SW-Test Case 8 (Barotropic Instability test by Galewsky, Tellus 2003)
!=======================================================================================================!
  real (kind=real_kind), private, parameter :: h0_gal    = 10000.0D0     ! mean fluid depth (meters)
  real (kind=real_kind), private, parameter :: hhat_gal  = 120.0D0       !  perturbation depth (meters)
  real (kind=real_kind), private, parameter :: gh0_gal   = g*h0_gal      ! mean geopotential
  real (kind=real_kind), private, parameter :: u0_gal    = 80.0D0        ! wind velocity
  real (kind=real_kind), private, parameter :: alpha_gal = 3.0D0         !  meridional scale factor in perturbation
  real (kind=real_kind), private, parameter :: beta_gal  = 15.0D0        ! zonal scale factor in perturbation
  real (kind=real_kind), private, parameter :: lat_gal   = DD_PI/4.0D0   !  latitude of perturbation
  real (kind=real_kind), private, parameter :: lon_gal   = DD_PI/4.0D0   !  longitude of perturbation
  real (kind=real_kind), private, parameter :: gamma_gal = DD_PI/18.0D0  ! zonal velocity field scale factor
  real (kind=real_kind), private, parameter :: eps_gal   = 1.0D-10       !  relative error of balance integral
  integer, parameter :: ngs_gal = 128   ! number of gauss points in balance integral
!=======================================================================================================! 
  public  :: sw1_init_state  ! Initialize test case 1: Cosine Bell
  public  :: sw1_init_pmean  ! Initialize pmean for test case 1
  public  :: sw1_velocity
  public  :: sw1_phi         ! test case 1 on unstaggered grid (nair) 
  public  :: sw1_errors
!=======================================================================================================!
  public  :: sw2_init_state  ! Initialize test case 2: Global steady state nonlinear geostrophic flow
  public  :: sw2_init_pmean
  public  :: sw2_coreolis_init
  public  :: sw2_errors
  public  :: sw2_phi
!=======================================================================================================!
  public  :: sw5_init_state  ! Initialize test case 5: Flow Over a Mountain
  public  :: sw5_init_pmean
  public  :: sw5_velocity
  public  :: sw5_coreolis_init
  public  :: sw5_errors
  public  :: sw5_invariants
  public  :: sw5_phi
!=======================================================================================================!
  public  :: sw1_init  
  public  :: sw2_init
  public  :: sw5_init      
!=======================================================================================================!
  public  :: galewsky_init_state   !routines for "Galewsky test" 
  public  :: galewsky_velocity
  public  :: galewsky_coriolis
  public  :: galewsky_lat_function
  public  :: galewsky_phi
  public  :: galewsky_init
!=======================================================================================================!
  real (kind=real_kind), public :: Imass0,Ienergy0,Ipenst0
  common /comintegral/Imass0,Ienergy0,Ipenst0

  real (kind=real_kind), public :: Ibalance_tot
  common /comtc8/Ibalance_tot
!=======================================================================================================!
 contains
!=======================================================================================================!
  ! ===========================================
  ! sw1_init_pmean
  ! ===========================================
  function sw1_init_pmean() result(pmean)
    implicit none
    real (kind=real_kind) :: pmean 

    pmean = gh0_sw1
  end function sw1_init_pmean
  ! ===========================================
  ! sw2_init_pmean
  ! ===========================================
  function sw2_init_pmean() result(pmean)
    implicit none
    real (kind=real_kind) :: pmean 

    pmean = gh0_sw2
  end function sw2_init_pmean
  ! ===========================================
  ! sw5_init_pmean
  ! ===========================================
  function sw5_init_pmean() result(pmean)
    implicit none
    real (kind=real_kind) :: pmean 

    pmean = gh0_sw5
  end function sw5_init_pmean

!=======================================================================================================!
!	sw1_init_state:
!	Initialize for shallow water testcase 1
!=======================================================================================================!
subroutine sw1_init_state(elem,nets,nete,pmean)
!=======================================================================================================!
    type(element_t) :: elem(:)
    integer, intent(in) :: nets,nete
    real (kind=real_kind):: pmean
    ! Local variables
    integer :: ie,k,nm1,n0,np1,nstep
!=======================================================================================================!
    nm1   = 1
    n0    = 2
    np1   = 3
    nstep = 0   

    pmean = sw1_init_pmean()       ! set pmean

    do ie=nets,nete
    do k=1,nlev
          elem(ie)%state%v(:,:,:,k,n0)= sw1_velocity(elem(ie)%spherep(:,:),elem(ie)%Dinv)	  
	  elem(ie)%state%ht(:,:,k) = sw1_phi(elem(ie)%spherep(:,:),pmean,latc_dg,lonc_dg)
	  elem(ie)%state%hs(:,:,k) = 0.0D0
!	   elem(ie)%state%psi(:,:,k)= elem(ie)%state%ht(:,:,k)
    end do
    end do

end subroutine sw1_init_state
!=======================================================================================================!
!	sw1_velocity: 											!
!	u=  u_0 ( cos(theta)*cos(alpha) + sin(theta)*cos(lambda)*sin(alpha) )				!
!	v= -u_0 ( sin(lambda)*sin(alpha) )								!
!=======================================================================================================!
function sw1_velocity(sphere,D) result(v)

    type (spherical_polar_t), intent(in) :: sphere(np,np)
    real (kind=real_kind),    intent(in) :: D(2,2,np,np)
    real (kind=real_kind)                :: v(np,np,2)

    ! Local variables

    real (kind=real_kind) :: snlon(np,np)
    real (kind=real_kind) :: cslon(np,np)
    real (kind=real_kind) :: snlat(np,np)
    real (kind=real_kind) :: cslat(np,np)
    real (kind=real_kind) :: csalpha
    real (kind=real_kind) :: snalpha
    real (kind=real_kind) :: V1,V2

    integer i,j,k

    csalpha = COS(alpha_dg)
    snalpha = SIN(alpha_dg)

    do j=1,np
       do i=1,np
          snlat(i,j) = SIN(sphere(i,j)%lat)
          cslat(i,j) = COS(sphere(i,j)%lat)
          snlon(i,j) = SIN(sphere(i,j)%lon)
          cslon(i,j) = COS(sphere(i,j)%lon)
       end do
    end do

    do k=1,nlev
       do j=1,np
          do i=1,np

             V1 =   u0*(cslat(i,j)*csalpha + snlat(i,j)*cslon(i,j)*snalpha)
             V2 =  -u0*(snlon(i,j)*snalpha)
             v(i,j,1)= V1
             v(i,j,2)= V2

             ! =====================================================
             ! map sphere velocities onto the contravariant cube velocities
             ! using the D^-T mapping matrix (see Loft notes for details)
             ! =====================================================

           !v(i,j,1)= V1*D(1,1,i,j) + V2*D(1,2,i,j)
           !v(i,j,2)= V1*D(2,1,i,j) + V2*D(2,2,i,j)
          end do
       end do
    end do

end function sw1_velocity
!=======================================================================================================!
!	sw1_height:
!	Initialize cosine bell centered at latc, lonc  on (np,np) points
!=======================================================================================================!
function sw1_phi(sphere,pmean,latc,lonc) result(ppsi)

    type (spherical_polar_t), intent(in) :: sphere(np,np) ! spherical coordinates of element
    real (kind=real_kind), intent(in)    :: pmean         ! mean geopotential
    real (kind=real_kind), intent(in)    :: latc          ! latitude of center of cosine bell
    real (kind=real_kind), intent(in)    :: lonc          ! longitude of center of cosine bell
    real (kind=real_kind)                :: ppsi(np,np)   ! geopotential

    real (kind=real_kind) :: snlat,cslat,lon, cslon, coef
    real (kind=real_kind) :: h, gh00, trm
    real (kind=real_kind) :: r,A, csalpha, snalpha
    integer i,j

    csalpha = COS(alpha_dg)
    snalpha = SIN(alpha_dg)

    do j=1,np
       do i=1,np
          snlat=SIN(sphere(i,j)%lat)
          cslat=COS(sphere(i,j)%lat)
          cslon=COS(sphere(i,j)%lon)
          lon  =sphere(i,j)%lon

          A = SIN(latc)*snlat + COS(latc)*cslat*COS(lon - lonc)   ! face 4 this == - cslat*snlon

          r=ACOS(A)
          if (r<rr) then
             h = (h0/2.0D0)*(1.0D0 + COS(DD_PI*r/rr))
          else
             h = 0.0D0
          end if
          ppsi(i,j)= h                                 !for sw1

       end do
    end do


  end function sw1_phi
!=======================================================================================================!
!	sw1_errors:
!=======================================================================================================!
subroutine sw1_errors(elem,iounit,tl,pmean,hybrid,nets,nete)
    type(element_t)                   :: elem(:)
    integer                           :: iounit
    type (TimeLevel_t)   , intent(in) :: tl         ! model time struct
    real (kind=real_kind), intent(in) :: pmean      ! mean geopotential
    type (hybrid_t)      , intent(in) :: hybrid     
    integer, intent(in)               :: nets
    integer, intent(in)               :: nete

    ! Local variables 

    real (kind=real_kind) :: latc,lonc
    real (kind=real_kind) :: lamdot,lon
    real (kind=real_kind) :: l1,l2,linf,time_tmp
    real (kind=real_kind) :: cslonc,snlonc
    real (kind=real_kind) :: cslatc,snlatc

    real (kind=real_kind) :: pt(np,np,nets:nete)
    real (kind=real_kind) :: p(np,np,nets:nete)

    integer ie
!=======================================================================================================!
    !$OMP BARRIER
    if (tl%nstep == 1) then
       if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
          open(iounit+0,file=TRIM(output_prefix)//"swtc1.l1.errors",form="formatted")
          open(iounit+1,file=TRIM(output_prefix)//"swtc1.l2.errors",form="formatted")
          open(iounit+2,file=TRIM(output_prefix)//"swtc1.linf.errors",form="formatted")
       end if
    end if

    ! Compute the expected lat,lon of the cosine bell:
    ! solid body rotations by Z axis and by X axis ..

    time_tmp   = Time_at(tl%nstep)
    lamdot = u0/rearth
    lon  = MODULO(lamdot*time_tmp,2.0D0*real(DD_PI,kind=real_kind))
    snlatc = SIN(lon)*SIN(alpha_dg)
    cslatc = SQRT(1.0D0 - snlatc*snlatc)

    cslonc = (COS(lon)*COS(lonc_dg) - SIN(lon)*SIN(lonc_dg)*COS(alpha_dg))/cslatc
    snlonc = (COS(lon)*SIN(lonc_dg) + SIN(lon)*COS(lonc_dg)*SIN(alpha_dg))/cslatc

    latc = ASIN(SIN(lon)*SIN(alpha_dg))
    lonc = ATAN2(snlonc,cslonc)

    if (lonc < 0.0D0)lonc = lonc + 2.0D0*DD_PI
!=======================================================================================================!
!    do ie=nets,nete
!       pt(:,:,ie)=sw1_phi(elem(ie)%spherep(:,:),pmean,latc,lonc)*g
!       pt(:,:,ie)=pt(:,:,ie) + pmean
!       p(:,:,ie) =elem(ie)%state%p(:,:,1,tl%n0) + pmean
!    end do
!=======================================================================================================!
    do ie=nets,nete
       pt(:,:,ie)=sw1_phi(elem(ie)%spherep(:,:),pmean,latc,lonc)
       p(:,:,ie) =elem(ie)%state%ht(:,:,1)
    end do
!=======================================================================================================!
    l1   = l1_snorm(elem,p(:,:,nets:nete),  pt(:,:,nets:nete),hybrid,np,nets,nete)
    l2   = l2_snorm(elem,p(:,:,nets:nete),  pt(:,:,nets:nete),hybrid,np,nets,nete)
    linf = linf_snorm(p(:,:,nets:nete),pt(:,:,nets:nete),hybrid,np,nets,nete)

    if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
        write(iounit+0,30)time_tmp/secpday,l1
        write(iounit+1,30)time_tmp/secpday,l2
        write(iounit+2,30)time_tmp/secpday,linf
        ! MNL: output error to stdout
        ! (Copied / pasted from tc1_errors in shallow_water_mod.F90)
        write(*,'(f6.2,a,3e15.7)') time_tmp/secpday,' days  l1,l2,linf=',&
                 l1,l2,linf
    end if

10  format("      time           l1              l2              linf           ")
20  format("====================================================================")
30  format(f11.6,4x,e13.6)

end subroutine sw1_errors
!=======================================================================================================!
!	sw2_init_state:
!	Initialize for shallow water testcase 2
!=======================================================================================================!
subroutine sw2_init_state(elem,nets,nete,pmean)
!=======================================================================================================!
    type(element_t) :: elem(:)
    integer, intent(in) :: nets,nete
    real (kind=real_kind):: pmean
    ! Local variables
    integer :: ie,k,nm1,n0,np1,nstep
!=======================================================================================================!
    nm1   = 1
    n0    = 2
    np1   = 3
    nstep = 0    
    pmean = sw2_init_pmean()     ! (m^2/s^2) set pmean value specified in Williamson, et. al.
    ! p 218, for test case 2.

    do ie=nets,nete
       elem(ie)%fcor=sw2_coreolis_init(elem(ie)%spherep)
       do k=1,nlev
         !elem(ie)%state%v(:,:,:,k,n0)= sw1_velocity(elem(ie)%spherep(:,:),elem(ie)%Dinv)  ! sw2 vel same as sw1
          elem(ie)%state%v(:,:,:,k,n0)= sw1_velocity(elem(ie)%spherep(:,:),elem(ie)%Dinv)  ! sw2 vel same as sw1
	  elem(ie)%state%ht(:,:,k) = sw2_phi(elem(ie)%spherep(:,:))/g
          elem(ie)%state%hs(:,:,k) = 0.0D0
!	  elem(ie)%state%psi(:,:,k)= elem(ie)%state%ht(:,:,k)
       end do
    end do

end subroutine sw2_init_state
!=======================================================================================================!
!	sw2_geopotential:										!
!=======================================================================================================!
function sw2_phi(sphere) result(p)

    type (spherical_polar_t), intent(in) :: sphere(np,np)
    real (kind=real_kind)                :: p(np,np)

    ! Local variables

    real (kind=real_kind) :: cslon
    real (kind=real_kind) :: snlat
    real (kind=real_kind) :: cslat
    real (kind=real_kind) :: csalpha
    real (kind=real_kind) :: snalpha
    real (kind=real_kind) :: coef

    integer i,j

    csalpha = COS(alpha_dg)
    snalpha = SIN(alpha_dg)
    coef = rearth*omega*u0 + (u0**2)/2.0D0

    do j=1,np
       do i=1,np
          snlat = SIN(sphere(i,j)%lat)
          cslat = COS(sphere(i,j)%lat)
          cslon = COS(sphere(i,j)%lon)
          p(i,j)= gh0_sw2 -coef*( -cslon*cslat*snalpha + snlat*csalpha )**2
       end do
    end do

end function sw2_phi
!=======================================================================================================!
!	sw2 coriolis parameter:										!
!=======================================================================================================!
function sw2_coreolis_init(sphere) result(fcor)

    type (spherical_polar_t), intent(in) :: sphere(np,np)
    real (kind=real_kind) :: fcor(np,np)

    ! Local variables

    real (kind=real_kind) :: cslon
    real (kind=real_kind) :: snlat
    real (kind=real_kind) :: cslat
    real (kind=real_kind) :: csalpha
    real (kind=real_kind) :: snalpha
    integer                  :: i,j

    csalpha = COS(alpha_dg)
    snalpha = SIN(alpha_dg)

    do j=1,np
       do i=1,np
          snlat = SIN(sphere(i,j)%lat)
          cslat = COS(sphere(i,j)%lat)
          cslon = COS(sphere(i,j)%lon)

          fcor(i,j) = 2.0D0*omega*(-cslon*cslat*snalpha + snlat*csalpha)
       end do
    end do

end function sw2_coreolis_init
!=======================================================================================================!
!	sw2_errors:											!
!=======================================================================================================!
subroutine sw2_errors(elem, iounit, tl, pmean, hybrid, nets, nete)
    type(element_t)                   :: elem(:)
    integer                           :: iounit
    type (TimeLevel_t)   , intent(in) :: tl         ! model time struct
    real (kind=real_kind), intent(in) :: pmean      ! mean geopotential
    type (hybrid_t)      , intent(in) :: hybrid     
    integer, intent(in)               :: nets
    integer, intent(in)               :: nete

    ! Local variables 

    real (kind=real_kind) :: l1,l2,linf,time_tmp
    real (kind=real_kind) :: pt(np,np,nets:nete)
    real (kind=real_kind) :: p(np,np,nets:nete)

    integer ie, npts
!=======================================================================================================!
    time_tmp= Time_at(tl%nstep)

    !$OMP BARRIER
    if (tl%nstep == 1) then
       if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
          open(iounit+0,file=TRIM(output_prefix)//"swtc2.l1.errors",form="formatted")
          open(iounit+1,file=TRIM(output_prefix)//"swtc2.l2.errors",form="formatted")
          open(iounit+2,file=TRIM(output_prefix)//"swtc2.linf.errors",form="formatted")
       end if
    end if
!=======================================================================================================!
!    do ie=nets,nete
!       pt(:,:,ie)=sw2_phi(elem(ie)%spherep(:,:))
!       pt(:,:,ie)=pt(:,:,ie) + pmean
!       p(:,:,ie) =elem(ie)%state%p(:,:,1,tl%n0) + pmean
!    end do
!=======================================================================================================!
    do ie=nets,nete
       pt(:,:,ie)=sw2_phi(elem(ie)%spherep(:,:))/g
       p(:,:,ie) =elem(ie)%state%ht(:,:,1)
    end do
!=======================================================================================================!
    npts = np
    l1   = l1_snorm(elem,p(:,:,nets:nete),  pt(:,:,nets:nete),hybrid,npts,nets,nete)
    l2   = l2_snorm(elem,p(:,:,nets:nete),  pt(:,:,nets:nete),hybrid,npts,nets,nete)
    linf = linf_snorm(p(:,:,nets:nete),pt(:,:,nets:nete),hybrid,npts,nets,nete)

    if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
        write(iounit+0,30)time_tmp/secpday,l1
        write(iounit+1,30)time_tmp/secpday,l2
        write(iounit+2,30)time_tmp/secpday,linf
        ! MNL: output error to stdout
        ! (Copied / pasted from sw1_errors)
        write(*,'(f6.2,a,3e15.7)') time_tmp/secpday,' days  l1,l2,linf=',&
                 l1,l2,linf
    end if

10  format("      time           l1              l2              linf           ")
20  format("====================================================================")
30  format(f11.6,4x,e13.6)

end subroutine sw2_errors
!=======================================================================================================!
! Barotropic Instability test (Galewsky test) 
!=======================================================================================================!
subroutine galewsky_init_state(elem,nets,nete,pmean,deriv)
!=======================================================================================================!
    type(element_t) :: elem(:)
    integer, intent(in) :: nets,nete
    real (kind=real_kind):: pmean
    type (derivative_t)  :: deriv
    type (quadrature_t)   :: gs

    ! Local variables
    integer :: ie,k,nm1,n0,np1,nstep
    real (kind=real_kind):: hts(np,np),hll(np,np)
!=======================================================================================================!
    nm1   = 1
    n0    = 2
    np1   = 3
    nstep = 0
    pmean = gh0_gal

    gs=gauss(ngs_gal)
    deallocate(gs%points)
    deallocate(gs%weights)


    do ie=nets,nete
       elem(ie)%fcor= galewsky_coriolis(elem(ie)%spherep)
       do k=1,nlev
          elem(ie)%state%v(:,:,:,k,n0)= galewsky_velocity(elem(ie)%spherep(:,:),elem(ie)%Dinv)
          elem(ie)%state%ht(:,:,k)= galewsky_phi(elem(ie)%spherep(:,:),gs) /  g
          elem(ie)%state%hs(:,:,k)= 0.0D0
       end do
    end do

  end subroutine galewsky_init_state
!=======================================================================================================!
  function galewsky_lat_integral(lat,gs) result(latint)
!======================================================================================================!
    real (kind=real_kind), intent(in):: lat
    type(quadrature_t), intent(in)   :: gs        ! gaussian points/wts
    real (kind=real_kind) :: latint

    real (kind=real_kind):: xlat(ngs_gal)
    real (kind=real_kind):: a,b, ab, sm , pi2
    integer :: i

    pi2  = 0.5D0*DD_PI

        a = -pi2
        b = lat
       ab = (b-a) * 0.5D0

         do i = 1, ngs_gal
           xlat(i) = (a + b + (b-a) * gs%points(i)) * 0.5D0
         enddo

         sm = 0.0D0
         do i = 1, ngs_gal
           sm = sm  !+   galewsky_lat_function(xlat(i)) * gs%weights(i)
         enddo

         latint = sm * ab 

  end function galewsky_lat_integral
!=======================================================================================================!
  function galewsky_lat_function(lat) result(val)
!=======================================================================================================!

    real (kind=real_kind), intent(in)    :: lat

    real (kind=real_kind) :: val
    real (kind=real_kind) :: umax,epsil,ulat, small 
    real (kind=real_kind) :: cf, lat0,lat1
    real (kind=real_kind) :: pi2

    pi2  = 0.5D0*DD_PI
    lat0 =  DD_PI / 7.0D0
    lat1 =  pi2 - lat0
    umax = 80.0D0
    epsil = umax * exp(4.0D0/(lat1-lat0)**2)
    cf = 2.0D0 * omega * sin(lat)
    small = 1.0D-12 

      if (abs(lat) == pi2) then
           val = 0.0D0
        else

             if ((lat >= lat0) .and. (lat <= lat1)) then
                ulat =  epsil * exp( 1.0D0 /((lat - lat0 +small)*(lat - lat1 +small)) )
              else
                ulat = 0.0D0
             endif

          val =  ulat * (cf *rearth + ulat*tan(lat))

        endif

  end function galewsky_lat_function
!=======================================================================================================!
  function galewsky_velocity(sphere,D) result(v)
!=======================================================================================================!
    type (spherical_polar_t), intent(in) :: sphere(np,np)
    real (kind=real_kind),    intent(in) :: D(2,2,np,np)
    real (kind=real_kind)                :: v(np,np,2)

    ! Local variables

    integer :: i,j
    real (kind=real_kind) :: umax,epsil,ulat, small 
    real (kind=real_kind) :: cf, lat0,lat1
    real (kind=real_kind) :: pi2
    real (kind=real_kind) :: lat
    real (kind=real_kind) :: V1,V2

    pi2  = 0.5D0*DD_PI
    lat0 =  DD_PI / 7.0D0
    lat1 =  pi2 - lat0
    umax = 80.0D0
    epsil = umax * exp(4.0D0/(lat1-lat0)**2)
    cf = 2.0D0 * omega * sin(lat)
    small = 1.0D-12 

    do j=1,np
       do i=1,np
          lat   = sphere(i,j)%lat

             if ((lat >= lat0) .and. (lat <= lat1)) then
                ulat =  epsil * exp( 1.0D0 /((lat - lat0 + small)*(lat - lat1 + small)) )
              else
                ulat = 0.0D0
             endif

             V1 =  ulat
             V2 =  0.0D0

          !v(i,j,1)= V1*D(1,1,i,j) + V2*D(1,2,i,j)   !contravariant vectors 
          !v(i,j,2)= V1*D(2,1,i,j) + V2*D(2,2,i,j)
          v(i,j,1)= V1
          v(i,j,2)= V2

       end do
    end do

  end function galewsky_velocity
!=======================================================================================================!
   function galewsky_coriolis(sphere) result(fcor)
!=======================================================================================================!
    type (spherical_polar_t), intent(in) :: sphere(np,np)
    real (kind=real_kind) :: fcor(np,np)

    ! Local variables

    real (kind=real_kind) :: snlat
    integer                  :: i,j

    do j=1,np
       do i=1,np
          snlat = SIN(sphere(i,j)%lat)
          fcor(i,j) = 2.0D0*omega* snlat
       end do
    end do

  end function galewsky_coriolis
!=======================================================================================================!
  function galewsky_phi(sphere,gs) result(p)
!=======================================================================================================!
    type (spherical_polar_t), intent(in) :: sphere(np,np)
    type(quadrature_t), intent(in)       :: gs        ! gaussian points/wts
    real (kind=real_kind)                :: p(np,np)

    ! Local variables

    real (kind=real_kind) :: lat
    real (kind=real_kind) :: lon
    real (kind=real_kind) :: alon, blat
    real (kind=real_kind) :: ghprime, galint 

    integer i,j

    do j=1,np
       do i=1,np

          lat   = sphere(i,j)%lat
          lon   = sphere(i,j)%lon

          if(lon.le.DD_PI) alon = alpha_gal*(lon)
          if(lon.gt.DD_PI) alon = alpha_gal*(2.0D0*DD_PI - lon)

          blat =  beta_gal*(lat_gal - lat)

          ghprime = cos(lat) * g*hhat_gal* exp(-alon*alon) * exp(-blat*blat)

          galint  =  galewsky_lat_integral(lat,gs) 
          p(i,j)  = gh0_gal + ghprime  - galint  

       end do
    end do

  end function galewsky_phi
!=======================================================================================================!
  subroutine galewsky_init(elem,tl,ie,k,pmean)
!=======================================================================================================!
    type(element_t) :: elem(:)
    type(TimeLevel_t), intent(in):: tl
    integer, intent(in)  :: ie,k
    type(quadrature_t)   :: gs        ! gaussian points/wts
    real (kind=real_kind):: pmean, a,b
! Local variables
    integer :: nm1,n0,np1,nstep
!=======================================================================================================!
    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep

    pmean = gh0_gal

    gs=gauss(ngs_gal)
    deallocate(gs%points)
    deallocate(gs%weights)

    elem(ie)%fcor = galewsky_coriolis(elem(ie)%spherep)
    elem(ie)%state%v(:,:,:,k,n0) = galewsky_velocity(elem(ie)%spherep(:,:),elem(ie)%Dinv)
    elem(ie)%state%ht(:,:,k) =  galewsky_phi(elem(ie)%spherep(:,:),gs)/g
    elem(ie)%state%hs(:,:,k) = 0.0D0

  end subroutine galewsky_init
!=======================================================================================================!
!	sw5_init_state:
!	Initialize for shallow water testcase 5
!=======================================================================================================!
subroutine sw5_init_state(elem,nets,nete,pmean,deriv)
!=======================================================================================================!
    type(element_t) :: elem(:)
    integer, intent(in) :: nets,nete
    real (kind=real_kind):: pmean
    type (derivative_t)  :: deriv
    ! Local variables
    integer :: ie,k,nm1,n0,np1,nstep
    real (kind=real_kind):: hts(np,np),hll(np,np)
!=======================================================================================================!
    nm1   = 1
    n0    = 2
    np1   = 3
    nstep = 0    
    pmean = sw5_init_pmean()     ! (m^2/s^2) set pmean value specified in Williamson, et. al.

    ! p 218, for test case 5
    do ie=nets,nete
       elem(ie)%fcor= sw5_coreolis_init(elem(ie)%spherep)
       do k=1,nlev
          elem(ie)%state%v(:,:,:,k,n0)= sw5_velocity(elem(ie)%spherep(:,:),elem(ie)%Dinv)	   
	  Call sw5_phi(elem(ie)%spherep(:,:),hts,hll)           
	  elem(ie)%state%ht(:,:,k)= hts(:,:)/g 
          elem(ie)%state%hs(:,:,k)= hll(:,:)/g
!	  elem(ie)%state%psi(:,:,k)= elem(ie)%state%ht(:,:,k)
       end do
    end do

end subroutine sw5_init_state
!=======================================================================================================!
!	sw5_velocity: test case 5									!
!	u=  u_0 ( cos(theta)*cos(0) + sin(theta)*cos(lambda)*sin(0) )					!
!	v= -u_0 ( sin(lambda)*sin(0) )									!
!=======================================================================================================!
function sw5_velocity(sphere,D) result(v)

    type (spherical_polar_t), intent(in) :: sphere(np,np)
    real (kind=real_kind),    intent(in) :: D(2,2,np,np)
    real (kind=real_kind)                :: v(np,np,2)

    ! Local variables

    real (kind=real_kind) :: snlon(np,np)
    real (kind=real_kind) :: cslon(np,np)
    real (kind=real_kind) :: snlat(np,np)
    real (kind=real_kind) :: cslat(np,np)
    real (kind=real_kind) :: csalpha
    real (kind=real_kind) :: snalpha
    real (kind=real_kind) :: V1,V2
    integer i,j

    csalpha = COS(0.0D0)
    snalpha = SIN(0.0D0)

    do j=1,np
       do i=1,np
          snlat(i,j) = SIN(sphere(i,j)%lat)
          cslat(i,j) = COS(sphere(i,j)%lat)
          snlon(i,j) = SIN(sphere(i,j)%lon)
          cslon(i,j) = COS(sphere(i,j)%lon)
       end do
    end do

    do j=1,np
       do i=1,np

          V1 =   u0_sw5*(cslat(i,j)*csalpha + snlat(i,j)*cslon(i,j)*snalpha)
          V2 =  -u0_sw5*(snlon(i,j)*snalpha)

          ! =====================================================
          ! map sphere velocities onto the contravariant cube velocities
          ! using the D mapping matrix (see Loft notes for details)
          ! =====================================================

           v(i,j,1)= V1
           v(i,j,2)= V2
          !v(i,j,1)= V1*D(1,1,i,j) + V2*D(1,2,i,j)
          !v(i,j,2)= V1*D(2,1,i,j) + V2*D(2,2,i,j)
       end do
    end do

end function sw5_velocity
!=======================================================================================================!
!	sw5 coriolis parameter:										!
!=======================================================================================================!
function sw5_coreolis_init(sphere) result(fcor)

    type (spherical_polar_t), intent(in) :: sphere(np,np)
    real (kind=real_kind) :: fcor(np,np)

    ! Local variables

    real (kind=real_kind) :: cslon
    real (kind=real_kind) :: snlat
    real (kind=real_kind) :: cslat
    real (kind=real_kind) :: csalpha
    real (kind=real_kind) :: snalpha
    integer               :: i,j

    csalpha = COS(0.0D0)
    snalpha = SIN(0.0D0)

    do j=1,np
       do i=1,np
          snlat = SIN(sphere(i,j)%lat)
          cslat = COS(sphere(i,j)%lat)
          cslon = COS(sphere(i,j)%lon)

          fcor(i,j) = 2.0D0*omega*(-cslon*cslat*snalpha + snlat*csalpha)
       end do
    end do

end function sw5_coreolis_init
!=======================================================================================================! 
!	sw5: ht, hs											!
!=======================================================================================================! 
subroutine sw5_phi(sphere,p,ps)

    type (spherical_polar_t), intent(in) :: sphere(np,np)
    real (kind=real_kind), intent(out)    :: ps(np,np)
    real (kind=real_kind), intent(out)    :: p(np,np)

    ! Local variables

    real (kind=real_kind) :: cslon
    real (kind=real_kind) :: snlat
    real (kind=real_kind) :: cslat
    real (kind=real_kind) :: csalpha
    real (kind=real_kind) :: snalpha
    real (kind=real_kind) :: coef
    real (kind=real_kind) :: lon
    real (kind=real_kind) :: lat
    real (kind=real_kind) :: rsq

    integer i,j,k

    csalpha = COS(0.0D0)
    snalpha = SIN(0.0D0)
    coef = rearth*omega*u0_sw5 + (u0_sw5**2)/2.0D0

    do j=1,np
       do i=1,np
          lat = sphere(i,j)%lat
          lon = sphere(i,j)%lon
          rsq=MIN((lat-lat_mtn)**2 + (lon-lon_mtn)**2,R_mtn**2)
          ps(i,j)=g*h_mtn*(1.0D0 - SQRT(rsq)/R_mtn)
       end do
    end do

    do j=1,np
       do i=1,np
          snlat = SIN(sphere(i,j)%lat)
          cslat = COS(sphere(i,j)%lat)
          cslon = COS(sphere(i,j)%lon)
          p(i,j)= gh0_sw5 - coef*( -cslon*cslat*snalpha + snlat*csalpha )**2  
       end do
    end do

end subroutine sw5_phi
!=======================================================================================================! 
!	sw5_errors:											!
!=======================================================================================================! 
subroutine sw5_errors(elem,iounit, tl, pmean, fstub, simday, hybrid, nets, nete)
    type(element_t)                   :: elem(:)
    integer              , intent(in) :: iounit
    type (TimeLevel_t)   , intent(in) :: tl         ! model time struct
    real (kind=real_kind), intent(in) :: pmean      ! mean geopotential
    character(len=*)     , intent(in) :: fstub      ! file path stub
    integer              , intent(in) :: simday     ! current day of simulation
    type (hybrid_t)      , intent(in) :: hybrid     
    integer              , intent(in) :: nets
    integer              , intent(in) :: nete

    ! Local variables 

    real (kind=real_kind) :: pt(np,np,nets:nete)
    real (kind=real_kind) :: p(np,np,nets:nete)

    real (kind=real_kind) :: vt(np,np,2,nets:nete)
    real (kind=real_kind) :: v(np,np,2,nets:nete)

    real (kind=real_kind) :: l1,l2,linf
    integer ie,npts
!=======================================================================================================!
    if (tl%nstep == 0) then
       if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
          open(iounit+0,file=TRIM(output_prefix)//"swtc5.l1.errors",form="formatted")
          open(iounit+1,file=TRIM(output_prefix)//"swtc5.l2.errors",form="formatted")
          open(iounit+2,file=TRIM(output_prefix)//"swtc5.linf.errors",form="formatted")
       end if
    end if

    ! ======================================================
    ! read in the reference state for this simulated day...
    ! ======================================================
!=======================================================================================================!
#ifdef _REFSOLN
    call ref_state_read(pt(:,:,nets:nete),vt(:,:,:,nets:nete),fstub,simday+1,nets,nete)
#endif
    ! Get the state variables from the simulation (layer 1 only)...

    do ie=nets,nete
       p(:,:,ie)   = elem(ie)%state%p(:,:,1,tl%n0) + pmean
       v(:,:,:,ie) = elem(ie)%state%v(:,:,:,1,tl%n0)
#ifndef _REFSOLN
       pt(:,:,ie)=p(:,:,ie)
       vt(:,:,:,ie)=v(:,:,:,ie)
#endif
    end do
!=======================================================================================================!
    npts=np

    l1   = l1_snorm(elem,p(:,:,nets:nete),  pt(:,:,nets:nete),hybrid,npts,nets,nete)
    l2   = l2_snorm(elem,p(:,:,nets:nete),  pt(:,:,nets:nete),hybrid,npts,nets,nete)
    linf = linf_snorm(p(:,:,nets:nete),pt(:,:,nets:nete),hybrid,npts,nets,nete)

    if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
        write(iounit+0,30)REAL(simday),l1
        write(iounit+1,30)REAL(simday),l2
        write(iounit+2,30)REAL(simday),linf
        ! MNL: it really seems to me like this whole routine should only run
        !      if _REFSOLN is defined, otherwise it will just output 0 for
        !      all the errors...
#ifdef _REFSOLN
        print *,simday, "L1=",l1
        print *,simday, "L2=",l2
        print *,simday, "Linf=",linf
#endif
    end if
    !$OMP BARRIER

30  format(f11.6,4x,e13.6)

end subroutine sw5_errors
!=======================================================================================================! 
!	sw5_invariants											!
!=======================================================================================================!
subroutine sw5_invariants(elem, iounit,tl,pmean,edge2,deriv,hybrid,nets,nete)
    type(element_t) :: elem(:)
    integer              , intent(in) :: iounit
    type (TimeLevel_t)  , intent(in)  :: tl    
    real (kind=real_kind), intent(in) :: pmean
    type (EdgeBuffer_t)               :: edge2
    type (derivative_t)               :: deriv
    type (hybrid_t)      , intent(in) :: hybrid
    integer              , intent(in) :: nets,nete

    ! Local variables

    integer :: ie,i,j,kptr

    real (kind=real_kind) :: vco(np,np,2)
    real (kind=real_kind) :: gv(np,np,2)

    real (kind=real_kind) :: zeta(np,np,nets:nete)
    real (kind=real_kind) :: E(np,np)

    real (kind=real_kind) :: mass(np,np,nets:nete)
    real (kind=real_kind) :: energy(np,np,nets:nete)
    real (kind=real_kind) :: penst(np,np,nets:nete)
    real (kind=real_kind) :: vort(np,np,nets:nete)
    real (kind=real_kind) :: div(np,np,nets:nete)
    real (kind=real_kind) :: divsq(np,np,nets:nete)

    real (kind=real_kind) :: v1,v2
    real (kind=real_kind) :: lenscale
    real (kind=real_kind) :: h,hstar,hs

    real (kind=real_kind) :: time
    real (kind=real_kind) :: Imass,Ienergy,Ipenst,Ivort,Idiv,Idivsq

    integer               :: npts
    integer               :: n0
!=======================================================================================================!
    n0 = tl%n0

    lenscale=rearth

    !$OMP BARRIER
    if (tl%nstep == 0) then
       if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
          open(iounit+0,file=TRIM(output_prefix)//"swtc5.mass",status="unknown",form="formatted")
          open(iounit+1,file=TRIM(output_prefix)//"swtc5.energy",status="unknown",form="formatted")
          open(iounit+2,file=TRIM(output_prefix)//"swtc5.penst",status="unknown",form="formatted")
          open(iounit+3,file=TRIM(output_prefix)//"swtc5.vort",status="unknown",form="formatted")
          open(iounit+4,file=TRIM(output_prefix)//"swtc5.div",status="unknown",form="formatted")
          open(iounit+5,file=TRIM(output_prefix)//"swtc5.rmsdiv",status="unknown",form="formatted")
       endif
    endif
!=======================================================================================================!
    do ie=nets,nete

       do j=1,np
          do i=1,np

             v1= elem(ie)%state%v(i,j,1,1,n0)
             v2= elem(ie)%state%v(i,j,2,1,n0)

             vco(i,j,1) = elem(ie)%met(1,1,i,j)*v1 + elem(ie)%met(1,2,i,j)*v2
             vco(i,j,2) = elem(ie)%met(2,1,i,j)*v1 + elem(ie)%met(2,2,i,j)*v2

             gv(i,j,1) = elem(ie)%metdet(i,j)*v1
             gv(i,j,2) = elem(ie)%metdet(i,j)*v2

          end do
       end do

       zeta(:,:,ie)  = vorticity(vco,deriv)/rearth
       div(:,:,ie)   = divergence(gv,deriv)/rearth

       do j=1,np
          do i=1,np
             zeta(i,j,ie)=elem(ie)%mp(i,j)*(zeta(i,j,ie)/elem(ie)%metdet(i,j))
             div(i,j,ie)=elem(ie)%mp(i,j)*(div(i,j,ie)/elem(ie)%metdet(i,j))
          end do
       end do

       kptr=0
       call edgeVpack(edge2, zeta(1,1,ie), 1, kptr,elem(ie)%desc)

       kptr=1
       call edgeVpack(edge2, div(1,1,ie), 1, kptr,elem(ie)%desc)

    end do
!=======================================================================================================!
    call bndry_exchangeV(hybrid,edge2)
!=======================================================================================================!
    do ie=nets,nete      

       kptr=0
       call edgeVunpack(edge2, zeta(1,1,ie), 1, kptr, elem(ie)%desc)

       kptr=1
       call edgeVunpack(edge2, div(1,1,ie), 1, kptr, elem(ie)%desc)

       do j=1,np
          do i=1,np
             zeta(i,j,ie) = elem(ie)%rmp(i,j)*zeta(i,j,ie)
             div(i,j,ie)  = elem(ie)%rmp(i,j)*div(i,j,ie)
             penst(i,j,ie)= 0.5D0*(zeta(i,j,ie) + elem(ie)%fcor(i,j))**2.0D0
          end do
       end do

       do j=1,np
          do i=1,np

             v1     = elem(ie)%state%v(i,j,1,1,n0)
             v2     = elem(ie)%state%v(i,j,2,1,n0)

             vco(i,j,1) = elem(ie)%met(1,1,i,j)*v1 + elem(ie)%met(1,2,i,j)*v2
             vco(i,j,2) = elem(ie)%met(2,1,i,j)*v1 + elem(ie)%met(2,2,i,j)*v2
             E(i,j) = 0.5D0*( vco(i,j,1)*v1 + vco(i,j,2)*v2 )

          end do
       end do

       do j=1,np
          do i=1,np
!             hstar         = (elem(ie)%state%p(i,j,1,n0) + pmean)/g
!             hs            = elem(ie)%state%ps(i,j)/g
             h	           = elem(ie)%state%ht(i,j,1)
             hs            = elem(ie)%state%hs(i,j,1)
             hstar         = h - hs
             mass(i,j,ie)  = hstar 
             energy(i,j,ie)= hstar*E(i,j) + 0.5D0*g*(h**2.0D0 - hs**2.0D0)
             penst(i,j,ie) = penst(i,j,ie)/hstar
             divsq(i,j,ie) = div(i,j,ie)*div(i,j,ie)  
             vort(i,j,ie)  = zeta(i,j,ie)
          end do
       end do

    end do
!=======================================================================================================!
    npts=np
    Imass   = global_integral(elem,mass(:,:,nets:nete),hybrid,npts,nets,nete)
    Ienergy = global_integral(elem,energy(:,:,nets:nete),hybrid,npts,nets,nete)
    Ipenst  = global_integral(elem,penst(:,:,nets:nete),hybrid,npts,nets,nete)
    Ivort   = global_integral(elem,vort(:,:,nets:nete),hybrid,npts,nets,nete)
    Idiv    = global_integral(elem,div(:,:,nets:nete),hybrid,npts,nets,nete)
    Idivsq  = global_integral(elem,divsq(:,:,nets:nete),hybrid,npts,nets,nete)

    time   = Time_at(tl%nstep)

    if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
       if (tl%nstep==0) then
          Imass0   = Imass
          Ienergy0 = Ienergy
          Ipenst0  = Ipenst
       end if

       write(iounit+0,*)time/secpday,(Imass-Imass0)/Imass0,Imass0
       write(iounit+1,*)time/secpday,(Ienergy-Ienergy0)/Ienergy0,Ienergy0
       write(iounit+2,*)time/secpday,(Ipenst-Ipenst0)/Ipenst0,Ipenst0
       write(iounit+3,*)time/secpday,Ivort
       write(iounit+4,*)time/secpday,Idiv
       write(iounit+5,*)time/secpday,SQRT(Idivsq)

       write(6,*)""
       write(6,*)time/secpday,"mass          =",(Imass-Imass0)/Imass0
       write(6,*)time/secpday,"total energy  =",(Ienergy-Ienergy0)/Ienergy0
       write(6,*)time/secpday,"pot enstrophy =",(Ipenst-Ipenst0)/Ipenst0
       write(6,*)time/secpday,"vorticity     =",Ivort
       write(6,*)time/secpday,"divergence    =",Idiv
       write(6,*)time/secpday,"RMS divergence=",SQRT(Idivsq)

    end if

  end subroutine sw5_invariants
!=======================================================================================================! 
subroutine sw1_init(elem,tl,ie,k,pmean)
!=======================================================================================================!
    type(element_t) :: elem(:)
    type(TimeLevel_t), intent(in):: tl
    integer, intent(in)  :: ie,k
    real (kind=real_kind):: pmean
! Local variables
    integer :: nm1,n0,np1,nstep
!=======================================================================================================!
    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep    

    elem(ie)%state%v(:,:,:,k,n0)= sw1_velocity(elem(ie)%spherep(:,:),elem(ie)%Dinv)	  
    elem(ie)%state%ht(:,:,k) = sw1_phi(elem(ie)%spherep(:,:),pmean,latc_dg,lonc_dg)
    elem(ie)%state%hs(:,:,k) = 0.0D0    
!    elem(ie)%state%psi(:,:,k)= elem(ie)%state%ht(:,:,k)

end subroutine sw1_init
!=======================================================================================================!
!=======================================================================================================! 
subroutine sw2_init(elem,tl,ie,k,pmean)
!=======================================================================================================!
    type(element_t) :: elem(:)
    type(TimeLevel_t), intent(in):: tl
    integer, intent(in)  :: ie,k
    real (kind=real_kind):: pmean
! Local variables
    integer :: nm1,n0,np1,nstep
!=======================================================================================================!
    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep  

    elem(ie)%fcor=sw2_coreolis_init(elem(ie)%spherep)
    elem(ie)%state%v(:,:,:,k,n0)=sw1_velocity(elem(ie)%spherep(:,:),elem(ie)%Dinv)
    elem(ie)%state%ht(:,:,k)= sw2_phi(elem(ie)%spherep(:,:))/g
    elem(ie)%state%hs(:,:,k)= 0.0D0
!    elem(ie)%state%psi(:,:,k)= elem(ie)%state%ht(:,:,k)	  

end subroutine sw2_init    
!=======================================================================================================!  
subroutine sw5_init(elem,tl,ie,k,pmean)
!=======================================================================================================!
    type(element_t) :: elem(:)
    type(TimeLevel_t), intent(in):: tl
    integer, intent(in)  :: ie,k
    real (kind=real_kind):: pmean
! Local variables
    integer :: nm1,n0,np1,nstep
    real (kind=real_kind):: hts(np,np),hll(np,np)
!=======================================================================================================!
    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep  
!=======================================================================================================!
    elem(ie)%fcor=sw5_coreolis_init(elem(ie)%spherep)
    elem(ie)%state%v(:,:,:,k,n0)= sw5_velocity(elem(ie)%spherep(:,:),elem(ie)%Dinv)	   
    Call sw5_phi(elem(ie)%spherep(:,:),hts,hll)           
    elem(ie)%state%ht(:,:,k)= hts(:,:)/g 
    elem(ie)%state%hs(:,:,k)= hll(:,:)/g
!    elem(ie)%state%psi(:,:,k)= elem(ie)%state%ht(:,:,k)

end subroutine sw5_init
!=======================================================================================================!  
end module dg_tests_mod
