#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module shallow_water_mod
  ! ------------------------
  use kinds, only : real_kind
  ! ------------------------
  use physical_constants, only : omega, g, rearth, rrearth, dd_pi
  ! ------------------------
  use dimensions_mod, only : nlev, np
  ! ------------------------
  use derivative_mod, only : derivative_t, vorticity_sphere, derivative_stag_t, gradient_wk, &
                             divergence, vorticity, laplace_sphere_wk, &
                             vlaplace_sphere_wk, divergence_sphere
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
  use time_mod, only : time_at, secpday, timelevel_t, nmax
  ! ------------------------
  use element_mod, only : element_t
  ! ------------------------
  use hybrid_mod, only : hybrid_t
  ! ------------------------
  use parallel_mod, only : parallel_t
  ! ------------------------
#ifdef _REFSOLN
  use ref_state_mod, only : ref_state_read, ref_state_write
#endif
  ! ------------------------
  use common_io_mod, only: output_prefix     ! Added to support output_prefix 
  ! ------------------------
  use viscosity_mod, only: biharmonic_wk, test_ibyp, neighbor_minmax, check_edge_flux
  ! ------------------------
  use control_mod, only : nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, limiter_option, integration, test_case, sub_case, kmass
  ! ------------------------

  implicit none
  private

  ! ============================================
  ! Shallow Water Test Case Parameters 
  ! for test case 1 and 2; Cosine Bell and 
  ! Steady State Geostrophic 
  ! ============================================
  real (kind=real_kind), private, parameter  :: h0     = 1000.0D0     ! height of cosine bell in meters
  real (kind=real_kind), private, parameter  :: gh0_tc1= g*h0         ! geopotential height of cosine bell
  real (kind=real_kind), private, parameter  :: rr     = 1.0/3.0      ! radius of cosine bell
  real (kind=real_kind), public,  parameter  :: latc0  = 0.0D0        ! init latitude of cosine bell center
  real (kind=real_kind), public,  parameter  :: lonc0  = 1.5D0*DD_PI  ! init longitude of cosine bell center
  real (kind=real_kind), private, parameter  :: u0     = 2.0D0*DD_PI*rearth/(12.0*86400.0) ! adv. wind velocity
  real (kind=real_kind), private, parameter  :: alpha  = 0.0D0        ! angle of wind advection rel. to equator
!  real (kind=real_kind), private, parameter  :: alpha  = DD_PI/2.0D0       ! angle of wind advection rel. to equator
!  real (kind=real_kind), private, parameter  :: alpha  = DD_PI/4.0D0  ! angle of wind advection rel. to equator
  real (kind=real_kind), private, parameter  :: gh0_tc2 = 2.94D4      ! m^2/s^2 (p. 218 of Williamson)

  ! ================================================
  ! Shallow Water Test Case Parameters 
  ! Parameters for test case 5, Flow over Mountain
  ! ================================================

  real (kind=real_kind), private, parameter  :: h0_tc5   = 5960.0D0  ! height field
  real (kind=real_kind), parameter, private  :: gh0_tc5  = g*h0_tc5  ! test case 5 geopotential
  real (kind=real_kind), private, parameter  :: u0_tc5   = 20.0D0    ! wind velocity

  real (kind=real_kind), private, parameter  :: h_mtn    = 2000.0D0     ! height of mountain
  real (kind=real_kind), private, parameter  :: R_mtn    = DD_PI/9.0D0  ! radius of mountain
  real (kind=real_kind), public,  parameter  :: lat_mtn  = DD_PI/6.0D0  ! latitude of mountain
  real (kind=real_kind), public,  parameter  :: lon_mtn  = 1.5D0*DD_PI  ! longitude of mountain

  ! ============================================
  ! Shallow Water Test Case Parameters 
  ! Parameters for test case 6, Rossby-Haurwitz
  ! ============================================

  real (kind=real_kind), parameter, private :: h0_tc6   = 8.0D3
  real (kind=real_kind), parameter, private :: gh0_tc6  = g*h0_tc6
  real (kind=real_kind), parameter, private :: omega_rh = 7.848D-6
  real (kind=real_kind), parameter, private :: K_rh     = omega_rh
  real (kind=real_kind), parameter, private :: R_rh     = 4 

  ! ============================================
  ! Shallow Water Test Case Parameters 
  ! Parameters for test case 8, 
  ! Barotropic Instability.
  ! ============================================

  real (kind=real_kind), parameter, private :: h0_tc8    = 10000.0D0     ! mean fluid depth (meters)
  real (kind=real_kind), parameter, private :: hhat_tc8  = 120.0D0       ! perturbation depth (meters)
  real (kind=real_kind), parameter, private :: gh0_tc8   = g*h0_tc8      ! mean geopotential
  real (kind=real_kind), private, parameter :: u0_tc8    = 80.0D0        ! wind velocity
  real (kind=real_kind), private, parameter :: alpha_tc8 = 3.0D0         ! meridional scale factor in perturbation
  real (kind=real_kind), private, parameter :: beta_tc8  = 15.0D0        ! zonal scale factor in perturbation
  real (kind=real_kind), public,  parameter :: lat_tc8   = DD_PI/4.0D0   ! latitude of perturbation
  real (kind=real_kind), public,  parameter :: lon_tc8   = DD_PI/4.0D0   ! longitude of perturbation
  real (kind=real_kind), public,  parameter :: gamma_tc8 = DD_PI/18.0D0  ! zonal velocity field scale factor
  real (kind=real_kind), public,  parameter :: eps_tc8   = 1.0D-10       ! relative error of balance integral
  integer, parameter :: Ngs_tc8 = 128   ! number of gauss points in balance integral

  ! ==========================================
  ! Vortex problem: Nair & Machenhaur, 2002
  ! ==========================================
  real (kind=real_kind), private, parameter :: alpha_vortex   =DD_PI*0.25D0 ! meridional scale factor in perturbation
  real (kind=real_kind), parameter, private :: a_vortex       = rearth
  real (kind=real_kind), private, parameter :: u0_vortex      = 2*DD_PI*a_vortex/10.368D5   ! wind velocity
  real (kind=real_kind), private, parameter :: omega_s_vortex = u0_vortex/a_vortex
  real (kind=real_kind), private, parameter :: lat0_vort1_vortex=0; 
  real (kind=real_kind), private, parameter :: lon0_vort1_vortex=DD_PI*0.5D0; !Vortex center at time t=0
  logical, private, parameter               :: static_vortex = .false.
  ! ============================================
  ! Shallow Water Test Case Parameters 
  ! Parameters for test case SJ1, 
  ! Strong jet northern hemisphere.
  ! ============================================
  real (kind=real_kind), parameter, private :: h0_sj1    = 10000.0D0     ! mean fluid depth (meters)
  real (kind=real_kind), parameter, private :: hhat_sj1  = 120.0D0       ! perturbation depth (meters)
  real (kind=real_kind),  private           :: gh0_sj1   = g*h0_sj1      ! mean geopotential
  real (kind=real_kind), private, parameter :: u0_sj1    = 80.0D0        ! wind velocity
  real (kind=real_kind), private, parameter :: alpha_sj1 = 3.0D0         ! meridional scale factor in perturbation
  real (kind=real_kind), private, parameter :: beta_sj1  = 15.0D0        ! zonal scale factor in perturbation
  real (kind=real_kind), public,  parameter :: lat_sj1   = DD_PI/4.0D0   ! latitude of perturbation
  real (kind=real_kind), public,  parameter :: lon_sj1   = DD_PI/4.0D0   ! longitude of perturbation
  real (kind=real_kind), public,  parameter :: gamma_sj1 = DD_PI/18.0D0  ! zonal velocity field scale factor
  real (kind=real_kind), public,  parameter :: eps_sj1   = 1.0D-10       ! relative error of balance integral
  integer, parameter :: Ngs_sj1 = 200   ! number of gauss points in balance integral
  real (kind=real_kind), public,  parameter :: lat0_sj1  = DD_PI/7.0D0       ! used for velocity
  real (kind=real_kind), public,  parameter :: lat1_sj1  = 5.0D0*DD_PI/14.0D0! used for velocity
  real (kind=real_kind), public,  parameter :: lat2_sj1  = DD_PI/4.0D0          ! used for velocity

  ! ==========================================
  ! Swirl problem: Nair & Lauritzen, 2010
  ! ==========================================
  integer, public,  parameter 	:: kmass_swirl=4 !at which level the const field sits
  real (kind=real_kind), private, parameter :: mean_pressure=0.0D0
  real (kind=real_kind), private, parameter :: bellradius=0.50D0
!gauss bell would use different radius
  real (kind=real_kind), private, parameter :: rotangle=0.0D0 ! IN RADIANS
  real (kind=real_kind), private, parameter :: tracer_lowest=0.1D0 
  real (kind=real_kind), private, parameter :: tracer_highest=0.9D0
  logical                                   :: add_pure_rotation

  real (kind=real_kind) :: lat1, lon1, lat2, lon2, Tperiod, Kcoef

  !sub_case = 1
  real (kind=real_kind), private, parameter :: lon1_case1=DD_PI 
  real (kind=real_kind), private, parameter :: lat1_case1=DD_PI/3.0D0; 
  real (kind=real_kind), private, parameter :: lon2_case1=DD_PI 
  real (kind=real_kind), private, parameter :: lat2_case1=-DD_PI/3.0D0; 

  !sub_case = 2
  real (kind=real_kind), private, parameter :: lon1_case2=5.0D0*DD_PI/6.0D0; 
  real (kind=real_kind), private, parameter :: lat1_case2=0.0D0; 
  real (kind=real_kind), private, parameter :: lon2_case2=7.0D0*DD_PI/6.0D0; 
  real (kind=real_kind), private, parameter :: lat2_case2=0.0D0;

  !sub_case = 3
  !real (kind=real_kind), private, parameter :: lon1_case3=3.0D0*DD_PI/4.0D0 
  !real (kind=real_kind), private, parameter :: lat1_case3=0.0D0; 
  !real (kind=real_kind), private, parameter :: lon2_case3=5.0D0*DD_PI/4.0D0 
  !real (kind=real_kind), private, parameter :: lat2_case3=0.0D0;

! i put subcase 2 center data here cause this is what pl uses for paper
  real (kind=real_kind), private, parameter :: lon1_case3=5.0D0*DD_PI/6.0D0; 
  real (kind=real_kind), private, parameter :: lat1_case3=0.0D0; 
  real (kind=real_kind), private, parameter :: lon2_case3=7.0D0*DD_PI/6.0D0; 
  real (kind=real_kind), private, parameter :: lat2_case3=0.0D0;

  !subcase 4 is the same as subcase 2 but with added zonal flow, see below
 

  public  :: sweq_invariants
  public  :: tc1_init_state  ! Initialize test case 1: Cosine Bell
  public  :: tc1_init_pmean  ! Initialize pmean for test case 1
  public  :: tc1_velocity
  public  :: tc1_geopotential
  public  :: tc1_phi         ! test case 1 on unstaggered grid (nair) 
  public  :: grad_cosine_bell
  public  :: tc1_errors

  public  :: tc2_init_state  ! Initialize test case 2: Global steady state nonlinear geostrophic flow
  public  :: tc2_init_pmean
  public  :: tc2_geopotential
  private :: tc2_coreolis_init
  public  :: tc2_errors
  public  :: tc2_phi

  public  :: tc5_init_state  ! Initialize test case 5: Flow Over a Mountain
  public  :: tc5_init_pmean
  public  :: tc5_velocity
  private :: tc5_coreolis_init
  public  :: tc5_mountain
  public  :: tc5_geopotential
  public  :: tc5_errors
  public  :: tc5_invariants
  public  :: tc5_phi
#if 0
  public  :: tc5_init_diag
#endif

  public  :: tc6_init_state  ! Initialize test case 6: Rossby-Haurwitz Wave
  public  :: tc6_init_pmean 
  private :: tc6_geopotential
  private :: tc6_velocity
  public  :: tc6_errors

  public  :: tc8_init_state  ! Initialize test case 8: Barotropic Instability
  public  :: tc8_init_pmean
  private :: tc8_balance
  private :: tc8_perturbation  
  private :: tc8_integrand
  private :: tc8_velocity

  public  :: vortex_init_state
  public  :: vortex_velocity
  public  :: vortex_errors
  public  :: vortex_exact

  public  :: swirl_velocity
  public  :: swirl_init_state
  public  :: swirl_errors


  public  :: sj1_init_state  ! Initialize strong jet case 1: Galewski et al.
  public  :: sj1_init_pmean
  private :: sj1_balance
  private :: sj1_perturbation  
  private :: sj1_integrand
  private :: sj1_velocity
  private :: sj1_velocity_cubedsphere
  !public  :: sj1_invariants

  real (kind=real_kind), public,save :: Imass0,Ienergy0,Ipenst0
!  common /comintegral/Imass0,Ienergy0,Ipenst0

  real (kind=real_kind), public,save :: Ibalance_tot
!  common /comtc8/Ibalance_tot


contains


  ! ======================================================
  ! shallow water equation invariants on a non-staggered grid.
  ! ======================================================
  subroutine sweq_invariants(elem, iounit,tl,pmean,edge3,deriv,hybrid,nets,nete)
    type(element_t), intent(inout), target :: elem(:)

    integer              , intent(in) :: iounit
    type (TimeLevel_t)  , intent(in) :: tl    
    real (kind=real_kind), intent(in) :: pmean
    type (EdgeBuffer_t)               :: edge3
    type (derivative_t)               :: deriv
    type (hybrid_t)      , intent(in) :: hybrid
    integer              , intent(in) :: nets,nete

    ! Local variables

    integer :: ie,i,j,kptr

    real (kind=real_kind) :: vco(np,np,2)
    real (kind=real_kind) :: gv(np,np,2)

    real (kind=real_kind) :: E(np,np)
    real (kind=real_kind) :: mass(np,np,nlev,nets:nete)
    real (kind=real_kind) :: kenergy(np,np,nets:nete)
    real (kind=real_kind) :: penergy(np,np,nets:nete)
    real (kind=real_kind) :: penst(np,np,nets:nete)
    real (kind=real_kind) :: pv(np,np,nets:nete)
    real (kind=real_kind) :: div(np,np,nets:nete)

    real (kind=real_kind) :: diss_p(np,np,nlev,nets:nete)
    real (kind=real_kind) :: diss_v(np,np,2,nlev,nets:nete)

    real (kind=real_kind) :: v1,v2,vlatlon(np,np,2)
    real (kind=real_kind) :: h,hstar,hs

    real (kind=real_kind) :: dm,dKE,dPE,dENS,dissE

    ! dont use Imass0,Ienergy0,Ipenst0 because they are being used by tc5_invariants
    real (kind=real_kind),save :: Imass_init(nlev),Ienergy_init,Ipenst_init,Ipv_init,Idiv_init
    real (kind=real_kind),save :: time_last,Imass_last(nlev),Ikenergy_last,Ipenergy_last,Ipenst_last,Ipv_last,Idiv_last
    real (kind=real_kind) :: time,Imass(nlev),Ikenergy,Ipenergy,Ipenst,Ienergy,Ipv,Idiv
    real (kind=real_kind), dimension(:,:,:), pointer :: v
    real (kind=real_kind), dimension(:,:),   pointer :: viscosity => NULL()

    integer               :: k,n0,k1

!    print *,'calling test_ibyp'
!    call test_ibyp(elem,hybrid,nets,nete)
!    call check_edge_flux(elem,deriv,nets,nete)
!    stop

    n0 = tl%n0
    k=1
#if (! defined ELEMENT_OPENMP)
    !$OMP BARRIER
#endif
    if (tl%nstep == 0) then
       if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
          open(iounit+0,file=TRIM(output_prefix)//"sweq.mass",status="unknown",form="formatted")
          open(iounit+1,file=TRIM(output_prefix)//"sweq.ke",status="unknown",form="formatted")
          open(iounit+2,file=TRIM(output_prefix)//"sweq.pe",status="unknown",form="formatted")
          open(iounit+3,file=TRIM(output_prefix)//"sweq.penst",status="unknown",form="formatted")
          open(iounit+4,file=TRIM(output_prefix)//"sweq.pv",status="unknown",form="formatted")
          open(iounit+5,file=TRIM(output_prefix)//"sweq.div",status="unknown",form="formatted")
       endif
    endif
    if (nu > 0 .or. nu_s > 0 ) then
       ! convert to lat-lon (needed to call biharmonic() )
       do ie=nets,nete
          v  => elem(ie)%state%v(:,:,:,k,n0)
          do j=1,np
             do i=1,np
                v1     = elem(ie)%state%v(i,j,1,k,n0)   ! contra
                v2     = elem(ie)%state%v(i,j,2,k,n0)   ! contra 
                v(i,j,1)=elem(ie)%D(1,1,i,j)*v1 + elem(ie)%D(1,2,i,j)*v2   ! contra->latlon
                v(i,j,2)=elem(ie)%D(2,1,i,j)*v1 + elem(ie)%D(2,2,i,j)*v2   ! contra->latlon
             enddo
          enddo
       enddo
       if (hypervis_order==2) then
          call biharmonic_wk(elem,diss_p,diss_v,deriv,edge3,hybrid,n0,nets,nete,nu_div/nu)
       else
          do ie=nets,nete
             v  => elem(ie)%state%v(:,:,:,k,n0)
             E(:,:)=elem(ie)%state%p(:,:,k,n0) + elem(ie)%state%ps(:,:)
             diss_p(:,:,k,ie)=laplace_sphere_wk(E,deriv,elem(ie),viscosity)
             diss_v(:,:,:,k,ie)=vlaplace_sphere_wk(v,deriv,elem(ie),viscosity)
          enddo
       endif
       ! convert lat-lon -> contra variant
       do ie=nets,nete
          v  => elem(ie)%state%v(:,:,:,k,n0)
          do j=1,np
             do i=1,np
                v1=v(i,j,1)
                v2=v(i,j,2)
                v(i,j,1) = elem(ie)%Dinv(1,1,i,j)*v1 + elem(ie)%Dinv(1,2,i,j)*v2
                v(i,j,2) = elem(ie)%Dinv(2,1,i,j)*v1 + elem(ie)%Dinv(2,2,i,j)*v2
             enddo
          enddo
       enddo
    else
       diss_p=0; diss_v=0;
    endif
    
    do ie=nets,nete      
       do j=1,np
          do i=1,np
             k1=1
             v1     = elem(ie)%state%v(i,j,1,k1,n0)   ! contra
             v2     = elem(ie)%state%v(i,j,2,k1,n0)   ! contra
             vlatlon(i,j,1) =elem(ie)%D(1,1,i,j)*v1 + elem(ie)%D(1,2,i,j)*v2   ! contra->latlon
             vlatlon(i,j,2) =elem(ie)%D(2,1,i,j)*v1 + elem(ie)%D(2,2,i,j)*v2   ! contra->latlon
             
             E(i,j) = 0.5D0*(vlatlon(i,j,1)**2 + vlatlon(i,j,2)**2) 
             hstar         = (elem(ie)%state%p(i,j,1,n0) + pmean)/g
             hs            = elem(ie)%state%ps(i,j)/g
             h             = hstar + hs



             kenergy(i,j,ie)  = hstar*E(i,j) 
             penergy(i,j,ie)  = 0.5D0*g*(h**2 - hs**2)
             diss_p(i,j,k1,ie) = nu_s*(E(i,j)+g*h)*diss_p(i,j,k1,ie)/g + &
                     nu*hstar*(vlatlon(i,j,1)*diss_v(i,j,1,k1,ie)+ &
                               vlatlon(i,j,2)*diss_v(i,j,2,k1,ie))
             ! weak dissipation operator has mass matrix "built in"
             ! integral will multiply by mass matrix again, so lets remove it:
             ! note: this is not the same as multiplying by rspheremp()! 
             diss_p(i,j,k1,ie)=diss_p(i,j,k1,ie)/elem(ie)%spheremp(i,j)  

!hs etc. is used only from the first level
!but mass conservation output will be done for all levels 
	     do k=1,nlev
		mass(i,j,k,ie)    = (elem(ie)%state%p(i,j,k,n0) + pmean)/g
	     enddo

          end do
       end do
       ! global integral < f*g > is computed correctly if either f or g is continious
       ! enstropy = < PV**2 >, so we need to make the PV continious before computing global integral
       div(:,:,ie) = divergence_sphere(vlatlon,deriv,elem(ie)) ! latlon vector -> scalar 
       penst(:,:,ie) = vorticity_sphere(vlatlon,deriv,elem(ie)) ! latlon vector -> scalar 
       penst(:,:,ie) = penst(:,:,ie)*elem(ie)%spheremp(:,:)
       kptr=0
       call edgeVpack(edge3, penst(1,1,ie), 1, kptr,elem(ie)%desc)
    end do
    call bndry_exchangeV(hybrid,edge3)
    do ie=nets,nete      
       kptr=0
       call edgeVunpack(edge3, penst(1,1,ie), 1, kptr, elem(ie)%desc)
       penst(:,:,ie) = penst(:,:,ie)*elem(ie)%rspheremp(:,:)
       do j=1,np
          do i=1,np
             pv(i,j,ie) = penst(i,j,ie) + elem(ie)%fcor(i,j)
             hstar         = (elem(ie)%state%p(i,j,1,n0) + pmean)/g
             penst(i,j,ie) = (0.5D0*(penst(i,j,ie) + elem(ie)%fcor(i,j))**2) / hstar
          end do
       end do
    end do


    do k=1,nlev
	Imass(k)   = global_integral(elem,mass(:,:,k,nets:nete),hybrid,np,nets,nete)
    enddo
    Ikenergy = global_integral(elem,kenergy(:,:,nets:nete),hybrid,np,nets,nete)
    Ipenergy = global_integral(elem,penergy(:,:,nets:nete),hybrid,np,nets,nete)
    Ienergy = Ikenergy + Ipenergy
    Ipenst  = global_integral(elem,penst(:,:,nets:nete),hybrid,np,nets,nete)
    Ipv  = global_integral(elem,pv(:,:,nets:nete),hybrid,np,nets,nete)
    Idiv  = global_integral(elem,div(:,:,nets:nete),hybrid,np,nets,nete)
    dissE  =  global_integral(elem,diss_p(:,:,k,nets:nete),hybrid,np,nets,nete)
    if (hypervis_order==2) dissE=-dissE  ! biharmonic is applied with a minus sign


    time   = Time_at(tl%nstep)

#if (! defined ELEMENT_OPENMP)
    !$OMP BARRIER
#endif
    if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
       if (tl%nstep==0) then
          Imass_init   = Imass
          Ienergy_init = Ienergy
          Ipenst_init  = Ipenst
          Ipv_init  = Ipv
          Idiv_init  = Idiv
       end if

       write(iounit+0,*)time/secpday,Imass(1)
       write(iounit+1,*)time/secpday,Ikenergy
       write(iounit+2,*)time/secpday,Ipenergy
       write(iounit+3,*)time/secpday,Ipenst
       write(iounit+4,*)time/secpday,Ipv
       write(iounit+5,*)time/secpday,Idiv

       !write(6,*)time/secpday,"mass          =",(Imass-Imass_init)/Imass_init
       !write(6,*)time/secpday,"total energy  =",(Ienergy-Ienergy_init)/Ienergy_init
       !write(6,*)time/secpday,"pot enstrophy =",(Ipenst-Ipenst_init)/Ipenst_init
       if (tl%nstep>0) then

          dKE = (Ikenergy-Ikenergy_last)/(time-time_last)/(Ikenergy+Ipenergy)
          dPE = (Ipenergy-Ipenergy_last)/(time-time_last)/(Ikenergy+Ipenergy)
          dENS = (Ipenst-Ipenst_last)/(time-time_last)/Ipenst
	  do k=1,nlev
	    dm = (Imass(k)-Imass_last(k))/(time-time_last)/Imass(k)
	    write(6,'(a,i3,a,e11.4,a,e13.6,a,e13.6)') "level =",k,"  M-M0/M0      =",(Imass(k)-Imass_init(k))/Imass_init(k), " dM/dt   /M  = ",dm
	  enddo
          if (test_case(1:5) /= "swtc1") then
             write(6,'(a,e11.4,a,e13.6,a,e13.6)') "ENS-ENS0/ENS0=",(Ipenst-Ipenst_init)/Ipenst_init, " dENS/dt /ENS= ",dENS
             write(6,'(a,e11.4,a,e13.6,a,e11.4,a,e11.4,a)') "E-E0/E0      =",(Ienergy-Ienergy_init)/Ienergy_init,&
                  " dE/dt   /E  = ",dKE+dPE,"  (",dKE,",",dPE,")" 
             write(6,'(a,e13.6,f10.3)') "                       Dissipation/E  = ",dissE/(Ikenergy+Ipenergy)
          endif
          dm = (Ipv-Ipv_last)/(time-time_last)
          write(6,'(a,e11.4,a,e13.6,a,e13.6)') "PV-PV0       =",(Ipv-Ipv_init), " dPV/dt      = ",dm
          dm = (Idiv-Idiv_last)/(time-time_last)
          write(6,'(a,e11.4,a,e13.6,a,e13.6)') "DV-DV0       =",(Idiv-Idiv_init), " dDV/dt      = ",dm
       endif

    end if

#if (! defined ELEMENT_OPENMP)
    !$OMP BARRIER
#endif
    Imass_last=Imass
    Ipenergy_last=Ipenergy
    Ikenergy_last=Ikenergy
    Ipenst_last=Ipenst
    Ipv_last=Ipv
    Idiv_last=Idiv
    time_last = time

  end subroutine sweq_invariants


  ! ===========================================
  ! tc1_init_pmean
  ! ===========================================
  function tc1_init_pmean() result(pmean)
    implicit none
    real (kind=real_kind) :: pmean 

    pmean = gh0_tc1
    ! pmean is needed for semi-implicit
    ! for explicit, it does not matter and it gives negative values
    ! with monotone limiter of order -1e-11 for some reason
    if (integration == "explicit" .and. limiter_option /= 0) pmean = 0
    if (integration == "runge_kutta" .and. limiter_option /= 0) pmean = 0
    if (integration == "full_imp" .and. limiter_option /= 0) pmean = 0
  end function tc1_init_pmean
  ! ===========================================
  ! tc2_init_pmean
  ! ===========================================
  function tc2_init_pmean() result(pmean)
    implicit none
    real (kind=real_kind) :: pmean 

    pmean = gh0_tc2
  end function tc2_init_pmean
  ! ===========================================
  ! tc5_init_pmean
  ! ===========================================
  function tc5_init_pmean() result(pmean)
    implicit none
    real (kind=real_kind) :: pmean 

    pmean = gh0_tc5
  end function tc5_init_pmean
  ! ===========================================
  ! tc6_init_pmean
  ! ===========================================
  function tc6_init_pmean() result(pmean)
    implicit none
    real (kind=real_kind) :: pmean 

    pmean = gh0_tc6
  end function tc6_init_pmean
  ! ===========================================
  ! tc8_init_pmean
  ! ===========================================
  function tc8_init_pmean() result(pmean)
    implicit none
    real (kind=real_kind) :: pmean 

    pmean = gh0_tc8
  end function tc8_init_pmean

  ! ===========================================
  ! sj1_init_pmean
  ! ===========================================
  function sj1_init_pmean() result(pmean)
    implicit none
    real (kind=real_kind) :: pmean 

    pmean = gh0_sj1
  end function sj1_init_pmean

  ! ===========================================
  !
  ! tc1_init_state:
  !
  ! Initialize for shallow water testcase 1
  !
  ! ===========================================
  subroutine tc1_init_state(elem, nets,nete,pmean)
    type(element_t), intent(inout) :: elem(:)
    integer, intent(in) :: nets
    integer, intent(in) :: nete
    real (kind=real_kind) :: pmean

! ------ where is intent(in) for pmean??? OG

    ! Local variables

    integer ie,k,case
    integer :: nm1 
    integer :: n0 
    integer :: np1

    nm1= 1
    n0 = 2
    np1= 3

    pmean = tc1_init_pmean()       ! set pmean
    do ie=nets,nete
       do k=1,nlev
          case = sub_case
          if (sub_case==1 .and. k>1) case=k  ! run one sub_case in each level
          elem(ie)%state%p(:,:,k,n0)=tc1_geopotential(elem(ie)%spherep(:,:),pmean,latc0,lonc0,case)
!          elem(ie)%state%p(:,:,k,n0)=0.0D0
          elem(ie)%state%p(:,:,k,nm1)=elem(ie)%state%p(:,:,k,n0)
          elem(ie)%state%p(:,:,k,np1)=0.0D0
          
          elem(ie)%state%v(:,:,:,k,n0)=tc1_velocity(elem(ie)%spherep(:,:),elem(ie)%Dinv)
          elem(ie)%state%v(:,:,:,k,nm1)=elem(ie)%state%v(:,:,:,k,n0)
          elem(ie)%state%v(:,:,:,k,np1)=0.0D0
       end do
       elem(ie)%state%ps(:,:)=0.0D0
       elem(ie)%state%gradps(:,:,:)=0.0D0
    end do

  end subroutine tc1_init_state

  ! ===========================================
  !
  ! tc1_velocity:
  !
  ! Reset velocities to initial values 
  ! for a single layer of one element for 
  ! test case 1, during the timestep.
  !
  ! ===========================================

  function tc1_velocity(sphere,D) result(v)

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

    csalpha = COS(alpha)
    snalpha = SIN(alpha)

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

             ! =====================================================
             ! map sphere velocities onto the contravariant cube velocities
             ! using the D^-T mapping matrix (see Loft notes for details)
             ! =====================================================

             v(i,j,1)= V1*D(1,1,i,j) + V2*D(1,2,i,j)
             v(i,j,2)= V1*D(2,1,i,j) + V2*D(2,2,i,j)
          end do
       end do
    end do

  end function tc1_velocity






!
!******************************************************************************
!
! Subroutines for rotation in most general settings, from 3d code
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
      real(kind=real_kind) :: pxreg,pyreg,&
                  pxrot,pyrot,&
                  pxcen,pycen
!
!-----------------------------------------------------------------------
!
      real(kind=real_kind) zsycen,zcycen,zxmxc,zsxmxc,zcxmxc,zsyreg,zcyreg, &
               zsyrot,zcyrot,zcxrot,zsxrot,zpi,zpih
      integer jy,jx

      zpih = DD_PI*0.5d0
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
         zsyrot = max(zsyrot,-1.0D0)
         zsyrot = min(zsyrot,+1.0D0)
         !
         pyrot = ASIN(zsyrot)
         !
         zcyrot = COS(pyrot)
         zcxrot = (zcycen*zcyreg*zcxmxc +zsycen*zsyreg)/zcyrot
         zcxrot = max(zcxrot,-1.0D0)
         zcxrot = min(zcxrot,+1.0D0)
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
         zsyreg = max(zsyreg,-1.0D0)
         zsyreg = min(zsyreg,+1.0D0)
         !
         pyreg = ASIN(zsyreg)
         !
         zcyreg = COS(pyreg)
         zcxmxc = (zcycen*zcyrot*zcxrot -&
              zsycen*zsyrot)/zcyreg
         zcxmxc = max(zcxmxc,-1.0D0)
         zcxmxc = min(zcxmxc,+1.0D0)
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
      real(kind=real_kind) puarg,pvarg,    &
               pures,pvres,    &
               pa,   pb,       &
               pc,   pd,       &
               pxreg,pyreg,    &
               pxrot,pyrot
      real(kind=real_kind) pxcen,pycen
!
!-----------------------------------------------------------------------
!
      integer jy,jx
      real(kind=real_kind) zpih,zsyc,zcyc,zsxreg,zcxreg,zsyreg,zcyreg,zxmxc,&
               zsxmxc,zcxmxc,zsxrot,zcxrot,zsyrot,zcyrot
!
!-----------------------------------------------------------------------
!
      IF (kcall.eq.1) then
         zpih = DD_PI*0.5d0
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
         zpih = DD_PI*0.5d0
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  ! ===========================================
  !
  ! tc1_geopotential:
  !
  ! Initialize cosine bell centered at latc, lonc
  !
  ! ===========================================

  function tc1_geopotential(sphere, pmean, latc, lonc,k) result(p)

    type (spherical_polar_t), intent(in) :: sphere(np,np) ! spherical coordinates of element
    real (kind=real_kind), intent(in)    :: pmean       ! mean geopotential
    real (kind=real_kind), intent(in)    :: latc        ! latitude of center of cosine bell
    real (kind=real_kind), intent(in)    :: lonc        ! longitude of center of cosine bell


    real (kind=real_kind)                :: p(np,np)    ! geopotential

    real (kind=real_kind) :: snlat,cslat,lon,lat
    real (kind=real_kind) :: h
    real (kind=real_kind) :: r,A
    real (kind=real_kind) :: h2,lon_min,lon_max
    integer i,j,k

    h2 = 0.25D0*rr
    lon_min = lonc - h2
    lon_max = lonc + h2


    do j=1,np
       do i=1,np
          snlat=SIN(sphere(i,j)%lat)
          cslat=COS(sphere(i,j)%lat)
          lon  =sphere(i,j)%lon
          lat  =sphere(i,j)%lat
          A    =SIN(latc)*snlat + COS(latc)*cslat*COS(lon - lonc)   ! face 4 this == - cslat*snlon
          r=ACOS(A)
          if (r<rr) then
             h = (h0/2.0d0)*(1.0D0 + COS(DD_PI*r/rr)) 
          else
             h = 0.0D0
          end if
          p(i,j)=(g*h) - pmean

          if (k==2) then
             ! gaussian, for testing with a smooth function
             ! at r=rr (radius of cosine bell) h=370
             p(i,j)=g*1000*exp(-(r/rr)**2)  - pmean
          endif
          if (k==3) then
             ! cylinder
             if (r<rr .and. (lon-lonc)<0 ) then
                p(i,j)=1000*g-pmean
             else
                p(i,j)=0 - pmean
             endif
          endif
          if (k==4) then
             ! slotted cylinder
             if (r<rr) then
                p(i,j)=1.0
             else
                p(i,j)=0
             endif

             if((lon > lon_min).and.(lon<lon_max) )then
                ! we are in the slot
                if((lat > -DD_PI).and.(lat<latc))then
                   p(i,j)=0.0D0
                else if((lat < DD_PI).and.(lat > 2.0D0*DD_PI + latc - rr))then
                   p(i,j)=0.0D0
                end if
             end if
             
             p(i,j)=g*p(i,j) - pmean

          end if
       end do
    end do

  end function tc1_geopotential

  ! ===========================================
  !
  ! tc1_phi:    (nair)   for unstaggered grids
  !
  ! Initialize cosine bell centered at latc, lonc  on (np,np) points
  !
  ! ===========================================

  function tc1_phi(sphere, pmean, latc, lonc) result(ppsi)

    type (spherical_polar_t), intent(in) :: sphere(np,np) ! spherical coordinates of element
    real (kind=real_kind), intent(in)    :: pmean         ! mean geopotential
    real (kind=real_kind), intent(in)    :: latc          ! latitude of center of cosine bell
    real (kind=real_kind), intent(in)    :: lonc          ! longitude of center of cosine bell
    real (kind=real_kind)                :: ppsi(np,np)   ! geopotential

    real (kind=real_kind) :: snlat,cslat,lon, cslon, coef
    real (kind=real_kind) :: h, gh00, trm
    real (kind=real_kind) :: r,A, csalpha, snalpha
    integer i,j

    csalpha = COS(alpha)
    snalpha = SIN(alpha)

    do j=1,np
       do i=1,np
          snlat=SIN(sphere(i,j)%lat)
          cslat=COS(sphere(i,j)%lat)
          cslon=COS(sphere(i,j)%lon)
          lon  =sphere(i,j)%lon

          A = SIN(latc)*snlat + COS(latc)*cslat*COS(lon - lonc)   ! face 4 this == - cslat*snlon

          r=ACOS(A)
          if (r<rr) then
             h = (h0/2.0d0)*(1.0D0 + COS(DD_PI*r/rr))
          else
             h = 0.0D0
          end if
          ppsi(i,j)= h                                 !for tc1

       end do
    end do


  end function tc1_phi


  ! ===========================================
  ! tc1_errors:
  !
  ! Compute Shallow Water Test Case 1 errors
  ! 
  ! ===========================================

  subroutine tc1_errors(elem, iounit, tl, pmean, hybrid, nets, nete)
    type(element_t), intent(inout) :: elem(:)

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

    integer ie,k

#if (! defined ELEMENT_OPENMP)
    !$OMP BARRIER
#endif

    if (tl%nstep == 0) then
       if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
          if(sub_case == 1)then

             open(iounit+0,file=TRIM(output_prefix)//"swtc1.cosi.l1.errors",form="formatted")
             open(iounit+1,file=TRIM(output_prefix)//"swtc1.cosi.l2.errors",form="formatted")
             open(iounit+2,file=TRIM(output_prefix)//"swtc1.cosi.linf.errors",form="formatted")

          else if(sub_case == 2)then

             open(iounit+0,file=TRIM(output_prefix)//"swtc1.gaus.l1.errors",form="formatted")
             open(iounit+1,file=TRIM(output_prefix)//"swtc1.gaus.l2.errors",form="formatted")
             open(iounit+2,file=TRIM(output_prefix)//"swtc1.gaus.linf.errors",form="formatted")

          else if(sub_case== 3)then

             open(iounit+0,file=TRIM(output_prefix)//"swtc1.cyli.l1.errors",form="formatted")
             open(iounit+1,file=TRIM(output_prefix)//"swtc1.cyli.l2.errors",form="formatted")
             open(iounit+2,file=TRIM(output_prefix)//"swtc1.cyli.linf.errors",form="formatted")

          else if(sub_case == 4)then

             open(iounit+0,file=TRIM(output_prefix)//"swtc1.slcy.l1.errors",form="formatted")
             open(iounit+1,file=TRIM(output_prefix)//"swtc1.slcy.l2.errors",form="formatted")
             open(iounit+2,file=TRIM(output_prefix)//"swtc1.slcy.linf.errors",form="formatted")

          end if
       end if
    end if

    ! Compute the expected lat,lon of the cosine bell:
    ! solid body rotations by Z axis and by X axis ..
    time_tmp   = Time_at(tl%nstep)
    lamdot = u0/rearth
    lon  = MODULO(lamdot*time_tmp,2.0D0*real(DD_PI,kind=real_kind))
    snlatc = SIN(lon)*SIN(alpha)
    cslatc = SQRT(1.0D0 - snlatc*snlatc)

    cslonc = (COS(lon)*COS(lonc0) - SIN(lon)*SIN(lonc0)*COS(alpha))/cslatc
    snlonc = (COS(lon)*SIN(lonc0) + SIN(lon)*COS(lonc0)*SIN(alpha))/cslatc

    latc = ASIN(SIN(lon)*SIN(alpha))
    lonc = ATAN2(snlonc,cslonc)

    if (lonc < 0.0D0) lonc = lonc + 2.0D0*DD_PI

    do k=1,nlev
       do ie=nets,nete
          pt(:,:,ie)=tc1_geopotential(elem(ie)%spherep(:,:),pmean,latc,lonc,sub_case)
          pt(:,:,ie)=pt(:,:,ie) + pmean
          p(:,:,ie) =elem(ie)%state%p(:,:,k,tl%n0) + pmean
       end do
       
       l1   = l1_snorm(elem, p(:,:,nets:nete),  pt(:,:,nets:nete),hybrid,np,nets,nete)
       l2   = l2_snorm(elem, p(:,:,nets:nete),  pt(:,:,nets:nete),hybrid,np,nets,nete)
       linf = linf_snorm(p(:,:,nets:nete),pt(:,:,nets:nete),hybrid,np,nets,nete)
       
       if (k==1) then
       if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
          write(iounit+0,30)time_tmp/secpday,l1
          write(iounit+1,30)time_tmp/secpday,l2
          write(iounit+2,30)time_tmp/secpday,linf
       end if
       end if

       if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
          write(*,'(a,i2,f6.2,a,3e15.7)') 'k=',k,time_tmp/secpday,' days  l1,l2,linf=',&
               l1,l2,linf
       end if
    enddo


#if (! defined ELEMENT_OPENMP)
    !$OMP BARRIER
#endif

10  format("      time           l1              l2              linf           ")
20  format("====================================================================")
30  format(f11.6,4x,e13.6)

  end subroutine tc1_errors


  ! ===========================================
  !
  ! grad_cosine_bell:
  !
  ! Initialize analytic gradient of 
  ! cosine bell centered at latc, lonc.
  !
  ! ===========================================

  function grad_cosine_bell(sphere,D,latc,lonc,npts) result(grad)

    integer, intent(in)                  :: npts
    type (spherical_polar_t), intent(in) :: sphere(npts,npts)
    real (kind=real_kind), intent(in)    :: D(2,2,npts,npts)
    real (kind=real_kind), intent(in)    :: latc,lonc
    real (kind=real_kind)                :: grad(npts,npts,2)

    real (kind=real_kind) :: snlat,cslat,snlon,cslon
    real (kind=real_kind) :: coef,A,r,grad1,grad2
    real (kind=real_kind) :: dA1,dA2,lon

    integer i,j

    do j=1,npts
       do i=1,npts

          snlat = SIN(sphere(i,j)%lat)
          cslat = COS(sphere(i,j)%lat)
          snlon = SIN(sphere(i,j)%lon)
          cslon = COS(sphere(i,j)%lon)
          lon  =sphere(i,j)%lon

          coef  = (g*h0*DD_PI/(2.0D0*rearth*rr))
          A     = SIN(latc)*snlat + COS(latc)*cslat*COS(lon - lonc)  ! face 4 this == - cslat*snlon

          r=ACOS(A)

          dA1 = COS(latc)*(-SIN(lon - lonc))                         ! face 4 this == -cslon
          dA2 = SIN(latc)*cslat - COS(latc)*snlat*COS(lon - lonc)    ! face 4 this == snlat*snlon

          if (r<rr) then
             if (ABS(ABS(A)-1.0D0) < 1.e-10) then
                grad(i,j,1)= 0.0D0
                grad(i,j,2)= 0.0D0
             else
                grad(i,j,1)= coef*SIN(DD_PI*r/rr)*(1.0/SQRT(1.0D0-A**2))*dA1
                grad(i,j,2)= coef*SIN(DD_PI*r/rr)*(1.0/SQRT(1.0D0-A**2))*dA2
             end if
          else
             grad(i,j,1) = 0.0D0
             grad(i,j,2) = 0.0D0
          end if

          grad1=D(1,1,i,j)*grad(i,j,1) + D(1,2,i,j)*grad(i,j,2)     ! assumes input is DV = v_i 
          grad2=D(2,1,i,j)*grad(i,j,1) + D(2,2,i,j)*grad(i,j,2)     !  "               "     "

          grad(i,j,1)=grad1
          grad(i,j,2)=grad2

       end do
    end do

  end function grad_cosine_bell

  ! ===========================================
  !
  ! tc2_init_state:
  !
  ! Initialize test case 2: Global steady state 
  ! nonlinear geostrophic flow
  !
  ! ===========================================

  subroutine tc2_init_state(elem,nets,nete,pmean)
    type(element_t), intent(inout) :: elem(:)

    integer, intent(in) :: nets
    integer, intent(in) :: nete
    real (kind=real_kind) :: pmean

    ! Local variables

    integer :: ie,k
    integer :: nm1 
    integer :: n0 
    integer :: np1

    nm1= 1
    n0 = 2
    np1= 3

    pmean = tc2_init_pmean()     ! (m^2/s^2) set pmean value specified in Williamson, et. al.
    ! p 218, for test case 2.

    do ie=nets,nete
       elem(ie)%fcor=tc2_coreolis_init(elem(ie)%spherep)
       do k=1,nlev
          elem(ie)%state%p(:,:,k,n0)=tc2_geopotential(elem(ie)%spherep(:,:))
          elem(ie)%state%p(:,:,k,nm1)=elem(ie)%state%p(:,:,k,n0)
          elem(ie)%state%p(:,:,k,np1)=0.0D0

          elem(ie)%state%v(:,:,:,k,n0)=tc1_velocity(elem(ie)%spherep(:,:),elem(ie)%Dinv)  ! tc2 vel same as tc1
          elem(ie)%state%v(:,:,:,k,nm1)=elem(ie)%state%v(:,:,:,k,n0)
          elem(ie)%state%v(:,:,:,k,np1)=0.0D0
       end do
       elem(ie)%state%ps(:,:)=0.0D0
       elem(ie)%state%gradps(:,:,:)=0.0D0
    end do

  end subroutine tc2_init_state

  ! ===========================================
  !
  ! tc2_geopotential:
  !
  ! Set geopotential to initial values 
  ! specified on page 218 of Williamson, et al.
  ! eq. 95
  !
  ! ===========================================

  function tc2_geopotential(sphere) result(p)

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

    csalpha = COS(alpha)
    snalpha = SIN(alpha)

    coef = rearth*omega*u0 + (u0**2)/2.0D0

    do j=1,np
       do i=1,np
          snlat = SIN(sphere(i,j)%lat)
          cslat = COS(sphere(i,j)%lat)
          cslon = COS(sphere(i,j)%lon)
          p(i,j)= -coef*( -cslon*cslat*snalpha + snlat*csalpha )**2
       end do
    end do

  end function tc2_geopotential

  ! ===========================================

  function tc2_phi(sphere) result(p)

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

    csalpha = COS(alpha)
    snalpha = SIN(alpha)
    coef = rearth*omega*u0 + (u0**2)/2.0D0

    do j=1,np
       do i=1,np
          snlat = SIN(sphere(i,j)%lat)
          cslat = COS(sphere(i,j)%lat)
          cslon = COS(sphere(i,j)%lon)
          !p(i,j)=  gh0_tc2 -coef*( -cslon*cslat*snalpha + snlat*csalpha )**2
          p(i,j)=  2.94D4  -coef*( -cslon*cslat*snalpha + snlat*csalpha )**2
       end do
    end do

  end function tc2_phi

  ! ========================================
  ! tc2_coreolis_init:
  !
  ! Initialize coreolis term for test case 2
  !
  ! ========================================

  function tc2_coreolis_init(sphere) result(fcor)

    type (spherical_polar_t), intent(in) :: sphere(np,np)
    real (kind=real_kind) :: fcor(np,np)

    ! Local variables

    real (kind=real_kind) :: cslon
    real (kind=real_kind) :: snlat
    real (kind=real_kind) :: cslat
    real (kind=real_kind) :: csalpha
    real (kind=real_kind) :: snalpha
    integer                  :: i,j

    csalpha = COS(alpha)
    snalpha = SIN(alpha)

    do j=1,np
       do i=1,np
          snlat = SIN(sphere(i,j)%lat)
          cslat = COS(sphere(i,j)%lat)
          cslon = COS(sphere(i,j)%lon)

          fcor(i,j) = 2.0D0*omega*(-cslon*cslat*snalpha + snlat*csalpha)
       end do
    end do

  end function tc2_coreolis_init


  ! ===========================================
  ! tc2_errors:
  !
  ! Compute Shallow Water Test Case 1 errors
  ! 
  ! ===========================================

  subroutine tc2_errors(elem,iounit, tl, pmean, hybrid, nets, nete)

    type(element_t), intent(inout) :: elem(:)
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

    time_tmp   = Time_at(tl%nstep)

#if (! defined ELEMENT_OPENMP)
    !$OMP BARRIER
#endif
    if ((tl%nstep == 1).or.(tl%nstep == 0)) then
       if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
          open(iounit+0,file=TRIM(output_prefix)//"swtc2.l1.errors",form="formatted")
          open(iounit+1,file=TRIM(output_prefix)//"swtc2.l2.errors",form="formatted")
          open(iounit+2,file=TRIM(output_prefix)//"swtc2.linf.errors",form="formatted")
       end if
    end if

    do ie=nets,nete
       pt(:,:,ie)=tc2_geopotential(elem(ie)%spherep(:,:))
       pt(:,:,ie)=pt(:,:,ie) + pmean
       p(:,:,ie) =elem(ie)%state%p(:,:,1,tl%n0) + pmean
    end do

    npts = np
    l1   = l1_snorm(elem, p(:,:,nets:nete),  pt(:,:,nets:nete),hybrid,npts,nets,nete)
    l2   = l2_snorm(elem, p(:,:,nets:nete),  pt(:,:,nets:nete),hybrid,npts,nets,nete)
    linf = linf_snorm(p(:,:,nets:nete),pt(:,:,nets:nete),hybrid,npts,nets,nete)

    if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
       write(iounit+0,30)time_tmp/secpday,l1
       write(iounit+1,30)time_tmp/secpday,l2
       write(iounit+2,30)time_tmp/secpday,linf
    end if
#if (! defined ELEMENT_OPENMP)
    !$OMP BARRIER
#endif

10  format("      time           l1              l2              linf           ")
20  format("====================================================================")
30  format(f11.6,4x,e13.6)

  end subroutine tc2_errors

  ! ===========================================
  !
  ! tc5_init_state:
  !
  ! Initialize test case 5: Flow over a mountain
  !
  ! ===========================================

  subroutine tc5_init_state(elem,nets,nete,pmean,deriv)
    type(element_t), intent(inout) :: elem(:)

    integer, intent(in) :: nets
    integer, intent(in) :: nete
    real (kind=real_kind) :: pmean
    type (derivative_t)   :: deriv

    ! Local variables

    integer ie,k
    integer :: nm1 
    integer :: n0 
    integer :: np1

    real (kind=real_kind) :: pmean_adjust

    nm1= 1
    n0 = 2
    np1= 3

    pmean = tc5_init_pmean()     ! (m^2/s^2) set pmean value specified in Williamson, et. al.
    ! p 218, for test case 2.

    ! default setup gives <hstar + pmean/g>=5619    but pmean = 5960
    pmean_adjust = 0  ! 340.284*g 

    do ie=nets,nete
       elem(ie)%fcor=tc2_coreolis_init(elem(ie)%spherep)
       elem(ie)%state%ps(:,:)=tc5_mountain(elem(ie)%spherep(:,:))

#ifdef _WK_GRAD
       elem(ie)%state%gradps(:,:,:)=gradient_wk(elem(ie)%state%ps,deriv)*rrearth
#else
       elem(ie)%state%gradps(:,:,:)=gradient(elem(ie)%state%ps,deriv)*rrearth
#endif
       
       do k=1,nlev
          elem(ie)%state%p(:,:,k,n0)=tc5_geopotential(elem(ie)%spherep(:,:)) + pmean_adjust
          elem(ie)%state%p(:,:,k,nm1)=elem(ie)%state%p(:,:,k,n0)
          elem(ie)%state%p(:,:,k,np1)=0.0D0

          elem(ie)%state%v(:,:,:,k,n0)=tc5_velocity(elem(ie)%spherep(:,:),elem(ie)%Dinv)  
          elem(ie)%state%v(:,:,:,k,nm1)=elem(ie)%state%v(:,:,:,k,n0)
          elem(ie)%state%v(:,:,:,k,np1)=0.0D0
       end do
    end do

    pmean = pmean - pmean_adjust


  end subroutine tc5_init_state

  ! ===========================================
  !
  ! tc5_velocity:
  !
  ! Reset velocities to initial values 
  ! for a single layer of one element for 
  ! test case 5, during the timestep.
  !
  ! ===========================================

  function tc5_velocity(sphere,D) result(v)

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

    csalpha = COS(alpha)
    snalpha = SIN(alpha)

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

          V1 =   u0_tc5*(cslat(i,j)*csalpha + snlat(i,j)*cslon(i,j)*snalpha)
          V2 =  -u0_tc5*(snlon(i,j)*snalpha)

          ! =====================================================
          ! map sphere velocities onto the contravariant cube velocities
          ! using the D mapping matrix (see Loft notes for details)
          ! =====================================================

          v(i,j,1)= V1*D(1,1,i,j) + V2*D(1,2,i,j)
          v(i,j,2)= V1*D(2,1,i,j) + V2*D(2,2,i,j)
       end do
    end do

  end function tc5_velocity

  ! ========================================
  ! tc5_coreolis_init:
  !
  ! Initialize coreolis term for test case 5
  !
  ! ========================================

  function tc5_coreolis_init(sphere) result(fcor)

    type (spherical_polar_t), intent(in) :: sphere(np,np)
    real (kind=real_kind) :: fcor(np,np)

    ! Local variables

    real (kind=real_kind) :: cslon
    real (kind=real_kind) :: snlat
    real (kind=real_kind) :: cslat
    real (kind=real_kind) :: csalpha
    real (kind=real_kind) :: snalpha
    integer               :: i,j

    csalpha = COS(alpha)
    snalpha = SIN(alpha)

    do j=1,np
       do i=1,np
          snlat = SIN(sphere(i,j)%lat)
          cslat = COS(sphere(i,j)%lat)
          cslon = COS(sphere(i,j)%lon)

          fcor(i,j) = 2.0D0*omega*(-cslon*cslat*snalpha + snlat*csalpha)
       end do
    end do

  end function tc5_coreolis_init
  ! ===========================================
  !
  ! tc5_mountain:
  !
  ! Set surface geopotential to initial values 
  ! specified on page 221 of Williamson, et al.
  ! eq. 134
  !
  ! ===========================================

  function tc5_mountain(sphere) result(ps)

    type (spherical_polar_t), intent(in) :: sphere(np,np)
    real (kind=real_kind)                :: ps(np,np)

    ! Local variables

    real (kind=real_kind) :: lon
    real (kind=real_kind) :: lat
    real (kind=real_kind) :: rsq

    integer i,j

    do j=1,np
       do i=1,np
          lat = sphere(i,j)%lat
          lon = sphere(i,j)%lon
          rsq=MIN((lat-lat_mtn)**2 + (lon-lon_mtn)**2,R_mtn**2)
          ps(i,j)=g*h_mtn*(1.0D0 - SQRT(rsq)/R_mtn)
       end do
    end do

  end function tc5_mountain

  !==================================================

  subroutine tc5_phi(sphere,p,ps)

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

    csalpha = COS(alpha)
    snalpha = SIN(alpha)

    coef = rearth*omega*u0_tc5 + (u0_tc5**2)/2.0D0

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
          p(i,j)= gh0_tc5 - coef*( -cslon*cslat*snalpha + snlat*csalpha )**2  
       end do
    end do

  end subroutine tc5_phi


  ! ===========================================
  !
  ! tc5_geopotential:
  !
  ! Set geopotential to initial values 
  ! specified on page 218 of Williamson, et al.
  ! eq. 95, with h0 = 5960 m and u0 = 20.0 m/s
  !
  ! ===========================================

  function tc5_geopotential(sphere) result(p)

    type (spherical_polar_t), intent(in) :: sphere(np,np)
    real (kind=real_kind)                :: p(np,np)

    ! Local variables

    real (kind=real_kind) :: ps(np,np)
    real (kind=real_kind) :: cslon
    real (kind=real_kind) :: snlat
    real (kind=real_kind) :: cslat
    real (kind=real_kind) :: csalpha
    real (kind=real_kind) :: snalpha
    real (kind=real_kind) :: coef

    integer i,j,k

    csalpha = COS(alpha)
    snalpha = SIN(alpha)

    coef = rearth*omega*u0_tc5 + (u0_tc5**2)/2.0D0

    ps(:,:)=tc5_mountain(sphere(:,:))
    do j=1,np
       do i=1,np
          snlat = SIN(sphere(i,j)%lat)
          cslat = COS(sphere(i,j)%lat)
          cslon = COS(sphere(i,j)%lon)
          p(i,j)= -coef*( -cslon*cslat*snalpha + snlat*csalpha )**2  - ps(i,j)
       end do
    end do

  end function tc5_geopotential

  ! ===========================================
  ! tc5_errors:
  !
  ! Compute Shallow Water Test Case 5 errors
  ! ===========================================

  subroutine tc5_errors(elem, iounit, tl, pmean, fstub, simday, hybrid, nets, nete, par)

    type(element_t), intent(in) :: elem(:)
    integer              , intent(in) :: iounit
    type (TimeLevel_t)   , intent(in) :: tl         ! model time struct
    real (kind=real_kind), intent(in) :: pmean      ! mean geopotential
    character(len=*)     , intent(in) :: fstub      ! file path stub
    integer              , intent(in) :: simday     ! current day of simulation
    type (hybrid_t)      , intent(in) :: hybrid     
    integer              , intent(in) :: nets
    integer              , intent(in) :: nete
    type(parallel_t)     , intent(in) :: par

    ! Local variables 

    real (kind=real_kind) :: pt(np,np,nets:nete)
    real (kind=real_kind) :: p(np,np,nets:nete)

    real (kind=real_kind) :: vt(np,np,2,nets:nete)
    real (kind=real_kind) :: v(np,np,2,nets:nete)

    real (kind=real_kind) :: l1,l2,linf
    integer ie,npts

#if (! defined ELEMENT_OPENMP)
    !$OMP BARRIER
#endif
    if ((tl%nstep == 1).or.(tl%nstep == 0)) then
       if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
          open(iounit+0,file=TRIM(output_prefix)//"swtc5.l1.errors",form="formatted")
          open(iounit+1,file=TRIM(output_prefix)//"swtc5.l2.errors",form="formatted")
          open(iounit+2,file=TRIM(output_prefix)//"swtc5.linf.errors",form="formatted")
       end if
    end if

    do ie=nets,nete
       p(:,:,ie)   = elem(ie)%state%p(:,:,1,tl%n0) + pmean
       v(:,:,:,ie) = elem(ie)%state%v(:,:,:,1,tl%n0)
    end do

    ! ======================================================
    ! write in the reference state for this simulated day...
    ! ======================================================

#ifdef _REFSOLN
!  Parallel version of ref_state, comment out if reading below
    call ref_state_write(p(:,:,nets:nete),v(:,:,:,nets:nete),fstub,simday,nets,nete,par)
    do ie=nets,nete
       pt(:,:,ie)=p(:,:,ie)
       vt(:,:,:,ie)=v(:,:,:,ie)
    end do
#endif

    ! ======================================================
    ! read in the reference state for this simulated day...
    ! ======================================================

#ifdef _REFSOLN
#if (! defined ELEMENT_OPENMP)
!    !$OMP BARRIER
#endif
!  Parallel version of ref_state, comment out if writing above
    call ref_state_read(pt(:,:,nets:nete),vt(:,:,:,nets:nete),fstub,simday,nets,nete,par)
#if (! defined ELEMENT_OPENMP)
    !$OMP BARRIER
#endif

    npts=np

    l1   = l1_snorm(elem,p(:,:,nets:nete),  pt(:,:,nets:nete),hybrid,npts,nets,nete)
    l2   = l2_snorm(elem,p(:,:,nets:nete),  pt(:,:,nets:nete),hybrid,npts,nets,nete)
    linf = linf_snorm(p(:,:,nets:nete),pt(:,:,nets:nete),hybrid,npts,nets,nete)

    if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
       write(iounit+0,30)REAL(simday),l1
       write(iounit+1,30)REAL(simday),l2
       write(iounit+2,30)REAL(simday),linf
       print *,simday, "L1=",l1
       print *,simday, "L2=",l2
       print *,simday, "Linf=",linf
    end if
#if (! defined ELEMENT_OPENMP)
    !$OMP BARRIER
#endif
#endif

30  format(f11.6,4x,e13.6)

  end subroutine tc5_errors
#if 0
  subroutine tc5_init_diag(iounit,t1,hybrid)
    implicit none

    type (hybrid_t) , intent(in) :: hybrid
    type (TimeLevel_t)  , intent(in) :: tl    
    integer ,         intent(in) :: iounit
    real(kind=real_kind)         :: time,time_io,tmp

    time   = Time_at(tl%nstep)/secpday

    if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
       open(iounit+0,file=TRIM(output_prefix)//"swtc5.mass",status="unknown",form="formatted")
       open(iounit+1,file=TRIM(output_prefix)//"swtc5.energy",status="unknown",form="formatted")
       open(iounit+2,file=TRIM(output_prefix)//"swtc5.penst",status="unknown",form="formatted")
       open(iounit+3,file=TRIM(output_prefix)//"swtc5.vort",status="unknown",form="formatted")
       open(iounit+4,file=TRIM(output_prefix)//"swtc5.div",status="unknown",form="formatted")
       open(iounit+5,file=TRIM(output_prefix)//"swtc5.rmsdiv",status="unknown",form="formatted")
       !==================================
       ! If we have a restart run
       !==================================
       if(tl%nstep .ne. 0) then 
          !========================================================
          ! Advance diagnostic output file to the proper time index 
          !========================================================
          do while ( time_io .lt. time ) 
             read(iounit+0,*)time_io,tmp,Imass0
             read(iounit+1,*)time_io,tmp,Ienergy0
             read(iounit+2,*)time_io,tmp,Ipenst0
             read(iounit+3,*)time_io,tmp
             read(iounit+4,*)time_io,tmp
             read(iounit+5,*)time_io,tmp
          enddo
       endif
    end if


  end subroutine tc5_init_diag
#endif


  ! ======================================================
  ! tc5_invariants_nonstag:
  ! compute invariants on page 221 of Williamson, et. al.
  ! on a non-staggered grid.
  ! ======================================================
  subroutine tc5_invariants(elem, iounit,tl,pmean,edge2,deriv,hybrid,nets,nete)
    type(element_t), intent(inout) :: elem(:)

    integer              , intent(in) :: iounit
    type (TimeLevel_t)  , intent(in) :: tl    
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
    real (kind=real_kind) :: h,hstar,hs

    real (kind=real_kind) :: time
    real (kind=real_kind) :: Imass,Ienergy,Ipenst,Ivort,Idiv,Idivsq

    integer               :: npts
    integer               :: n0

    n0 = tl%n0

#if (! defined ELEMENT_OPENMP)
    !$OMP BARRIER
#endif
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

    do ie=nets,nete

       do j=1,np
          do i=1,np

             v1     = elem(ie)%state%v(i,j,1,1,n0)
             v2     = elem(ie)%state%v(i,j,2,1,n0)

             vco(i,j,1) = elem(ie)%met(1,1,i,j)*v1 + elem(ie)%met(1,2,i,j)*v2
             vco(i,j,2) = elem(ie)%met(2,1,i,j)*v1 + elem(ie)%met(2,2,i,j)*v2

             gv(i,j,1) = elem(ie)%metdet(i,j)*v1
             gv(i,j,2) = elem(ie)%metdet(i,j)*v2

          end do
       end do

       zeta(:,:,ie)  = vorticity(vco,deriv)*rrearth
       div(:,:,ie)   = divergence(gv,deriv)*rrearth

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

    call bndry_exchangeV(hybrid,edge2)

    do ie=nets,nete      

       kptr=0
       call edgeVunpack(edge2, zeta(1,1,ie), 1, kptr, elem(ie)%desc)

       kptr=1
       call edgeVunpack(edge2, div(1,1,ie), 1, kptr, elem(ie)%desc)

       do j=1,np
          do i=1,np
             zeta(i,j,ie) = elem(ie)%rmp(i,j)*zeta(i,j,ie)
             div(i,j,ie) = elem(ie)%rmp(i,j)*div(i,j,ie)
             penst(i,j,ie) = 0.5D0*(zeta(i,j,ie) + elem(ie)%fcor(i,j))**2
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
             hstar         = (elem(ie)%state%p(i,j,1,n0) + pmean)/g
             hs            = elem(ie)%state%ps(i,j)/g
             h             = hstar + hs
             mass(i,j,ie)    = hstar
             energy(i,j,ie)  = hstar*E(i,j) + 0.5D0*g*(h**2 - hs**2)
             penst(i,j,ie)   = penst(i,j,ie)/hstar
             divsq(i,j,ie)   = div(i,j,ie)*div(i,j,ie)  
             vort(i,j,ie)    = zeta(i,j,ie)
          end do
       end do

    end do

    npts=np
    Imass   = global_integral(elem,mass(:,:,nets:nete),hybrid,npts,nets,nete)
    Ienergy = global_integral(elem,energy(:,:,nets:nete),hybrid,npts,nets,nete)
    Ipenst  = global_integral(elem,penst(:,:,nets:nete),hybrid,npts,nets,nete)
    Ivort   = global_integral(elem,vort(:,:,nets:nete),hybrid,npts,nets,nete)
    Idiv    = global_integral(elem,div(:,:,nets:nete),hybrid,npts,nets,nete)
    Idivsq  = global_integral(elem,divsq(:,:,nets:nete),hybrid,npts,nets,nete)

    time   = Time_at(tl%nstep)

#if (! defined ELEMENT_OPENMP)
    !$OMP BARRIER
#endif
    if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
       !print *,'mass, pmean/g = ',Imass,pmean/g
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

#if (! defined ELEMENT_OPENMP)
    !$OMP BARRIER
#endif

  end subroutine tc5_invariants

  ! ===========================================
  !
  ! tc6_init_state:
  !
  ! Initialize test case 6: Rossby-Haurwitz Wave
  !
  ! ===========================================

  subroutine tc6_init_state(elem,nets,nete,pmean)
    type(element_t), intent(inout) :: elem(:)

    integer, intent(in) :: nets
    integer, intent(in) :: nete
    real (kind=real_kind) :: pmean

    ! Local variables

    integer :: ie,k
    integer :: nm1 
    integer :: n0 
    integer :: np1

    nm1= 1
    n0 = 2
    np1= 3

    pmean = tc6_init_pmean()     ! (m^2/s^2) set pmean value specified in Williamson, et. al.
    ! p 222, for test case 6.

    do ie=nets,nete
       do k=1,nlev
          elem(ie)%state%p(:,:,k,n0)=tc6_geopotential(elem(ie)%spherep(:,:),np)
          elem(ie)%state%p(:,:,k,nm1)=elem(ie)%state%p(:,:,k,n0)
          elem(ie)%state%p(:,:,k,np1)=0.0D0

          elem(ie)%state%v(:,:,:,k,n0)=tc6_velocity(elem(ie)%spherep(:,:),elem(ie)%Dinv)
          elem(ie)%state%v(:,:,:,k,nm1)=elem(ie)%state%v(:,:,:,k,n0)
          elem(ie)%state%v(:,:,:,k,np1)=0.0D0
       end do
       elem(ie)%state%ps(:,:)=0.0D0
       elem(ie)%state%gradps(:,:,:)=0.0D0
    end do

  end subroutine tc6_init_state

  ! ===========================================
  !
  ! tc6_geopotential:
  !
  ! Set geopotential to initial values 
  ! specified on page 222 of Williamson, et al.
  ! eqs. 146-149
  !
  ! ===========================================

  function tc6_geopotential(sphere,npts) result(p)

    integer, intent(in)                  :: npts
    type (spherical_polar_t), intent(in) :: sphere(npts,npts)
    real (kind=real_kind)                :: p(npts,npts)

    ! Local variables

    real (kind=real_kind) :: lat
    real (kind=real_kind) :: lon
    real (kind=real_kind) :: cslat
    real (kind=real_kind) :: A,B,C
    real (kind=real_kind) :: Afac,Bfac

    integer i,j

    do j=1,npts
       do i=1,npts
          lat   = sphere(i,j)%lat
          lon   = sphere(i,j)%lon
          cslat = COS(lat)
          Afac  = (R_rh+1)*(cslat**(2*(R_rh+1)))                &
               +  (2*(R_rh**2) - R_rh - 2)*(cslat**(2*R_rh))   &
               -   2*(R_rh**2)*(cslat**(2*(R_rh-1)))
          A     = 0.5D0*omega_rh*(2*omega +omega_rh)*(cslat**2) + 0.25D0*(K_rh**2)*Afac
          Bfac  = 2.0D0*(omega + omega_rh)*K_rh/((R_rh+1)*(R_rh+2))
          B     = Bfac*(cslat**R_rh)*((R_rh**2 + 2*R_rh + 2) - ((R_rh+1)**2)*(cslat**2))
          C     = 0.25d0*(K_rh**2)*(cslat**(2*R_rh))*((R_rh+1)*(cslat**2) - (R_rh+2))
          p(i,j)= rearth*rearth*(A + B*COS(R_rh*lon) + C*COS(2*R_rh*lon))
       end do
    end do

  end function tc6_geopotential

  ! ===========================================
  !
  ! tc6_velocity:
  !
  ! Set geopotential to initial values 
  ! specified on page 222 of Williamson, et al.
  ! eqs. 143-144
  !
  ! ===========================================

  function tc6_velocity(sphere,D) result(v) 

    type (spherical_polar_t), intent(in) :: sphere(np,np)
    real (kind=real_kind),    intent(in) :: D(2,2,np,np)
    real (kind=real_kind)                :: v(np,np,2)

    ! Local variables

    real (kind=real_kind) :: lat
    real (kind=real_kind) :: lon
    real (kind=real_kind) :: cslat
    real (kind=real_kind) :: snlat
    real (kind=real_kind) :: V1,V2

    integer i,j

    do j=1,np
       do i=1,np
          lat   = sphere(i,j)%lat
          lon   = sphere(i,j)%lon
          cslat = COS(lat)
          snlat = SIN(lat)

          V1 =  rearth*(omega_rh*cslat + K_rh*(cslat**(R_rh-1))*(R_rh*(snlat**2) - (cslat**2))*COS(R_rh*lon))
          V2 = -rearth*K_rh*R_rh*(cslat**(R_rh-1))*snlat*SIN(R_rh*lon)

          ! =====================================================
          ! map sphere velocities onto the contravariant cube velocities
          ! using the D mapping matrix (see Loft notes for details)
          ! =====================================================

          v(i,j,1)= V1*D(1,1,i,j) + V2*D(1,2,i,j)
          v(i,j,2)= V1*D(2,1,i,j) + V2*D(2,2,i,j)

       end do
    end do

  end function tc6_velocity

  ! ===========================================
  ! tc6_errors:
  !
  ! Compute Shallow Water Test Case 6 errors
  ! 
  ! ===========================================

  subroutine tc6_errors(elem, iounit, tl, pmean, fstub, simday, hybrid, nets, nete, par)
    type(element_t), intent(in) :: elem(:)

    integer              , intent(in) :: iounit
    type (TimeLevel_t)   , intent(in) :: tl         ! model time struct
    real (kind=real_kind), intent(in) :: pmean      ! mean geopotential
    character(len=*)     , intent(in) :: fstub      ! file path stub
    integer              , intent(in) :: simday     ! current day of simulation
    type (hybrid_t)      , intent(in) :: hybrid     
    integer              , intent(in) :: nets
    integer              , intent(in) :: nete
    type(parallel_t)     , intent(in) :: par

    ! Local variables 

    real (kind=real_kind) :: pt(np,np,nets:nete)
    real (kind=real_kind) :: p(np,np,nets:nete)

    real (kind=real_kind) :: vt(np,np,2,nets:nete)
    real (kind=real_kind) :: v(np,np,2,nets:nete)

    real (kind=real_kind) :: l1,l2,linf
    integer ie,npts

#if (! defined ELEMENT_OPENMP)
    !$OMP BARRIER
#endif
    if (tl%nstep == 0) then
       if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
          open(iounit+0,file=TRIM(output_prefix)//"swtc6.l1.errors",form="formatted")
          open(iounit+1,file=TRIM(output_prefix)//"swtc6.l2.errors",form="formatted")
          open(iounit+2,file=TRIM(output_prefix)//"swtc6.linf.errors",form="formatted")
       end if
    end if

    do ie=nets,nete
       p(:,:,ie)   = elem(ie)%state%p(:,:,1,tl%n0) + pmean
       v(:,:,:,ie) = elem(ie)%state%v(:,:,:,1,tl%n0)
    end do

    ! ======================================================
    ! read in the reference state for this simulated day...
    ! ======================================================

#ifdef _REFSOLN
!  Parallel version of ref_state, comment out read below if writing here
    call ref_state_write(p(:,:,nets:nete),v(:,:,:,nets:nete),fstub,simday,nets,nete,par)
    do ie=nets,nete
       pt(:,:,ie)=p(:,:,ie)
       vt(:,:,:,ie)=v(:,:,:,ie)
    end do
#endif

    ! ======================================================
    ! read in the reference state for this simulated day...
    ! ======================================================
#ifdef _REFSOLN
#if (! defined ELEMENT_OPENMP)
!    !$OMP BARRIER
#endif
!  Parallel version of ref_state, comment out if writing above
    call ref_state_read(pt(:,:,nets:nete),vt(:,:,:,nets:nete),fstub,simday,nets,nete,par)
#if (! defined ELEMENT_OPENMP)
    !$OMP BARRIER
#endif

    ! Get the state variables from the simulation (layer 1 only)...

    npts = np

    l1   = l1_snorm(elem,p(:,:,nets:nete),  pt(:,:,nets:nete),hybrid,npts,nets,nete)
    l2   = l2_snorm(elem,p(:,:,nets:nete),  pt(:,:,nets:nete),hybrid,npts,nets,nete)
    linf = linf_snorm(p(:,:,nets:nete),pt(:,:,nets:nete),hybrid,npts,nets,nete)

    if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
       write(iounit+0,30)REAL(simday),l1
       write(iounit+1,30)REAL(simday),l2
       write(iounit+2,30)REAL(simday),linf
       print *,simday, "L1=",l1
       print *,simday, "L2=",l2
       print *,simday, "Linf=",linf
    end if
#endif
#if (! defined ELEMENT_OPENMP)
    !$OMP BARRIER
#endif

10  format("      time           l1              l2              linf           ")
20  format("====================================================================")
30  format(f11.6,4x,e13.6)

  end subroutine tc6_errors

  ! ===================================================================
  !
  ! tc8_init_state:
  !
  ! Initialize test case 8: Barotropic Instability Test Case of Polvani
  !
  ! ===================================================================

  subroutine tc8_init_state(elem, nets,nete,hybrid,pmean)
    type(element_t), intent(inout) :: elem(:)

    integer, intent(in)   :: nets
    integer, intent(in)   :: nete
    type (hybrid_t)       :: hybrid
    real (kind=real_kind) :: pmean

    ! Local variables

    real (kind=real_kind) :: a,b
    real (kind=real_kind) :: mean_balance
    real (kind=real_kind) :: pert(np,np)
    real (kind=real_kind) :: phi_balance(np,np,nets:nete)
    type (quadrature_t)   :: gs

    integer :: ie,k
    integer :: nm1 
    integer :: n0 
    integer :: np1

    nm1= 1
    n0 = 2
    np1= 3

    ! ================================================
    ! Initialize the total balance integral... RDL
    ! ================================================

    a = lat_tc8 - gamma_tc8
    b = lat_tc8 + gamma_tc8

!!!RDL Ibalance_tot = trapezoid(tc8_integrand,a,b,eps_tc8)
!!!RDL Ibalance_tot = simpsons(tc8_integrand,a,b,eps_tc8)

    gs=gauss(Ngs_tc8)
    Ibalance_tot = gaussian_int(tc8_integrand,a,b,gs)    
    print *,"Ibalance_tot=",Ibalance_tot
    deallocate(gs%points)
    deallocate(gs%weights)
    print *,"done..."
    do ie=nets,nete
       phi_balance(:,:,ie) = tc8_balance(elem(ie)%spherep(:,:),np)
    end do

    pmean         = tc8_init_pmean()
    mean_balance  = global_integral(elem,phi_balance(:,:,nets:nete),hybrid,np,nets,nete)
    print *,"mean_balance=",mean_balance," pmean=",pmean

    do ie=nets,nete
       pert(:,:) = tc8_perturbation(elem(ie)%spherep(:,:),np)
       elem(ie)%state%ps(:,:) = phi_balance(:,:,ie) - mean_balance + pert(:,:)
       elem(ie)%state%ps(:,:) = phi_balance(:,:,ie)
       do k=1,nlev
          elem(ie)%state%p(:,:,k,n0)=elem(ie)%state%ps(:,:)
          elem(ie)%state%p(:,:,k,nm1)=elem(ie)%state%p(:,:,k,n0)
          elem(ie)%state%p(:,:,k,np1)=0.0D0

          elem(ie)%state%v(:,:,:,k,n0)=tc8_velocity(elem(ie)%spherep(:,:),elem(ie)%Dinv)
          elem(ie)%state%v(:,:,:,k,nm1)=elem(ie)%state%v(:,:,:,k,n0)
          elem(ie)%state%v(:,:,:,k,np1)=0.0D0
       end do
       elem(ie)%state%ps(:,:)=0.0D0
       elem(ie)%state%gradps(:,:,:)=0.0D0
    end do

  end subroutine tc8_init_state

  ! =======================================================
  !
  ! tc8_integrand:
  ! 
  ! Scaled integrand function needed to integrate the 
  ! geostrophic balance equation using the trapezoid rule.
  ! Integrand in units of u0_tc8*f*rearth (m^2/s^2)
  !
  ! =======================================================

  function tc8_integrand(lat) result(Ival)

    real (kind=real_kind), intent(in) :: lat
    real (kind=real_kind) :: Ival

    real (kind=real_kind) :: cs_tc8
    real (kind=real_kind) :: ulat    ! dimensionless velocity
    real (kind=real_kind) :: epsilon ! unitless factor in integral

    real (kind=real_kind) :: a
    real (kind=real_kind) :: pi2

    pi2    = 0.5D0*DD_PI
    a      = pi2/gamma_tc8
    cs_tc8 = COS(a*(lat - lat_tc8))
    ulat   =  cs_tc8*cs_tc8
    epsilon = u0_tc8/(rearth*2.0D0*omega)    


    Ival = ulat*(sin(lat) + epsilon*(TAN(lat)*ulat))


  end function tc8_integrand

  ! ==========================================================
  !
  ! tc8_balance:
  !
  ! Compute geostrophic integral balance equation for the
  ! Barotropic Instability Test Case of Polvani.
  !
  ! ===========================================================

  function tc8_balance(sphere,npts) result(p)

    integer, intent(in)                  :: npts
    type (spherical_polar_t), intent(in) :: sphere(npts,npts)
    real (kind=real_kind)                :: p(npts,npts)

    ! Local variables

    real (kind=real_kind) :: lat
    real (kind=real_kind) :: lon

    real (kind=real_kind) :: a, b
    real (kind=real_kind) :: Ibalance
    real (kind=real_kind) :: f
    real (kind=real_kind) :: pi2  ! pi/2

    integer i,j
    type (quadrature_t) :: gs


    gs=gauss(Ngs_tc8)

    do j=1,npts
       do i=1,npts

          lat   = sphere(i,j)%lat
          lon   = sphere(i,j)%lon
          f     = 2.0D0*omega  ! coreolis coefficient
          pi2   = 0.5D0*DD_PI

          if (lat <= (lat_tc8-gamma_tc8) ) then
             Ibalance = -Ibalance_tot
          else
             if (lat >= (lat_tc8 + gamma_tc8)) then
                Ibalance = 0.0D0
             else
                a = lat_tc8 + gamma_tc8
                b = lat
!!!DBG            Ibalance = trapezoid(tc8_integrand,a,b,eps_tc8)
!!!DBG            Ibalance = simpsons(tc8_integrand,a,b,eps_tc8)
                Ibalance = gaussian_int(tc8_integrand,a,b,gs)    
             end if

          end if

          p(i,j)  = - (rearth*f*u0_tc8)*Ibalance ! may be wrong.. RDL

       end do
    end do

    deallocate(gs%weights)
    deallocate(gs%points)

  end function tc8_balance

  ! ===========================================
  !
  ! tc8_perturbation:
  !
  ! Compute the geopotential perturbation 
  ! from geostrophic balance for the Barotropic 
  ! Instability Test Case of Polvani.
  !
  ! ===========================================

  function tc8_perturbation(sphere,npts) result(p)

    integer, intent(in)                  :: npts
    type (spherical_polar_t), intent(in) :: sphere(npts,npts)
    real (kind=real_kind)                :: p(npts,npts)

    ! Local variables

    real (kind=real_kind) :: lat
    real (kind=real_kind) :: lon
    real (kind=real_kind) :: alon, blat
    real (kind=real_kind) :: secha,sechb
    real (kind=real_kind) :: ghprime

    integer i,j

    do j=1,npts
       do i=1,npts

          lat   = sphere(i,j)%lat
          lon   = sphere(i,j)%lon

          if(lon.le.DD_PI) alon = alpha_tc8*(lon)
          if(lon.gt.DD_PI) alon = alpha_tc8*(2.0D0*DD_PI - lon)

          blat  = beta_tc8*(lat - lat_tc8) 
          secha = 2.0D0/(EXP(alon) + EXP(-alon))
          sechb = 2.0D0/(EXP(blat) + EXP(-blat))

          ghprime = g*hhat_tc8*(secha*sechb)**2
          p(i,j)  = ghprime
       end do
    end do

  end function tc8_perturbation

  ! ===========================================
  !
  ! tc8_velocity:
  !
  ! Set initial zonal velocity field for barotropic 
  ! instability test case of Polvani.
  !
  ! ===========================================

  function tc8_velocity(sphere,D) result(v) 

    type (spherical_polar_t), intent(in) :: sphere(np,np)
    real (kind=real_kind),    intent(in) :: D(2,2,np,np)
    real (kind=real_kind)                :: v(np,np,2)

    ! Local variables

    real (kind=real_kind) :: lat
    real (kind=real_kind) :: cs_tc8
    real (kind=real_kind) :: V1,V2

    integer i,j

    do j=1,np
       do i=1,np
          lat   = sphere(i,j)%lat

          if (ABS(lat-lat_tc8) <= gamma_tc8) then 
             cs_tc8 = COS((DD_PI/2.0D0)*((lat - lat_tc8)/gamma_tc8))
             V1 =  u0_tc8*(cs_tc8**2)
          else
             V1 = 0.0D0
          end if
          V2 =  0.0D0

          ! =====================================================
          ! map sphere velocities onto the contravariant cube velocities
          ! using the D^-T mapping matrix (see Loft notes for details)
          ! =====================================================

          v(i,j,1)= V1*D(1,1,i,j) + V2*D(1,2,i,j)
          v(i,j,2)= V1*D(2,1,i,j) + V2*D(2,2,i,j)

       end do
    end do

  end function tc8_velocity
  
  ! ===========================================
  !
  ! vortex_init_state:
  !
  ! Initializes vortex
  !
  ! ===========================================
  subroutine vortex_init_state(elem, nets,nete,pmean)
    type(element_t), intent(inout) :: elem(:)
    integer, intent(in) :: nets
    integer, intent(in) :: nete
    real (kind=real_kind) :: pmean, phiuv(3)

    ! Local variables

    integer ie,k, i, j
    integer :: nm1 
    integer :: n0 
    integer :: np1

    nm1= 1
    n0 = 2
    np1= 3

    pmean = 0.0D0

    do ie=nets,nete
       do k=1,nlev
          do j=1,np
             do i=1,np
                phiuv = dynamic_velocity(elem(ie)%spherep(i,j)%lon,elem(ie)%spherep(i,j)%lat,0.0D0)
                elem(ie)%state%p(i,j,k,n0)=phiuv(1)
             end do
          end do

          elem(ie)%state%v(:,:,:,k,n0)=vortex_velocity(0.0D0,elem(ie)%spherep(:,:),elem(ie)%Dinv)
          elem(ie)%state%v(:,:,:,k,nm1)=elem(ie)%state%v(:,:,:,k,n0)
          elem(ie)%state%v(:,:,:,k,np1)=elem(ie)%state%v(:,:,:,k,n0)

          elem(ie)%state%p(:,:,k,nm1)=elem(ie)%state%p(:,:,k,n0)
          elem(ie)%state%p(:,:,k,np1)=elem(ie)%state%p(:,:,k,n0)

       end do
       elem(ie)%state%ps(:,:)=0.0D0
       elem(ie)%state%gradps(:,:,:)=0.0D0
    end do

  end subroutine vortex_init_state

!--------------------------------------------------------------------------

  function vortex_exact(t,sphere) result(p)

    type (spherical_polar_t), intent(in) :: sphere(np,np)
    real (kind=real_kind), intent(in) :: t
    real (kind=real_kind) :: phiuv(3)
    real (kind=real_kind) :: p(np,np)

    ! Local variables

    integer i, j

    do j=1,np
       do i=1,np
          !phiuv = dynamic_vortex(sphere(i,j)%lon,sphere(i,j)%lat,t,u0_vortex,rearth,lon0_vort1_vortex,lat0_vort1_vortex)
          phiuv = dynamic_velocity(sphere(i,j)%lon,sphere(i,j)%lat,t)
          p(i,j)=phiuv(1)
       end do
    end do
    
  end function vortex_exact

!----------------------------------------------------------------------

  function vortex_velocity(t,sphere,D) result(v)

    type (spherical_polar_t), intent(in) :: sphere(np,np)
    real (kind=real_kind),    intent(in) :: t,D(2,2,np,np)
    real (kind=real_kind)                :: v(np,np,2)

    ! Local variables

    real (kind=real_kind) :: snlon(np,np)
    real (kind=real_kind) :: cslon(np,np)
    real (kind=real_kind) :: snlat(np,np)
    real (kind=real_kind) :: cslat(np,np)
    real (kind=real_kind) :: csalpha
    real (kind=real_kind) :: snalpha
    real (kind=real_kind) :: V1,V2
    real (kind=real_kind) :: phiuv(3)

    integer i,j,k

    csalpha = COS(alpha_vortex)
    snalpha = SIN(alpha_vortex)

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

             !phiuv = dynamic_vortex(sphere(i,j)%lon,sphere(i,j)%lat,t,u0_vortex,rearth,lon0_vort1_vortex,lat0_vort1_vortex)
             phiuv = dynamic_velocity(sphere(i,j)%lon,sphere(i,j)%lat,t)

             V1 =   u0_vortex*(cslat(i,j)*csalpha + snlat(i,j)*cslon(i,j)*snalpha)
             V2 =  -u0_vortex*(snlon(i,j)*snalpha)

             ! =====================================================
             ! map sphere velocities onto the contravariant cube velocities
             ! using the D^-T mapping matrix (see Loft notes for details)
             ! =====================================================

             !phiuv=0.0*phiuv
             if(static_vortex)then
                V1=0.0D0
                V2=0.0D0
             endif
             v(i,j,1)= (V1+phiuv(2))*D(1,1,i,j) + (V2+phiuv(3))*D(1,2,i,j)
             v(i,j,2)= (V1+phiuv(2))*D(2,1,i,j) + (V2+phiuv(3))*D(2,2,i,j)
          end do
       end do
    end do

  end function vortex_velocity

!--------------------------------------------------------------------------------

  subroutine vortex_rotatedsphere(lc,tc,la,th,rol,rot)
    ! Rotate to new North Pole at (lc,tc)
    ! (rol,rot) are the rotated coordinates coorsponding to (la,th) in 
    ! the regular sphere
    implicit none
    real (kind=real_kind), intent(in)  :: lc,tc,la,th
    real (kind=real_kind), intent(out) :: rol,rot
    
    real (kind=real_kind) :: cost,sint,sinc,cosc
    real (kind=real_kind) ::  trm, trm1,trm2,trm3
    real (kind=real_kind) ::  pi2
    
    pi2=2.0D0*DD_PI

    sinc = sin(tc)
    cosc = cos(tc)
    cost = cos(th)
    sint = sin(th)
    
    trm  = cost * cos(la- lc)
    trm1 = cost * sin(la- lc)
    trm2 = sinc * trm  - cosc * sint
    trm3 = sinc * sint + cosc * trm
    
    rol = atan2(trm1,trm2)
    !if (rol > pi2  ) rol = rol - pi2
    !if (rol < 0.0D0 ) rol = rol + pi2

    rot = asin(trm3)      

  end subroutine vortex_rotatedsphere
  
!----------------------------------------------------------------------

  subroutine vortex_rotatedsphereback(lmc,thc,lam,the,rla,rth)
    
    implicit none
    real (kind=real_kind), intent(in)  :: lmc,thc,lam,the
    real (kind=real_kind), intent(out) :: rla,rth
    !
    real (kind=real_kind) :: cost,sint,cosp,sinp,clam,slam 
    real (kind=real_kind) :: trm, t1,t2,t3
    real (kind=real_kind) ::  pi2
    
    pi2=2.0D0*DD_PI

    !
    !* Back to unrotated system

    cost = cos(the)
    sint = sin(the)
    clam = cos(lam)
    slam = sin(lam)
    cosp = cos(thc)
    sinp = sin(thc)

    t1 = slam * cost
    t2 = sint*cosp + cost*clam*sinp
    t3 = sint*sinp - cost*clam*cosp
    rla =  lmc + atan2(t1,t2)

    !if (rla < 0.0D0 )  rla = rla + pi2
    !if (rla > pi2 )   rla = rla - pi2

    rth =  asin(t3)
    
  end subroutine vortex_rotatedsphereback

!----------------------------------------------------------------------------

  function dynamic_vortex(lm_in,th_in,itn,unot,a,lm0,th0) result(fld)

    implicit none
    real (kind=real_kind), intent(in) :: lm_in,th_in,itn,unot,a,lm0,th0

    ! itn is the absolute time
 
    real (kind=real_kind)   :: fld(3), tow
    
    real (kind=real_kind)   :: lm,th,plat, plon, ang
    real (kind=real_kind)   :: rot,rol, ldt,tdt
    
    real (kind=real_kind)   :: lmc,thc 
    real (kind=real_kind)   :: pl,pt
    
    real (kind=real_kind)   :: dd,rd, rl,rt, trm
    real (kind=real_kind)   :: tt,tv,st,ct, omg, vmax,rho,lam,the
    
    lm = lm_in
    th = th_in

    !
    !    Find the position of the vortex center as a function of
    !    alpha and dt
    !

    if(.NOT.static_vortex)then
       lmc =  DD_PI
       thc =  0.5D0*DD_PI - alpha_vortex
       !
       !  Go to the rotated system
       !
       call vortex_rotatedsphere(lmc,thc,lm,th,rl,rt)
       lam = rl  + omega_s_vortex*itn  
       the = rt

       !
       !  Back to unrotated system 
       !
       call vortex_rotatedsphereback(lmc,thc,lam,the,rl,rt)
       lm =  rl
       th =  rt

    end if

    lmc = lm0             !center of the vortex at t=0
    thc = th0
    
    !    Note:  (lmc,thc) is the new north pole position of the rotated coordinates
    !     relative to the given regular (la,th) coordinates.
    pl = lmc
    pt = thc 

    ! Vortex parameters
    dd = 5.0D0
    rd = 3.0D0
    vmax = 1.5D0 *SQRT(3.0D0)*u0_vortex

    !  Analytic Vortex field  at any given time (iteration)

    !Rotated coordinates w.r.t (pl,pt)
    call vortex_rotatedsphere(pl,pt,lm,th,rl,rt)

    rho = rd * COS(rt)
    tt = TANH(rho)
    ! eqn 13
    tv = (1.0D0 - tt**2) * tt * vmax
    
    if(ABS(rho)<1.0D-13)then
       omg = 0.0D0
    else
       omg = tv / (rho*a)
    end if

    tow   = itn 
    !Exact soln.
    trm =  rho*sin(rl - tow*omg) / dd

    !print *,"a*omg = ", a*omg

    ! exact sol
    fld(1) =  (1.0D0 - tanh(trm))
    ! Velo just rotation no translation
    !fld(2) = a*omg*(sin(th)*cos(th_in) - cos(th)*cos(lm_in-lm)*sin(th_in))
    !fld(3) = a*omg*cos(th)*sin(lm_in-lm)
    !fld(2) = a*omg*(sin(the)*cos(th_in) - cos(the)*cos(lm_in-lam)*sin(th_in))
    !fld(3) = a*omg*cos(the)*sin(lm_in-lam)
#if 1
    if(.NOT.static_vortex)then
       if(abs(alpha_vortex)<1.0D-13)then
       else
          lmc =  0.0D0
          thc =  0.5D0*DD_PI - alpha_vortex
          !
          !  Go to the rotated system
          !
          call vortex_rotatedsphere(lmc,thc,lm0,th0,rl,rt)
          lam = rl  + omega_s_vortex*itn  
          the = rt
          
          !
          !  Back to unrotated system 
          !
          call vortex_rotatedsphereback(lmc,thc,lam,the,rl,rt)
          lmc =  rl
          thc =  rt
       end if
    end if

    fld(2) = a*omg*(sin(thc)*cos(th_in) - cos(thc)*cos(lm_in-lmc)*sin(th_in))
    fld(3) = a*omg*cos(thc)*sin(lm_in-lmc)
    !fld(2) = a*omg*(sin(the)*cos(th_in) - cos(the)*cos(lm_in-lam)*sin(th_in))
    !fld(3) = a*omg*cos(the)*sin(lm_in-lam)

#endif

  end function dynamic_vortex

!------------------------------------------------------------------------------

  function  dynamic_velocity(lm,th,itn) result(vout)      
    
    implicit none
    
    real (kind=real_kind), intent(in) :: itn
    real (kind=real_kind), intent(in) :: lm,th
    real (kind=real_kind)             :: vout(3)

    real (kind=real_kind) :: pl, pt, lm0, th0, lmc,thc,rol,rot,ldt,tdt
    real (kind=real_kind) :: pi2, t1,t2,t3, dd,rd, half, one
    real (kind=real_kind) :: sint,cost,tant, tanl,su,sv, rla,rth, trm, c1,c2
    real (kind=real_kind) :: sinp,cosp,slam,clam, x,y , omg, tow, scale
    real (kind=real_kind) :: tt,tv,st,ct, vmax,rho , three , unot, us,vs ,ditn
    
    integer ::   i,j, n1,n2 ,l,k, ip
    
    pi2 = DD_PI * 0.5D0
    
    !  alpha:  Flow direction and  Rotational parameters (given)    
    
    unot = omega_s_vortex
 
    !intial center position of the bell
    lm0 = lon0_vort1_vortex
    th0 = lat0_vort1_vortex
    lmc = lm0
    thc = th0

    c1 = cos(alpha_vortex)
    c2 = sin(alpha_vortex)
    
   !To find the  central position of the distribution as a function of
   !alpha and dt

    ditn = itn
    
    if (abs(alpha_vortex)< 1.0D-13) then
       
       ldt = lm0 + unot * ditn
       tdt = th0 + 0.0D0
       rol  = ldt
       rot  = tdt

    else

       lmc = 0.0D0
       thc = pi2 - alpha_vortex
       ldt = lm0 + unot * ditn
       tdt = th0 + 0.0D0
       call vortex_rotatedsphere(lmc,thc,ldt,tdt,rol,rot)
      
    endif
    
    pl = rol       !center of the solid-body dbn.
    pt = rot
    
    ! Velocity Components in Global (spu,spv) spherical coordinates
    !  & Vortex parameters (see Nair et al. 2002)
    
    ct = cos(pt)
    st = sin(pt)
    
    c1 = cos(alpha_vortex)
    c2 = sin(alpha_vortex)
    
    dd  = 5.0D0      
    rd  = 3.0D0 
    tow = 0.0D0     !Time units

    scale = a_vortex   !vortex rotation rate
    unot  = u0_vortex

    vmax = unot*1.5D0 *sqrt(3.0D0)       !Analytically determined
    
    cost = cos(th)
    sint = sin(th)
    clam = cos(lm)
    slam = sin(lm)
    cosp = cos(lm - pl)
    sinp = sin(lm - pl)
    
        
    !!Rotated coordinates w.r.t (pl,pt)
    t1 = sinp * cost
    t2 = st * cosp*cost - ct * sint
    t3 = st * sint + ct * cost*cosp
    rth = asin(t3)    
    rla = atan2(t1,t2)
    if (rla <  0.0D0) rla = rla + 2.0D0*DD_PI
    
    !! or Call rotated_sphere(pl,pt,lm,th,rla,rth)
    call vortex_rotatedsphere(pl,pt,lm,th,rla,rth)

    rho = rd * cos(rth) 
    tt = tanh(rho)
    tv = (1.0D0  - tt**2) * tt * vmax   !Normalized tangential velocity
    
    if (abs(rho) <1.0D-13 ) then
       omg = 0.0D0
    else
       omg = tv / (rho*scale) 
    endif
    

    ! Combined velocity (keep us=vs=0.0 for static vortex)
    
    vout(2) =  scale *omg * (cost*st - cosp*sint*ct)  
    vout(3) =  scale *omg * sinp*ct 
    
    !call vortex_rotatedsphere(pl,pt,lm,th,rol,rot)
    ! pl,pt,rol,rot
    call vortex_rotatedsphere(pl,pt,lm,th,rla,rth)

    rho = rd * cos(rth) 
    tt  = tanh(rho)
    tv  = (1.0D0  - tt**2) * tt * vmax   !Normalized tangential velocity
    
    if (abs(rho) <1.0D-13 ) then
       omg = 0.0D0
    else
       omg = tv / (rho*scale) 
    endif
    
    trm =  rho*sin(rla - ditn*omg) / dd

    vout(1) =  (1.0D0 - tanh(trm))

  end function dynamic_velocity


!-----------------------------------------------------------------------------

  subroutine vortex_errors(elem, iounit, tl, hybrid, nets, nete)

    type(element_t), intent(inout) :: elem(:)

    integer                           :: iounit
    type (TimeLevel_t)   , intent(in) :: tl         ! model time struct
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

    integer ie,k

#if (! defined ELEMENT_OPENMP)
    !$OMP BARRIER
#endif
    if (tl%nstep == 0) then
       if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
          open(iounit+0,file=TRIM(output_prefix)//"vortex.l1.errors",form="formatted")
          open(iounit+1,file=TRIM(output_prefix)//"vortex.l2.errors",form="formatted")
          open(iounit+2,file=TRIM(output_prefix)//"vortex.linf.errors",form="formatted")
       end if
    end if

    time_tmp   = Time_at(tl%nstep)

    if (lonc < 0.0D0) lonc = lonc + 2.0D0*DD_PI

    do k=1,nlev
       do ie=nets,nete
          pt(:,:,ie)=vortex_exact(time_tmp,elem(ie)%spherep)
          p(:,:,ie) =elem(ie)%state%p(:,:,k,tl%n0)
       end do
       
       l1   = l1_snorm(elem, p(:,:,nets:nete),  pt(:,:,nets:nete),hybrid,np,nets,nete)
       l2   = l2_snorm(elem, p(:,:,nets:nete),  pt(:,:,nets:nete),hybrid,np,nets,nete)
       linf = linf_snorm(p(:,:,nets:nete),pt(:,:,nets:nete),hybrid,np,nets,nete)
       
       if (k==1) then
       if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
          write(iounit+0,30)time_tmp/secpday,l1
          write(iounit+1,30)time_tmp/secpday,l2
          write(iounit+2,30)time_tmp/secpday,linf
       end if
       end if

       if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
          write(*,'(a,i2,f6.2,a,3e15.7)') 'k=',k,time_tmp/secpday,' days  l1,l2,linf=',&
               l1,l2,linf
       end if
    enddo


#if (! defined ELEMENT_OPENMP)
    !$OMP BARRIER
#endif

10  format("      time           l1              l2              linf           ")
20  format("====================================================================")
30  format(f11.6,4x,e13.6)

  end subroutine vortex_errors

!-----------------------------------------------------------------------------

  subroutine swirl_errors(elem, iounit, tl, hybrid, nets, nete)

    type(element_t), intent(inout) :: elem(:)

    integer                           :: iounit
    type (TimeLevel_t)   , intent(in) :: tl         
    type (hybrid_t)      , intent(in) :: hybrid     
    integer, intent(in)               :: nets
    integer, intent(in)               :: nete

    ! Local variables 

    real (kind=real_kind) :: latc,lonc
    real (kind=real_kind) :: lamdot,lon
    real (kind=real_kind) :: l1,l2,linf,time_tmp
    real (kind=real_kind) :: cslonc,snlonc
    real (kind=real_kind) :: cslatc,snlatc

    real (kind=real_kind) :: pt(np,np,nets:nete) ! true tracer value
    real (kind=real_kind) :: p(np,np,nets:nete)  ! obtained tracer value

    integer ie,k,i,j

#if (! defined ELEMENT_OPENMP)
    !$OMP BARRIER
#endif
    if (tl%nstep == 0) then
       if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then

          open(iounit+0,file=TRIM(output_prefix)//"swirl.l1.errors",form="formatted")
          open(iounit+1,file=TRIM(output_prefix)//"swirl.l2.errors",form="formatted")
          open(iounit+2,file=TRIM(output_prefix)//"swirl.linf.errors",form="formatted")
       end if
    end if

    time_tmp   = Time_at(tl%nstep)

    if (lonc < 0.0D0) lonc = lonc + 2.0D0*DD_PI

    do k=1,nlev
       if((k.ne.kmass).and.(kmass>0))then
         do ie=nets,nete
            pt(:,:,ie)=swirl_init_tracer(elem(ie)%spherep(:,:),k)/&
                       swirl_init_tracer(elem(ie)%spherep(:,:),kmass)
            p(:,:,ie) =elem(ie)%state%p(:,:,k,tl%n0)/elem(ie)%state%p(:,:,kmass,tl%n0)
         end do
       else
         do ie=nets,nete
            pt(:,:,ie)=swirl_init_tracer(elem(ie)%spherep(:,:),k)
            p(:,:,ie) =elem(ie)%state%p(:,:,k,tl%n0)
         end do
       endif
     
       l1   = l1_snorm(elem, p(:,:,nets:nete),  pt(:,:,nets:nete),hybrid,np,nets,nete)

       l2   = l2_snorm(elem, p(:,:,nets:nete),  pt(:,:,nets:nete),hybrid,np,nets,nete)

       linf = linf_snorm(p(:,:,nets:nete),pt(:,:,nets:nete),hybrid,np,nets,nete)
       
       if (k==1) then
       if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
          write(iounit+0,30)time_tmp/secpday,l1
          write(iounit+1,30)time_tmp/secpday,l2
          write(iounit+2,30)time_tmp/secpday,linf
       end if
       end if

       if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
          write(*,'(a,i2,f6.2,a,3e15.7)') 'k=',k,time_tmp/secpday,' days  l1,l2,linf=',&
               l1,l2,linf
       end if
    enddo


#if (! defined ELEMENT_OPENMP)
    !$OMP BARRIER
#endif

10  format("      time           l1              l2              linf           ")
20  format("====================================================================")
30  format(f11.6,4x,e13.6)

  end subroutine swirl_errors


  ! ===========================================================================
  !
  ! swirl_init_state:
  !
  ! Initializes swirl
  ! modified by OG
  !
  ! =============================================================================
  subroutine swirl_init_state(elem,nets,nete,pmean,hybrid,edge3)
    type(element_t), intent(inout) :: elem(:)
    real (kind=real_kind), intent(inout) :: pmean
    integer, intent(in) :: nets
    integer, intent(in) :: nete
    real (kind=real_kind) :: phiuv(3)
!I am passing edge3 and hybrid because i want to test limiters
!this is temporary change though
    type (EdgeBuffer_t)  , intent(inout),optional :: edge3 
    type (hybrid_t)      , intent(in),optional :: hybrid  

    ! Local variables
    integer ie,k, i, j
    integer :: nm1 
    integer :: n0 
    integer :: np1
    real (kind=real_kind) :: pmin(nlev,nets:nete),pmax(nlev,nets:nete)


!temporary code to test limiters
!i still need it, OG
    real (kind=real_kind) :: min_var(nlev,nets:nete),max_var(nlev,nets:nete)
    real (kind=real_kind) :: local_var,maxvar_neighb,minvar_neighb,value
    real (kind=real_kind), parameter :: eps=1e-8
!!!!!!!!!!!!!!!!!!11


    pmean=mean_pressure

    if (sub_case==1)then
       lon1=lon1_case1
       lat1=lat1_case1
       lon2=lon2_case1
       lat2=lat2_case1
       Kcoef=2.4   
       add_pure_rotation = .false.
    elseif(sub_case==2)then 
       lon1=lon1_case2
       lat1=lat1_case2
       lon2=lon2_case2
       lat2=lat2_case2 
       Kcoef=2.0   
       add_pure_rotation = .false.
    elseif(sub_case==3)then
       lon1=lon1_case3
       lat1=lat1_case3
       lon2=lon2_case3
       lat2=lat2_case3
       Kcoef=1.0   
       add_pure_rotation = .true.
!as subcase 2 but with rotation
    elseif(sub_case==4)then
       lon1=lon1_case2
       lat1=lat1_case2
       lon2=lon2_case2
       lat2=lat2_case2
       Kcoef=2.0   
       add_pure_rotation = .true.
    endif

    Tperiod=5*24*60*60

    Kcoef=Kcoef*rearth/24/60/60

    nm1= 1
    n0 = 2
    np1= 3

    do ie=nets,nete
       do k=1,nlev

          elem(ie)%state%p(:,:,k,n0)=swirl_init_tracer(elem(ie)%spherep(:,:),k)

          elem(ie)%state%v(:,:,:,k,n0)=swirl_velocity(0.0D0,elem(ie)%spherep(:,:),elem(ie)%Dinv)
          elem(ie)%state%v(:,:,:,k,nm1)=elem(ie)%state%v(:,:,:,k,n0)
          elem(ie)%state%v(:,:,:,k,np1)=elem(ie)%state%v(:,:,:,k,n0)
 
          elem(ie)%state%p(:,:,k,nm1)=elem(ie)%state%p(:,:,k,n0)
          elem(ie)%state%p(:,:,k,np1)=elem(ie)%state%p(:,:,k,n0)

       end do

       elem(ie)%state%ps(:,:)=0.0D0
       elem(ie)%state%gradps(:,:,:)=0.0D0
    end do


!temporary code to test limiters
!i still need it, OG
!     if (present(hybrid)) then
!       if(nlev>13)then
!       !initialize layers with smoothness metric
!       !this is a duplicate of whats in advance_nonstag_rk
! 	call neighbor_minmax(elem,hybrid,edge3,nets,nete,n0,pmin,pmax,min_var,max_var)
! 
! 	do ie=nets,nete
! 	    do k=1,3
! 	      local_var=maxval(elem(ie)%state%p(:,:,k,n0))-minval(elem(ie)%state%p(:,:,k,n0))
! 	      maxvar_neighb=max_var(k,ie)
! 	      minvar_neighb=min_var(k,ie)
! 	      value=(max(maxvar_neighb,local_var))/&
! 		  (min(minvar_neighb,local_var)+eps)
! 	      elem(ie)%state%p(:,:,k+7,:)=min(1000.0d0,value)
! 	    enddo
! 	enddo
!       endif
!     endif
!!!!!!!!!!!!!!!!!!!!!

  end subroutine swirl_init_state

  ! ===========================================
  !
  ! swirl_init_tracer:
  !
  ! Initializes many tracers
  !
  ! ===========================================

  recursive function swirl_init_tracer(sphere,level) result(p)

    type (spherical_polar_t), intent(in) :: sphere(np,np) ! spherical coordinates of element
    integer, intent(in) :: level
    real (kind=real_kind)                :: p(np,np)    ! geopotential
    real (kind=real_kind) :: snlat,cslat,lon,lat, lon_min, lon_max
    real (kind=real_kind) :: h1,h2
    real (kind=real_kind) :: r1,r2,A1,A2, bellrad, latc, lonc,xc,yc,zc,snlon,cslon,x1,y1,z1
    real (kind=real_kind) :: lon_1,lat_1,lon_2,lat_2
    integer i,j,k

! layers in swirl:
! lev=1: COS BELLS
! lev=2: GAUSS BELLS
! lev=3: SLOTTED CYLINDERS
! lev=4: CONST=1

!------------------------ levels below are more or less temporary
! lev=5: FLIPPED COS BELLS
! lev=6: TRACER 6 = -0.8 Tracer1^2+0.9, to test mixing by Thuburn&P paper
! lev=7: Tracer 6 + Tracer 2, to test if linear relations are preserved
! lev=8: metric for level 1
! lev=9: metric for level 2
! lev=10:  metric for level 3
! lev=11: temporary (right now is used to mark elements)
! lev=12: is filter applied to level 1 and where
! lev=13: is filter applied to level 2 and where
! lev=14: is filter applied to level 3 and where

    if(level==1)then
	bellrad=bellradius !bellradius=1/2
        !these values for true swirl test cases
	lon_1=lon1; lat_1=lat1;lon_2=lon2; lat_2=lat2;
        !these values below make bells apart by more than 1 element if NE=20
        !lon_1=3.0*DD_PI/4.0; lat_1=0.0;lon_2=5.0*DD_PI/4.0; lat_2=0.0;

	do j=1,np
	  do i=1,np
	    snlat=SIN(sphere(i,j)%lat)
	    cslat=COS(sphere(i,j)%lat)
	    lon  =sphere(i,j)%lon
	    lat  =sphere(i,j)%lat

	    A1=SIN(lat_1)*snlat + COS(lat_1)*cslat*COS(lon - lon_1)   
	    A2=SIN(lat_2)*snlat + COS(lat_2)*cslat*COS(lon - lon_2)
	    r1=ACOS(A1)
	    r2=ACOS(A2)

	    if (r1<bellrad) then
!temp code to fix lim8 for badly scaled
!	      h1 = 0.0d0+1e9*(1.0d0 + COS(DD_PI*r1/bellrad))
	      h1 = tracer_lowest+tracer_highest*(1.0d0 + COS(DD_PI*r1/bellrad)) 
	    else
!	      h1 = 0.0d0
	      h1 = tracer_lowest
	    endif
	    if (r2<bellrad) then
!	      h2 = 0.0d0+1e9*(1.0d0 + COS(DD_PI*r2/bellrad))
	      h2 = tracer_lowest+tracer_highest*(1.0d0 + COS(DD_PI*r2/bellrad)) 
	    else
!	      h2 = 0.0d0
	      h2 = tracer_lowest
	    endif

	    p(i,j)=(h1+h2)/2.0D0 

	  end do
	end do

    ! initializing tracer, GAUSS BELL (not like in swtc1(too small?) but like in N&P paper)
     elseif(level==2)then
	bellrad=5.0D0
	do j=1,np
	  do i=1,np

	      lon_1=lon1; lat_1=lat1;lon_2=lon2; lat_2=lat2;

	      snlat=SIN(sphere(i,j)%lat)
	      cslat=COS(sphere(i,j)%lat)
	      snlon=SIN(sphere(i,j)%lon)
	      cslon=COS(sphere(i,j)%lon)

              xc = cos(lat_1)*cos(lon_1)
              yc = cos(lat_1)*sin(lon_1)
              zc = sin(lat_1)
              x1 = cslat*cslon
              y1 = cslat*snlon
              z1 = snlat
              r1 = (x1-xc)**2 + (y1-yc)**2 + (z1-zc)**2
              h1 = 2.0D0*(0.95d0)*exp(-bellrad*r1)

              xc = cos(lat_2)*cos(lon_2)
              yc = cos(lat_2)*sin(lon_2)
              zc = sin(lat_2)
              r2 = (x1-xc)**2 + (y1-yc)**2 + (z1-zc)**2
              h2 = 2.0D0*(0.95d0)*exp(-bellrad*r2)

	      p(i,j)=(h1+h2)/2.0D0
	  end do
	end do

    !slotted cylinders
    elseif(level==3)then

        bellrad=0.5d0
	p(:,:)=tracer_lowest
	h2 = bellrad/6.0d0

	latc  = lat1
	lonc  = lon1  

	do j=1,np
	  do i=1,np
	    snlat=SIN(sphere(i,j)%lat)
	    cslat=COS(sphere(i,j)%lat)
	    lon  =sphere(i,j)%lon
	    lat  =sphere(i,j)%lat

	    A1    =SIN(latc)*snlat + COS(latc)*cslat*COS(lon - lonc)   

	    r1=ACOS(A1)

	    if ((r1<bellrad).and.(abs(lon-lonc)>=h2)) then
		p(i,j)=tracer_highest+tracer_lowest
	    endif
	    if ((r1<bellrad).and.(abs(lon-lonc)<h2).and.(lat-latc<-5.d0/12.d0*bellrad)) then
		p(i,j)=tracer_highest+tracer_lowest
	    endif
	  enddo
	enddo

	latc  = lat2
	lonc  = lon2 

	do j=1,np
	  do i=1,np
	    snlat=SIN(sphere(i,j)%lat)
	    cslat=COS(sphere(i,j)%lat)
	    lon  =sphere(i,j)%lon
	    lat  =sphere(i,j)%lat

	    A1    =SIN(latc)*snlat + COS(latc)*cslat*COS(lon - lonc)   

	    r1=ACOS(A1)

	    if ((r1<bellrad).and.(abs(lon-lonc)>=h2)) then
		p(i,j)=tracer_highest+tracer_lowest
	    endif
	    if ((r1<bellrad).and.(abs(lon-lonc)<h2).and.(lat-latc>5.d0/12.d0*bellrad)) then
		p(i,j)=tracer_highest+tracer_lowest
	    endif
	  enddo
	enddo

    !initializing density
    elseif(level==4)then
        p(:,:)=1.0D0
    elseif(level==5)then
	p(:,:)=1.0d0
	!flipped tracer
	!p(:,:)=tracer_highest+2.0d0*tracer_lowest-swirl_init_tracer(sphere,1)
    !tracer made of level1 to test mixing
    elseif(level==6)then
	p(:,:)=swirl_init_tracer(sphere,1)**2*(-0.8d0)+0.9d0
    !tracer to test if scheme preserve linear relationship
    elseif(level==7)then
	p(:,:)=swirl_init_tracer(sphere,6)+swirl_init_tracer(sphere,2)
    !this level is for observing metric of level 1
    elseif(level==8)then
	p(:,:)=0.0d0 
    !this level is for observing metric of level 2
    elseif(level==9)then
	p(:,:)=0.0d0 
    !this level is for observing metric of level 3
    elseif(level==10)then
	p(:,:)=0.0d0 
    !this is another level for debugging or for smoothness measure
    !or for int smooth cyllinder
    elseif(level==11)then
	do j=1,np
	  do i=1,np
	      !lon_1=lon1; lat_1=lat1;lon_2=lon2; lat_2=lat2;
	      lon_1=3.0*DD_PI/4.0; lat_1=0.0;lon_2=5.0*DD_PI/4.0; lat_2=0.0;
	      
	      snlat=SIN(sphere(i,j)%lat)
	      cslat=COS(sphere(i,j)%lat)
	      snlon=SIN(sphere(i,j)%lon)
	      cslon=COS(sphere(i,j)%lon)

              xc = cos(lat_2)*cos(lon_2)
              yc = cos(lat_2)*sin(lon_2)
              zc = sin(lat_2)
              x1 = cslat*cslon
              y1 = cslat*snlon
              z1 = snlat
              r1 = (x1-xc)**2 + (y1-yc)**2 + (z1-zc)**2
              h1 = 0.5d0+0.5d0*tanh(8.0d0*r1+3.d0)*tanh(3.d0-8.0d0*r1)

	      p(i,j)=h1
	  end do
	end do
    !this level is for observing limiter applied to level 1
    elseif(level==12)then
	p(:,:)=0.0d0 
    !this level is for observing limiter applied to level 2
    elseif(level==13)then
	p(:,:)=0.0d0 
    !this level is for observing limiter applied to level 3
    elseif(level==14)then
	p(:,:)=0.0d0 
    endif
  end function swirl_init_tracer


!--------------------------------------------------------------------
! velocities for test SWIRL
!
! modif by OG
!--------------------------------------------------------------------

  function swirl_velocity(t,sphere,D) result(v)

    type (spherical_polar_t), intent(in) :: sphere(np,np)
    real (kind=real_kind),    intent(in) :: t,D(2,2,np,np)
    real (kind=real_kind)                :: v(np,np,2)

    ! Local variables

    real (kind=real_kind) :: snhalflon(np,np)
    real (kind=real_kind) :: snlon(np,np)
    real (kind=real_kind) :: sndoublelat(np,np)
    real (kind=real_kind) :: sndoublelon(np,np)
    real (kind=real_kind) :: cslat(np,np)
    real (kind=real_kind) :: csperiod, PI
    real (kind=real_kind) :: V1,V2
    integer i,j,k

    real (kind=real_kind) :: sphlon,sphlat
    real (kind=real_kind) :: v1_old,v2_old,v1_new,v2_new,lonnew,latnew,shift_for_rotation
    real (kind=real_kind) :: vel_pure_rotation_marker

    if(add_pure_rotation)then
      shift_for_rotation=t*DD_PI/Tperiod
      vel_pure_rotation_marker=1.0D0
    else
      shift_for_rotation=0.0D0
      vel_pure_rotation_marker=0.0D0
    endif

    csperiod=COS(t*DD_PI/Tperiod)

    if(sub_case==1)then
      do j=1,np
	do i=1,np

	  !peter's code for adding zonal flow
	  !            tt = 5
	  !            ck = 2
	  !            slon  = sin(lon-(t)*pi2/tt)!solid-body rotation added
	  !            slon2 = sin(two*(lon-(t)*pi2/tt))!solid-body rotation added
	  !            clon  = cos(lon-(t)*pi2/tt)!solid-body rotation added
	  !            clon2 = cos(two*(lon-(t)*pi2/tt))!solid-body rotation added
	  !            sslm  = (sin(half*(lon-(t)*pi2/tt)))**2
	  ! 
	  !            uexact =  ck*slon*slon*SIN(two*lat)*COS(t*omega) + clat*pi2/(tt)
	  !            vexact =  ck*slon2*clat*COS(t*omega)
	  !            udc    =  two*ck*slon*slon*slat*COS(t*omega) + pi2/tt

	  sphlon=sphere(i,j)%lon-shift_for_rotation*2.0d0
	  sphlat=sphere(i,j)%lat
	  CALL regrot(sphlon,sphlat,lonnew,latnew,0.0d0,-0.5d0*DD_PI+rotangle,1)
	  snlon(i,j) = SIN(lonnew)
	  snhalflon(i,j) = SIN(lonnew/2.0D0)
	  cslat(i,j) = COS(latnew)
	  sndoublelat(i,j) = SIN(latnew*2.0D0)
	  sndoublelon(i,j) = SIN(lonnew*2.0D0)
	  do k=1,nlev
	      V1 =  Kcoef*(snhalflon(i,j))**2*sndoublelat(i,j)*csperiod &
		+cslat(i,j)*2.0D0*DD_PI*rearth/Tperiod*vel_pure_rotation_marker

	      V2 =  Kcoef/2.0D0*snlon(i,j)*cslat(i,j)*csperiod
	      v1_old=V1
	      v2_old=V2
	      CALL turnwi(v1_old,v2_old,&
	      v1_new,v2_new,sphlon,sphlat,lonnew,latnew,0.0d0,-0.5d0*DD_PI+rotangle,-1)
	      IF (ABS(DD_PI*0.5d0-sphlat)<1.0E-8.OR.ABS(DD_PI*0.5d0+sphlat)<1.0E-8) THEN
		v2_new = 0.0d0
	      endif
	      IF (ABS(v1_new)<1.0E-10) v1_new=0.0d0
	      V1=v1_new
	      V2=v2_new
	      v(i,j,1)= V1*D(1,1,i,j) + V2*D(1,2,i,j)
	      v(i,j,2)= V1*D(2,1,i,j) + V2*D(2,2,i,j)
	  end do
	end do
    end do

  elseif(sub_case==2)then
    do j=1,np
      do i=1,np
	sphlon=sphere(i,j)%lon-shift_for_rotation*2
	sphlat=sphere(i,j)%lat
	CALL regrot(sphlon,sphlat,lonnew,latnew,0.0d0,-0.5d0*DD_PI+rotangle,1)
	snlon(i,j) = SIN(lonnew)
	snhalflon(i,j) = SIN(lonnew/2.0D0)
	cslat(i,j) = COS(latnew)
	sndoublelat(i,j) = SIN(latnew*2.0D0)
	sndoublelon(i,j) = SIN(lonnew*2.0D0)
	do k=1,nlev
	    V1 = Kcoef*(snlon(i,j))**2*sndoublelat(i,j)*csperiod +&
	      cslat(i,j)*2.0D0*DD_PI*rearth/Tperiod*vel_pure_rotation_marker
	    V2 = Kcoef*sndoublelon(i,j)*cslat(i,j)*csperiod

	    v1_old=V1
	    v2_old=V2
	    CALL turnwi(v1_old,v2_old,&
	    v1_new,v2_new,sphlon,sphlat,lonnew,latnew,0.0d0,-0.5d0*DD_PI+rotangle,-1)
	    IF (ABS(DD_PI*0.5d0-sphlat)<1.0E-8.OR.ABS(DD_PI*0.5d0+sphlat)<1.0E-8) THEN
	      v2_new = 0.0d0
	    endif
	    IF (ABS(v1_new)<1.0E-10) v1_new=0.0d0
	    V1=v1_new
	    V2=v2_new
	    v(i,j,1)= V1*D(1,1,i,j) + V2*D(1,2,i,j)
	    v(i,j,2)= V1*D(2,1,i,j) + V2*D(2,2,i,j)
	  end do
	end do
      end do

    elseif(sub_case==3)then
      do j=1,np
	do i=1,np
	  sphlon=sphere(i,j)%lon-shift_for_rotation*2.0d0
	  sphlat=sphere(i,j)%lat
	  CALL regrot(sphlon,sphlat,lonnew,latnew,0.0d0,-0.5d0*DD_PI+rotangle,1)
	  snlon(i,j) = SIN(lonnew)
	  snhalflon(i,j) = SIN(lonnew/2.0D0)
	  cslat(i,j) = COS(latnew)
	  sndoublelat(i,j) = SIN(latnew*2.0D0)
	  sndoublelon(i,j) = SIN(lonnew*2.0D0)
	  do k=1,nlev
	    V1 =  -Kcoef*(snhalflon(i,j))**2*sndoublelat(i,j)*cslat(i,j)**2*csperiod &
		+cslat(i,j)*2.0D0*DD_PI*rearth/Tperiod*vel_pure_rotation_marker
	    V2 =  Kcoef/2.0D0*snlon(i,j)*cslat(i,j)**3*csperiod
	    v1_old=V1
	    v2_old=V2
	    CALL turnwi(v1_old,v2_old,&
	    v1_new,v2_new,sphlon,sphlat,lonnew,latnew,0.0d0,-0.5d0*DD_PI+rotangle,-1)
	    IF (ABS(DD_PI*0.5d0-sphlat)<1.0E-8.OR.ABS(DD_PI*0.5d0+sphlat)<1.0E-8) THEN
		v2_new = 0.0d0
	    endif
	    IF (ABS(v1_new)<1.0E-10) v1_new=0.0d0
	    V1=v1_new
	    V2=v2_new
	    v(i,j,1)= V1*D(1,1,i,j) + V2*D(1,2,i,j)
	    v(i,j,2)= V1*D(2,1,i,j) + V2*D(2,2,i,j)
	  end do
	end do
      end do

    elseif(sub_case==4)then
      do j=1,np
        do i=1,np
          sphlon=sphere(i,j)%lon-shift_for_rotation*2.0d0
          sphlat=sphere(i,j)%lat
          CALL regrot(sphlon,sphlat,lonnew,latnew,0.0d0,-0.5d0*DD_PI+rotangle,1)
          snlon(i,j) = SIN(lonnew)
          snhalflon(i,j) = SIN(lonnew/2.0D0)
          cslat(i,j) = COS(latnew)
          sndoublelat(i,j) = SIN(latnew*2.0D0)
          sndoublelon(i,j) = SIN(lonnew*2.0D0)
          do k=1,nlev
!true case 4
            V1 =  Kcoef*(snlon(i,j))**2*sndoublelat(i,j)*csperiod &
               +cslat(i,j)*2.0D0*DD_PI*rearth/Tperiod*vel_pure_rotation_marker
            V2 =  Kcoef*sndoublelon(i,j)*cslat(i,j)*csperiod

!experimental case 4 with flow 'away' from poles
!            V1 =  0.5*Kcoef*(snlon(i,j))**2*sndoublelat(i,j)*csperiod &
!               +cslat(i,j)*2.0D0*DD_PI*rearth/Tperiod*vel_pure_rotation_marker
!            V2 =  0.5*Kcoef*sndoublelon(i,j)*cslat(i,j)*csperiod

!only rotation to test limiters
!              V1 =  cslat(i,j)*2.0D0*DD_PI*rearth/Tperiod*vel_pure_rotation_marker
!              V2 =  0.0d0

             v1_old=V1
             v2_old=V2
             CALL turnwi(v1_old,v2_old,&
             v1_new,v2_new,sphlon,sphlat,lonnew,latnew,0.0d0,-0.5d0*DD_PI+rotangle,-1)
             IF (ABS(DD_PI*0.5d0-sphlat)<1.0E-8.OR.ABS(DD_PI*0.5d0+sphlat)<1.0E-8) THEN
                v2_new = 0.0d0
             endif
             IF (ABS(v1_new)<1.0E-10) v1_new=0.0d0
             V1=v1_new
             V2=v2_new
             v(i,j,1)= V1*D(1,1,i,j) + V2*D(1,2,i,j)
             v(i,j,2)= V1*D(2,1,i,j) + V2*D(2,2,i,j)
          end do
       end do
      end do
    endif

  end function swirl_velocity

  ! ===================================================================
  !
  ! sj1_init_state:
  !
  ! Initialize strong jet case 1: Galewski et al.
  !
  ! ===================================================================
 
  subroutine sj1_init_state(elem,nets,nete,hybrid,pmean,deriv)

    type(element_t), target,intent(inout) :: elem(:)
    integer, intent(in)   :: nets
    integer, intent(in)   :: nete
    type (hybrid_t)       :: hybrid
    real (kind=real_kind) :: pmean
    type (derivative_t)   :: deriv

    ! Local variables
  
    real (kind=real_kind) :: a,b
    real (kind=real_kind) :: mean_balance
    real (kind=real_kind) :: pert(np,np)
    real (kind=real_kind) :: phi_balance(np,np,nets:nete)
    type (quadrature_t)   :: gs

    integer :: ie,k
    integer :: nm1 
    integer :: n0 
    integer :: np1
    type (element_t), pointer :: pElem,sElem

    logical :: Debug=.FALSE.

    nm1= 1
    n0 = 2
    np1= 3

    a =-DD_PI/2.0D0
    b = DD_PI/2.0D0

    gs=gauss(Ngs_sj1)
    Ibalance_tot = gaussian_int(sj1_integrand,a,b,gs)    
    print *,"Ibalance_tot=",Ibalance_tot
    deallocate(gs%points)
    deallocate(gs%weights)
    print *,"done..."

    do ie=nets,nete
       pElem => elem(ie)
       phi_balance(:,:,ie) = sj1_balance(pElem%spherep(:,:),np)
    end do
    pmean    = g*h0_sj1

    gh0_sj1  = pmean

    mean_balance  = global_integral(elem,phi_balance(:,:,nets:nete),hybrid,np,nets,nete)

    if(Debug)print *,"BEFORE mean height=",mean_balance/g
    if(Debug)print *,"BEFORE pmean=",pmean/g

    pmean = gh0_sj1 + mean_balance
    gh0_sj1 = pmean

    mean_balance  = pmean-global_integral(elem,phi_balance(:,:,nets:nete),hybrid,np,nets,nete)

    if(Debug)print *,"AFTER mean height=",mean_balance/g
    if(Debug)print *,"AFTER pmean=",pmean/g

    do ie=nets,nete
       pElem => elem(ie)

       pert(:,:)           = sj1_perturbation(pElem%spherep(:,:),np)
       pElem%state%ps(:,:) = 0.0D0
       do k=1,nlev
          pElem%state%p(:,:,k,n0) =pmean - phi_balance(:,:,ie)+ g*pert(:,:)
          pElem%state%p(:,:,k,nm1)=pElem%state%p(:,:,k,n0)
          pElem%state%p(:,:,k,np1)=pElem%state%p(:,:,k,n0)

          pElem%state%v(:,:,:,k,n0)=sj1_velocity_cubedsphere(pElem%spherep(:,:),pElem%Dinv)
          pElem%state%v(:,:,:,k,nm1)=pElem%state%v(:,:,:,k,n0)
          pElem%state%v(:,:,:,k,np1)=pElem%state%v(:,:,:,k,n0)
       end do
       pElem%state%gradps(:,:,:) = 0.0D0
    end do
    pmean=0.0D0
  end subroutine sj1_init_state

  function sj1_velocity(lat) result(ulat)
    
    real (kind=real_kind),intent(in) :: lat
    real (kind=real_kind)            :: ulat,en

    en    = exp(-((2.0D0/(lat1_sj1 - lat0_sj1))**2.0D0))
        
    if((lat0_sj1.lt.lat).AND.(lat.lt.lat1_sj1))then
       ulat  = (u0_sj1/en)*exp(1.0D0/((lat-lat0_sj1)*(lat-lat1_sj1)))
    else
       ulat  = 0.0_real_kind
    endif
    
  end function sj1_velocity

! =======================================================
!
! sj1_integrand:
! 
! Scaled integrand function needed to integrate the 
! geostrophic balance equation using the trapezoid rule.
! Integrand in units of u0_sj1*f*rearth (m^2/s^2)
!
! =======================================================

  function sj1_integrand(lat) result(Ival)

    real (kind=real_kind), intent(in) :: lat
    real (kind=real_kind) :: Ival
       
    real (kind=real_kind) :: ulat    ! dimensionless velocity

    ulat = sj1_velocity(lat)
    Ival = rearth*ulat*(2.0D0*omega*sin(lat) + tan(lat)*ulat/rearth)

  end function sj1_integrand

! ==========================================================
!
! sj1_balance:
!
! Compute geostrophic integral balance equation for the
! Barotropic Instability Test Case of Polvani.
!
! ===========================================================
 
  function sj1_balance(sphere,npts) result(p)

    integer, intent(in)                  :: npts
    type (spherical_polar_t), intent(in) :: sphere(npts,npts)
    real (kind=real_kind)                :: p(npts,npts)

    ! Local variables

    real (kind=real_kind) :: lat
    real (kind=real_kind) :: lon

    real (kind=real_kind) :: a, b
    real (kind=real_kind) :: Ibalance
    real (kind=real_kind) :: f
    real (kind=real_kind) :: pi2  ! pi/2

    integer i,j
    type (quadrature_t) :: gs

    gs=gauss(Ngs_sj1)

    do j=1,npts
       do i=1,npts

          lat   = sphere(i,j)%lat
          lon   = sphere(i,j)%lon
          pi2   = 0.5D0*DD_PI

          b = lat
          a = -pi2
          Ibalance = gaussian_int(sj1_integrand,a,b,gs)    
          p(i,j)  =  Ibalance

       end do
    end do
    
    deallocate(gs%weights)
    deallocate(gs%points)

  end function sj1_balance

! ===========================================
!
! sj1_perturbation:
!
! Compute the geopotential perturbation 
! from geostrophic balance for the Barotropic 
! Instability Test Case of Polvani.
!
! ===========================================
 
  function sj1_perturbation(sphere,npts) result(p)

    integer, intent(in)                  :: npts
    type (spherical_polar_t), intent(in) :: sphere(npts,npts)
    real (kind=real_kind)                :: p(npts,npts)

    ! Local variables

    real (kind=real_kind) :: lat
    real (kind=real_kind) :: lon
    real (kind=real_kind) :: alon, blat
    real (kind=real_kind) :: secha,sechb
    real (kind=real_kind) :: ghprime

    integer i,j

    do j=1,npts
       do i=1,npts

          lat   = sphere(i,j)%lat
          lon   = sphere(i,j)%lon
          if(lon.ge.DD_PI)then
             lon = lon -2.0D0*DD_PI
          end if

          p(i,j) = hhat_sj1*cos(lat)*exp(-((lon*alpha_sj1)**2.0D0 + ((lat2_sj1 - lat)*beta_sj1)**2.0D0 ))
       end do
    end do

  end function sj1_perturbation

! ===========================================
!
! sj1_velocity:
!
! Set initial zonal velocity field for barotropic 
! instability test case of Polvani.
!
! ===========================================

  function sj1_velocity_cubedsphere(sphere,D) result(v) 

    type (spherical_polar_t), intent(in) :: sphere(np,np)
    real (kind=real_kind),    intent(in) :: D(2,2,np,np)
    real (kind=real_kind)                :: v(np,np,2)

    ! Local variables

    real (kind=real_kind) :: lat
    real (kind=real_kind) :: V1
        
    integer i,j

    do j=1,np
       do i=1,np
          lat = sphere(i,j)%lat        
          V1  = sj1_velocity(lat)
          
          ! Map onto contravariant velocities
          v(i,j,1)= V1*D(1,1,i,j)
          v(i,j,2)= V1*D(2,1,i,j)
       end do
    end do

  end function sj1_velocity_cubedsphere


end module shallow_water_mod
