



#ifndef CAM
#include "config.h"

module dcmip16_wrapper

! Implementation of the dcmip2012 dycore tests for the preqx dynamics target

use dcmip12_wrapper,      only: pressure_thickness, set_tracers, get_evenly_spaced_z, set_hybrid_coefficients
use control_mod,          only: test_case, dcmip4_moist, dcmip4_X
use baroclinic_wave,      only: baroclinic_wave_test
use supercell,            only: supercell_init, supercell_test, supercell_z
use tropical_cyclone,     only: tropical_cyclone_test
use derivative_mod,       only: derivative_t, gradient_sphere
use dimensions_mod,       only: np, nlev, nlevp , qsize, qsize_d, nelemd
use element_mod,          only: element_t
use element_state,        only: nt=>timelevels
use hybrid_mod,           only: hybrid_t
use hybvcoord_mod,        only: hvcoord_t, set_layer_locations
use kinds,                only: rl=>real_kind, iulog
use parallel_mod,         only: abortmp
use element_ops,          only: set_state, get_state, tests_finalize, set_forcing_rayleigh_friction, set_thermostate
use physical_constants,   only: p0, g, Rgas, kappa, Cp, Rwater_vapor

implicit none

real(rl), dimension(:,:,:,:), allocatable :: u0, v0                     ! storage for dcmip2-x sponge layer
real(rl):: zi(nlevp), zm(nlev)                                          ! z coordinates
real(rl):: ddn_hyai(nlevp), ddn_hybi(nlevp)                             ! vertical derivativess of hybrid coefficients
real(rl):: tau
real(rl):: ztop

real(rl), parameter ::rh2o    = 461.5d0,            & ! Gas constant for water vapor (J/kg/K)
                      Mvap    = (Rwater_vapor/Rgas) - 1.d0    ! Constant for virtual temp. calc. (~0.608)
contains

!_____________________________________________________________________
subroutine dcmip2016_test1(elem,hybrid,hvcoord,nets,nete)

  ! moist baroclinic wave test

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index

  integer,  parameter :: zcoords = 0                                    ! use vertical pressure coordinates
  integer,  parameter :: deep    = 0                                    ! use shallow atmosphere approximation
  integer,  parameter :: pertt   = 0                                    ! use exponential perturbation type
  integer,  parameter :: moist   = 1                                    ! use moist version
  real(rl), parameter :: dcmip_X = 1.0_rl                               ! full scale planet

  integer :: i,j,k,ie                                                   ! loop indices
  real(rl):: lon,lat                                                    ! pointwise coordiantes
  real(rl):: p,z,phis,u,v,w,T,thetav, phis_ps,ps,rho,rho_mean,q(1),dp    ! pointwise field values

  if (hybrid%masterthread) write(iulog,*) 'initializing dcmip2016 test 1: moist baroclinic wave'

  ! set initial conditions
  do ie = nets,nete
    do k=1,nlev;

      ! no surface topography
      p  =  p0*hvcoord%etam(k)
      dp = (hvcoord%etai(k+1)-hvcoord%etai(k))*p0

      do j=1,np; do i=1,np
        lon = elem(ie)%spherep(i,j)%lon
        lat = elem(ie)%spherep(i,j)%lat

        call baroclinic_wave_test(deep,moist,pertt,dcmip_X,lon,lat,p,z,zcoords,u,v,T,thetav,phis,ps,rho,q(1))

        call set_state(u,v,w,T,ps,phis,p,dp,zm(k),g, i,j,k,elem(ie),1,nt)
        call set_tracers(q,1, dp,i,j,k,lat,lon,elem(ie))

      enddo; enddo
    enddo
  enddo

end subroutine

!_____________________________________________________________________
subroutine dcmip2016_test2(elem,hybrid,hvcoord,nets,nete)

  ! tropical cyclone test

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index

  real(rl), parameter :: ztop    = 30000_rl                             ! top of model at 30km

  integer :: i,j,k,ie , ierr                                                  ! loop indices
  real(rl):: lon,lat,ntop                                               ! pointwise coordiantes
  real(rl):: p,z,u,v,w,T,thetav,phis,ps,rho,q(4),dp                     ! pointwise field values

  if (hybrid%masterthread) write(iulog,*) 'initializing dcmip2016 test 2: tropical cyclone'
  !use vertical levels specificed in cam30 file

  ! set initial conditions
  do ie = nets,nete
     do k=1,nlev; do j=1,np; do i=1,np

        ! get surface pressure
        lon = elem(ie)%spherep(i,j)%lon
        lat = elem(ie)%spherep(i,j)%lat
        z=0; call tropical_cyclone_test(lon,lat,p,z,1,u,v,T,thetav,phis,ps,rho,q(1))

        ! get pressure at level midpoints
        p = p0*hvcoord%hyam(k) + ps*hvcoord%hybm(k)

        ! get initial conditions at pressure level p
        call tropical_cyclone_test(lon,lat,p,z,0,u,v,T,thetav,phis,ps,rho,q(1))

        dp = pressure_thickness(ps,k,hvcoord)
        w  = 0
        q(2)=0
        q(3)=0
        q(4)=0

        call set_state(u,v,w,T,ps,phis,p,dp,z,g,i,j,k,elem(ie),1,nt)
        call set_tracers(q,4,dp,i,j,k,lat,lon,elem(ie))
     enddo; enddo; enddo;

    call tests_finalize(elem(ie),hvcoord,1,nt)
  enddo

end subroutine

!_____________________________________________________________________
subroutine dcmip2016_test3(elem,hybrid,hvcoord,nets,nete)

  ! supercell storm test case

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index

  integer,  parameter :: zcoords  = 0                                   ! 0 -> use p coords
  integer,  parameter :: pert     = 1                                   ! 1 -> add thermal perturbation
  real(rl), parameter :: ztop     = 20000_rl                            ! top of model at 20km

  integer :: i,j,k,ie                                                   ! loop indices
  real(rl):: lon,lat                                                    ! pointwise coordiantes
  real(rl):: p,dp, z,phis,u,v,w,T,thetav,phis_ps,ps,rho,rhom,q(qsize_d)            ! pointwise field values
  real(rl):: p_i1,p_i2

  if (hybrid%masterthread) write(iulog,*) 'initializing dcmip2016 test 3: supercell storm'

  ! initialize hydrostatic state
  call supercell_init()

  ! get evenly spaced z levels from 0 to 20km
  call get_evenly_spaced_z(zi,zm, 0.0_rl, ztop)                          ! get evenly spaced z levels

  ! get eta levels matching z at equator
  do k=1,nlevp
    lon =0; lat = 0;
    call supercell_z(lon, lat, zi(k), p, thetav, rho, q(1), pert)
    hvcoord%etai(k) = p/p0
    if (hybrid%masterthread) print *,"k=",k," z=",zi(k)," p=",p,"etai=",hvcoord%etai(k)
  enddo
  call set_hybrid_coefficients(hvcoord,hybrid,  hvcoord%etai(1), 1.0_rl)
  call set_layer_locations(hvcoord, .true., hybrid%masterthread)

  ! set initial conditions
  do ie = nets,nete
  if (hybrid%masterthread) write(*,"(A,I5,A)",advance="NO") " ie=",ie,achar(13)

    do k=1,nlev

      do j=1,np; do i=1,np
        lon = elem(ie)%spherep(i,j)%lon
        lat = elem(ie)%spherep(i,j)%lat

        ! get surface pressure at lat, lon, z=0
        z=0; call supercell_test(lon,lat,p,z,1,u,v,T,thetav,ps,rho,q(1),pert)

        ! get hydrostatic pressure at level k
        p  =  p0*hvcoord%hyam(k) + ps*hvcoord%hybm(k)
        dp = pressure_thickness(ps,k,hvcoord)

        call supercell_test(lon,lat,p,z,zcoords,u,v,T,thetav,ps,rhom,q(1),pert)

        w    = 0 ! no vertical motion
        phis = 0 ! no topography
        q(2) = 0 ! no initial clouds
        q(3) = 0 ! no initial rain
        q(4) = 0 ! precip

        ! convert virtual temp to dry temp (?)
        !T = T * (1.0d0+q(1)) ! * (1.0d0+Mvap*q(1))

        call set_state(u,v,w,T,ps,phis,p,dp,z,g,i,j,k,elem(ie),1,nt)
        call set_tracers(q,qsize_d, dp,i,j,k,lat,lon,elem(ie))
      enddo; enddo
    enddo

    ! set density (dphi) to ensure hydrostatic balance
    call tests_finalize(elem(ie),hvcoord,1,nt)

  enddo

end subroutine


!_______________________________________________________________________
subroutine dcmip2016_forcing(elem,hybrid,hvcoord,nets,nete,nt,ntQ,dt,test)

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(in)            :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  integer,            intent(in)            :: nt, ntQ                  ! time level index
  integer,            intent(in)            :: test                     ! dcmip2016 test number
  real(rl),           intent(in)            :: dt                       ! time-step size

  integer   :: pbl_type=-1, prec_type=0

  integer :: i,j,k,ie                                                     ! loop indices
  real(rl):: lat
  real(rl), dimension(np,np,nlev) :: u,v,w,T,theta,exner,p,dp,cp_star,rho,z,qv,qc,qr
  real(rl), dimension(np,np,nlev) :: T0,qv0,qc0,qr0
  real(rl), dimension(np,np)      :: precl
  real(rl), dimension(np,np,nlev) :: p_kess, theta_kess, exner_kess
  real(rl), dimension(np,np,nlev) :: theta_inv,qv_inv,qc_inv,qr_inv,rho_inv,exner_inv,z_inv ! inverted columns

  do ie = nets,nete

    ! get current element state
    call get_state(u,v,w,T,theta,exner,p,dp,cp_star,rho,z,g,i,j,elem(ie),hvcoord,nt,ntQ)

    ! get mixing ratios
    qv  = elem(ie)%state%Qdp(:,:,:,1,ntQ)/dp
    qc  = elem(ie)%state%Qdp(:,:,:,2,ntQ)/dp
    qr  = elem(ie)%state%Qdp(:,:,:,3,ntQ)/dp

    ! save un-forced prognostics
    T0=T; qv0=qv; qc0=qc; qr0=qr

    ! given rho,T,qv: get quantities consistent with Kessler's equation of state
    !p_kess      = rho * Rgas * T*(1.0_rl + Mvap * qv)
    !exner_kess  = (p_kess/p0)**(Rgas/Cp)
    !theta_kess = T/exner_kess

    ! invert columns (increasing z)
    !theta_inv= theta_kess(:,:,nlev:1:-1)
theta_inv= theta(:,:,nlev:1:-1)

    qv_inv   = qv        (:,:,nlev:1:-1)
    qc_inv   = qc        (:,:,nlev:1:-1)
    qr_inv   = qr        (:,:,nlev:1:-1)
    rho_inv  = rho       (:,:,nlev:1:-1)
!    exner_inv= exner_kess(:,:,nlev:1:-1)
exner_inv= exner(:,:,nlev:1:-1)

    z_inv    = z         (:,:,nlev:1:-1)

    ! apply forcing to columns
    do j=1,np; do i=1,np

      ! apply physics to (u,v,p, qv,qc,qr). (rho held constant)
      CALL KESSLER(       &
        theta_inv(i,j,:), &
        qv_inv(i,j,:),    &
        qc_inv(i,j,:),    &
        qr_inv(i,j,:),    &
        rho_inv(i,j,:),   &
        exner_inv(i,j,:), &
        dt,               &
        z_inv(i,j,:),     &
        nlev,             &
        precl(i,j))

    enddo; enddo;

    ! revert columns (increasing eta)
!    theta_kess= theta_inv(:,:,nlev:1:-1)
theta= theta_inv(:,:,nlev:1:-1)

    qv        = qv_inv   (:,:,nlev:1:-1)
    qc        = qc_inv   (:,:,nlev:1:-1)
    qr        = qr_inv   (:,:,nlev:1:-1)

    ! convert theta back to T using Kessler's EOS
    !exner_kess = (rho*Rgas*theta_kess*(1.0_rl + Mvap * qv)/p0)**(Rgas/(cp-Rgas))
    !T = theta_kess*exner_kess


elem(ie)%state%theta_dp_cp(:,:,:,nt) = theta*(Cp_star*dp)

    ! set dynamics forcing
    elem(ie)%derived%FM(:,:,1,:) = 0
    elem(ie)%derived%FM(:,:,2,:) = 0
    !elem(ie)%derived%FT(:,:,:)   = (T - T0)/dt

    ! set tracer-mass forcing
    elem(ie)%derived%FQ(:,:,:,1) = dp*(qv-qv0)/dt
    elem(ie)%derived%FQ(:,:,:,2) = dp*(qc-qc0)/dt
    elem(ie)%derived%FQ(:,:,:,3) = dp*(qr-qr0)/dt

    !elem(ie)%state%Qdp(:,:,1,4,ntQ) = precl*dp(:,:,1) ! store precl in level 1 of tracer #4

  enddo

end subroutine

#if 0
!_______________________________________________________________________
subroutine dcmip2016_forcing(elem,hybrid,hvcoord,nets,nete,nt,ntQ,dt,test)

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(in)            :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  integer,            intent(in)            :: nt, ntQ                  ! time level index
  integer,            intent(in)            :: test                     ! dcmip2016 test number
  real(rl),           intent(in)            :: dt                       ! time-step size

  integer   :: pbl_type=-1, prec_type=0

  integer :: i,j,k,ie                                                     ! loop indices
  real(rl):: lat
  real(rl), dimension(np,np,nlev) :: u,v,w,T,theta,thetav,exner,p,dp,cp_star,qv,qc,qr,rho,z
  real(rl), dimension(np,np,nlevp) :: zi

  real(rl), dimension(np,np,nlev) :: u0,v0,T0,qv0,qc0,qr0,  pm
  real(rl), dimension(np,np)      :: precl

  real(rl), dimension(nlev ) :: u1,v1,p1,qv1,qc1,qr1,rho1,z1
  real(rl), dimension(nlevp) :: zi1

  ! set boundary layer and precitipation types
  select case(test)
    case(1);      pbl_type= 0; prec_type=0;
    case(2);      pbl_type= 0; prec_type=0;
    case(3);      pbl_type=-1; prec_type=0;
    case default; stop
  end select

  do ie = nets,nete

    ! get current element state
    call get_state(u,v,w,T,theta,exner,p,dp,cp_star,z,g,i,j,elem(ie),hvcoord,nt,ntQ)

    ! get mixing ratios
    qv  = elem(ie)%state%Qdp(:,:,:,1,ntQ)/dp
    qc  = elem(ie)%state%Qdp(:,:,:,2,ntQ)/dp
    qr  = elem(ie)%state%Qdp(:,:,:,3,ntQ)/dp

    ! get rho from p and T
    rho = p/(T*Rgas)

    ! save un-forced fields
    u0=u; v0=v; T0=T; qv0=qv; qc0=qc; qr0=qr

    ! get layer interfaces from layer midpoints
    zi(:,:,nlevp) = 0.0
    zi(:,:,2:nlev)= ( z(:,:,1:nlev-1)+z(:,:,2:nlev) )/2
    zi(:,:,1)     = zi(:,:,2) + (z(:,:,1)-z(:,:,2))

!pm = p*(1.0d0+qv)*(1.0d0+Mvap*qv)

    ! apply dcmip forcing to columns
    do j=1,np; do i=1,np

      lat = elem(ie)%spherep(i,j)%lat

      ! invert columns (increasing z)
      u1    = u   (i,j,nlev:1:-1)
      v1    = v   (i,j,nlev:1:-1)
      p1    = p   (i,j,nlev:1:-1)
!p1    = pm   (i,j,nlev:1:-1)

      qv1   = qv  (i,j,nlev:1:-1)
      qc1   = qc  (i,j,nlev:1:-1)
      qr1   = qr  (i,j,nlev:1:-1)
      rho1  = rho (i,j,nlev:1:-1)
      z1    = z   (i,j,nlev:1:-1)
      zi1   = zi  (i,j,nlevp:1:-1)

      ! apply physics to (u,v,p, qv,qc,qr). (rho held constant)
      call DCMIP2016_PHYSICS(test, u1,v1,p1, qv1,qc1,qr1, rho1,dt,z1,zi1,lat,nlev,precl,pbl_type,prec_type)

      ! revert columns (increasing eta)
      u   (i,j,:) = u1 (nlev:1:-1)
      v   (i,j,:) = v1 (nlev:1:-1)
      p   (i,j,:) = p1 (nlev:1:-1)
!pm   (i,j,:) = p1 (nlev:1:-1)
      qv  (i,j,:) = qv1(nlev:1:-1)
      qc  (i,j,:) = qc1(nlev:1:-1)
      qr  (i,j,:) = qr1(nlev:1:-1)

    enddo; enddo;

!p = pm/((1.0d0+qv)*(1.0d0+Mvap*qv))

    ! get T from p and rho
    T = p/(rho*Rgas)

    ! get dynamics forcing
    elem(ie)%derived%FM(:,:,1,:) = (u - u0)/dt
    elem(ie)%derived%FM(:,:,2,:) = (v - v0)/dt
    elem(ie)%derived%FT(:,:,:)   = (T - T0)/dt

    ! get tracer-mass forcing
    elem(ie)%derived%FQ(:,:,:,1) = dp*(qv-qv0)/dt
    elem(ie)%derived%FQ(:,:,:,2) = dp*(qc-qc0)/dt
    elem(ie)%derived%FQ(:,:,:,3) = dp*(qr-qr0)/dt

    !elem(ie)%state%Qdp(:,:,1,4,ntQ) = precl*dp(:,:,1) ! store precl in level 1 of tracer #4

  enddo

end subroutine
#endif

end module dcmip16_wrapper
#endif
