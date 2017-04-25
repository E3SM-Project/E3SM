



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
use physical_constants,   only: p0, g, Rgas

implicit none

real(rl), dimension(:,:,:,:), allocatable :: u0, v0                     ! storage for dcmip2-x sponge layer
real(rl):: zi(nlevp), zm(nlev)                                          ! z coordinates
real(rl):: ddn_hyai(nlevp), ddn_hybi(nlevp)                             ! vertical derivativess of hybrid coefficients
real(rl):: tau
real(rl):: ztop
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

  integer,  parameter :: zcoords = 0                                    ! we are not using z coords
  real(rl), parameter :: ztop    = 30000_rl                             ! top of model at 30km

  integer :: i,j,k,ie                                                   ! loop indices
  real(rl):: lon,lat,ntop                                               ! pointwise coordiantes
  real(rl):: p,z,u,v,w,T,thetav,phis,ps,rho,q(1),dp                     ! pointwise field values

  if (hybrid%masterthread) write(iulog,*) 'initializing dcmip2016 test 2: tropical cyclone'

  ! get pressure at ztop
  lon=0; lat=0; z=ztop
  call tropical_cyclone_test(lon,lat,p,z,1,u,v,t,thetav,phis,ps,rho,q(1))
  ntop = p/p0
  if (hybrid%masterthread) print *,"ntop = ",ntop

  ! get evenly spaced eta levels from 1 to ntop
  forall(k=1:nlevp) hvcoord%etai(k) = ntop - (k-1)*(ntop-1)/(nlevp-1)
  call set_hybrid_coefficients(hvcoord,hybrid,  hvcoord%etai(1), 1.0_rl)
  call set_layer_locations(hvcoord, .true., hybrid%masterthread)

  ! set initial conditions
  do ie = nets,nete
     do k=1,nlev; do j=1,np; do i=1,np

        lon = elem(ie)%spherep(i,j)%lon
        lat = elem(ie)%spherep(i,j)%lat
        p   = p0*hvcoord%hyam(k) + ps*hvcoord%hybm(k)

        call tropical_cyclone_test(lon,lat,p,z,zcoords,u,v,T,thetav,phis,ps,rho,q(1))
        dp = pressure_thickness(ps,k,hvcoord)
        w  = 0

        call set_state(u,v,w,T,ps,phis,p,dp,z,g,i,j,k,elem(ie),1,nt)
        call set_tracers(q,1, dp,i,j,k,lat,lon,elem(ie))
     enddo; enddo; enddo;

  !call tests_finalize(elem(ie),hvcoord,1,nt)
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
  real(rl):: p,dp, z,phis,u,v,w,T,thetav,phis_ps,ps,rho,q(qsize_d)            ! pointwise field values
  real(rl):: p_i1,p_i2

  if (hybrid%masterthread) write(iulog,*) 'initializing dcmip2016 test 3: supercell storm'

  ! initialize hydrostatic state
  call supercell_init()

  ! get evenly spaced z coords from 0 to 20km
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

    do k=1,nlev;

      do j=1,np; do i=1,np
        lon = elem(ie)%spherep(i,j)%lon
        lat = elem(ie)%spherep(i,j)%lat


        if(zcoords==0) then
          ! set initial conditions at const eta levels

          ! get surface pressure at lat, lon, z=0
          z=0; call supercell_test(lon,lat,p,z,1,u,v,T,thetav,ps,rho,q(1),pert)

          ! get hydrostatic pressure level
          p  =  p0*hvcoord%hyam(k) + ps*hvcoord%hybm(k)
          dp = pressure_thickness(ps,k,hvcoord)

        endif

        if(zcoords==1) then
          ! set initial conditions at const z levels
          call supercell_z(lon, lat, zi(k)  , p_i1, thetav, rho, q(1), pert)
          call supercell_z(lon, lat, zi(k+1), p_i2, thetav, rho, q(1), pert)
          dp = p_i2-p_i1
          z = zm(k)
        endif

        call supercell_test(lon,lat,p,z,zcoords,u,v,T,thetav,ps,rho,q(1),pert)

        w    = 0 ! no vertical motion
        phis = 0 ! no topography
        q(2) = 0 ! no initial clouds
        q(3) = 0 ! no initial rain
        q(4) = 0 ! precip

        ! store fields in q to debug ICs
        ! q(2) = thetav; q(3) = T; q(4) = p; q(5) = u

        call set_state(u,v,w,T,ps,phis,p,dp,z,g,i,j,k,elem(ie),1,nt)
        call set_tracers(q,qsize_d, dp,i,j,k,lat,lon,elem(ie))
      enddo; enddo
    enddo
    call tests_finalize(elem(ie),hvcoord,1,nt)

  enddo

end subroutine

!_______________________________________________________________________
subroutine dcmip2016_test3_forcing(elem,hybrid,hvcoord,nets,nete,nt,ntQ,dt,test,prec_type,pbl_type)

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(in)            :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  integer,            intent(in)            :: nt, ntQ                  ! time level index
  real(rl),           intent(in)            :: dt                       ! time-step size
  integer,            intent(in)            :: test                     ! dcmip16 test number
  integer,            intent(in)            :: prec_type                ! precipitation type
  integer,            intent(in)            :: pbl_type                 ! planetary boundary layer type

  integer :: i,j,k,ie                                                     ! loop indices
  real(rl):: lat
  real(rl), dimension(np,np,nlev) :: u,v,w,T,theta,exner,p,dp,qv,qc,qr,rho,z,zi
  real(rl), dimension(np,np,nlev) :: u0,v0,T0,qv0,qc0,qr0
  real(rl), dimension(np,np)      :: ps,phis,precl

  do ie = nets,nete

    call get_state(u,v,w,T,theta,exner,p,dp,z,g,i,j,elem(ie),hvcoord,nt,ntQ)
    rho = p/(Rgas*T)
    zi  = z ! todo: get z levels at interfaces?
    qv  = elem(ie)%state%Qdp(:,:,:,1,ntQ)/dp
    qc  = elem(ie)%state%Qdp(:,:,:,2,ntQ)/dp
    qr  = elem(ie)%state%Qdp(:,:,:,3,ntQ)/dp

    u0=u; v0=v; T0=T; qv0=qv; qc0=qc; qr0=qr

    do j=1,np; do i=1,np

      lat = elem(ie)%spherep(i,j)%lat

      !call DCMIP2016_PHYSICS(test, u(i,j,:), v(i,j,:), p(i,j,:),&
      !  qv(i,j,:), qc(i,j,:), qr(i,j,:), rho(i,j,:), dt, z(i,j,:), zi(i,j,:), &
      !  lat, nlev, precl, pbl_type, prec_type)

      CALL KESSLER(        &
      theta(i,j,:),        &
      qv(i,j,:),           &
      qc(i,j,:),           &
      qr(i,j,:),           &
      rho(i,j,:),          &
      exner(i,j,:),        &
      dt,                  &
      z(i,j,:),            &
      nlev,                  &
      precl(i,j))

    enddo; enddo;

    T = theta*exner

    !print *,"u=",u(1)," v=",v(1), " w= ",w(1), "T=",t(1)," p=",p(1)," dp=",dp(1)," z=",z(1)," qv=",qv(1)," qc=",qc(1)," qr=",qr(1)," precl=",precl

    elem(ie)%derived%FM(:,:,1,:) = (u - u0)/dt
    elem(ie)%derived%FM(:,:,2,:) = (v - v0)/dt
    elem(ie)%derived%FM(:,:,3,:) = 0
    elem(ie)%derived%FT(:,:,:)   = (T - T0)/dt
    elem(ie)%derived%FQ(:,:,:,1) = (qv-qv0)/dt
    elem(ie)%derived%FQ(:,:,:,2) = (qc-qc0)/dt
    elem(ie)%derived%FQ(:,:,:,3) = (qr-qr0)/dt
!print *,"T-T0= ",T-T0
    elem(ie)%state%Qdp(:,:,1,4,ntQ) = precl*dp(:,:,1) ! store precl in level 1 of tracer #4
  enddo

end subroutine


end module dcmip16_wrapper
#endif
