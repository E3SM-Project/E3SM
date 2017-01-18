!
! DCMIP Tests for the PESE Target
!_______________________________________________________________________
#ifndef CAM

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module dcmip_tests

! Implementation of the dcmip2012 dycore tests for the PESE dynamics target

use control_mod,          only: test_case
use dcmip2012_test1_2_3,  only: test1_advection_deformation, test1_advection_hadley, test1_advection_orography, &
                                test2_steady_state_mountain, test2_schaer_mountain,test3_gravity_wave
!use tropical_cyclone,     only: tropical_cyclone_test
use derivative_mod,       only: derivative_t, gradient_sphere
use dimensions_mod,       only: np, nlev, nlevp, qsize, qsize_d, nelemd
use element_mod,          only: element_t
use element_state,        only: elem_state_t, derived_state_t, nt=>timelevels
use hybrid_mod,           only: hybrid_t
use hybvcoord_mod,        only: hvcoord_t, set_layer_locations
use kinds,                only: rl=>real_kind, iulog
use parallel_mod,         only: abortmp
use vertical_se,          only: init_vertical_interp, ddn_hyai,ddn_hybi
implicit none

! physical constants used by dcmip2012 test 3.1
real(rl), parameter ::              &
  g       = 9.80616,                & ! grav const
  a       = 6371229.0,              & ! earth radius in meters
  Rd      = 287.0,                  & ! dry gas const
  cp      = 1004.5,                 & ! heat capacity const pressure
  kappa   = Rd/cp,                  &
  pi      = 3.141592654,            &
  p0      = 100000.0                  ! reference pressure

real(rl), dimension(:,:,:,:), allocatable :: u0, v0                     ! storage for dcmip2-x sponge layer
real(rl):: zi(nlevp), zm(nlev)                                          ! z coordinates

contains

!_____________________________________________________________________
subroutine dcmip2012_test1_1(elem,hybrid,hvcoord,nets,nete,time,n0,n1)

  ! 3d deformational flow

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  real(rl),           intent(in)            :: time                     ! current time
  integer,            intent(in)            :: n0,n1                    ! time level indices

  logical ::  initialized = .false.

  integer,  parameter :: zcoords = 0                                    ! we are not using z coords
  logical,  parameter :: use_eta = .true.                               ! we are using hybrid eta coords
  real(rl), parameter ::      &
      T0      = 300.d0,       &                                         ! temperature (K)
      ztop    = 12000.d0,     &                                         ! model top (m)
      H       = Rd * T0 / g                                             ! scale height

  integer :: i,j,k,ie                                                   ! loop indices
  real(rl):: lon,lat                                                    ! pointwise coordiantes
  real(rl):: p,z,phis,u,v,w,T,phis_ps,ps,rho,q(4),dp,eta_dot,dp_dn      ! pointwise field values
  real(rl):: etai(nlevp),etam(nlev)                                           ! vertical interpolation levels

  ! set vertical interpolation coordinates
  if(.not. initialized) then
    initialized = .true.
    if (hybrid%masterthread) write(iulog,*) 'initializing dcmip2012 test 1-1: 3d deformational flow'
    call get_evenly_spaced_z(zi,zm, 0.0_rl,ztop)                        ! get evenly spaced z levels
    etai  = exp(-zi/H)                                                  ! set interpolated eta levels from z
    etam  = 0.5_rl*(etai(1:nlev)+etai(2:nlevp))
    call init_vertical_interp(etam,nlev)
  endif

  ! set prescribed state at level midpoints
  do ie = nets,nete; do k=1,nlev; do j=1,np; do i=1,np
      lon  = elem(ie)%spherep(i,j)%lon; lat  = elem(ie)%spherep(i,j)%lat
      z = H * log(1.0d0/hvcoord%etam(k))
      p = p0 * hvcoord%etam(k)
      call test1_advection_deformation(time,lon,lat,p,z,zcoords,u,v,w,T,phis,ps,rho,q(1),q(2),q(3),q(4))

      dp = pressure_thickness(ps,k,hvcoord)
      call set_state(u,v,T,ps,phis,dp,zm(k), i,j,k,elem(ie),n0,n1)
      if(time==0) then
        call set_tracers(q,qsize,dp,i,j,k,lat,lon,elem(ie))
      endif

  enddo; enddo; enddo; enddo

  ! set prescribed state at level interfaces
  do ie = nets,nete; do k=1,nlevp; do j=1,np; do i=1,np
      lon  = elem(ie)%spherep(i,j)%lon; lat  = elem(ie)%spherep(i,j)%lat
      z = H  * log(1.0d0/hvcoord%etai(k))
      p = p0 * hvcoord%etai(k)
      call test1_advection_deformation(time,lon,lat,p,z,zcoords,u,v,w,T,phis,ps,rho,q(1),q(2),q(3),q(4))

      ! get vertical derivative of p at point i,j,k
      dp_dn = ddn_hyai(k)*p0 + ddn_hybi(k)*ps

      ! get vertical eta velocity at point i,j,k
      eta_dot = -g*rho*w/p0

      ! store vertical mass flux
      elem(ie)%derived%eta_dot_dpdn(i,j,k) = eta_dot * dp_dn

  enddo; enddo; enddo; enddo

end subroutine

!_____________________________________________________________________
subroutine dcmip2012_test1_2(elem,hybrid,hvcoord,nets,nete,time,n0,n1)

  !  Hadley-like Meridional Circulation

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  real(rl),           intent(in)            :: time                     ! current time
  integer,            intent(in)            :: n0,n1                    ! time level indices

  logical ::  initialized = .false.

  integer,  parameter :: zcoords = 0                                    ! we are not using z coords
  logical,  parameter :: use_eta = .true.                               ! we are using hybrid eta coords
  real(rl), parameter ::      &
      T0      = 300.d0,       &                                         ! temperature (K)
      ztop    = 12000.d0,     &                                         ! model top (m)
      H       = Rd * T0 / g                                             ! scale height

  integer :: i,j,k,ie                                                   ! loop indices
  real(rl):: lon,lat                                                    ! pointwise coordiantes
  real(rl):: p,z,phis,u,v,w,T,phis_ps,ps,rho,q(2),dp,eta_dot,dp_dn      ! pointwise field values
  real(rl):: etai(nlevp),etam(nlev)                                     ! vertical interpolation levels

  ! set vertical interpolation coordinates
  if(.not. initialized) then
    initialized = .true.
    if (hybrid%masterthread) write(iulog,*) 'initializing dcmip2012 test 1-2: Hadley-like Meridional Circulation'
    call get_evenly_spaced_z(zi,zm, 0.0_rl,ztop)                        ! get evenly spaced z levels
    etai  = exp(-zi/H)                                                  ! set interpolated eta levels from z
    etam  = 0.5_rl*(etai(1:nlev)+etai(2:nlevp))
    call init_vertical_interp(etam,nlev)
  endif

  ! set prescribed state at level midpoints
  do ie = nets,nete; do k=1,nlev; do j=1,np; do i=1,np
      lon  = elem(ie)%spherep(i,j)%lon; lat  = elem(ie)%spherep(i,j)%lat
      z = H * log(1.0d0/hvcoord%etam(k))
      p = p0 * hvcoord%etam(k)
      call test1_advection_hadley(time,lon,lat,p,z,zcoords,u,v,w,t,phis,ps,rho,q(1),q(2))
      dp = pressure_thickness(ps,k,hvcoord)
      call set_state(u,v,T,ps,phis,dp,zm(k), i,j,k,elem(ie),n0,n1)
      if(time==0) call set_tracers(q,qsize,dp,i,j,k,lat,lon,elem(ie))

  enddo; enddo; enddo; enddo

  ! set prescribed state at level interfaces
  do ie = nets,nete; do k=1,nlevp; do j=1,np; do i=1,np
      lon  = elem(ie)%spherep(i,j)%lon; lat  = elem(ie)%spherep(i,j)%lat
      z = H  * log(1.0d0/hvcoord%etai(k))
      p = p0 * hvcoord%etai(k)
      call test1_advection_hadley(time,lon,lat,p,z,zcoords,u,v,w,T,phis,ps,rho,q(1),q(2))

      ! get vertical derivative of p at point i,j,k
      dp_dn = ddn_hyai(k)*p0 + ddn_hybi(k)*ps

      ! get vertical eta velocity at point i,j,k
      eta_dot = -g*rho*w/p0

      ! store vertical mass flux
      elem(ie)%derived%eta_dot_dpdn(i,j,k) = eta_dot * dp_dn

  enddo; enddo; enddo; enddo

end subroutine

!_____________________________________________________________________
subroutine dcmip2012_test1_3(elem,hybrid,hvcoord,nets,nete,time,n0,n1,deriv)

  !  Horizontal advection of thin cloud-like tracers over orography

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type (derivative_t),intent(in)            :: deriv
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  real(rl),           intent(in)            :: time                     ! current time
  integer,            intent(in)            :: n0,n1                    ! time level indices

  logical ::  initialized = .false.

  integer,  parameter :: cfv     = 0                                    ! h-vel is not coordinate following
  integer,  parameter :: zcoords = 0                                    ! we are not using z coords
  logical,  parameter :: use_eta = .true.                               ! we are using hybrid eta coords
  real(rl), parameter ::      &
      T0      = 300.d0,       &                                         ! temperature (K)
      ztop    = 12000.d0,     &                                         ! model top (m)
      H       = Rd * T0 / g                                             ! scale height

  integer :: i,j,k,ie                                                   ! loop indices
  real(rl):: lon,lat,hyam,hybm,hyai,hybi                                ! pointwise coordiantes
  real(rl):: p,z,phis,u,v,w,T,phis_ps,ps,rho,q(4),dp,gc                 ! pointwise field values
  real(rl):: grad_p(np,np,2),p_i(np,np),u_i(np,np),v_i(np,np)
  real(rl):: etai(nlevp),etam(nlev)                                     ! vertical interpolation levels

  ! set vertical interpolation coordinates
  if(.not. initialized) then
    initialized = .true.
    if (hybrid%masterthread) write(iulog,*) 'initializing dcmip2012 test 1-1: 3d deformational flow'
    call get_evenly_spaced_z(zi,zm, 0.0_rl,ztop)                        ! get evenly spaced z levels
    etai  = exp(-zi/H)                                                  ! set interpolated eta levels from z
    etam  = 0.5_rl*(etai(1:nlev)+etai(2:nlevp))
    call init_vertical_interp(etam,nlev)
  endif

  ! set prescribed state at level midpoints
  do ie = nets,nete; do k=1,nlev; do j=1,np; do i=1,np
      hyam=hvcoord%hyam(k); hybm=hvcoord%hybm(k)
      lon  = elem(ie)%spherep(i,j)%lon; lat  = elem(ie)%spherep(i,j)%lat
      call test1_advection_orography(lon,lat,p,z,zcoords,cfv,use_eta,hyam,hybm,gc,u,v,w,t,phis,ps,rho,q(1),q(2),q(3),q(4))
      dp = pressure_thickness(ps,k,hvcoord)
      call set_state(u,v,T,ps,phis,dp,zm(k), i,j,k,elem(ie),n0,n1)
      if(time==0) call set_tracers(q,qsize,dp,i,j,k,lat,lon,elem(ie))

  enddo; enddo; enddo; enddo

  ! set prescribed state at level interfaces
  do ie = nets,nete;
    do k=1,nlevp;
      do j=1,np; do i=1,np
        hyai=hvcoord%hyai(k); hybi=hvcoord%hybi(k)
        lon  = elem(ie)%spherep(i,j)%lon; lat  = elem(ie)%spherep(i,j)%lat
        call test1_advection_orography (lon,lat,p,z,zcoords,cfv,use_eta,hyai,hybi,gc,u,v,w,t,phis,ps,rho,q(1),q(2),q(3),q(4))
        p_i(i,j) = p
        u_i(i,j) = u
        v_i(i,j) = v
      enddo; enddo

      ! get vertical mass flux
      grad_p = gradient_sphere(p_i,deriv,elem(ie)%Dinv)
      elem(ie)%derived%eta_dot_dpdn(:,:,k) = -u_i*grad_p(:,:,1) - v_i*grad_p(:,:,2)
    enddo;
    elem(ie)%derived%eta_dot_dpdn(:,:,1)     = 0
    elem(ie)%derived%eta_dot_dpdn(:,:,nlevp) = 0
  enddo;

end subroutine

!_____________________________________________________________________
subroutine dcmip2012_test2_0(elem,hybrid,hvcoord,nets,nete)

  ! steady state atmosphere with orography

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index

  integer,  parameter :: zcoords = 0                                    ! we are not using z coords
  logical,  parameter :: use_eta = .true.                               ! we are using hybrid eta coords
  real(rl), parameter ::      &
      T0      = 300.d0,       &                                         ! temperature (K)
      gamma   = 0.0065d0,     &                                         ! temperature lapse rate (K/m)
      ztop    = 12000.d0,     &                                         ! model top (m)
      exponent= g/(Rd*gamma)

  integer :: i,j,k,ie                                                   ! loop indices
  real(rl):: lon,lat,hyam,hybm                                          ! pointwise coordiantes
  real(rl):: p,z,phis,u,v,w,T,phis_ps,ps,rho,q(1),dp    ! pointwise field values
  real(rl):: etai(nlevp),etam(nlev)                                     ! vertical interpolation levels

  if (hybrid%masterthread) write(iulog,*) 'initializing dcmip2012 test 2-0: steady state atmosphere with orography'

  ! set vertical interpolation coordinates
  call get_evenly_spaced_z(zi,zm, 0.0_rl,ztop)                          ! get evenly spaced z levels
  etai  = (1.d0 - gamma/T0*zi)**exponent                                ! set eta levels from z in orography-free region
  etam  = 0.5_rl*(etai(1:nlev)+etai(2:nlevp))
  call init_vertical_interp(etam,nlev)

  ! set initial conditions
  do ie = nets,nete; do k=1,nlev; do j=1,np; do i=1,np
    call get_coordinates(lat,lon,hyam,hybm, i,j,k,elem(ie),hvcoord)
    call test2_steady_state_mountain(lon,lat,p,z,zcoords,use_eta,hyam,hybm,u,v,w,T,phis,ps,rho,q(1))
    dp = pressure_thickness(ps,k,hvcoord)
    call set_state(u,v,T,ps,phis,dp,zm(k), i,j,k,elem(ie),1,nt)
    call set_tracers(q,1,dp,i,j,k,lat,lon,elem(ie))
  enddo; enddo; enddo; enddo

end subroutine

!_____________________________________________________________________
subroutine dcmip2012_test2_x(elem,hybrid,hvcoord,nets,nete,shear)

  ! nonhydrostatic orographic waves (with or without shear)

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  integer,            intent(in)            :: shear                    ! flag: 1=shear 0=no shear

  integer,  parameter :: zcoords = 0                                    ! we are not using z coords
  logical,  parameter :: use_eta = .true.                               ! we are using hybrid eta coords
  real(rl), parameter ::   &
      Teq     = 300.d0,    &                                            ! temperature at equator
      ztop    = 30000.d0,	 &                                            ! model top (m)
      H       = Rd*Teq/g                                                ! characteristic height scale

  integer :: i,j,k,ie                                                   ! loop indices
  real(rl):: lon,lat,hyam,hybm                                          ! pointwise coordiantes
  real(rl):: p,z,phis,u,v,w,T,phis_ps,ps,rho,q(1),dp    ! pointwise field values
  real(rl):: etai(nlevp),etam(nlev)                                     ! vertical interpolation levels

  if (hybrid%masterthread) write(iulog,*) 'initializing dcmip2012 test 2-x: steady state atmosphere with orography'

  ! set vertical interpolation coordinates
  call get_evenly_spaced_z(zi,zm, 0.0_rl,ztop)                        ! get evenly spaced z levels
  etai  = exp(-zi/H)                                                  ! set interpolated eta levels from z
  etam  = 0.5_rl*(etai(1:nlev)+etai(2:nlevp))
  call init_vertical_interp(etam,nlev)

  ! set initial conditions
  do ie = nets,nete; do k=1,nlev; do j=1,np; do i=1,np
    call get_coordinates(lat,lon,hyam,hybm, i,j,k,elem(ie),hvcoord)
    call test2_schaer_mountain(lon,lat,p,z,zcoords,use_eta,hyam,hybm,shear,u,v,w,T,phis,ps,rho,q(1))
    dp = pressure_thickness(ps,k,hvcoord)
    call set_state(u,v,T,ps,phis,dp,zm(k), i,j,k,elem(ie),1,nt)
    call set_tracers(q,1,dp,i,j,k,lat,lon,elem(ie))
  enddo; enddo; enddo; enddo

  ! store initial velocity fields for use in sponge layer
  allocate( u0(np,np,nlev,nelemd) )
  allocate( v0(np,np,nlev,nelemd) )

  do ie = nets,nete
    u0(:,:,:,ie) = elem(ie)%state%v(:,:,1,:,1)
    v0(:,:,:,ie) = elem(ie)%state%v(:,:,2,:,1)
  enddo

end subroutine

!_____________________________________________________________________
subroutine dcmip2012_test2_x_forcing(elem,hybrid,hvcoord,nets,nete,n,dt)

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(in)            :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  integer,            intent(in)            :: n                        ! time level index
  real(rl),           intent(in)            :: dt                       ! time-step size

  integer :: ie, k

  real(rl), parameter ::  &
    tau     = 25.0,       &   ! rayleigh function relaxation time
    ztop    = 30000.d0,		&   ! model top
    zh      = 20000.d0        ! sponge-layer cutoff height

  real(rl) :: f_d(nlev)                                                 ! damping function
  real(rl) :: z(np,np,nlev)

  ! Compute damping as a function of layer-midpoint height
  where(zm .ge. zh)
    f_d = sin(pi/2 *(zm - zh)/(ztop - zh))**2
  elsewhere
    f_d = 0.0d0
  end where

  ! apply sponge layer forcing to momentum terms
  do ie=nets,nete
    do k=1,nlev
      elem(ie)%derived%FM(:,:,1,k,n) = -f_d(k)/tau * ( elem(ie)%state%v(:,:,1,k,n) - u0(:,:,k,ie) )
      elem(ie)%derived%FM(:,:,2,k,n) = -f_d(k)/tau * ( elem(ie)%state%v(:,:,2,k,n) - v0(:,:,k,ie) )
    enddo
  enddo

end subroutine

!_____________________________________________________________________
subroutine dcmip2012_test3(elem,hybrid,hvcoord,nets,nete)

  ! nonhydrostatic gravity waves

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index

  integer,  parameter :: zcoords = 0                                    ! we are not using z coords
  logical,  parameter :: use_eta = .true.                               ! we are using hybrid eta coords

  real(rl), parameter ::    &                                           ! parameters needed to get eta from z
    T0      = 300.d0,       &	! temperature (k)
    ztop    = 10000.d0,     & ! model top (m)
    N       = 0.01d0,       & ! Brunt-Vaisala frequency
    bigG    = (g*g)/(N*N*Cp)  ! temperature, isothermal

  integer :: i,j,k,ie                                                   ! loop indices
  real(rl):: lon,lat,hyam,hybm                                          ! pointwise coordiantes
  real(rl):: p,z,phis,u,v,w,T,T_mean,phis_ps,ps,rho,rho_mean,q(1),dp    ! pointwise field values
  real(rl):: etai(nlevp),etam(nlev)                                     ! vertical interpolation levels

  if (hybrid%masterthread) write(iulog,*) 'initializing dcmip2012 test 3-0: nonhydrostatic gravity waves'

  ! set  vertical interpolation coordinates
  call get_evenly_spaced_z(zi,zm, 0.0_rl,ztop)                                   ! get evenly spaced z levels
  etai  = ( (bigG/T0)*(exp(-zi*N*N/g) -1 )+1 ) **(1.0/kappa)    ! set eta levels from z at equator
  etam  = 0.5_rl*(etai(1:nlev)+etai(2:nlevp))
  call init_vertical_interp(etam,nlev)

  ! set initial conditions
  do ie = nets,nete; do k=1,nlev; do j=1,np; do i=1,np
    call get_coordinates(lat,lon,hyam,hybm, i,j,k,elem(ie),hvcoord)
    call test3_gravity_wave(lon,lat,p,z,zcoords,use_eta,hyam,hybm,u,v,w,T,T_mean,phis,ps,rho,rho_mean,q(1))
    dp = pressure_thickness(ps,k,hvcoord)
    call set_state(u,v,T,ps,phis,dp,zm(k), i,j,k,elem(ie),1,nt)
    call set_tracers(q,1, dp,i,j,k,lat,lon,elem(ie))
  enddo; enddo; enddo; enddo

end subroutine

#if 0

!_____________________________________________________________________
subroutine dcmip2016_test1(elem,hybrid,hvcoord,nets,nete)

  ! moist baroclininc wave

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index

  integer,  parameter :: zcoords = 0                                    ! we are not using z coords
  logical,  parameter :: use_eta = .true.                               ! we are using hybrid eta coords

  real(rl), parameter ::    &                                           ! parameters needed to get eta from z
    T0      = 300.d0,       &	! temperature (k)
    ztop    = 10000.d0,     & ! model top (m)
    N       = 0.01d0,       & ! Brunt-Vaisala frequency
    bigG    = (g*g)/(N*N*Cp)  ! temperature, isothermal

  integer :: i,j,k,ie                                                   ! loop indices
  real(rl):: lon,lat,hyam,hybm                                          ! pointwise coordiantes
  real(rl):: p,z,phis,u,v,w,T,T_mean,phis_ps,ps,rho,rho_mean,q(1),dp    ! pointwise field values

  if (hybrid%masterthread) write(iulog,*) 'initializing dcmip2012 test 3-0: nonhydrostatic gravity waves'

  ! set analytic vertical coordinates
  call get_evenly_spaced_z(zi,zm, 0.0_rl,ztop)                                   ! get evenly spaced z levels
  hvcoord%etai  = ( (bigG/T0)*(exp(-zi*N*N/g) -1 )+1 ) **(1.0/kappa)    ! set eta levels from z at equator
  call set_hybrid_coefficients(hvcoord,hybrid,  hvcoord%etai(1), 1.0_rl)! set hybrid A and B from eta levels
  call set_layer_locations(hvcoord, .true., hybrid%masterthread)

  ! set initial conditions
  do ie = nets,nete; do k=1,nlev; do j=1,np; do i=1,np
    call get_coordinates(lat,lon,hyam,hybm, i,j,k,elem(ie),hvcoord)
    !call test3_gravity_wave(lon,lat,p,z,zcoords,use_eta,hyam,hybm,u,v,w,T,T_mean,phis,ps,rho,rho_mean,q(1))
    dp = pressure_thickness(ps,k,hvcoord)
    call set_state(u,v,T,ps,phis,dp,zm(k), i,j,k,elem(ie),1,nt)
    call set_tracers(q,1, dp,i,j,k,lat,lon,elem(ie))
  enddo; enddo; enddo; enddo

end subroutine

!_____________________________________________________________________
subroutine dcmip2016_test2(elem,hybrid,hvcoord,nets,nete)

  ! tropical cyclone

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index

  logical,  parameter :: use_eta = .true.                               ! we are using hybrid eta coords

  real(rl), parameter :: ztop = 15000.d0                                ! model top (m)

  integer :: i,j,k,ie,zcoords                                           ! loop indices
  real(rl):: lon,lat                                                    ! pointwise coordiantes
  real(rl):: p,z,phis,u,v,w,T,thetav,ps,rho,q(3),dp,eta                 ! pointwise field values
  real(rl):: A,B,etai(nlevp),etam(nlev)

  if (hybrid%masterthread) write(iulog,*) 'initializing dcmip2016 test 2: tropical cyclone'

  ! compute interpolation points at evenly spaced z
  call get_evenly_spaced_z(zi,zm, 0.0_rl,ztop)                          ! get evenly spaced z levels

  zcoords=1; lon =0; lat = 0
  do k=1,nlevp
    B = 1.0 - zi(k)/zi(1)
    call tropical_cyclone_test(lon,lat,p,zi(k),1,u,v,t,thetav,phis,ps,rho,q(1))
    A = (p - B*ps)/p0
    etai(k) = A+B
  enddo

  etam  = 0.5_rl*(etai(1:nlev)+etai(2:nlevp))
  call init_vertical_interp(etam,nlev)

  ! set initial conditions
  do ie = nets,nete; do k=1,nlev; do j=1,np; do i=1,np
    lon  = elem(ie)%spherep(i,j)%lon
    lat  = elem(ie)%spherep(i,j)%lat

    ! call routine just to get ps value
    z=0; call tropical_cyclone_test(lon,lat,p,z,1,u,v,t,thetav,phis,ps,rho,q(1))

    ! get pressure from ps and level
    p = hvcoord%hyam(k)*p0 + hvcoord%hybm(k)*ps

    ! call routine again using vertical pressure on eta level
    call tropical_cyclone_test(lon,lat,p,z,0,u,v,t,thetav,phis,ps,rho,q(1))
    q(2)=0
    q(3)=0
    dp = pressure_thickness(ps,k,hvcoord)
    call set_state(u,v,T,ps,phis,dp,zm(k), i,j,k,elem(ie),1,nt)
    call set_tracers(q,3, dp,i,j,k,lat,lon,elem(ie))
  enddo; enddo; enddo; enddo

end subroutine

!_____________________________________________________________________
subroutine dcmip2016_test2_forcing(elem,hybrid,hvcoord,nets,nete,n,n_q,dt)

  use dcmip2016_physics_z, only: dcmip2016_physics
  use vertical_se, only: compress, decompress, n_unique

  ! get forcing terms for tropical cyclone tests

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(in)            :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  integer,            intent(in)            :: n, n_q                   ! dynamics and tracer time indices
  real(rl),           intent(in)            :: dt                       ! time-step size

  real(rl), dimension(n_unique) :: u0,v0,T0,qv0,qc0,qr0,Am,Bm
  real(rl), dimension(n_unique) :: u,v,p,qv,qc,qr,rho,T,dp,z_m
  real(rl) :: lat, precl, z_i(n_unique+1), ps

  integer :: ie, i,j,k
  integer, parameter :: test     = 2
  integer, parameter :: pbl_type = 0 ! 0=Reed_jablonowski, 1=Bryan boundary layer
  integer, parameter :: prec_type= 0 ! 0=Kessler physics,  1=Reed-Jabolonowski microphysics

  Am = compress(hvcoord%hyam)
  Bm = compress(hvcoord%hybm)

  ! call physics routine one column at a time
  do ie = nets,nete; do j=1,np; do i=1,np

    ps = elem(ie)%state%ps_v(i,j,n)
    forall(k=1:n_unique) p(k)= Am(k)*p0 + Bm(k)*ps

    T   = compress( elem(ie)%state%T    (i,j,:,n) )
    z_m = compress( elem(ie)%derived%phi(i,j,:)/g )
    z_i(1:n_unique)=z_m; z_i(n_unique+1)=0        ! (TODO: get accurate levels for precl from Reed-Jablonowski)
    z_i(2:n_unique) = 0.5*( z_m(1:n_unique-1) + z_m(2:n_unique))
    z_i(1) = z_m(1)
    z_i(n_unique+1) = z_m(n_unique)

    ! get dry air density from ideal gas law
    rho = p/(Rd*T)

    ! get each state vector at column i,j
    lat = elem(ie)%spherep(i,j)%lat
    u   = compress( elem(ie)%state%v(i,j,1,:,n) )
    v   = compress( elem(ie)%state%v(i,j,2,:,n) )
    qv  = compress( elem(ie)%state%Q(i,j,:,1)   ) ! vapor mixing ratio        (gm/gm)
    qc  = compress( elem(ie)%state%Q(i,j,:,2)   ) ! cloud water mixing ratio  (gm/gm)
    qr  = compress( elem(ie)%state%Q(i,j,:,3)   ) ! rain water mixing ratio   (gm/gm)

    precl = 0 ! precipitation rate (output from physics routine)

    ! store current column values to compute differences after update
    u0=u;   v0=v;   T0=T;   qv0=qv;   qc0=qc;   qr0=qr

    ! apply KESSLER/Reed-Jabolonowski microphysics and Reed-Jablonoswki/Bryan boundary layer
    call dcmip2016_physics(test,u,v,p,qv,qc,qr,rho,dt,z_m,z_i,lat,n_unique,precl,pbl_type,prec_type)

    ! compute new temperatures from ideal gas law applied to new pressures
    T = p/(Rd*rho)

    ! get forcing terms from difference between old and new column values
    elem(ie)%derived%FM(i,j,1,:,n  )  = decompress( u  - u0 )
    elem(ie)%derived%FM(i,j,2,:,n  )  = decompress( v  - v0 )
    elem(ie)%derived%FT(i,j,  :,n  )  = decompress( T  - T0 )
    elem(ie)%derived%FQ(i,j,:,1,n_q)  = decompress( qv - qv0)
    elem(ie)%derived%FQ(i,j,:,2,n_q)  = decompress( qc - qc0)
    elem(ie)%derived%FQ(i,j,:,3,n_q)  = decompress( qr - qr0)

if(any(isnan(u  - u0)))   call abortmp("nan in u")
if(any(isnan(v  - v0)))   call abortmp("nan in v")
if(any(isnan(T  - T0)))   call abortmp("nan in T")
if(any(isnan(qv  - qv0))) call abortmp("nan in qv")
if(any(isnan(qc  - qc0))) call abortmp("nan in qc")
if(any(isnan(qr  - qr0))) call abortmp("nan in qr")

    enddo; enddo;

  enddo

end subroutine

!_____________________________________________________________________
subroutine dcmip2016_test3(elem,hybrid,hvcoord,nets,nete)

  ! supercell storm

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index

  integer,  parameter :: zcoords = 0                                    ! we are not using z coords
  logical,  parameter :: use_eta = .true.                               ! we are using hybrid eta coords

  real(rl), parameter ::    &                                           ! parameters needed to get eta from z
    T0      = 300.d0,       &	! temperature (k)
    ztop    = 10000.d0,     & ! model top (m)
    N       = 0.01d0,       & ! Brunt-Vaisala frequency
    bigG    = (g*g)/(N*N*Cp)  ! temperature, isothermal

  integer :: i,j,k,ie                                                   ! loop indices
  real(rl):: lon,lat,hyam,hybm                                          ! pointwise coordiantes
  real(rl):: p,z,phis,u,v,w,T,T_mean,phis_ps,ps,rho,rho_mean,q(1),dp    ! pointwise field values

  if (hybrid%masterthread) write(iulog,*) 'initializing dcmip2012 test 3-0: nonhydrostatic gravity waves'

  ! set analytic vertical coordinates
  call get_evenly_spaced_z(zi,zm, 0.0_rl,ztop)                                   ! get evenly spaced z levels
  hvcoord%etai  = ( (bigG/T0)*(exp(-zi*N*N/g) -1 )+1 ) **(1.0/kappa)    ! set eta levels from z at equator
  call set_hybrid_coefficients(hvcoord,hybrid,  hvcoord%etai(1), 1.0_rl)! set hybrid A and B from eta levels
  call set_layer_locations(hvcoord, .true., hybrid%masterthread)

  ! set initial conditions
  do ie = nets,nete; do k=1,nlev; do j=1,np; do i=1,np
    call get_coordinates(lat,lon,hyam,hybm, i,j,k,elem(ie),hvcoord)
    !call test3_gravity_wave(lon,lat,p,z,zcoords,use_eta,hyam,hybm,u,v,w,T,T_mean,phis,ps,rho,rho_mean,q(1))
    dp = pressure_thickness(ps,k,hvcoord)
    call set_state(u,v,T,ps,phis,dp,zm(k), i,j,k,elem(ie),1,nt)
    call set_tracers(q,1, dp,i,j,k,lat,lon,elem(ie))
  enddo; enddo; enddo; enddo

end subroutine

#endif

!_____________________________________________________________________
subroutine get_evenly_spaced_z(zi,zm, zb,zt)

  real(rl), intent(in)    :: zb,zt      ! top and bottom coordinates
  real(rl), intent(inout) :: zi(nlevp)  ! z at interfaces
  real(rl), intent(inout) :: zm(nlev)   ! z at midpoints
  integer :: k

  forall(k=1:nlevp) zi(k) = zt-(k-1)*(zt-zb)/(nlevp-1)
  zm = 0.5_rl*( zi(2:nlevp) + zi(1:nlev) )

end subroutine

!_____________________________________________________________________
subroutine set_hybrid_coefficients(hv, hybrid, eta_t, c)

  ! create an analytical set of A,B coefficients, given known eta levels

  type(hvcoord_t),    intent(inout) :: hv       ! hybrid vertical coordinate stucture
  type(hybrid_t),     intent(in)    :: hybrid   ! hybrid parallal structure
  real(rl),           intent(in)    :: eta_t    ! top eta level
  real(rl),           intent(in)    :: c        ! exponent

  real(rl)  :: eta_c, tmp
  integer   :: k

  ! place cutoff halfway between bottom and top eta coordiantes
  eta_c = hv%etai(nlev/2)

  do k=1,nlevp
    ! get values of hybrid coefficients
    tmp        = max( (hv%etai(k)-eta_c)/(1.0-eta_c), 0.0_rl)
    hv%hybi(k) = tmp**c
    hv%hyai(k) = hv%etai(k) - hv%hybi(k)
    if(hybrid%par%masterproc) print *,k,': etai = ',hv%etai(k),' Ai = ',hv%hyai(k),' Bi = ',hv%hybi(k);

    ! get derivatives of hybrid coefficients
    ddn_hybi(k) = c*tmp**(c-1)
    if(hv%etai(k)>eta_c) ddn_hybi(k)=0.0d0
    ddn_hyai(k) = 1.0d0 - ddn_hybi(k)
  enddo

  hv%hyam = 0.5_rl *(hv%hyai(2:nlev+1) + hv%hyai(1:nlev))
  hv%hybm = 0.5_rl *(hv%hybi(2:nlev+1) + hv%hybi(1:nlev))
  hv%etam = hv%hyam + hv%hybm

end subroutine

!_____________________________________________________________________
subroutine get_coordinates(lat,lon,hyam,hybm, i,j,k,elem,hvcoord)

  ! get lat,lon, vertical coords at node(i,j,k)

  real(rl),         intent(out):: lon,lat,hyam,hybm
  integer,          intent(in) :: i,j,k
  type(element_t),  intent(in) :: elem
  type(hvcoord_t),  intent(in) :: hvcoord

  ! get horizontal coordinates at column i,j
  lon  = elem%spherep(i,j)%lon
  lat  = elem%spherep(i,j)%lat

  ! get hybrid coeffiecients at midpoint of vertical level k
  hyam = hvcoord%hyam(k)
  hybm = hvcoord%hybm(k)

end subroutine

!_____________________________________________________________________
real(rl) function pressure_thickness(ps,k,hv)

  real(rl),         intent(in) :: ps
  integer,          intent(in) :: k
  type(hvcoord_t),  intent(in) :: hv
  pressure_thickness = (hv%hyai(k+1)-hv%hyai(k))*p0 + (hv%hybi(k+1)-hv%hybi(k))*ps

end function

!_____________________________________________________________________
subroutine set_state(u,v,T,ps,phis,dp,zm, i,j,k,elem,n0,n1)

  ! set state variables at node(i,j,k) at layer midpoints

  real(rl),         intent(in)    :: u,v,T,ps,phis,dp,zm
  integer,          intent(in)    :: i,j,k,n0,n1
  type(element_t),  intent(inout) :: elem

  ! set prognostic state variables at level midpoints
  elem%state%v   (i,j,1,k,n0:n1) = u
  elem%state%v   (i,j,2,k,n0:n1) = v
  elem%state%T   (i,j,k,n0:n1)   = T
  elem%state%dp3d(i,j,k,n0:n1)   = dp
  elem%state%ps_v(i,j,n0:n1)     = ps
  elem%state%phis(i,j)           = phis

  ! set some diagnostic variables
  elem%derived%dp(i,j,k)         = dp
  elem%derived%phi(i,j,k)        = g*zm

end subroutine

!_____________________________________________________________________
subroutine set_tracers(q,nq, dp,i,j,k,lat,lon,elem)

  ! set tracer values at node(i,j,k)

  real(rl),         intent(in)    :: q(nq), dp, lat, lon
  integer,          intent(in)    :: i,j,k,nq
  type(element_t),  intent(inout) :: elem

  real(rl), parameter :: wl = 1.0 ! checkerboard wavelength in dg
  integer :: qi

  ! set known tracers to q and the rest to a checkerboard pattern
  elem%state%Q(i,j,k,:)    = 0.5d0*(1.0d0+sign(1.0d0,sin(lat*360.0/wl)*sin(lon*360.0/wl)))
  elem%state%Q(i,j,k,1:nq) = q

  ! compute tracer mass qdp from mixing ratio q
  do qi = 1,nq
    elem%state%Qdp (i,j,k,qi,:) = q(qi)*dp
  enddo

end subroutine
end module dcmip_tests

#endif
