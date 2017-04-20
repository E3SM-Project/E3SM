



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
use dimensions_mod,       only: np, nlev, nlevp, qsize, qsize_d, nelemd
use element_mod,          only: element_t
use element_state,        only: nt=>timelevels
use hybrid_mod,           only: hybrid_t
use hybvcoord_mod,        only: hvcoord_t, set_layer_locations
use kinds,                only: rl=>real_kind, iulog
use parallel_mod,         only: abortmp
use element_ops,          only: set_state, copy_state, tests_finalize, set_forcing_rayleigh_friction
use physical_constants,   only: p0, g

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
  logical,  parameter :: use_eta = .true.                               ! we are using hybrid eta coords

  integer :: i,j,k,ie                                                   ! loop indices
  real(rl):: lon,lat,hyam,hybm                                          ! pointwise coordiantes
  real(rl):: p,z,phis,u,v,w,T,T_mean,phis_ps,ps,rho,rho_mean,q(1),dp    ! pointwise field values

  if (hybrid%masterthread) write(iulog,*) 'initializing dcmip2016 test 2: tropical cyclone'
  if (hybrid%masterthread) call abortmp('dcmip2016 test 2 not yet implemented for this solver')

  ! set analytic vertical coordinates
  !call get_evenly_spaced_z(zi,zm, 0.0_rl,ztop)                                   ! get evenly spaced z levels
  !hvcoord%etai  = ( (bigG/T0)*(exp(-zi*N*N/g) -1 )+1 ) **(1.0/kappa)    ! set eta levels from z at equator
  !call set_hybrid_coefficients(hvcoord,hybrid,  hvcoord%etai(1), 1.0_rl)! set hybrid A and B from eta levels
  !call set_layer_locations(hvcoord, .true., hybrid%masterthread)

  ! set initial conditions
  !do ie = nets,nete
  !   do k=1,nlev; do j=1,np; do i=1,np
  !      call get_coordinates(lat,lon,hyam,hybm, i,j,k,elem(ie),hvcoord)
  !      call test3_gravity_wave(lon,lat,p,z,zcoords,use_eta,hyam,hybm,u,v,w,T,T_mean,phis,ps,rho,rho_mean,q(1))
  !      dp = pressure_thickness(ps,k,hvcoord)
  !      call set_state(u,v,w,T,ps,phis,p,dp,zm(k),g, i,j,k,elem(ie),1,nt)
  !      call set_tracers(q,1, dp,i,j,k,lat,lon,elem(ie))
  !   enddo; enddo; enddo;
  !   call tests_finalize(elem(ie),hvcoord,1,nt)
  !enddo

end subroutine

!_____________________________________________________________________
subroutine dcmip2016_test3(elem,hybrid,hvcoord,nets,nete)

  ! supercell storm test case

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index

  integer,  parameter :: zcoords  = 1                                   ! 0 -> use p coords
  integer,  parameter :: pert     = 1                                   ! 1 -> add thermal perturbation
  real(rl), parameter :: ztop     = 20000_rl                            ! top of model at 20km

  integer :: i,j,k,ie                                                   ! loop indices
  real(rl):: lon,lat                                                    ! pointwise coordiantes
  real(rl):: p,dp, z,phis,u,v,w,T,thetav,phis_ps,ps,rho,q(3)            ! pointwise field values
  real(rl):: p_i1,p_i2,p_save            ! pointwise field values

  if (hybrid%masterthread) write(iulog,*) 'initializing dcmip2016 test 3: supercell storm'

  ! compute hydrostatic initial state for this test case
  call supercell_init()

  ! get evenly spaced z coords from 0 to 20km
  call get_evenly_spaced_z(zi,zm, 0.0_rl, ztop)                          ! get evenly spaced z levels

  ! get eta levels matching z levels, at equator
  do k=1,nlevp
    lon =0; lat = 0;
    call supercell_z(lon, lat, zi(k), p, thetav, rho, q(1), pert)
    hvcoord%etai(k) = p/p0
    if (hybrid%masterthread) print *,"k=",k," z=",zm(k)," p=",p,"etai=",hvcoord%etai(k)
  enddo
  call set_hybrid_coefficients(hvcoord,hybrid,  hvcoord%etai(1), 1.0_rl)
  call set_layer_locations(hvcoord, .true., hybrid%masterthread)

  ! set initial conditions
  do ie = nets,nete
  if (hybrid%masterthread) print *,"ie=",ie

    do k=1,nlev;

      do j=1,np; do i=1,np
        lon = elem(ie)%spherep(i,j)%lon
        lat = elem(ie)%spherep(i,j)%lat

        ! set initial conditions on at constant eta levels
        if(zcoords==0) then
          ! get surface pressure at lat, lon, z=0
          z=0; call supercell_test(lon,lat,p,z,1,u,v,T,thetav,ps,rho,q(1),pert)

          ! get hydrostatic pressure level
          p  =  p0*hvcoord%hyam(k) + ps*hvcoord%hybm(k)
          dp = (hvcoord%hyai(k+1)-hvcoord%hyai(k))*p0 + (hvcoord%hybi(k+1)-hvcoord%hybi(k))*ps
        endif

        ! set initial conditions as if on constant z surfaces
        if(zcoords==1) then
          call supercell_z(lon, lat, zi(k)  , p_i1, thetav, rho, q(1), pert)
          call supercell_z(lon, lat, zi(k+1), p_i2, thetav, rho, q(1), pert)
          dp = p_i2-p_i1

          z = zm(k)
        endif

        call supercell_test(lon,lat,p,z,zcoords,u,v,T,thetav,ps,rho,q(1),pert)

        w    = 0 ! no vertical motion
        phis = 0 ! no topography
        q(2) = thetav !0 ! no initial clouds
        q(3) = T !0 ! no initial rain

        !print *," u=",u," T=",T," ps=",ps," z=",z," p=",p," dp=",dp, " q=",q(1)," dp=",dp

        call set_state(u,v,w,T,ps,phis,p,dp,z,g,i,j,k,elem(ie),1,nt)
        call set_tracers(q,3, dp,i,j,k,lat,lon,elem(ie))

      enddo; enddo
    enddo

  enddo

end subroutine







! get eta levels corresponding to midpoint z levels, at the equator
!hvcoord%etai(nlevp)=1
!do k=nlev,1,-1
!  lon =0; lat = 0;
!  call supercell_z(lon, lat, zm(k), p, thetav, rho, q(1), pert=0)
!  hvcoord%etam(k) = p/p0
!  hvcoord%etai(k) = 2.0*hvcoord%etam(k) - hvcoord%etai(k+1)

!  if (hybrid%masterthread) print *,"k=",k," z=",zm(k)," p=",p,"etam=",hvcoord%etam(k), "etai=",hvcoord%etai(k)
!enddo

end module dcmip16_wrapper
#endif
