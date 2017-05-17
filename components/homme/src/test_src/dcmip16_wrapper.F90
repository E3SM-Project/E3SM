



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
use element_ops,          only: set_state, set_elem_state, get_state, tests_finalize, set_forcing_rayleigh_friction, set_thermostate
use physical_constants,   only: p0, g, Rgas, kappa, Cp, Rwater_vapor
use reduction_mod,        only: parallelmax
use time_mod,             only: time_at, TimeLevel_t

implicit none

real(rl), dimension(:,:,:,:), allocatable :: u0, v0                     ! storage for dcmip2-x sponge layer
real(rl):: zi(nlevp), zm(nlev)                                          ! z coordinates
real(rl):: ddn_hyai(nlevp), ddn_hybi(nlevp)                             ! vertical derivativess of hybrid coefficients
real(rl):: tau
real(rl), parameter :: ztop3  = 20000_rl !20000_rl                            ! top of model at 20km

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

  real(rl), parameter :: ztop2    = 30000_rl                             ! top of model at 30km

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

  integer :: i,j,k,ie                                                   ! loop indices
  real(rl):: lon,lat                                                    ! pointwise coordiantes

  real(rl), dimension(np,np,nlev,qsize_d) :: q
  real(rl), dimension(np,np,nlev) :: p,dp,z,u,v,w,T,thetav,rho,rhom
  real(rl), dimension(np,np) :: phis, ps

  real(rl) :: p1,thetav1,rho1,q1

  if (hybrid%masterthread) write(iulog,*) 'initializing dcmip2016 test 3: supercell storm'

  ! initialize hydrostatic state
  call supercell_init()

  ! get evenly spaced z levels from 0 to 20km
  call get_evenly_spaced_z(zi,zm, 0.0_rl, ztop3)                          ! get evenly spaced z levels

  ! get eta levels matching z at equator
  do k=1,nlevp
    lon =0; lat = 0;
    call supercell_z(lon, lat, zi(k), p1, thetav1, rho1, q1, pert)
    hvcoord%etai(k) = p1/p0
    if (hybrid%masterthread) print *,"k=",k," z=",zi(k)," p=",p1,"etai=",hvcoord%etai(k)
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

        ! get surface pressure ps(i,j) at lat, lon, z=0
        z(i,j,k)=0; call supercell_test(lon,lat,p(i,j,k),z(i,j,k),1,u(i,j,k),v(i,j,k),T(i,j,k),thetav(i,j,k),ps(i,j),rho(i,j,k),q(i,j,k,1),pert)

        ! get hydrostatic pressure at level k
        p(i,j,k)  = p0*hvcoord%hyam(k) + ps(i,j)*hvcoord%hybm(k)
        dp(i,j,k) = pressure_thickness(ps(i,j),k,hvcoord)

        call supercell_test(lon,lat,p(i,j,k),z(i,j,k),zcoords,u(i,j,k),v(i,j,k),T(i,j,k),thetav(i,j,k),ps(i,j),rho(i,j,k),q(i,j,k,1),pert)

        w   (i,j,k)  = 0 ! no vertical motion
        phis(i,j)    = 0 ! no topography
        q   (i,j,k,2)= 0 ! no initial clouds
        q   (i,j,k,3)= 0 ! no initial rain
        q   (i,j,k,4)= 0 ! precip

        call set_tracers(q(i,j,k,:),qsize_d,dp(i,j,k),i,j,k,lat,lon,elem(ie))

      enddo; enddo

      call set_elem_state(u,v,w,T,ps,phis,p,dp,z,g,elem(ie),1,nt,ntQ=1)

    enddo

    ! set density to ensure hydrostatic balance and save initial state
    call tests_finalize(elem(ie),hvcoord,1,nt,ie)

  enddo

end subroutine

!_______________________________________________________________________
subroutine dcmip2016_append_measurements(max_w,max_precl,tl,hybrid)
  real(rl),           intent(in) :: max_w, max_precl
  type(TimeLevel_t),  intent(in) :: tl
  type(hybrid_t),     intent(in) :: hybrid                   ! hybrid parallel structure

  real(rl),         parameter :: sample_period  = 60.0_rl
  character(len=*), parameter :: w_filename     = "measurement_wmax.txt"
  character(len=*), parameter :: precl_filename = "measurement_prect_rate.txt"
  character(len=*), parameter :: time_filename  = "measurement_time.txt"

  real(rl) :: pmax_w, pmax_precl, time
  real(rl) :: next_sample_time = 0.0

  time = time_at(tl%nstep)

  ! initialize output file
  if(next_sample_time == 0.0) then
    open(unit=10,file=w_filename,    form="formatted",status="replace")
    close(10)

    open(unit=11,file=precl_filename,form="formatted",status="replace")
    close(11)

    open(unit=12,file=time_filename, form="formatted",status="replace")
    close(12)

  endif

  ! append measurements at regular intervals
  if(time .ge. next_sample_time) then

    next_sample_time = next_sample_time + sample_period
    pmax_w     = parallelMax(max_w,    hybrid)
    pmax_precl = parallelMax(max_precl,hybrid)

    if (hybrid%masterthread) then
      print *,"time=",time_at(tl%nstep)," pmax_w=",pmax_w," pmax_precl=",pmax_precl

      open(unit=10,file=w_filename,form="formatted",position="append")
        write(10,'(99E24.15)') pmax_w
      close(10)

      open(unit=11,file=precl_filename,form="formatted",position="append")
        write(11,'(99E24.15)') pmax_precl
      close(11)

      open(unit=12,file=time_filename,form="formatted",position="append")
        write(12,'(99E24.15)') time
      close(12)
    endif
  endif

  end subroutine

!_______________________________________________________________________
subroutine dcmip2016_forcing(elem,hybrid,hvcoord,nets,nete,nt,ntQ,dt,tl, test)

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(in)            :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  integer,            intent(in)            :: nt, ntQ                  ! time level index
  real(rl),           intent(in)            :: dt                       ! time-step size
  type(TimeLevel_t),  intent(in)            :: tl                       ! time level structure
  integer,            intent(in)            :: test                     ! dcmip2016 test number

  real(rl), parameter :: zc   = 15000_rl                                ! sponge layer cutoff at 13km
  real(rl), parameter :: tau  = 120_rl                                  ! rayleigh damping time-scale (s)

  integer :: i,j,k,ie                                                     ! loop indices
  real(rl):: lat
  real(rl), dimension(np,np,nlev) :: u,v,w,T,theta,exner,exner_kess,p,dp,cp_star,rho,z,qv,qc,qr
  real(rl), dimension(np,np,nlev) :: T0,qv0,qc0,qr0
  real(rl), dimension(np,np)      :: precl
  real(rl), dimension(np,np,nlev) :: theta_inv,qv_inv,qc_inv,qr_inv,rho_inv,exner_inv,z_inv ! inverted columns
  real(rl) :: max_w, pmax_w, max_precl, pmax_precl

  max_w     = -huge(rl)
  max_precl = -huge(rl)

  do ie = nets,nete

    ! get current element state
    call get_state(u,v,w,T,theta,exner,p,dp,cp_star,rho,z,g,i,j,elem(ie),hvcoord,nt,ntQ)

    ! get mixing ratios
    qv  = elem(ie)%state%Qdp(:,:,:,1,ntQ)/dp
    qc  = elem(ie)%state%Qdp(:,:,:,2,ntQ)/dp
    qr  = elem(ie)%state%Qdp(:,:,:,3,ntQ)/dp

    ! save un-forced prognostics
    T0=T; qv0=qv; qc0=qc; qr0=qr

    ! compute form of exner pressure expected by Kessler physics
    exner_kess = (p/p0)**(Rgas/Cp)

    ! invert columns (increasing z)
    theta_inv= theta(:,:,nlev:1:-1)
    qv_inv   = qv   (:,:,nlev:1:-1)
    qc_inv   = qc   (:,:,nlev:1:-1)
    qr_inv   = qr   (:,:,nlev:1:-1)
    rho_inv  = rho  (:,:,nlev:1:-1)
    exner_inv= exner_kess(:,:,nlev:1:-1)
    z_inv    = z    (:,:,nlev:1:-1)

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
    theta = theta_inv(:,:,nlev:1:-1)
    qv    = qv_inv   (:,:,nlev:1:-1)
    qc    = qc_inv   (:,:,nlev:1:-1)
    qr    = qr_inv   (:,:,nlev:1:-1)
    T     = theta*exner

    ! add sponge layer at top of model?
    ! call set_forcing_rayleigh_friction(elem(ie),z,ztop,zc,tau,u0(:,:,:,ie),v0(:,:,:,ie),nt)

    ! set dynamics forcing
    elem(ie)%derived%FM(:,:,1,:) = 0
    elem(ie)%derived%FM(:,:,2,:) = 0
    elem(ie)%derived%FT(:,:,:)   = (T - T0)/dt

    ! set tracer-mass forcing
    elem(ie)%derived%FQ(:,:,:,1) = dp*(qv-qv0)/dt
    elem(ie)%derived%FQ(:,:,:,2) = dp*(qc-qc0)/dt
    elem(ie)%derived%FQ(:,:,:,3) = dp*(qr-qr0)/dt

    ! debug: set theta directly, instead of FT
    ! elem(ie)%state%theta_dp_cp(:,:,:,nt) = theta*(Cp_star*dp)

    !elem(ie)%state%Qdp(:,:,1,4,ntQ) = precl*dp(:,:,1) ! store precl in level 1 of tracer #4

    ! perform measurements of max w, and max prect
    max_w     = max( max_w    , maxval(w    ) )
    max_precl = max( max_precl, maxval(precl) )

  enddo

  call dcmip2016_append_measurements(max_w,max_precl,tl,hybrid)

end subroutine

end module dcmip16_wrapper
#endif
