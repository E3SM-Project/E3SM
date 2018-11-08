



#ifndef CAM
#include "config.h"

module dcmip16_wrapper

! Implementation of the dcmip2012 dycore tests for the preqx dynamics target

use dcmip12_wrapper,      only: pressure_thickness, set_tracers, get_evenly_spaced_z, set_hybrid_coefficients
use control_mod,          only: test_case, dcmip16_pbl_type, dcmip16_prec_type, use_moisture, theta_hydrostatic_mode
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
use parallel_mod,         only: abortmp,iam
use element_ops,          only: set_state, set_state_i, set_elem_state, get_state, tests_finalize,&
     set_forcing_rayleigh_friction, set_thermostate
use physical_constants,   only: p0, g, Rgas, kappa, Cp, Rwater_vapor, pi=>dd_pi
use reduction_mod,        only: parallelmax, parallelmin
use terminator,           only: initial_value_terminator, tendency_terminator
use time_mod,             only: time_at, TimeLevel_t

implicit none

! this cannot be made stack variable, nelemd is not compile option
real(rl),dimension(:,:,:), allocatable :: precl ! storage for column precip
real(rl):: zi(nlevp), zm(nlev)                                          ! z coordinates
real(rl):: ddn_hyai(nlevp), ddn_hybi(nlevp)                             ! vertical derivativess of hybrid coefficients
real(rl):: tau
real(rl), parameter :: rh2o    = 461.5d0,            &                  ! Gas constant for water vapor (J/kg/K)
                       Mvap    = (Rwater_vapor/Rgas) - 1.d0             ! Constant for virtual temp. calc. (~0.608)

real(rl) :: sample_period  = 60.0_rl
real(rl) :: rad2dg = 180.0_rl/pi
contains

!---------------------------------------------------------------------
!init routine to call before any dcmip16 inin routines, including restart runs
subroutine dcmip2016_init()
  implicit none
!$OMP BARRIER
!$OMP MASTER
  if (.not.allocated(precl)) then
    allocate(precl(np,np,nelemd))
    precl(:,:,:) = 0.0
  else
    call abortmp('ERROR: in dcmip2016_init() precl has already been allocated') 
  endif
!$OMP END MASTER
!$OMP BARRIER

end subroutine dcmip2016_init

!_____________________________________________________________________
subroutine dcmip2016_test1(elem,hybrid,hvcoord,nets,nete)

  ! moist baroclinic wave test

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index

  integer,  parameter :: use_zcoords  = 0                               ! use vertical pressure coordinates
  integer,  parameter :: is_deep      = 0                               ! use shallow atmosphere approximation
  integer,  parameter :: pertt        = 0                               ! use exponential perturbation type
  real(rl), parameter :: dcmip_X      = 1.0_rl                          ! full scale planet
  integer :: moist                                                      ! use moist version
  integer :: i,j,k,ie                                                   ! loop indices
  real(rl):: lon,lat                                                    ! pointwise coordiantes

  real(rl), dimension(np,np,nlev):: p,z,u,v,w,T,thetav,rho,dp           ! field values
  real(rl), dimension(np,np,nlevp):: p_i,w_i,z_i
  real(rl), dimension(np,np):: ps, phis
  real(rl), dimension(np,np,nlev,6):: q

  real(rl) :: min_thetav, max_thetav
  min_thetav = +huge(rl)
  max_thetav = -huge(rl)

  moist = 0
  if (use_moisture) moist=1

  if (hybrid%masterthread) write(iulog,*) 'initializing dcmip2016 test 1: moist baroclinic wave'

  if (qsize<5) call abortmp('ERROR: test requires qsize>=5')

  ! set initial conditions
  do ie = nets,nete
    do k=1,nlevp; do j=1,np; do i=1,np

        ! no surface topography
        p_i(i,j,k)  = p0*hvcoord%etai(k)

        lon = elem(ie)%spherep(i,j)%lon
        lat = elem(ie)%spherep(i,j)%lat

        w_i(i,j,k)   = 0.0d0
        ! call this only to compute z_i, will ignore all other output
        call baroclinic_wave_test(is_deep,moist,pertt,dcmip_X,lon,lat,p_i(i,j,k),&
            z_i(i,j,k),use_zcoords,u(i,j,1),v(i,j,1),T(i,j,1),thetav(i,j,1),phis(i,j),ps(i,j),rho(i,j,1),q(i,j,1,1))
    enddo; enddo; enddo
    do k=1,nlev; do j=1,np; do i=1,np

        ! no surface topography
        p(i,j,k)  = p0*hvcoord%etam(k)
        dp(i,j,k) = (hvcoord%etai(k+1)-hvcoord%etai(k))*p0

        lon = elem(ie)%spherep(i,j)%lon
        lat = elem(ie)%spherep(i,j)%lat

        q(i,j,k,1:5) = 0.0d0
        q(i,j,k,6) = 1
        w(i,j,k)   = 0.0d0

        call baroclinic_wave_test(is_deep,moist,pertt,dcmip_X,lon,lat,p(i,j,k),&
            z(i,j,k),use_zcoords,u(i,j,k),v(i,j,k),T(i,j,k),thetav(i,j,k),phis(i,j),ps(i,j),rho(i,j,k),q(i,j,k,1))

        ! initialize tracer chemistry
        call initial_value_terminator( lat*rad2dg, lon*rad2dg, q(i,j,k,4), q(i,j,k,5) )
        call set_tracers(q(i,j,k,1:6),6,dp(i,j,k),i,j,k,lat,lon,elem(ie))

        min_thetav =  min( min_thetav,   thetav(i,j,k) )
        max_thetav =  max( max_thetav,   thetav(i,j,k) )

    enddo; enddo; enddo

    call set_elem_state(u,v,w,w_i,T,ps,phis,p,dp,z,z_i,g,elem(ie),1,nt,ntQ=1)
    call tests_finalize(elem(ie),hvcoord,1,nt)

  enddo
  sample_period = 1800.0 ! sec
  !print *,"min thetav = ",min_thetav, "max thetav=",max_thetav


end subroutine

!_____________________________________________________________________
subroutine dcmip2016_test2(elem,hybrid,hvcoord,nets,nete)

  ! tropical cyclone test

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index

  integer :: i,j,k,ie , ierr                                                  ! loop indices
  real(rl):: lon,lat,ntop                                               ! pointwise coordiantes
  real(rl), dimension(np,np,nlev):: p,z,u,v,w,T,rho,dp                     ! pointwise field values
  real(rl), dimension(np,np):: ps,phis
  real(rl), dimension(np,np,nlevp):: w_i,z_i,p_i
  real(rl) :: thetav,q(3)

  if (hybrid%masterthread) write(iulog,*) 'initializing dcmip2016 test 2: tropical cyclone'
  !use vertical levels specificed in cam30 file

  ! set initial conditions
  do ie = nets,nete
     do k=1,nlevp; do j=1,np; do i=1,np

        ! get surface pressure ps(i,j) at lat, lon, z=0, ignore all other output
        lon = elem(ie)%spherep(i,j)%lon
        lat = elem(ie)%spherep(i,j)%lat
        z_i(i,j,k)=0; call tropical_cyclone_test(lon,lat,p(i,j,1),z_i(i,j,k),1,u(i,j,1),&
             v(i,j,1),T(i,j,1),thetav,phis(i,j),ps(i,j),rho(i,j,1),q(1))

        ! get hydrostatic pressure at level k
        p_i(i,j,k)  = p0*hvcoord%hyai(k) + ps(i,j)*hvcoord%hybi(k)

        ! call this only to compute z_i, will ignore all other output
        call tropical_cyclone_test(lon,lat,p_i(i,j,k),z_i(i,j,k),0,u(i,j,1),v(i,j,1),&
             T(i,j,1),thetav,phis(i,j),ps(i,j),rho(i,j,1),q(1))

        w_i(i,j,k)  = 0

     enddo; enddo; enddo;

     do k=1,nlev; do j=1,np; do i=1,np

        ! get surface pressure
        lon = elem(ie)%spherep(i,j)%lon
        lat = elem(ie)%spherep(i,j)%lat
        z=0; call tropical_cyclone_test(lon,lat,p(i,j,k),z(i,j,k),1,u(i,j,k),v(i,j,k),&
             T(i,j,k),thetav,phis(i,j),ps(i,j),rho(i,j,k),q(1))

        ! get pressure at level midpoints
        p(i,j,k) = p0*hvcoord%hyam(k) + ps(i,j)*hvcoord%hybm(k)

        ! get initial conditions at pressure level p
        call tropical_cyclone_test(lon,lat,p(i,j,k),z(i,j,k),0,u(i,j,k),v(i,j,k),&
             T(i,j,k),thetav,phis(i,j),ps(i,j),rho(i,j,k),q(1))

        dp(i,j,k) = pressure_thickness(ps(i,j),k,hvcoord)
        w(i,j,k)  = 0
        q(2)=0
        q(3)=0

        call set_tracers(q(:),3,dp(i,j,k),i,j,k,lat,lon,elem(ie))
     enddo; enddo; enddo;

    call set_elem_state(u,v,w,w_i,T,ps,phis,p,dp,z,z_i,g,elem(ie),1,nt,ntQ=1)
    call tests_finalize(elem(ie),hvcoord,1,nt)
  enddo

  sample_period = 1800.0 ! sec

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

  real(rl), parameter :: ztop3  = 20000_rl                              ! top of model at 20km

  real(rl), dimension(np,np,nlev,3) :: q
  real(rl), dimension(np,np,nlev) :: p,dp,z,u,v,w,T,thetav,rho,rhom
  real(rl), dimension(np,np,nlevp) :: p_i,z_i,w_i
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

    do k=1,nlevp

      do j=1,np; do i=1,np

        lon = elem(ie)%spherep(i,j)%lon
        lat = elem(ie)%spherep(i,j)%lat

        ! get surface pressure ps(i,j) at lat, lon, z=0, ignore all other output
        z_i(i,j,k)=0; call supercell_test(lon,lat,p(i,j,1),z_i(i,j,k),1,u(i,j,1),v(i,j,1),&
             T(i,j,1),thetav(i,j,1),ps(i,j),rho(i,j,1),q(i,j,1,1),pert)

        ! get hydrostatic pressure at level k
        p_i(i,j,k)  = p0*hvcoord%hyai(k) + ps(i,j)*hvcoord%hybi(k)

        ! call this only to compute z_i, will ignore all other output
        call supercell_test(lon,lat,p_i(i,j,k),z_i(i,j,k),zcoords,u(i,j,1),v(i,j,1),T(i,j,1),&
             thetav(i,j,1),ps(i,j),rho(i,j,1),q(i,j,1,1),pert)

        w_i (i,j,k)  = 0 ! no vertical motion

      enddo; enddo
    enddo
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

        call set_tracers(q(i,j,k,:),3,dp(i,j,k),i,j,k,lat,lon,elem(ie))

      enddo; enddo
    enddo
    call set_elem_state(u,v,w,w_i,T,ps,phis,p,dp,z,z_i,g,elem(ie),1,nt,ntQ=1)

    ! set density to ensure hydrostatic balance and save initial state
    call tests_finalize(elem(ie),hvcoord,1,nt,ie)

  enddo

  sample_period = 60.0 ! sec
end subroutine

!_______________________________________________________________________
subroutine dcmip2016_append_measurements(max_w,max_precl,min_ps,tl,hybrid)
  real(rl),           intent(in) :: max_w, max_precl, min_ps
  type(TimeLevel_t),  intent(in) :: tl
  type(hybrid_t),     intent(in) :: hybrid                   ! hybrid parallel structure

  character(len=*), parameter :: w_filename     = "measurement_wmax.txt"
  character(len=*), parameter :: precl_filename = "measurement_prect_rate.txt"
  character(len=*), parameter :: time_filename  = "measurement_time.txt"
  character(len=*), parameter :: ps_filename    = "measurement_psmin.txt"

  real(rl) :: pmax_w, pmax_precl, pmin_ps, time
  real(rl) :: next_sample_time = 0.0

  time = time_at(tl%nstep)
  ! initialize output file
  if(next_sample_time == 0.0 .and. hybrid%masterthread) then
    open(unit=10,file=w_filename,    form="formatted",status="UNKNOWN")
    close(10)

    open(unit=11,file=precl_filename,form="formatted",status="UNKNOWN")
    close(11)

    open(unit=12,file=time_filename, form="formatted",status="UNKNOWN")
    close(12)

    open(unit=12,file=ps_filename, form="formatted",status="UNKNOWN")
    close(13)
  endif
  ! append measurements at regular intervals
  if(time .ge. next_sample_time) then
!$OMP BARRIER
!$OMP MASTER
    next_sample_time = next_sample_time + sample_period
!$OMP END MASTER
!$OMP BARRIER
    pmax_w     = parallelMax(max_w,    hybrid)
    pmax_precl = parallelMax(max_precl,hybrid)
    pmin_ps    = parallelMin(min_ps,   hybrid)

    if (hybrid%masterthread) then
      print *,"time=",time_at(tl%nstep)," pmax_w (m/s)=",pmax_w," pmax_precl (mm/day)=",pmax_precl*(1000.0)*(24.0*3600)," pmin_ps (Pa)=",pmin_ps

      open(unit=10,file=w_filename,form="formatted",position="append")
        write(10,'(99E24.15)') pmax_w
      close(10)

      open(unit=11,file=precl_filename,form="formatted",position="append")
        write(11,'(99E24.15)') pmax_precl
      close(11)

      open(unit=12,file=time_filename,form="formatted",position="append")
        write(12,'(99E24.15)') time
      close(12)

      open(unit=13,file=ps_filename,form="formatted",position="append")
        write(13,'(99E24.15)') pmin_ps
      close(13)

    endif
  endif

  end subroutine

!_______________________________________________________________________
subroutine dcmip2016_test1_forcing(elem,hybrid,hvcoord,nets,nete,nt,ntQ,dt,tl)

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(in)            :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  integer,            intent(in)            :: nt, ntQ                  ! time level index
  real(rl),           intent(in)            :: dt                       ! time-step size
  type(TimeLevel_t),  intent(in)            :: tl                       ! time level structure

  integer :: i,j,k,ie                                                     ! loop indices
  real(rl), dimension(np,np,nlev) :: u,v,w,T,exner_kess,theta_kess,p,dp,rho,z,qv,qc,qr
  real(rl), dimension(np,np,nlev) :: u0,v0,T0,qv0,qc0,qr0,cl,cl2,ddt_cl,ddt_cl2
  real(rl), dimension(np,np,nlev) :: rho_dry,rho_new,Rstar,p_pk
  real(rl), dimension(nlev)       :: u_c,v_c,p_c,qv_c,qc_c,qr_c,rho_c,z_c, th_c
  real(rl) :: max_w, max_precl, min_ps
  real(rl) :: lat, lon, dz_top(np,np), zi(np,np,nlevp),zi_c(nlevp), ps(np,np)

  integer :: pbl_type, prec_type, qi
  integer, parameter :: test = 1

  prec_type = dcmip16_prec_type
  pbl_type  = dcmip16_pbl_type

  max_w     = -huge(rl)
  max_precl = -huge(rl)
  min_ps    = +huge(rl)

  do ie = nets,nete

    precl(:,:,ie) = -1.0d0

    ! get current element state
    call get_state(u,v,w,T,p,dp,ps,rho,z,zi,g,elem(ie),hvcoord,nt,ntQ)

    ! compute form of exner pressure expected by Kessler physics
    exner_kess = (p/p0)**(Rgas/Cp)
    theta_kess = T/exner_kess

    ! get mixing ratios
    ! use qi to avoid compiler warnings when qsize_d<5
    qi=1;  qv  = elem(ie)%state%Qdp(:,:,:,qi,ntQ)/dp
    qi=2;  qc  = elem(ie)%state%Qdp(:,:,:,qi,ntQ)/dp
    qi=3;  qr  = elem(ie)%state%Qdp(:,:,:,qi,ntQ)/dp
    qi=4;  cl  = elem(ie)%state%Qdp(:,:,:,qi,ntQ)/dp
    qi=5;  cl2 = elem(ie)%state%Qdp(:,:,:,qi,ntQ)/dp

    ! ensure positivity
    where(qv<0); qv=0; endwhere
    where(qc<0); qc=0; endwhere
    where(qr<0); qr=0; endwhere


    rho_dry = (1-qv)*rho  ! convert to dry density using wet mixing ratio

    ! convert to dry mixing ratios
    qv  = qv*rho/rho_dry
    qc  = qc*rho/rho_dry
    qr  = qr*rho/rho_dry

    ! save un-forced prognostics
    u0=u; v0=v; T0=T; qv0=qv; qc0=qc; qr0=qr

    ! apply forcing to columns
    do j=1,np; do i=1,np

      ! invert column
      u_c  = u  (i,j,nlev:1:-1)
      v_c  = v  (i,j,nlev:1:-1)
      qv_c = qv (i,j,nlev:1:-1)
      qc_c = qc (i,j,nlev:1:-1)
      qr_c = qr (i,j,nlev:1:-1)
      p_c  = p  (i,j,nlev:1:-1)
      rho_c= rho_dry(i,j,nlev:1:-1)
      z_c  = z  (i,j,nlev:1:-1)
      zi_c = zi (i,j,nlevp:1:-1)
      th_c = theta_kess(i,j,nlev:1:-1)

      ! get forced versions of u,v,p,qv,qc,qr. rho is constant
      call DCMIP2016_PHYSICS(test, u_c, v_c, p_c, th_c, qv_c, qc_c, qr_c, rho_c, dt, z_c, zi_c, lat, nlev, &
                             precl(i,j,ie), pbl_type, prec_type)

      ! revert column
      u(i,j,:)  = u_c(nlev:1:-1)
      v(i,j,:)  = v_c(nlev:1:-1)
      p(i,j,:)  = p_c(nlev:1:-1)
      qv(i,j,:) = qv_c(nlev:1:-1)
      qc(i,j,:) = qc_c(nlev:1:-1)
      qr(i,j,:) = qr_c(nlev:1:-1)
      theta_kess(i,j,:) = th_c(nlev:1:-1)

      lon = elem(ie)%spherep(i,j)%lon
      lat = elem(ie)%spherep(i,j)%lat

      do k=1,nlev
        call tendency_terminator( lat*rad2dg, lon*rad2dg, cl(i,j,k), cl2(i,j,k), dt, ddt_cl(i,j,k), ddt_cl2(i,j,k))
      enddo

    enddo; enddo;

    if (theta_hydrostatic_mode) then
       ! hydrostatic model assumes physics does not change pressure
       ! so assume T,PHI change, with P held fixed
       ! ps_v will be adjusted after physics to conserve dry mass
    else
       rho_new = rho_dry*(1+qv)
       Rstar = (Rgas+(Rwater_vapor-Rgas)*qv*rho_dry/rho_new)
       p_pk = rho_new*Rstar*theta_kess
       exner_kess = ( p_pk / p0)**( (Rgas/Cp) / ( 1 - (Rgas/Cp)))
    endif
    T = exner_kess*theta_kess

    ! set dynamics forcing
    elem(ie)%derived%FM(:,:,1,:) = (u - u0)/dt
    elem(ie)%derived%FM(:,:,2,:) = (v - v0)/dt
    elem(ie)%derived%FT(:,:,:)   = (T - T0)/dt

    ! set tracer-mass forcing. conserve tracer mass
    elem(ie)%derived%FQ(:,:,:,1) = (rho_dry/rho)*dp*(qv-qv0)/dt
    elem(ie)%derived%FQ(:,:,:,2) = (rho_dry/rho)*dp*(qc-qc0)/dt
    elem(ie)%derived%FQ(:,:,:,3) = (rho_dry/rho)*dp*(qr-qr0)/dt

    qi=4; elem(ie)%derived%FQ(:,:,:,qi) = dp*ddt_cl
    qi=5; elem(ie)%derived%FQ(:,:,:,qi) = dp*ddt_cl2

    ! perform measurements of max w, and max prect
    max_w     = max( max_w    , maxval(w    ) )
    max_precl = max( max_precl, maxval(precl(:,:,ie)) )
    min_ps    = min( min_ps,    minval(ps) )

  enddo

  call dcmip2016_append_measurements(max_w,max_precl,min_ps,tl,hybrid)

end subroutine

!_______________________________________________________________________
subroutine dcmip2016_test2_forcing(elem,hybrid,hvcoord,nets,nete,nt,ntQ,dt,tl, test)

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(in)            :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  integer,            intent(in)            :: nt, ntQ                  ! time level index
  real(rl),           intent(in)            :: dt                       ! time-step size
  type(TimeLevel_t),  intent(in)            :: tl                       ! time level structure
  integer,            intent(in)            :: test                     ! dcmip2016 test number

  integer :: i,j,k,ie                                                     ! loop indices
  real(rl), dimension(np,np,nlev) :: u,v,w,T,exner_kess,theta_kess,p,dp,rho,z,qv,qc,qr
  real(rl), dimension(np,np,nlev) :: u0,v0,T0,qv0,qc0,qr0
  real(rl), dimension(np,np,nlev) :: rho_dry,rho_new,Rstar,p_pk
  real(rl), dimension(nlev)       :: u_c,v_c,p_c,qv_c,qc_c,qr_c,rho_c,z_c, th_c
  real(rl) :: max_w, max_precl, min_ps
  real(rl) :: lat, dz_top(np,np), zi(np,np,nlevp),zi_c(nlevp), ps(np,np)

  integer :: pbl_type, prec_type

  prec_type = dcmip16_prec_type
  pbl_type  = dcmip16_pbl_type

  if(test==3) prec_type=0 ! kessler only for test 3

  max_w     = -huge(rl)
  max_precl = -huge(rl)
  min_ps    = +huge(rl)

  do ie = nets,nete

    precl(:,:,ie) = -1.0d0

    ! get current element state
    call get_state(u,v,w,T,p,dp,ps,rho,z,zi,g,elem(ie),hvcoord,nt,ntQ)

    ! compute form of exner pressure expected by Kessler physics
    exner_kess = (p/p0)**(Rgas/Cp)
    theta_kess = T/exner_kess

    ! get wet mixing ratios
    qv  = elem(ie)%state%Qdp(:,:,:,1,ntQ)/dp
    qc  = elem(ie)%state%Qdp(:,:,:,2,ntQ)/dp
    qr  = elem(ie)%state%Qdp(:,:,:,3,ntQ)/dp

    ! ensure positivity
    where(qv<0); qv=0; endwhere
    where(qc<0); qc=0; endwhere
    where(qr<0); qr=0; endwhere

    rho_dry = (1-qv)*rho  ! convert to dry density using wet mixing ratio

    ! convert to dry mixing ratios
    qv  = qv*rho/rho_dry
    qc  = qc*rho/rho_dry
    qr  = qr*rho/rho_dry


    ! save un-forced prognostics (DRY)
    u0=u; v0=v; T0=T; qv0=qv; qc0=qc; qr0=qr


    ! apply forcing to columns
    do j=1,np; do i=1,np

      ! invert column
      u_c  = u  (i,j,nlev:1:-1)
      v_c  = v  (i,j,nlev:1:-1)
      qv_c = qv (i,j,nlev:1:-1)
      qc_c = qc (i,j,nlev:1:-1)
      qr_c = qr (i,j,nlev:1:-1)
      p_c  = p  (i,j,nlev:1:-1)
      rho_c= rho_dry(i,j,nlev:1:-1)
      z_c  = z  (i,j,nlev:1:-1)
      zi_c = zi (i,j,nlevp:1:-1)
      th_c = theta_kess(i,j,nlev:1:-1)

      ! get forced versions of u,v,p,qv,qc,qr. rho is constant
      call DCMIP2016_PHYSICS(test, u_c, v_c, p_c, th_c, qv_c, qc_c, qr_c, rho_c, dt, z_c, zi_c, lat, nlev, &
                             precl(i,j,ie), pbl_type, prec_type)

      ! revert column
      u(i,j,:)  = u_c(nlev:1:-1)
      v(i,j,:)  = v_c(nlev:1:-1)
      p(i,j,:)  = p_c(nlev:1:-1)
      qv(i,j,:) = qv_c(nlev:1:-1)
      qc(i,j,:) = qc_c(nlev:1:-1)
      qr(i,j,:) = qr_c(nlev:1:-1)
      theta_kess(i,j,:) = th_c(nlev:1:-1)

    enddo; enddo;
    if (theta_hydrostatic_mode) then
       ! hydrostatic model assumes physics does not change pressure
       ! so assume T,PHI change, with P held fixed
       ! ps_v will be adjusted after physics to conserve dry mass
    else
       rho_new = rho_dry*(1+qv)
       Rstar = (Rgas+(Rwater_vapor-Rgas)*qv*rho_dry/rho_new)
       p_pk = rho_new*Rstar*theta_kess
       exner_kess = ( p_pk / p0)**( (Rgas/Cp) / ( 1 - (Rgas/Cp)))
    endif
    T = exner_kess*theta_kess


    ! set dynamics forcing
    elem(ie)%derived%FM(:,:,1,:) = (u - u0)/dt
    elem(ie)%derived%FM(:,:,2,:) = (v - v0)/dt
    elem(ie)%derived%FT(:,:,:)   = (T-T0)/dt


    ! set tracer-mass forcing. conserve tracer mass
    elem(ie)%derived%FQ(:,:,:,1) = (rho_dry/rho)*dp*(qv-qv0)/dt
    elem(ie)%derived%FQ(:,:,:,2) = (rho_dry/rho)*dp*(qc-qc0)/dt
    elem(ie)%derived%FQ(:,:,:,3) = (rho_dry/rho)*dp*(qr-qr0)/dt


    ! perform measurements of max w, and max prect
    max_w     = max( max_w    , maxval(w    ) )
    max_precl = max( max_precl, maxval(precl(:,:,ie)) )
    min_ps    = min( min_ps,    minval(ps) )

  enddo

  call dcmip2016_append_measurements(max_w,max_precl,min_ps,tl,hybrid)

end subroutine dcmip2016_test2_forcing

!_______________________________________________________________________
subroutine dcmip2016_test3_forcing(elem,hybrid,hvcoord,nets,nete,nt,ntQ,dt,tl)

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(in)            :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  integer,            intent(in)            :: nt, ntQ                  ! time level index
  real(rl),           intent(in)            :: dt                       ! time-step size
  type(TimeLevel_t),  intent(in)            :: tl                       ! time level structure

  integer :: i,j,k,ie                                                     ! loop indices
  real(rl):: lat
  real(rl), dimension(np,np,nlev) :: u,v,w,T,theta_kess,exner_kess,p,dp,rho,z,qv,qc,qr
  real(rl), dimension(np,np,nlev) :: T0,qv0,qc0,qr0
  real(rl), dimension(np,np,nlev) :: rho_dry,rho_new,Rstar,p_pk
  real(rl), dimension(np,np,nlev) :: theta_inv,qv_inv,qc_inv,qr_inv,rho_inv,exner_inv,z_inv ! inverted columns
  real(rl), dimension(np,np,nlevp):: zi
  real(rl), dimension(np,np)      :: ps
  real(rl) :: max_w, max_precl, min_ps

  max_w     = -huge(rl)
  max_precl = -huge(rl)
  min_ps    = +huge(rl)

  do ie = nets,nete

    precl(:,:,ie) = 0.0

    ! get current element state
    call get_state(u,v,w,T,p,dp,ps,rho,z,zi,g,elem(ie),hvcoord,nt,ntQ)

    ! get mixing ratios
    qv  = elem(ie)%state%Qdp(:,:,:,1,ntQ)/dp
    qc  = elem(ie)%state%Qdp(:,:,:,2,ntQ)/dp
    qr  = elem(ie)%state%Qdp(:,:,:,3,ntQ)/dp

    ! ensure positivity
    where(qv<0); qv=0; endwhere
    where(qc<0); qc=0; endwhere
    where(qr<0); qr=0; endwhere

    rho_dry = (1-qv)*rho  ! convert to dry density using wet mixing ratio

    ! convert to dry mixing ratios
    qv  = qv*rho/rho_dry
    qc  = qc*rho/rho_dry
    qr  = qr*rho/rho_dry


    ! compute form of exner pressure expected by Kessler physics
    exner_kess = (p/p0)**(Rgas/Cp)
    theta_kess = T/exner_kess

    ! save un-forced prognostics
    T0=T; qv0=qv; qc0=qc; qr0=qr


    ! invert columns (increasing z)
    theta_inv= theta_kess(:,:,nlev:1:-1)
    qv_inv   = qv   (:,:,nlev:1:-1)
    qc_inv   = qc   (:,:,nlev:1:-1)
    qr_inv   = qr   (:,:,nlev:1:-1)
    rho_inv  = rho_dry  (:,:,nlev:1:-1)
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
        precl(i,j,ie))

    enddo; enddo;

    ! revert columns (increasing eta)
    theta_kess = theta_inv(:,:,nlev:1:-1)
    qv    = qv_inv   (:,:,nlev:1:-1)
    qc    = qc_inv   (:,:,nlev:1:-1)
    qr    = qr_inv   (:,:,nlev:1:-1)


    if (theta_hydrostatic_mode) then
       ! hydrostatic model assumes physics does not change pressure
       ! so assume T,PHI change, with P held fixed
       ! ps_v will be adjusted after physics to conserve dry mass
    else
       rho_new = rho_dry*(1+qv)
       Rstar = (Rgas+(Rwater_vapor-Rgas)*qv*rho_dry/rho_new)
       p_pk = rho_new*Rstar*theta_kess
       exner_kess = ( p_pk / p0)**( (Rgas/Cp) / ( 1 - (Rgas/Cp)))
    endif
    T = exner_kess*theta_kess

    ! set dynamics forcing
    elem(ie)%derived%FM(:,:,1,:) = 0
    elem(ie)%derived%FM(:,:,2,:) = 0
    elem(ie)%derived%FT(:,:,:)   = (T-T0)/dt

    ! set tracer-mass forcing. conserve tracer mass
    elem(ie)%derived%FQ(:,:,:,1) = (rho_dry/rho)*dp*(qv-qv0)/dt
    elem(ie)%derived%FQ(:,:,:,2) = (rho_dry/rho)*dp*(qc-qc0)/dt
    elem(ie)%derived%FQ(:,:,:,3) = (rho_dry/rho)*dp*(qr-qr0)/dt

    ! perform measurements of max w, and max prect
    max_w     = max( max_w    , maxval(w    ) )
    max_precl = max( max_precl, maxval(precl(:,:,ie)) )
    min_ps    = min( min_ps,    minval(ps) )

  enddo

  call dcmip2016_append_measurements(max_w,max_precl,min_ps,tl,hybrid)

end subroutine
end module dcmip16_wrapper
#endif
