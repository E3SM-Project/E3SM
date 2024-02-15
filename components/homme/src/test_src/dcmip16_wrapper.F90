



#ifndef CAM
#include "config.h"

module dcmip16_wrapper

! Implementation of the dcmip2012 dycore tests for the preqx dynamics target

use dcmip12_wrapper,      only: pressure_thickness, set_tracers, get_evenly_spaced_z, set_hybrid_coefficients
use control_mod,          only: test_case, dcmip16_pbl_type, dcmip16_prec_type, use_moisture, theta_hydrostatic_mode,&
     sub_case, case_planar_bubble, bubble_prec_type
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
     set_forcing_rayleigh_friction
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

real(rl) :: sample_period  = 2.0_rl
real(rl) :: rad2dg = 180.0_rl/pi

type :: PhysgridData_t
   integer :: nphys
   real(rl), allocatable :: ps(:,:), zs(:,:), T(:,:,:), uv(:,:,:,:), omega_p(:,:,:), q(:,:,:,:)
end type PhysgridData_t

type (PhysgridData_t) :: pg_data

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

  if (qsize<6) call abortmp('ERROR: test requires qsize>=6')

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
    call tests_finalize(elem(ie),hvcoord)

  enddo
  sample_period = 1800.0 ! sec
  !print *,"min thetav = ",min_thetav, "max thetav=",max_thetav


end subroutine

subroutine dcmip2016_pg_init(elem,hybrid,hvcoord,nets,nete,nphys)
  use gllfvremap_mod, only: gfr_init

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  integer,            intent(in)            :: nphys                    ! pgN, N parameter, for physgrid

  integer :: ncol

  if (hybrid%ithr == 0) then
     ncol = nphys*nphys
     pg_data%nphys = nphys
     call gfr_init(hybrid%par, elem, nphys, boost_pg1=.true.)
     allocate(pg_data%ps(ncol,nelemd), pg_data%zs(ncol,nelemd), pg_data%T(ncol,nlev,nelemd), &
          pg_data%omega_p(ncol,nlev,nelemd), pg_data%uv(ncol,2,nlev,nelemd), &
          pg_data%q(ncol,nlev,qsize,nelemd))
  end if
  !$omp barrier
end subroutine dcmip2016_pg_init

subroutine dcmip2016_test1_pg(elem,hybrid,hvcoord,nets,nete,nphys)
  use gllfvremap_mod, only: gfr_init

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  integer,            intent(in)            :: nphys                    ! pgN, N parameter, for physgrid

  call dcmip2016_pg_init(elem,hybrid,hvcoord,nets,nete,nphys)
  call dcmip2016_test1(elem,hybrid,hvcoord,nets,nete)
  sample_period = 3600*24
end subroutine dcmip2016_test1_pg

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
    call tests_finalize(elem(ie),hvcoord)
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

  integer :: i,j,k,ie,imod                                              ! loop indices
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
     imod=max(1,(nete-nets)/50) ! limit output to 50 lines. 
     !if (hybrid%masterthread) write(*,"(A,I5,A)",advance="NO") " ie=",ie,achar(13)
     if (hybrid%masterthread .and. mod(ie,imod)==0) &
          write(*,"(A,2I5)") " ie=",ie,nete

    do k=1,nlevp

      do j=1,np; do i=1,np

        lon = elem(ie)%spherep(i,j)%lon
        lat = elem(ie)%spherep(i,j)%lat
        if (sub_case==2) lon=mod(lon+pi,2*pi)  ! shift initial condition

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
        if (sub_case==2) lon=mod(lon+pi,2*pi)  ! shift initial condition
        
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
    call tests_finalize(elem(ie),hvcoord,ie)

  enddo

  sample_period = 1 ! 60 orig sec
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
      print *,"time=",time_at(tl%nstep)," pmax_w (m/s)=",pmax_w
      print *,"time=",time_at(tl%nstep)," pmax_precl (mm/day)=",pmax_precl*(1000.0)*(24.0*3600)
      print *,"time=",time_at(tl%nstep)," pmin_ps (Pa)=",pmin_ps

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
  real(rl), dimension(np,np)      :: delta_ps(np,np)
  real(rl) :: max_w, max_precl, min_ps
  real(rl) :: lat, lon, dz_top(np,np), zi(np,np,nlevp),zi_c(nlevp), ps(np,np)

  integer :: pbl_type, prec_type, qi
  integer, parameter :: test = 1
  logical :: toy_chemistry_on

  if (case_planar_bubble) then
    toy_chemistry_on = .false.
    prec_type = bubble_prec_type
    if (qsize .ne. 3) call abortmp('ERROR: moist bubble test requires qsize=3')
  else
    toy_chemistry_on = .true.
    prec_type = dcmip16_prec_type
  endif

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
    if (toy_chemistry_on) then
      qi=4;  cl  = elem(ie)%state%Qdp(:,:,:,qi,ntQ)/dp
      qi=5;  cl2 = elem(ie)%state%Qdp(:,:,:,qi,ntQ)/dp
    endif

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
      if (case_planar_bubble) then
         lat = 0
      else
         lat = elem(ie)%spherep(i,j)%lat
      end if
      call DCMIP2016_PHYSICS(test, u_c, v_c, p_c, th_c, qv_c, qc_c, qr_c, rho_c, dt, z_c, zi_c, lat, nlev, &
                             precl(i,j,ie), pbl_type, prec_type)

      ! revert column
      u(i,j,:)  = u_c(nlev:1:-1)
      v(i,j,:)  = v_c(nlev:1:-1)
      qv(i,j,:) = qv_c(nlev:1:-1)
      qc(i,j,:) = qc_c(nlev:1:-1)
      qr(i,j,:) = qr_c(nlev:1:-1)
      theta_kess(i,j,:) = th_c(nlev:1:-1)

      if (toy_chemistry_on) then 
        lon = elem(ie)%spherep(i,j)%lon
        lat = elem(ie)%spherep(i,j)%lat

        do k=1,nlev
          call tendency_terminator( lat*rad2dg, lon*rad2dg, cl(i,j,k), cl2(i,j,k), dt, ddt_cl(i,j,k), ddt_cl2(i,j,k))
        enddo
      endif


    enddo; enddo;

    ! convert from theta to T w.r.t. new model state
    ! assume hydrostatic pressure pi changed by qv forcing
    ! assume NH pressure perturbation unchanged
    delta_ps = sum( (rho_dry/rho)*dp*(qv-qv0) , 3 )
    do k=1,nlev
       p(:,:,k) = p(:,:,k) + hvcoord%hybm(k)*delta_ps(:,:)
    enddo
    exner_kess = (p/p0)**(Rgas/Cp)
    T = exner_kess*theta_kess

    ! set dynamics forcing
    elem(ie)%derived%FM(:,:,1,:) = (u - u0)/dt
    elem(ie)%derived%FM(:,:,2,:) = (v - v0)/dt
    elem(ie)%derived%FT(:,:,:)   = (T - T0)/dt

    ! set tracer-mass forcing. conserve tracer mass
    ! rho_dry*(qv-qv0)*dz = FQ deta, dz/deta = -dp/(g*rho)
    elem(ie)%derived%FQ(:,:,:,1) = (rho_dry/rho)*dp*(qv-qv0)/dt
    elem(ie)%derived%FQ(:,:,:,2) = (rho_dry/rho)*dp*(qc-qc0)/dt
    elem(ie)%derived%FQ(:,:,:,3) = (rho_dry/rho)*dp*(qr-qr0)/dt


    if (toy_chemistry_on) then 
    qi=4; elem(ie)%derived%FQ(:,:,:,qi) = dp*ddt_cl
    qi=5; elem(ie)%derived%FQ(:,:,:,qi) = dp*ddt_cl2
    endif


    ! perform measurements of max w, and max prect
    max_w     = max( max_w    , maxval(w    ) )
    max_precl = max( max_precl, maxval(precl(:,:,ie)) )
    min_ps    = min( min_ps,    minval(ps) )

  enddo

  call dcmip2016_append_measurements(max_w,max_precl,min_ps,tl,hybrid)

end subroutine

subroutine toy_init(rcd)
  real(rl), intent(inout) :: rcd(6)
  rcd(1) = 1; rcd(2) = -1; rcd(3) = 1; rcd(4) = -1; rcd(5) = 1; rcd(6) = -1
end subroutine toy_init

subroutine toy_rcd(q, rcd)
  real(rl), intent(in) :: q(:,:,:,:)
  real(rl), intent(inout) :: rcd(6)

  rcd(1) = min(rcd(1), minval(q(:,:,:,1)))
  rcd(2) = max(rcd(2), maxval(q(:,:,:,1)))
  rcd(3) = min(rcd(3), minval(q(:,:,:,2)))
  rcd(4) = max(rcd(4), maxval(q(:,:,:,2)))
  rcd(5) = min(rcd(5), minval(q(:,:,:,1) + 2*q(:,:,:,2)))
  rcd(6) = max(rcd(6), maxval(q(:,:,:,1) + 2*q(:,:,:,2)))
end subroutine toy_rcd

subroutine toy_print(hybrid, nstep, rcd)
  use reduction_mod, only: ParallelMin, ParallelMax

  type(hybrid_t), intent(in) :: hybrid
  integer, intent(in) :: nstep
  real(rl), intent(inout) :: rcd(6)

  if (modulo(nstep,50) == 0) then
     rcd(1) = ParallelMin(rcd(1), hybrid)
     rcd(2) = ParallelMax(rcd(2), hybrid)
     rcd(3) = ParallelMin(rcd(3), hybrid)
     rcd(4) = ParallelMax(rcd(4), hybrid)
     rcd(5) = ParallelMin(rcd(5), hybrid)
     rcd(6) = ParallelMax(rcd(6), hybrid)
     if (hybrid%masterthread) &
          write(*,'(a,i5,es11.3,es11.3,es11.3,es11.3,es11.3,es11.3)') 'toy>', nstep, rcd
  end if
end subroutine toy_print

subroutine dcmip2016_test1_pg_forcing(elem,hybrid,hvcoord,nets,nete,nt,ntQ,dt,tl)
  use element_ops, only: get_field
  use gllfvremap_mod
  use perf_mod, only: t_startf, t_stopf
  use control_mod, only: ftype

  ! to DSS precl
  use edge_mod, only: edgevpack_nlyr, edgevunpack_nlyr, edge_g
  use bndry_mod, only: bndry_exchangev

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(in)            :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  integer,            intent(in)            :: nt, ntQ                  ! time level index
  real(rl),           intent(in)            :: dt                       ! time-step size
  type(TimeLevel_t),  intent(in)            :: tl                       ! time level structure

  integer, parameter :: iqv = 1
  integer, parameter :: test = 1

  integer :: i,j,k,ie,qi
  real(rl), dimension(np,np,nlev) :: w,dp
  real(rl), dimension(nlev)       :: u_c,v_c,p_c,qv_c,qc_c,qr_c,rho_c,z_c,th_c
  real(rl) :: max_w, max_precl, min_ps
  real(rl) :: lat, lon, zi_c(nlevp), wrk3(np,np,nlev), wrk4(np,np,nlev,2), wf(np*np,1)

  integer :: nf, ncol
  real(rl), dimension(np,np,nlev) :: dp_fv, p_fv, u_fv, v_fv, T_fv, exner_kess_fv, &
       theta_kess_fv, Rstar, rho_fv, rho_dry_fv, u0, v0, T0, z_fv, ddt_cl, ddt_cl2
  real(rl), dimension(np,np,nlev,qsize) :: Q_fv, Q0_fv
  real(rl), dimension(np,np,nlevp) :: phi_i, zi_fv
  real(rl), dimension(np,np) :: zs_fv, ps_fv, delta_ps
  real(rl) :: precl_fv(np,np,1), rcd(6)
  real(rl), allocatable :: qmin(:,:,:), qmax(:,:,:)
  integer :: pbl_type, prec_type
  logical :: toy_chemistry_on

  nf = pg_data%nphys
  ncol = nf*nf

  if (case_planar_bubble) then
     toy_chemistry_on = .false.
     prec_type = bubble_prec_type
     if (qsize .ne. 3) call abortmp('ERROR: moist bubble test requires qsize=3')
  else
     toy_chemistry_on = .true.
     prec_type = dcmip16_prec_type
  endif

  pbl_type  = dcmip16_pbl_type

  max_w     = -huge(rl)
  max_precl = -huge(rl)
  min_ps    = +huge(rl)

  call t_startf('gfr_dyn_to_fv_phys')
  call gfr_dyn_to_fv_phys(hybrid, nt, hvcoord, elem, nets, nete, &
       pg_data%ps, pg_data%zs, pg_data%T, pg_data%uv, pg_data%omega_p, pg_data%q)
  call t_stopf('gfr_dyn_to_fv_phys')

  do ie = nets,nete
     ! for max_w
     call get_field(elem(ie), 'w', w, hvcoord, nt, ntQ)

     T_fv(:nf,:nf,:) = reshape(pg_data%T(:,:,ie), (/nf,nf,nlev/))
     u_fv(:nf,:nf,:) = reshape(pg_data%uv(:,1,:,ie), (/nf,nf,nlev/))
     v_fv(:nf,:nf,:) = reshape(pg_data%uv(:,2,:,ie), (/nf,nf,nlev/))
     Q_fv(:nf,:nf,:,:) = reshape(pg_data%q(:,:,:,ie), (/nf,nf,nlev,qsize/))
     ps_fv(:nf,:nf) = reshape(pg_data%ps(:,ie), (/nf,nf/))
     zs_fv(:nf,:nf) = reshape(pg_data%zs(:,ie), (/nf,nf/))
     zs_fv(:nf,:nf) = zs_fv(:nf,:nf)/g
     do k = 1,nlev
        p_fv(:nf,:nf,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*ps_fv(:nf,:nf)
        dp_fv(:nf,:nf,k) = (hvcoord%hyai(k+1) - hvcoord%hyai(k))*hvcoord%ps0 + &
                           (hvcoord%hybi(k+1) - hvcoord%hybi(k))*ps_fv(:nf,:nf)
     end do

     ! Rederive the remaining vars so they are self-consistent; use
     ! hydrostatic assumption.
     if (use_moisture) then
        Rstar(:nf,:nf,:) = Rgas + (Rwater_vapor - Rgas)*Q_fv(:nf,:nf,:,iqv)
     else
        Rstar(:nf,:nf,:) = Rgas
     end if
     rho_fv(:nf,:nf,:) = p_fv(:nf,:nf,:)/(Rstar(:nf,:nf,:)*T_fv(:nf,:nf,:))
     phi_i(:nf,:nf,nlevp) = g*zs_fv(:nf,:nf)
     do k = nlev,1,-1
        phi_i(:nf,:nf,k) = phi_i(:nf,:nf,k+1) + &
             (Rstar(:nf,:nf,k)*(dp_fv(:nf,:nf,k)*T_fv(:nf,:nf,k)))/p_fv(:nf,:nf,k)
     end do
     do k=1,nlev
        z_fv(:nf,:nf,k) = (phi_i(:nf,:nf,k)+phi_i(:nf,:nf,k+1))/(2*g)
     end do
     do k=1,nlevp
        zi_fv(:nf,:nf,k) = phi_i(:nf,:nf,k)/g
     end do

     rho_dry_fv(:nf,:nf,:) = (1 - Q_fv(:nf,:nf,:,iqv))*rho_fv(:nf,:nf,:)

     ! Compute form of exner pressure expected by Kessler physics.
     exner_kess_fv(:nf,:nf,:) = (p_fv(:nf,:nf,:)/p0)**(Rgas/Cp)
     theta_kess_fv(:nf,:nf,:) = T_fv(:nf,:nf,:)/exner_kess_fv(:nf,:nf,:)

     u0(:nf,:nf,:) = u_fv(:nf,:nf,:); v0(:nf,:nf,:) = v_fv(:nf,:nf,:); T0(:nf,:nf,:) = T_fv(:nf,:nf,:)
     Q0_fv(:nf,:nf,:,:) = Q_fv(:nf,:nf,:,:);

     ! Convert to dry mixing ratios.
     do i = 1,3
        Q_fv(:nf,:nf,:,i) = (rho_fv(:nf,:nf,:)/rho_dry_fv(:nf,:nf,:))*Q_fv(:nf,:nf,:,i)
     end do

     do j = 1,nf
        do i = 1,nf
           u_c  = u_fv(i,j,nlev:1:-1)
           v_c  = v_fv(i,j,nlev:1:-1)
           qv_c = Q_fv(i,j,nlev:1:-1,1)
           qc_c = Q_fv(i,j,nlev:1:-1,2)
           qr_c = Q_fv(i,j,nlev:1:-1,3)
           p_c  = p_fv(i,j,nlev:1:-1)
           rho_c= rho_dry_fv(i,j,nlev:1:-1)
           z_c  = z_fv(i,j,nlev:1:-1)
           zi_c = zi_fv(i,j,nlevp:1:-1)
           th_c = theta_kess_fv(i,j,nlev:1:-1)

           ! Get forced versions of u,v,p,qv,qc,qr. rho is constant.
           if (case_planar_bubble) then
              lat = 0
           else
              call gfr_f_get_latlon(ie, i, j, lat, lon)
           end if
           call DCMIP2016_PHYSICS(test, u_c, v_c, p_c, th_c, qv_c, qc_c, qr_c, rho_c, dt, &
                z_c, zi_c, lat, nlev, precl_fv(i,j,1), pbl_type, prec_type)

           u_fv(i,j,:)   = u_c(nlev:1:-1)
           v_fv(i,j,:)   = v_c(nlev:1:-1)
           Q_fv(i,j,:,1) = qv_c(nlev:1:-1)
           Q_fv(i,j,:,2) = qc_c(nlev:1:-1)
           Q_fv(i,j,:,3) = qr_c(nlev:1:-1)
           theta_kess_fv(i,j,:) = th_c(nlev:1:-1)

           if (toy_chemistry_on) then
              call gfr_f_get_latlon(ie, i, j, lat, lon)
              qi = 4
              do k=1,nlev
                 call tendency_terminator(lat*rad2dg, lon*rad2dg, Q_fv(i,j,k,qi), Q_fv(i,j,k,qi+1), &
                      dt, ddt_cl(i,j,k), ddt_cl2(i,j,k))
              end do
           end if
        enddo
     enddo

     do i = 1,3
        Q_fv(:nf,:nf,:,i) = (rho_dry_fv(:nf,:nf,:)/rho_fv(:nf,:nf,:))*Q_fv(:nf,:nf,:,i)
     end do
     if (toy_chemistry_on) then
        qi = 4
        Q_fv(:nf,:nf,:,qi  ) = Q_fv(:nf,:nf,:,qi  ) + dt*ddt_cl (:nf,:nf,:)
        Q_fv(:nf,:nf,:,qi+1) = Q_fv(:nf,:nf,:,qi+1) + dt*ddt_cl2(:nf,:nf,:)
     end if

     ! Convert from theta to T w.r.t. new model state.
     ! Assume hydrostatic pressure pi changed by qv forcing.
     ! Assume NH pressure perturbation is unchanged.
     delta_ps(:nf,:nf) = sum(dp_fv(:nf,:nf,:)*(Q_fv(:nf,:nf,:,iqv) - Q0_fv(:nf,:nf,:,iqv)), 3)
     do k = 1,nlev
        p_fv(:nf,:nf,k) = p_fv(:nf,:nf,k) + hvcoord%hybm(k)*delta_ps(:nf,:nf)
     enddo
     exner_kess_fv(:nf,:nf,:) = (p_fv(:nf,:nf,:)/p0)**(Rgas/Cp)
     T_fv(:nf,:nf,:) = exner_kess_fv(:nf,:nf,:)*theta_kess_fv(:nf,:nf,:)

     ! These gfr calls are special to this routine, to handle
     ! DCMIP-specific precl.
     wf(:ncol,1) = reshape(precl_fv(:nf,:nf,1), (/ncol/))
     call gfr_f2g_scalar(ie, elem(ie)%metdet, wf(:,:1), wrk3(:,:,:1))
     call gfr_g_make_nonnegative(elem(ie)%metdet, wrk3(:,:,:1))
     precl(:,:,ie) = wrk3(:,:,1)

     ! T, uv tendencies
     pg_data%T(:,:,ie) = reshape((T_fv(:nf,:nf,:) - T0(:nf,:nf,:))/dt, (/ncol,nlev/))
     pg_data%uv(:,1,:,ie) = reshape((u_fv(:nf,:nf,:) - u0(:nf,:nf,:))/dt, (/ncol,nlev/))
     pg_data%uv(:,2,:,ie) = reshape((v_fv(:nf,:nf,:) - v0(:nf,:nf,:))/dt, (/ncol,nlev/))
     ! q state
     do i = 1,qsize
        pg_data%q(:,:,i,ie) = reshape(Q_fv(:nf,:nf,:,i), (/ncol,nlev/))
     end do
     
     ! Measure max w and max prect. w is not used in the physics, so
     ! just look at the GLL values.
     max_w = max(max_w, maxval(w))
     ! ps isn't updated by the physics, so just look at the GLL values.
     min_ps = min(min_ps, minval(elem(ie)%state%ps_v(:,:,nt)))
  enddo

  call t_startf('gfr_fv_phys_to_dyn')
  call gfr_fv_phys_to_dyn(hybrid, nt, hvcoord, elem, nets, nete, &
       pg_data%T, pg_data%uv, pg_data%q)
  call t_stopf('gfr_fv_phys_to_dyn')
  ! dp_coupling doesn't do the DSS; stepon does. Thus, this DCMIP test
  ! also needs to do its own DSS.
  call gfr_f2g_dss(hybrid, elem, nets, nete)
  call gfr_pg1_reconstruct(hybrid, nt, hvcoord, elem, nets, nete)

  if (toy_chemistry_on) then
     call toy_init(rcd)
     do ie = nets,nete
        do i = 1,2
           wrk4(:,:,:,i) = elem(ie)%state%Q(:,:,:,i+3)
        end do
        call toy_rcd(wrk4, rcd)
     end do
     call toy_print(hybrid, tl%nstep, rcd)
  end if

  ! In standalone Homme, for all ftype values, FQ is tendency.
  if (.true.) then !(ftype == 0) then
     ! Convert FQ from state to Qdp tendency.
     do ie = nets,nete
        do k = 1,nlev
           dp(:,:,k) = (hvcoord%hyai(k+1) - hvcoord%hyai(k))*hvcoord%ps0 + &
                       (hvcoord%hybi(k+1) - hvcoord%hybi(k))*elem(ie)%state%ps_v(:,:,nt)
        end do
        do i = 1,qsize
           elem(ie)%derived%FQ(:,:,:,i) = &
                dp*(elem(ie)%derived%FQ(:,:,:,i) - elem(ie)%state%Q(:,:,:,i))/dt
        end do
     end do
  end if

  ! DSS precl
  do ie = nets,nete
     precl(:,:,ie) = precl(:,:,ie)*elem(ie)%spheremp
     call edgeVpack_nlyr(edge_g, elem(ie)%desc, precl(:,:,ie), 1, 0, 1)
  end do
  call bndry_exchangeV(hybrid, edge_g)
  do ie = nets,nete
     call edgeVunpack_nlyr(edge_g, elem(ie)%desc, precl(:,:,ie), 1, 0, 1)
     precl(:,:,ie) = precl(:,:,ie)*elem(ie)%rspheremp
     max_precl = max( max_precl, maxval(precl(:,:,ie)) )
  end do

  call dcmip2016_append_measurements(max_w,max_precl,min_ps,tl,hybrid)
end subroutine dcmip2016_test1_pg_forcing

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
  real(rl), dimension(np,np)      :: delta_ps(np,np)
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

      lat=0.0 ! unused in test 2
      ! get forced versions of u,v,p,qv,qc,qr. rho is constant
      call DCMIP2016_PHYSICS(test, u_c, v_c, p_c, th_c, qv_c, qc_c, qr_c, rho_c, dt, z_c, zi_c, lat, nlev, &
                             precl(i,j,ie), pbl_type, prec_type)

      ! revert column
      u(i,j,:)  = u_c(nlev:1:-1)
      v(i,j,:)  = v_c(nlev:1:-1)
      qv(i,j,:) = qv_c(nlev:1:-1)
      qc(i,j,:) = qc_c(nlev:1:-1)
      qr(i,j,:) = qr_c(nlev:1:-1)
      theta_kess(i,j,:) = th_c(nlev:1:-1)

    enddo; enddo;
    ! convert from theta to T w.r.t. new model state
    ! assume hydrostatic pressure pi changed by qv forcing
    ! assume NH pressure perturbation unchanged
    delta_ps = sum( (rho_dry/rho)*dp*(qv-qv0) , 3 )
    do k=1,nlev
       p(:,:,k) = p(:,:,k) + hvcoord%hybm(k)*delta_ps(:,:)
    enddo
    exner_kess = (p/p0)**(Rgas/Cp)
    T = exner_kess*theta_kess


    ! set dynamics forcing
    elem(ie)%derived%FM(:,:,1,:) = (u - u0)/dt
    elem(ie)%derived%FM(:,:,2,:) = (v - v0)/dt
    elem(ie)%derived%FT(:,:,:)   = (T-T0)/dt


    ! set tracer-mass forcing. conserve tracer mass
    ! rho_dry*(qv-qv0)*dz = FQ deta, dz/deta = -dp/(g*rho)
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
  real(rl), dimension(np,np)      :: ps,delta_ps(np,np)
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


    ! convert from theta to T w.r.t. new model state
    ! assume hydrostatic pressure pi changed by qv forcing
    ! assume NH pressure perturbation unchanged
    delta_ps = sum( (rho_dry/rho)*dp*(qv-qv0) , 3 )
    do k=1,nlev
       p(:,:,k) = p(:,:,k) + hvcoord%hybm(k)*delta_ps(:,:)
    enddo
    exner_kess = (p/p0)**(Rgas/Cp)
    T = exner_kess*theta_kess

    ! set dynamics forcing
    elem(ie)%derived%FM(:,:,1,:) = 0
    elem(ie)%derived%FM(:,:,2,:) = 0
    elem(ie)%derived%FT(:,:,:)   = (T-T0)/dt

    ! set tracer-mass forcing. conserve tracer mass.  
    ! rho_dry*(qv-qv0)*dz = FQ deta, dz/deta = -dp/(g*rho)
    elem(ie)%derived%FQ(:,:,:,1) = (rho_dry/rho)*dp*(qv-qv0)/dt
    elem(ie)%derived%FQ(:,:,:,2) = (rho_dry/rho)*dp*(qc-qc0)/dt
    elem(ie)%derived%FQ(:,:,:,3) = (rho_dry/rho)*dp*(qr-qr0)/dt
    ! after forcings above are applied:
    ! Qdp_new = (rho_dry/rho)*dp*qv   


    ! perform measurements of max w, and max prect
    max_w     = max( max_w    , maxval(w    ) )
    max_precl = max( max_precl, maxval(precl(:,:,ie)) )
    min_ps    = min( min_ps,    minval(ps) )

  enddo

  call dcmip2016_append_measurements(max_w,max_precl,min_ps,tl,hybrid)

end subroutine
end module dcmip16_wrapper
#endif
