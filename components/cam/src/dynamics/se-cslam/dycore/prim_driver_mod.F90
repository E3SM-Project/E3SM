!#define _DBG_ print *,"file: ",__FILE__," line: ",__LINE__," ithr: ",hybrid%ithr
#define _DBG_
module prim_driver_mod
  use shr_kind_mod,           only: r8=>shr_kind_r8
  use cam_logfile,            only: iulog
  use dimensions_mod,         only: np, nlev, nelem, nelemd, GlobalUniqueCols, qsize, nc,nhc
  use hybrid_mod,             only: hybrid_t, config_thread_region, PrintHybrid
  use derivative_mod,         only: derivative_t
  use fvm_control_volume_mod, only: fvm_struct

  use element_mod,            only: element_t, timelevels, allocate_element_desc
  use thread_mod ,            only: horz_num_threads, vert_num_threads, tracer_num_threads
  use perf_mod,               only: t_startf, t_stopf
  use prim_init,              only: gp, fvm_corners, fvm_points

  implicit none
  private
  public :: prim_init2, prim_run_subcycle, prim_finalize
  public :: prim_set_dry_mass

contains

!=============================================================================!

  subroutine prim_init2(elem, fvm, hybrid, nets, nete, tl, hvcoord)
    use dimensions_mod,         only: irecons_tracer
    use dimensions_mod,         only: fv_nphys, ntrac, nc
    use cam_abortutils,         only: endrun
    use parallel_mod,           only: syncmp
    use time_mod,               only: timelevel_t, tstep, phys_tscale, nsplit, TimeLevel_Qdp
    use prim_state_mod,         only: prim_printstate
    use control_mod,            only: runtype, &
         topology, rsplit, qsplit, rk_stage_user,         &
         nu, nu_q, nu_div, hypervis_subcycle, hypervis_subcycle_q
    use fvm_control_volume_mod, only: fvm_supercycling,n0_fvm
    use fvm_mod,                only: fill_halo_fvm,ghostBufQnhc
    use thread_mod,             only: omp_get_thread_num
    use global_norms_mod,       only: test_global_integral, print_cfl
    use hybvcoord_mod,          only: hvcoord_t
    use prim_advection_mod,     only: prim_advec_init2,deriv
    use prim_advance_mod,       only: prim_advance_init, compute_omega

    type (element_t), intent(inout) :: elem(:)
    type (fvm_struct), intent(inout)    :: fvm(:)
    type (hybrid_t), intent(in) :: hybrid

    type (TimeLevel_t), intent(inout)    :: tl              ! time level struct
    type (hvcoord_t), intent(inout)      :: hvcoord         ! hybrid vertical coordinate struct

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)


    ! ==================================
    ! Local variables
    ! ==================================

!   variables used to calculate CFL
    real (kind=r8) :: dtnu            ! timestep*viscosity parameter
    real (kind=r8) :: dt_dyn_vis      ! viscosity timestep used in dynamics
    real (kind=r8) :: dt_tracer_vis      ! viscosity timestep used in tracers

    real (kind=r8) :: dp

    integer :: i,j,k,ie,t,q
    integer :: n0,n0_qdp


    do ie=nets,nete
      elem(ie)%derived%FM=0.0_r8
      elem(ie)%derived%FT=0.0_r8
      elem(ie)%derived%FQ=0.0_r8
    end do

    ! ==========================
    ! begin executable code
    ! ==========================
    call prim_advance_init(hybrid%par,elem)

    if (topology == "cube") then
       call test_global_integral(elem, hybrid,nets,nete)
    end if


    ! compute most restrictive dt*nu for use by variable res viscosity:
    ! compute timestep seen by viscosity operator:
    dt_dyn_vis = tstep
    dt_tracer_vis=tstep*qsplit

    ! compute most restrictive condition:
    ! note: dtnu ignores subcycling
    dtnu=max(dt_dyn_vis*max(nu,nu_div), dt_tracer_vis*nu_q)
    ! compute actual viscosity timesteps with subcycling
    dt_tracer_vis = dt_tracer_vis/hypervis_subcycle_q
    dt_dyn_vis = dt_dyn_vis/hypervis_subcycle

    ! ==================================
    ! Initialize derivative structure
    ! ==================================
    call Prim_Advec_Init2(fvm_corners, fvm_points)
    if (fv_nphys>0.and.nc.ne.fv_nphys) then
      !
      ! need to fill halo for dp_coupling for fvm2phys mapping
      !
      call fill_halo_fvm(ghostBufQnhc,elem,fvm,hybrid,nets,nete,n0_fvm,nhc,1,nlev)
    end if
    !$OMP BARRIER
    if (hybrid%ithr==0) then
       call syncmp(hybrid%par)
    end if
    !$OMP BARRIER

    if (topology /= "cube") then
       call endrun('Error: only cube topology supported for primaitve equations')
    endif

    ! timesteps to use for advective stability:  tstep*qsplit and tstep
    call print_cfl(elem,hybrid,nets,nete,dtnu)

    if (hybrid%masterthread) then
       ! CAM has set tstep based on dtime before calling prim_init2(),
       ! so only now does HOMME learn the timstep.  print them out:
       write(iulog,'(a,2f9.2)') "dt_remap: (0=disabled)   ",tstep*qsplit*rsplit

       if (ntrac>0) then
          write(iulog,'(a,2f9.2)') "dt_tracer (fvm)          ",tstep*qsplit*fvm_supercycling
       end if
       if (qsize>0) then
          write(iulog,'(a,2f9.2)') "dt_tracer (SE), per RK stage: ",tstep*qsplit,(tstep*qsplit)/(rk_stage_user-1)
       end if
       write(iulog,'(a,2f9.2)')    "dt_dyn:                  ",tstep
       write(iulog,'(a,2f9.2)')    "dt_dyn (viscosity):      ",dt_dyn_vis
       write(iulog,'(a,2f9.2)')    "dt_tracer (viscosity):   ",dt_tracer_vis


       if (phys_tscale/=0) then
          write(iulog,'(a,2f9.2)') "CAM physics timescale:       ",phys_tscale
       endif
       write(iulog,'(a,2f9.2)') "CAM dtime (dt_phys):         ",tstep*nsplit*qsplit*max(rsplit,1)

       write(iulog,*) "CAM-SE uses dry-mass vertical coordinates"
     end if

     n0=tl%n0
     call TimeLevel_Qdp( tl, qsplit, n0_qdp)
     call compute_omega(hybrid,n0,n0_qdp,elem,deriv,nets,nete,tstep,hvcoord)

    if (hybrid%masterthread) write(iulog,*) "initial state:"
    call prim_printstate(elem, tl, hybrid,nets,nete, fvm)

  end subroutine prim_init2

!=======================================================================================================!


  subroutine prim_run_subcycle(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord,nsubstep)
!
!   advance all variables (u,v,T,ps,Q,C) from time t to t + dt_q
!
! for the RK schemes:
!   input:
!       tl%nm1   not used
!       tl%n0    data at time t
!       tl%np1   new values at t+dt_q
!
!   then we update timelevel pointers:
!       tl%nm1 = tl%n0
!       tl%n0  = tl%np1
!   so that:
!       tl%nm1   tracers:  t    dynamics:  t+(qsplit-1)*dt
!       tl%n0    time t + dt_q
!
! for the implicit schemes:
!
!   input:
!       tl%nm1   variables at t-1 level are stored fro BDF2 scheme
!       tl%n0    data at time t
!       tl%np1   new values at t+dt_q
!       generally dt_q = t for BDF2, so its t+1
!
!   then we update timelevel pointers:
!       tl%nm1 = tl%n0
!       tl%n0  = tl%np1
!   so that:
!       tl%nm1   tracers:  t    dynamics:  t+(qsplit-1)*dt
!       tl%n0    time t + dt_q
!
!
    use hybvcoord_mod, only : hvcoord_t
    use time_mod,               only: TimeLevel_t, timelevel_update, timelevel_qdp, nsplit
    use control_mod,            only: statefreq,disable_diagnostics,qsplit, rsplit
    use prim_advance_mod,       only: applycamforcing
    use prim_advance_mod,       only: calc_tot_energy_dynamics,compute_omega
    use prim_state_mod,         only: prim_printstate
    use prim_advection_mod,     only: vertical_remap, deriv
    use fvm_control_volume_mod, only: n0_fvm
    use thread_mod,             only: omp_get_thread_num
    use perf_mod   ,            only: t_startf, t_stopf
    use fvm_mod    ,            only: fill_halo_fvm, ghostBufQnhc
    use dimensions_mod,         only: ntrac,fv_nphys

    type (element_t) , intent(inout) :: elem(:)
    type(fvm_struct), intent(inout)  :: fvm(:)
    type (hybrid_t), intent(in)      :: hybrid  ! distributed parallel structure (shared)
    type (hvcoord_t), intent(in)     :: hvcoord         ! hybrid vertical coordinate struct
    integer, intent(in)              :: nets  ! starting thread element number (private)
    integer, intent(in)              :: nete  ! ending thread element number   (private)
    real(kind=r8), intent(in)        :: dt  ! "timestep dependent" timestep
    type (TimeLevel_t), intent(inout):: tl
    integer, intent(in)              :: nsubstep  ! nsubstep = 1 .. nsplit

    real(kind=r8)   :: dt_q, dt_remap
    integer         :: ie, q,k,n0_qdp,np1_qdp,r, nstep_end,region_num_threads
    real (kind=r8)  :: dp_np1(np,np)
    logical         :: compute_diagnostics
!    type (hybrid_t) :: vybrid

    ! ===================================
    ! Main timestepping loop
    ! ===================================
    dt_q = dt*qsplit
    dt_remap = dt_q
    nstep_end = tl%nstep + qsplit
    dt_remap=dt_q*rsplit
    nstep_end = tl%nstep + qsplit*rsplit  ! nstep at end of this routine

    ! compute diagnostics for STDOUT
    compute_diagnostics=.false.

    if (statefreq>0) then
      if (MODULO(nstep_end,statefreq)==0 .or. nstep_end==tl%nstep0) then
        compute_diagnostics=.true.
      endif
    end if

    if(disable_diagnostics) compute_diagnostics=.false.


    call TimeLevel_Qdp( tl, qsplit, n0_qdp)

    call calc_tot_energy_dynamics(elem,fvm,nets,nete,tl%n0,n0_qdp,n0_fvm,'dAF')
    call ApplyCAMForcing(elem,fvm,tl%n0,n0_qdp,dt_remap,nets,nete,nsubstep)


    ! loop over rsplit vertically lagrangian timesteps
    call prim_step(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord,1)
    do r=2,rsplit
       call TimeLevel_update(tl,"leapfrog")
       call prim_step(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord,r)
    enddo
    ! defer final timelevel update until after remap and diagnostics


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !  apply vertical remap
    !  always for tracers
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !compute timelevels for tracers (no longer the same as dynamics)
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp)
    ! note: time level update for fvm tracers takes place in fvm_mod

    call calc_tot_energy_dynamics(elem,fvm,nets,nete,tl%np1,np1_qdp,n0_fvm,'dAD')

    call t_startf('vertical_remap')
    call vertical_remap(hybrid,elem,fvm,hvcoord,dt_remap,tl%np1,np1_qdp,n0_fvm,nets,nete)
    call t_stopf('vertical_remap')



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! time step is complete.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call calc_tot_energy_dynamics(elem,fvm,nets,nete,tl%np1,np1_qdp,n0_fvm,'dAR')

    if (nsubstep==nsplit) &
         call compute_omega(hybrid,tl%np1,np1_qdp,elem,deriv,nets,nete,dt,hvcoord)

    ! now we have:
    !   u(nm1)   dynamics at  t+dt_remap - 2*dt
    !   u(n0)    dynamics at  t+dt_remap - dt
    !   u(np1)   dynamics at  t+dt_remap
    !
    !   Q(1)   Q at t+dt_remap



    ! =================================
    ! update dynamics time level pointers
    ! =================================
    call TimeLevel_update(tl,"leapfrog")
    ! note: time level update for fvm tracers takes place in fvm_mod

    ! now we have:
    !   u(nm1)   dynamics at  t+dt_remap - dt       (Robert-filtered)
    !   u(n0)    dynamics at  t+dt_remap
    !   u(np1)   undefined


    ! ============================================================
    ! Print some diagnostic information
    ! ============================================================
    if (compute_diagnostics) then
       call prim_printstate(elem, tl, hybrid,nets,nete, fvm)
    end if

    if (ntrac>0.and.nsubstep==nsplit.and.nc.ne.fv_nphys) then
      !
      ! fill the fvm halo for mapping in d_p_coupling if
      ! physics grid resolution is different than fvm resolution
      !
      call fill_halo_fvm(ghostBufQnhc, elem,fvm,hybrid,nets,nete,n0_fvm,nhc,1,nlev)
    end if

  end subroutine prim_run_subcycle


  subroutine prim_step(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord, rstep)
    !
    !   Take qsplit dynamics steps and one tracer step
    !   for vertically lagrangian option, this subroutine does only the horizontal step
    !
    !   input:
    !       tl%nm1   not used
    !       tl%n0    data at time t
    !       tl%np1   new values at t+dt_q
    !
    !   then we update timelevel pointers:
    !       tl%nm1 = tl%n0
    !       tl%n0  = tl%np1
    !   so that:
    !       tl%nm1   tracers:  t    dynamics:  t+(qsplit-1)*dt
    !       tl%n0    time t + dt_q
    !
    use hybvcoord_mod,          only: hvcoord_t
    use time_mod,               only: TimeLevel_t, timelevel_update
    use control_mod,            only: statefreq, qsplit, nu_p
    use control_mod,            only: TRACERTRANSPORT_CONSISTENT_SE_FVM, tracer_transport_type
    use thread_mod,             only: omp_get_thread_num
    use prim_advance_mod,       only: prim_advance_exp
    use prim_advection_mod,     only: prim_advec_tracers_remap, prim_advec_tracers_fvm, deriv
    use derivative_mod,         only: subcell_integration
    use fvm_control_volume_mod, only: fvm_supercycling
    use hybrid_mod,             only: set_region_num_threads, config_thread_region
    use dimensions_mod,         only: ntrac

    type (element_t) ,  intent(inout) :: elem(:)
    type(fvm_struct),   intent(inout) :: fvm(:)
    type (hybrid_t),    intent(in)    :: hybrid  ! distributed parallel structure (shared)
    type (hvcoord_t),   intent(in)    :: hvcoord         ! hybrid vertical coordinate struct
    integer,            intent(in)    :: nets  ! starting thread element number (private)
    integer,            intent(in)    :: nete  ! ending thread element number   (private)
    real(kind=r8),      intent(in)    :: dt  ! "timestep dependent" timestep
    type (TimeLevel_t), intent(inout) :: tl
    integer, intent(in)               :: rstep ! vertical remap subcycling step

    type (hybrid_t):: hybridnew
    real(kind=r8)  :: st, st1, dp, dt_q
    integer        :: ie,t,q,k,i,j,n, n_Q
    integer        :: ithr
    integer        :: region_num_threads

    real (kind=r8) :: tempdp3d(np,np), x
    real (kind=r8) :: tempmass(nc,nc)
    real (kind=r8) :: tempflux(nc,nc,4)

    real (kind=r8) :: dp_np1(np,np)

    dt_q = dt*qsplit
    if (ntrac>0.and.rstep==1) then
      do ie=nets,nete
        elem(ie)%sub_elem_mass_flux=0
      end do
    end if

    ! ===============
    ! initialize mean flux accumulation variables and save some variables at n0
    ! for use by advection
    ! ===============
    do ie=nets,nete
      elem(ie)%derived%vn0=0              ! mean horizontal mass flux
      elem(ie)%derived%omega=0
      if (nu_p>0) then
         elem(ie)%derived%dpdiss_ave=0
         elem(ie)%derived%dpdiss_biharmonic=0
      endif

      ! dp at time t:  use floating lagrangian levels:
      elem(ie)%derived%dp(:,:,:)=elem(ie)%state%dp3d(:,:,:,tl%n0)
    enddo

    ! ===============
    ! Dynamical Step
    ! ===============
    n_Q = tl%n0  ! n_Q = timelevel of FV tracers at time t.  need to save this
                 ! FV tracers still carry 3 timelevels
                 ! SE tracers only carry 2 timelevels

    call t_startf('prim_advance_exp')
!    ithr   = 0 ! omp_get_thread_num()
!    vybrid = hybrid_create(hybrid%par,ithr)

    call prim_advance_exp(elem, fvm, deriv, hvcoord,   &
         hybrid, dt, tl, nets, nete)

    call t_stopf('prim_advance_exp')

    do n=2,qsplit
       call TimeLevel_update(tl,"leapfrog")

       call t_startf('prim_advance_exp')

       call prim_advance_exp(elem, fvm, deriv, hvcoord,   &
            hybrid, dt, tl, nets, nete)

    call t_stopf('prim_advance_exp')

       ! defer final timelevel update until after Q update.
    enddo
#ifdef HOMME_TEST_SUB_ELEMENT_MASS_FLUX
    if (ntrac>0.and.rstep==1) then
      do ie=nets,nete
      do k=1,nlev
        tempdp3d = elem(ie)%state%dp3d(:,:,k,tl%np1) - &
                   elem(ie)%derived%dp(:,:,k)
        call subcell_integration(tempdp3d, np, nc, elem(ie)%metdet,tempmass)
        tempflux = dt_q*elem(ie)%sub_elem_mass_flux(:,:,:,k)
        do i=1,nc
        do j=1,nc
          x = SUM(tempflux(i,j,:))
          if (ABS(tempmass(i,j)).lt.1e-11_r8 .and. 1e-11_r8.lt.ABS(x)) then
            print *,__FILE__,__LINE__,"**********",ie,k,i,j,tempmass(i,j),x
          elseif (1e-5_r8.lt.ABS((tempmass(i,j)-x)/tempmass(i,j))) then
            print *,__FILE__,__LINE__,"**********",ie,k,i,j,tempmass(i,j),x,&
                   ABS((tempmass(i,j)-x)/tempmass(i,j))
          endif
        end do
        end do
      end do
      end do
    end if
#endif

    ! current dynamics state variables:
    !    derived%dp              =  dp at start of timestep
    !    derived%vn0             =  mean horiz. flux:   U*dp
    ! rsplit>0
    !        state%v(:,:,:,np1)      = velocity on lagrangian levels
    !        state%dp3d(:,:,:,np1)   = dp3d
    !


    ! ===============
    ! Tracer Advection.
    ! in addition, this routine will apply the DSS to:
    !        derived%omega           =
    ! Tracers are always vertically lagrangian.
    ! ===============
    ! Advect tracers if their count is > 0.
    ! special case in CAM: if CSLAM tracers are turned on , qsize=1 but this tracer should
    ! not be advected.  This will be cleaned up when the physgrid is merged into CAM trunk
    ! Currently advecting all species
    if (qsize > 0) then

      call t_startf('prim_advec_tracers_remap')
      region_num_threads = tracer_num_threads
#ifdef _OPENMP
      call omp_set_nested(.true.)
#endif
!JMD      !$OMP PARALLEL NUM_THREADS(region_num_threads), DEFAULT(SHARED), PRIVATE(hybridnew)
!JMD     hybridnew = config_thread_region(hybrid,'tracer')
      call Prim_Advec_Tracers_remap(elem, deriv,hvcoord,hybrid,dt_q,tl,nets,nete)
!JMD      !$OMP END PARALLEL
#ifdef _OPENMP
      call omp_set_nested(.false.)
#endif
      call t_stopf('prim_advec_tracers_remap')
    end if
    !
    ! only run fvm transport every fvm_supercycling rstep
    !
    if (ntrac>0 .and. (mod(rstep,fvm_supercycling) == 0)) then
       !
       ! FVM transport
       !
      if (tracer_transport_type == TRACERTRANSPORT_CONSISTENT_SE_FVM) &
      call Prim_Advec_Tracers_fvm(elem,fvm,hvcoord,hybrid,&
           dt_q,tl,nets,nete)
    endif

  end subroutine prim_step


!=======================================================================================================!


  subroutine prim_finalize(hybrid)
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)

    ! ==========================
    ! end of the hybrid program
    ! ==========================
  end subroutine prim_finalize

!=========================================================================================

    subroutine prim_set_dry_mass(elem, hvcoord,initial_global_ave_dry_ps,q)
      use element_mod,      only: element_t
      use hybvcoord_mod ,   only: hvcoord_t
      use dimensions_mod,   only: nelemd, nlev, np
      use constituents,     only: cnst_type, qmin, pcnst
      use cam_logfile,      only: iulog
      use spmd_utils,       only: masterproc

      type (element_t)     , intent(inout):: elem(:)
      type (hvcoord_t)     , intent(in)   :: hvcoord
      real (kind=r8), intent(in)   :: initial_global_ave_dry_ps
      real (kind=r8), intent(inout):: q(np,np,nlev,nelemd,pcnst)

      ! local
      real (kind=r8)               :: global_ave_ps_inic,dp_tmp, factor(np,np,nlev)
      integer                      :: ie, i, j ,k, m_cnst

      if (initial_global_ave_dry_ps == 0) return;

      call get_global_ave_surface_pressure(elem, global_ave_ps_inic)

      do ie=1,nelemd
        elem(ie)%state%psdry(:,:)=elem(ie)%state%psdry(:,:)*(initial_global_ave_dry_ps/global_ave_ps_inic)
        do k=1,nlev
          do j = 1,np
            do i = 1,np
              dp_tmp =  ((hvcoord%hyai(k+1) - hvcoord%hyai(k))*hvcoord%ps0)+&
                   ((hvcoord%hybi(k+1) - hvcoord%hybi(k))*elem(ie)%state%psdry(i,j))
              factor(i,j,k) = elem(ie)%state%dp3d(i,j,k,1)/dp_tmp
              elem(ie)%state%dp3d(i,j,k,:) = dp_tmp
            end do
          end do
        end do
        !
        ! conserve initial condition mass of 'wet' tracers (following dryairm.F90 for FV dycore)
        ! and conserve mixing ratio (not mass) of 'dry' tracers
        !
        do  m_cnst=1,pcnst
          if (cnst_type(m_cnst).ne.'dry') then
            do k=1,nlev
              do j = 1,np
                do i = 1,np
                  q(i,j,k,ie,m_cnst) = q(i,j,k,ie,m_cnst)*factor(i,j,k)
                  q(i,j,k,ie,m_cnst) = max(qmin(m_cnst),q(i,j,k,ie,m_cnst))
                end do
              end do
            end do
          end if
        end do
      end do
      if (masterproc) then
        write (iulog,*) "------ info from prim_set_dry_mass -----------------------------------------------------------"
        write (iulog,*) "Scaling dry surface pressure to global average of = ",&
             initial_global_ave_dry_ps/100.0_r8,"hPa"
        write (iulog,*) "Average dry surface pressure in initial condition = ",global_ave_ps_inic/100.0_r8,"hPa"
        write (iulog,*) "Average dry surface pressure change               = ",&
             initial_global_ave_dry_ps-global_ave_ps_inic,"Pa"
        write (iulog,*) "Mixing ratios that are wet have been scaled so that total mass of tracer is conserved"
        write (iulog,*) "Mixing ratios that are dry have not been changed (mass not conserved in scaling process)"
        write (iulog,*) "------ end info from prim_set_dry_mass -------------------------------------------------------"
      endif
    end subroutine prim_set_dry_mass

    subroutine get_global_ave_surface_pressure(elem, global_ave_ps_inic)
      use element_mod       , only : element_t
      use dimensions_mod    , only : np
      use global_norms_mod  , only : global_integral
      use hybrid_mod        , only : config_thread_region, get_loop_ranges, hybrid_t
      use parallel_mod      , only : par

      type (element_t)     , intent(in)   :: elem(:)
      real (kind=r8), intent(out)  :: global_ave_ps_inic

      ! local
      real (kind=r8), allocatable  :: tmp(:,:,:)
      type (hybrid_t)                     :: hybrid
      integer                             :: ie, nets, nete

      !JMD $OMP PARALLEL NUM_THREADS(horz_num_threads), DEFAULT(SHARED), PRIVATE(hybrid,nets,nete,n)
      !JMD        hybrid = config_thread_region(par,'horizontal')
      hybrid = config_thread_region(par,'serial')
      call get_loop_ranges(hybrid,ibeg=nets,iend=nete)
      allocate(tmp(np,np,nets:nete))

      do ie=nets,nete
        tmp(:,:,ie)=elem(ie)%state%psdry(:,:)
      enddo
      global_ave_ps_inic = global_integral(elem, tmp(:,:,nets:nete),hybrid,np,nets,nete)
      deallocate(tmp)
    end subroutine get_global_ave_surface_pressure
end module prim_driver_mod
