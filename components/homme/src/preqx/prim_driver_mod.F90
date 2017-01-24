module prim_driver_mod

  use prim_driver_base, only:&
    deriv1,&
    prim_init1,&
    prim_init2 ,&
    prim_finalize,&
    smooth_topo_datasets

  use control_mod,        only: energy_fixer, statefreq, ftype, qsplit, rsplit
  use dimensions_mod,     only: np, nlev, nlevp, nelem, nelemd, qsize
  use element_mod,        only: element_t, allocate_element_desc
  use element_state,      only: timelevels
  use hybrid_mod,         only: hybrid_t
  use hybvcoord_mod,      only: hvcoord_t
  use kinds,              only: real_kind, iulog
  use parallel_mod,       only: abortmp
  use perf_mod,           only: t_startf, t_stopf
  use prim_advection_mod, only: prim_advec_tracers_remap
  use prim_state_mod,     only: prim_printstate, prim_diag_scalars, prim_energy_halftimes
  use time_mod,           only: timeLevel_t, timelevel_update, timelevel_qdp

#ifndef CAM
  use test_mod,           only: compute_test_forcing
#endif

  implicit none
  contains

  subroutine prim_run_subcycle(elem, hybrid,nets,nete, dt, tl, hvcoord,nsubstep)

    !   advance dynamic variables and tracers (u,v,T,ps,Q,C) from time t to t + dt_q
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

    use control_mod,        only: statefreq, ftype, qsplit, rsplit, disable_diagnostics
    use hybvcoord_mod,      only: hvcoord_t
    use parallel_mod,       only: abortmp
    use prim_advance_mod,   only: applycamforcing, applycamforcing_dynamics
    use prim_state_mod,     only: prim_printstate, prim_diag_scalars, prim_energy_halftimes
    use vertremap_mod,      only: vertical_remap
    use reduction_mod,      only: parallelmax
    use time_mod,           only: TimeLevel_t, timelevel_update, timelevel_qdp, nsplit

#if USE_OPENACC
    use openacc_utils_mod,  only: copy_qdp_h2d, copy_qdp_d2h
#endif

    type (element_t) ,    intent(inout) :: elem(:)
    type (hybrid_t),      intent(in)    :: hybrid                       ! distributed parallel structure (shared)
    type (hvcoord_t),     intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct
    integer,              intent(in)    :: nets                         ! starting thread element number (private)
    integer,              intent(in)    :: nete                         ! ending thread element number   (private)
    real(kind=real_kind), intent(in)    :: dt                           ! "timestep dependent" timestep
    type (TimeLevel_t),   intent(inout) :: tl
    integer,              intent(in)    :: nsubstep                     ! nsubstep = 1 .. nsplit

    real(kind=real_kind) :: dp, dt_q, dt_remap
    real(kind=real_kind) :: dp_np1(np,np)
    integer :: ie,i,j,k,n,q,t
    integer :: n0_qdp,np1_qdp,r,nstep_end
    logical :: compute_diagnostics, compute_energy

    ! compute timesteps for tracer transport and vertical remap
    dt_q      = dt*qsplit
    dt_remap  = dt_q
    nstep_end = tl%nstep + qsplit
    if (rsplit>0) then
       dt_remap  = dt_q*rsplit   ! rsplit=0 means use eulerian code, not vert. lagrange
       nstep_end = tl%nstep + qsplit*rsplit  ! nstep at end of this routine
    endif

    ! activate diagnostics periodically for display to stdout
    compute_energy = .false.
    compute_diagnostics   = .false.
    if (MODULO(nstep_end,statefreq)==0 .or. nstep_end==tl%nstep0) then
       compute_diagnostics= .true.
       compute_energy     = .true.
    endif
    if(disable_diagnostics) compute_diagnostics= .false.

    ! compute scalar diagnostics if currently active
    if (compute_diagnostics) then
      call t_startf("prim_diag_scalars")
      call prim_diag_scalars(elem,hvcoord,tl,4,.true.,nets,nete)
      call t_stopf("prim_diag_scalars")
    endif

#ifndef CAM
    ! Get HOMME test case forcing
    call compute_test_forcing(elem,hybrid,hvcoord,tl%n0,n0_qdp,dt_remap,nets,nete)
#endif

    ! Apply CAM Physics forcing

    !   ftype= 2: Q was adjusted by physics, but apply u,T forcing here
    !   ftype= 1: forcing was applied time-split in CAM coupling layer
    !   ftype= 0: apply all forcing here
    !   ftype=-1: do not apply forcing

    call TimeLevel_Qdp(tl, qsplit, n0_qdp)

    if (ftype==0) then
      call t_startf("ApplyCAMForcing")
      call ApplyCAMForcing(elem, hvcoord,tl%n0,n0_qdp, dt_remap,nets,nete)
      call t_stopf("ApplyCAMForcing")

    elseif (ftype==2) then
      call t_startf("ApplyCAMForcing_dynamics")
      call ApplyCAMForcing_dynamics(elem, hvcoord,tl%n0,dt_remap,nets,nete)
      call t_stopf("ApplyCAMForcing_dynamics")
    endif

    ! E(1) Energy after CAM forcing
    if (compute_energy) then
      call t_startf("prim_energy_halftimes")
      call prim_energy_halftimes(elem,hvcoord,tl,1,.true.,nets,nete)
      call t_stopf("prim_energy_halftimes")
    endif

    ! qmass and variance, using Q(n0),Qdp(n0)
    if (compute_diagnostics) then
      call t_startf("prim_diag_scalars")
      call prim_diag_scalars(elem,hvcoord,tl,1,.true.,nets,nete)
      call t_stopf("prim_diag_scalars")
    endif

    ! initialize dp3d from ps
    do ie=nets,nete
       do k=1,nlev
          elem(ie)%state%dp3d(:,:,k,tl%n0)=&
               ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%n0)
       enddo
    enddo

#if (USE_OPENACC)
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp)
    call t_startf("copy_qdp_h2d")
    call copy_qdp_h2d( elem , n0_qdp )
    call t_stopf("copy_qdp_h2d")
#endif

    ! Loop over rsplit vertically lagrangian timesteps
    call t_startf("prim_step_rX")
    call prim_step(elem, hybrid,nets,nete, dt, tl, hvcoord,compute_diagnostics,1)
    call t_stopf("prim_step_rX")

    do r=2,rsplit
       call TimeLevel_update(tl,"leapfrog")
       call t_startf("prim_step_rX")
       call prim_step(elem, hybrid,nets,nete, dt, tl, hvcoord,.false.,r)
       call t_stopf("prim_step_rX")
    enddo
    ! defer final timelevel update until after remap and diagnostics

#if (USE_OPENACC)
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp)
    call t_startf("copy_qdp_h2d")
    call copy_qdp_d2h( elem , np1_qdp )
    call t_stopf("copy_qdp_h2d")
#endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  apply vertical remap
    !  always for tracers
    !  if rsplit>0:  also remap dynamics and compute reference level ps_v
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !compute timelevels for tracers (no longer the same as dynamics)
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp)
    call vertical_remap(hybrid,elem,hvcoord,dt_remap,tl%np1,np1_qdp,nets,nete)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! time step is complete.  update some diagnostic variables:
    ! Q    (mixing ratio)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call t_startf("prim_run_subcyle_diags")
    do ie=nets,nete
#if (defined COLUMN_OPENMP)
       !$omp parallel do default(shared), private(k,q,dp_np1)
#endif
       do k=1,nlev    !  Loop inversion (AAM)
          dp_np1(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%np1)
          !dir$ simd
          do q=1,qsize
             elem(ie)%state%Q(:,:,k,q)=elem(ie)%state%Qdp(:,:,k,q,np1_qdp)/dp_np1(:,:)
          enddo
       enddo
    enddo
    call t_stopf("prim_run_subcyle_diags")

    ! now we have:
    !   u(nm1)   dynamics at  t+dt_remap - 2*dt
    !   u(n0)    dynamics at  t+dt_remap - dt
    !   u(np1)   dynamics at  t+dt_remap
    !
    !   Q(1)   Q at t+dt_remap
    if (compute_diagnostics) then
      call t_startf("prim_diag_scalars")
      call prim_diag_scalars(elem,hvcoord,tl,2,.false.,nets,nete)
      call t_stopf("prim_diag_scalars")
    endif
    if (compute_energy) then
      call t_startf("prim_energy_halftimes")
      call prim_energy_halftimes(elem,hvcoord,tl,2,.false.,nets,nete)
      call t_stopf("prim_energy_halftimes")
    endif

    if (compute_diagnostics) then
      call t_startf("prim_diag_scalars")
      call prim_diag_scalars(elem,hvcoord,tl,3,.false.,nets,nete)
      call t_stopf("prim_diag_scalars")

      call t_startf("prim_energy_halftimes")
      call prim_energy_halftimes(elem,hvcoord,tl,3,.false.,nets,nete)
      call t_stopf("prim_energy_halftimes")
    endif

    ! =================================
    ! update dynamics time level pointers
    ! =================================
    call TimeLevel_update(tl,"leapfrog")
    ! now we have:
    !   u(nm1)   dynamics at  t+dt_remap - dt       
    !   u(n0)    dynamics at  t+dt_remap
    !   u(np1)   undefined

    ! ============================================================
    ! Print some diagnostic information
    ! ============================================================
    if (compute_diagnostics) then
       call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete)
    end if

  end subroutine prim_run_subcycle

  !_____________________________________________________________________
  subroutine prim_step(elem, hybrid,nets,nete, dt, tl, hvcoord, compute_diagnostics,rstep)

    !  Take qsplit dynamics steps and one tracer step
    !  for vertically lagrangian option, this subroutine does only the horizontal step
    !
    !  input:
    !      tl%nm1   not used
    !      tl%n0    data at time t
    !      tl%np1   new values at t+dt_q
    !
    !  then we update timelevel pointers:
    !      tl%nm1 = tl%n0
    !      tl%n0  = tl%np1
    !  so that:
    !      tl%nm1   tracers:  t    dynamics:  t+(qsplit-1)*dt
    !      tl%n0    time t + dt_q

    use control_mod,        only: statefreq, integration, ftype, qsplit, nu_p, rsplit
    use control_mod,        only: use_semi_lagrange_transport
    use hybvcoord_mod,      only : hvcoord_t
    use parallel_mod,       only: abortmp
    use prim_advance_mod,   only: prim_advance_exp
    use prim_advection_mod, only: prim_advec_tracers_remap
    use reduction_mod,      only: parallelmax
    use time_mod,           only: time_at,TimeLevel_t, timelevel_update, nsplit

    type(element_t),      intent(inout) :: elem(:)
    type(hybrid_t),       intent(in)    :: hybrid   ! distributed parallel structure (shared)
    type(hvcoord_t),      intent(in)    :: hvcoord  ! hybrid vertical coordinate struct
    integer,              intent(in)    :: nets     ! starting thread element number (private)
    integer,              intent(in)    :: nete     ! ending thread element number   (private)
    real(kind=real_kind), intent(in)    :: dt       ! "timestep dependent" timestep
    type(TimeLevel_t),    intent(inout) :: tl
    integer,              intent(in)    :: rstep    ! vertical remap subcycling step

    real(kind=real_kind) :: st, st1, dp, dt_q
    real(kind=real_kind) :: maxcflx, maxcfly
    real(kind=real_kind) :: dp_np1(np,np)
    logical :: compute_diagnostics
    integer :: ie, t, q,k,i,j,n, n_Q

    dt_q = dt*qsplit
 
    ! ===============
    ! initialize mean flux accumulation variables and save some variables at n0
    ! for use by advection
    ! ===============
    do ie=nets,nete
      elem(ie)%derived%eta_dot_dpdn=0     ! mean vertical mass flux
      elem(ie)%derived%vn0=0              ! mean horizontal mass flux
      elem(ie)%derived%omega_p=0
      if (nu_p>0) then
         elem(ie)%derived%dpdiss_ave=0
         elem(ie)%derived%dpdiss_biharmonic=0
      endif
      if (use_semi_lagrange_transport) then
        elem(ie)%derived%vstar=elem(ie)%state%v(:,:,:,:,tl%n0)
      end if
      elem(ie)%derived%dp(:,:,:)=elem(ie)%state%dp3d(:,:,:,tl%n0)
    enddo

    ! ===============
    ! Dynamical Step
    ! ===============
    call t_startf("prim_step_dyn")
    n_Q = tl%n0  ! n_Q = timelevel of FV tracers at time t.  need to save this
                 ! FV tracers still carry 3 timelevels
                 ! SE tracers only carry 2 timelevels
    call prim_advance_exp(elem, deriv1, hvcoord,   &
         hybrid, dt, tl, nets, nete, compute_diagnostics)
    do n=2,qsplit
       call TimeLevel_update(tl,"leapfrog")
       call prim_advance_exp(elem, deriv1, hvcoord,hybrid, dt, tl, nets, nete, .false.)
       ! defer final timelevel update until after Q update.
    enddo
    call t_stopf("prim_step_dyn")

    ! current dynamics state variables:
    !    derived%dp              =  dp at start of timestep
    !    derived%vstar           =  velocity at start of tracer timestep
    !    derived%vn0             =  mean horiz. flux:   U*dp
    !    state%dp3d(:,:,:,np1)   = dp3d
    ! rsplit=0
    !        state%v(:,:,:,np1)      = velocity on reference levels
    ! rsplit>0
    !        state%v(:,:,:,np1)      = velocity on lagrangian levels 
    !        
    ! Tracer Advection.  
    ! in addition, this routine will apply the DSS to:
    !        derived%eta_dot_dpdn    =  mean vertical velocity (used for remap below)
    !        derived%omega           =
    ! Tracers are always vertically lagrangian.  
    ! For rsplit=0: 
    !   if tracer scheme needs v on lagrangian levels it has to vertically interpolate
    !   if tracer scheme needs dp3d, it needs to derive it from ps_v

    call t_startf("prim_step_advec")
    if (qsize > 0) then
      call t_startf("PAT_remap")
      call Prim_Advec_Tracers_remap(elem, deriv1,hvcoord,hybrid,dt_q,tl,nets,nete)
      call t_stopf("PAT_remap")
    end if
    call t_stopf("prim_step_advec")

  end subroutine prim_step

end module
