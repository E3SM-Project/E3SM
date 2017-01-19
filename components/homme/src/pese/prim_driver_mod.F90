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
  use prim_advance_mod,   only: compute_and_apply_rhs
  use prim_advection_mod, only: prim_advect_tracers
  use prim_state_mod,     only: prim_printstate, prim_diag_scalars, prim_energy_halftimes
  use time_mod,           only: timeLevel_t, timelevel_update, timelevel_qdp
  use test_mod,           only: compute_test_forcing

  implicit none
  contains

! customized: prim_run_subcycle , prim_step

!_____________________________________________________________________
  subroutine apply_forcing(elem,hybrid,hvcoord,tl,dt_remap,nets,nete, compute_diagnostics, compute_energy)

    ! apply forcing from CAM or HOMME stand-alone tests

    use prim_advance_mod, only: ApplyCAMForcing, ApplyCAMForcing_dynamics

    implicit none
    type(element_t),      intent(inout) :: elem(:)
    type(hybrid_t),       intent(in)    :: hybrid
    type(hvcoord_t),      intent(in)    :: hvcoord
    type (TimeLevel_t),   intent(inout) :: tl
    real(kind=real_kind), intent(in)    :: dt_remap
    integer,              intent(in)    :: nets,nete
    logical,              intent(in)    :: compute_diagnostics
    logical,              intent(in)    :: compute_energy

    integer :: n0_qdp

    call TimeLevel_Qdp(tl, qsplit, n0_qdp)

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

  end subroutine

!_______________________________________________________________________
subroutine prim_run_subcycle(elem,hybrid,nets,nete,dt,tl,hvcoord,nsubstep)

    type (element_t) ,    intent(inout) :: elem(:)
    type (hybrid_t),      intent(in)    :: hybrid                       ! distributed parallel structure (shared)
    type (hvcoord_t),     intent(inout) :: hvcoord                      ! hybrid vertical coordinate struct
    integer,              intent(in)    :: nets                         ! starting thread element number (private)
    integer,              intent(in)    :: nete                         ! ending thread element number   (private)
    real(kind=real_kind), intent(in)    :: dt                           ! "timestep dependent" timestep
    type (TimeLevel_t),   intent(inout) :: tl
    integer,              intent(in)    :: nsubstep                     ! nsubstep = 1 .. nsplit

    logical :: compute_diagnostics, compute_energy
    real(kind=real_kind) :: dt_q, dt_remap
    integer :: n0_qdp,np1_qdp,r,nstep_end


    ! compute dt_q, dt_remp, and nstep_end
    call get_subcycle_stepsize(dt_q,dt_remap,nstep_end,dt,tl)

    ! enable/disable diagnositics for this time-step
    call set_diagnostic_state(compute_diagnostics,compute_energy,tl,nstep_end)

    ! compute scalar diagnostics if currently active
    if (compute_diagnostics) then
      call t_startf("prim_diag_scalars")
      call prim_diag_scalars(elem,hvcoord,tl,3,.true.,nets,nete)
      call t_stopf("prim_diag_scalars")

      call t_startf("prim_energy_halftimes")
      call prim_energy_halftimes(elem,hvcoord,tl,3,.true.,nets,nete)
      call t_stopf("prim_energy_halftimes")
    endif

    ! apply cam forcing or stand-alone test forcing to state variables
    call apply_forcing(elem,hybrid,hvcoord,tl,dt_remap,nets,nete, compute_diagnostics, compute_energy)

    ! get E(1): energy after CAM forcing
    if (compute_energy) call prim_energy_halftimes(elem,hvcoord,tl,1,.true.,nets,nete)

    ! get qmass and variance, using Q(n0),Qdp(n0)
    if (compute_diagnostics) call prim_diag_scalars(elem,hvcoord,tl,1,.true.,nets,nete)

    ! loop over rsplit vertically lagrangian timesteps
    call prim_step(elem, hybrid,nets,nete, dt, tl, hvcoord,compute_diagnostics,1)
    do r=2,rsplit
       call TimeLevel_update(tl,"leapfrog")
       call prim_step(elem, hybrid,nets,nete,dt,tl,hvcoord,.false.,r)
    enddo

    ! compute scalar diagnostics if currently active
    if (compute_diagnostics) then
      ! E(1): Energy after CAM forcing
      call t_startf("prim_energy_halftimes")
      call prim_energy_halftimes(elem,hvcoord,tl,1,.true.,nets,nete)
      call t_stopf("prim_energy_halftimes")

      ! qmass and variance, using Q(n0),Qdp(n0)
      call t_startf("prim_diag_scalars")
      call prim_diag_scalars(elem,hvcoord,tl,1,.true.,nets,nete)
      call t_stopf("prim_diag_scalars")
    endif

    ! update dynamics time level pointers
    call TimeLevel_update(tl,"leapfrog")

    ! print some diagnostic information
    if (compute_diagnostics) call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete)

  end subroutine prim_run_subcycle

!_______________________________________________________________________
subroutine prim_step(elem,hybrid,nets,nete,dt,tl,hvcoord,compute_diagnostics,rstep)

    use control_mod,        only: statefreq, integration, ftype, qsplit, nu_p, rsplit
    use control_mod,        only: use_semi_lagrange_transport
    use hybvcoord_mod,      only: hvcoord_t
    use parallel_mod,       only: abortmp
    use prim_advance_mod,   only: prim_advance_exp
    use time_mod,           only: timelevel_t, timelevel_update

    type(element_t),      intent(inout) :: elem(:)
    type(hybrid_t),       intent(in)    :: hybrid   ! distributed parallel structure (shared)
    type(hvcoord_t),      intent(in)    :: hvcoord  ! hybrid vertical coordinate struct
    integer,              intent(in)    :: nets     ! starting thread element number (private)
    integer,              intent(in)    :: nete     ! ending thread element number   (private)
    real(kind=real_kind), intent(in)    :: dt       ! "timestep dependent" timestep
    type(TimeLevel_t),    intent(inout) :: tl
    logical,              intent(in)    :: compute_diagnostics
    integer,              intent(in)    :: rstep    ! vertical remap subcycling step

    real(kind=real_kind) :: st, st1, dp, dt_q
    integer :: ie, t, q,k,i,j,n, n_Q
    real (kind=real_kind) :: maxcflx, maxcfly
    real (kind=real_kind) :: dp_np1(np,np)

    dt_q = dt*qsplit

    ! clear derived quantities

    call t_startf("prim_step_init")
    dt_q = dt*qsplit

    do ie=nets,nete
      if (nu_p>0) then
         elem(ie)%derived%dpdiss_ave        = 0
         elem(ie)%derived%dpdiss_biharmonic = 0
      endif
    enddo
    call t_stopf("prim_step_init")

    ! take qsplit dynamics steps, followed by one tracer step
    do n=1,qsplit
      call prim_advance_exp(elem, deriv1, hvcoord, hybrid, dt, tl, nets, nete, compute_diagnostics)
      if(n<qsplit) call timelevel_update(tl,"leapfrog")
    enddo

    ! take one tracer timestep
    if (qsize > 0) call prim_advect_tracers(elem, deriv1, hybrid, nets, nete, hvcoord, dt_q, tl)

  end subroutine prim_step

  !_____________________________________________________________________
  subroutine get_subcycle_stepsize(dt_q,dt_remap,nstep_end,dt,tl)

    ! compute timesteps for tracer transport and vertical remap

    use control_mod, only: qsplit, rsplit

    real(kind=real_kind), intent(out) :: dt_q, dt_remap
    integer,              intent(out) :: nstep_end
    real(kind=real_kind), intent(in)  :: dt
    type (TimeLevel_t),   intent(in)  :: tl

    dt_q      = dt*qsplit
    dt_remap  = dt_q
    nstep_end = tl%nstep + qsplit

    if (rsplit>0) then
       dt_remap  = dt_q*rsplit   ! rsplit=0 means use eulerian code, not vert. lagrange
       nstep_end = tl%nstep + qsplit*rsplit  ! nstep at end of this routine
    endif

  end subroutine

  !_____________________________________________________________________
  subroutine set_diagnostic_state(compute_diagnostics,compute_energy,tl,nstep_end)

    ! enable or disable scalar and energy diagnostics for this timestep

    use control_mod, only: energy_fixer, disable_diagnostics

    logical,              intent(out) :: compute_diagnostics, compute_energy
    type (TimeLevel_t),   intent(in)  :: tl
    integer,              intent(in)  :: nstep_end

    ! activate energy diagnostics if using an energy fixer
    compute_energy = energy_fixer > 0

    ! activate diagnostics periodically for display to stdout
    compute_diagnostics   = .false.
    if (MODULO(nstep_end,statefreq)==0 .or. nstep_end==tl%nstep0) then
       compute_diagnostics= .true.
       compute_energy     = .true.
    endif
    if(disable_diagnostics) compute_diagnostics= .false.

  end subroutine

end module
