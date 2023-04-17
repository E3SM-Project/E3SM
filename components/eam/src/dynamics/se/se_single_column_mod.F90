module se_single_column_mod
!--------------------------------------------------------
! 
! Module for the SE single column model

use element_mod, only: element_t
use scamMod
use constituents, only: cnst_get_ind, pcnst
use dimensions_mod, only: nelemd, np
use time_manager, only: get_nstep, dtime, is_first_step, &
  is_first_restart_step, is_last_step
use ppgrid, only: begchunk
use pmgrid, only: plev, plon
use parallel_mod,            only: par
#ifdef MODEL_THETA_L
use element_ops, only: get_R_star, get_pottemp
use eos, only: pnh_and_exner_from_eos
#endif
use element_ops, only: get_temperature
use spmd_utils, only: masterproc

implicit none

public scm_setinitial
public scm_setfield
public scm_broadcast
public apply_SC_forcing

!=========================================================================
contains
!=========================================================================

subroutine scm_setinitial(elem)

  use kinds, only : real_kind
  use constituents, only: qmin

  implicit none

  type(element_t), intent(inout) :: elem(:)

  integer i, j, k, cix, ie, thelev
  integer inumliq, inumice, icldliq, icldice

  if (.not. use_replay .and. get_nstep() .eq. 0 .and. par%dynproc) then
    call cnst_get_ind('NUMLIQ', inumliq, abrtf=.false.)
    call cnst_get_ind('NUMICE', inumice, abrtf=.false.)
    call cnst_get_ind('CLDLIQ', icldliq)
    call cnst_get_ind('CLDICE', icldice)
    
    do ie=1,nelemd
      do j=1,np
        do i=1,np

          ! Find level where tobs is no longer zero
          thelev=1
          do k=1, PLEV
            if (tobs(k) .ne. 0) then
              thelev=k
              go to 1000
            endif
          enddo

1000 continue

          if (get_nstep() .le. 1) then
            do k=1,thelev-1
#ifdef MODEL_THETA_L
              tobs(k)=elem(ie)%derived%FT(i,j,k)
#else
              tobs(k)=elem(ie)%state%T(i,j,k,1)
#endif
              qobs(k)=elem(ie)%state%Q(i,j,k,1)
            enddo
          else
#ifdef MODEL_THETA_L
            tobs(:)=elem(ie)%derived%FT(i,j,:)
#else
            tobs(:)=elem(ie)%state%T(i,j,:,1)
#endif
            qobs(:)=elem(ie)%state%Q(i,j,:,1)
          endif

          if (get_nstep() .eq. 0) then
            do cix = 1, pcnst
               if (scm_zero_non_iop_tracers) elem(ie)%state%Q(i,j,:,cix) = qmin(cix)
            end do
            do k=thelev, PLEV
#ifdef MODEL_THETA_L
              if (have_t) elem(ie)%derived%FT(i,j,k)=tobs(k)
#else
              if (have_t) elem(ie)%state%T(i,j,k,1)=tobs(k)
#endif
              if (have_q) elem(ie)%state%Q(i,j,k,1)=qobs(k)
            enddo

            do k=1,PLEV
              if (have_ps) elem(ie)%state%ps_v(i,j,1) = psobs
              if (have_u) elem(ie)%state%v(i,j,1,k,1) = uobs(k)
              if (have_v) elem(ie)%state%v(i,j,2,k,1) = vobs(k)
              if (have_numliq) elem(ie)%state%Q(i,j,k,inumliq) = numliqobs(k)
              if (have_cldliq) elem(ie)%state%Q(i,j,k,icldliq) = cldliqobs(k)
              if (have_numice) elem(ie)%state%Q(i,j,k,inumice) = numiceobs(k)
              if (have_cldice) elem(ie)%state%Q(i,j,k,icldice) = cldiceobs(k)
              !  If DP-CRM mode we do NOT want to write over the dy-core vertical
              !    velocity with the large-scale one.  wfld is used in forecast.F90
              !    for the compuation of the large-scale subsidence.
              if (have_omega .and. .not. dp_crm) elem(ie)%derived%omega_p(i,j,k) = wfld(k)
              if (dp_crm) elem(ie)%derived%omega_p(i,j,k) = 0.0_real_kind
            enddo

          endif

        enddo
      enddo
    enddo
  endif

  ! If DP-CRM mode then SHOC/CLUBB needs to know about grid
  !   length size.  The calculations of this based on a sphere in the
  !   SHOC and CLUBB interefaces are not valid for a planar grid, thus
  !   save the grid length from the dycore. Note that planar dycore
  !   only supports uniform grids, thus we only save one value.
  ! Set this if it is the first time step or the first restart step
  if ((get_nstep() .eq. 0 .or. is_first_restart_step()) .and. dp_crm .and. par%dynproc) then
    do ie=1,nelemd
        dyn_dx_size = elem(ie)%dx_short * 1000.0_real_kind
    enddo
  endif

end subroutine scm_setinitial

subroutine scm_broadcast()

  use mpishorthand
  
#ifdef SPMD  

  call mpibcast(have_ps,1,mpilog,0,mpicom)
  call mpibcast(have_tg,1,mpilog,0,mpicom)
  call mpibcast(have_lhflx,1,mpilog,0,mpicom)
  call mpibcast(have_shflx,1,mpilog,0,mpicom)
  call mpibcast(have_t,1,mpilog,0,mpicom)
  call mpibcast(have_q,1,mpilog,0,mpicom)
  call mpibcast(have_u,1,mpilog,0,mpicom)
  call mpibcast(have_v,1,mpilog,0,mpicom)
  call mpibcast(have_omega,1,mpilog,0,mpicom)
  call mpibcast(have_cldliq,1,mpilog,0,mpicom)
  call mpibcast(have_divt,1,mpilog,0,mpicom)
  call mpibcast(have_divq,1,mpilog,0,mpicom)
  call mpibcast(have_divt3d,1,mpilog,0,mpicom)
  call mpibcast(have_divq3d,1,mpilog,0,mpicom)
  call mpibcast(use_3dfrc,1,mpilog,0,mpicom)

  call mpibcast(psobs,1,mpir8,0,mpicom)
  call mpibcast(tground,1,mpir8,0,mpicom)
  call mpibcast(lhflxobs,1,mpir8,0,mpicom)
  call mpibcast(shflxobs,1,mpir8,0,mpicom)
!  
  call mpibcast(tobs,plev,mpir8,0,mpicom)
  call mpibcast(qobs,plev,mpir8,0,mpicom)
  call mpibcast(uobs,plev,mpir8,0,mpicom)
  call mpibcast(vobs,plev,mpir8,0,mpicom)
  call mpibcast(cldliqobs,plev,mpir8,0,mpicom)
  call mpibcast(wfld,plev,mpir8,0,mpicom) 
  
  call mpibcast(divt,plev,mpir8,0,mpicom)
  call mpibcast(divq,plev,mpir8,0,mpicom)
  call mpibcast(divt3d,plev,mpir8,0,mpicom)
  call mpibcast(divq3d,plev,mpir8,0,mpicom)
  
#endif

end subroutine scm_broadcast

subroutine scm_setfield(elem,iop_update_phase1)

  implicit none

  logical, intent(in) :: iop_update_phase1
  type(element_t), intent(inout) :: elem(:)

  integer i, j, k, ie

  do ie=1,nelemd
    if (have_ps .and. use_replay .and. .not. iop_update_phase1) elem(ie)%state%ps_v(:,:,:) = psobs
    if (have_ps .and. .not. use_replay) elem(ie)%state%ps_v(:,:,:) = psobs
    do i=1, PLEV
      ! If DP CRM mode do NOT write over dycore vertical velocity
      if ((have_omega .and. iop_update_phase1) .and. .not. dp_crm) elem(ie)%derived%omega_p(:,:,i)=wfld(i)  !     set t to tobs at first
    end do
  end do

end subroutine scm_setfield

subroutine apply_SC_forcing(elem,hvcoord,hybrid,tl,n,t_before_advance,nets,nete)
! 
    use scamMod
    use kinds, only : real_kind
    use dimensions_mod, only : np, np, nlev, npsq, nelem
    use control_mod, only : use_cpstar, qsplit, qsplit
    use hybvcoord_mod, only : hvcoord_t
    use element_mod, only : element_t
    use physical_constants, only : Cp, Rgas, cpwater_vapor
    use time_mod
    use hybrid_mod, only : hybrid_t
    use constituents, only: pcnst
    use time_manager, only: get_nstep
    use cam_history, only: outfld
    use shr_const_mod, only: SHR_CONST_PI

    integer :: t1,t2,n,nets,nete,pp
    type (element_t)     , intent(inout), target :: elem(:)
    type (hvcoord_t)                  :: hvcoord
    type (TimeLevel_t), intent(in)       :: tl
    type(hybrid_t),             intent(in) :: hybrid    
    logical :: t_before_advance, do_column_scm
    real(kind=real_kind), parameter :: rad2deg = 180.0_real_kind / SHR_CONST_PI

    integer :: ie,k,i,j,t,m,nm_f,tlQdp
    integer :: nelemd_todo, np_todo
    real (kind=real_kind), dimension(np,np,nlev)  :: dpt1,dpt2   ! delta pressure
    real (kind=real_kind), dimension(np,np)  :: E
    real (kind=real_kind), dimension(np,np)  :: suml,suml2,v1,v2
    real (kind=real_kind), dimension(np,np,nlev)  :: sumlk, suml2k
    real (kind=real_kind), dimension(np,np,nlev)  :: p,T_v,phi, pnh
    real (kind=real_kind), dimension(np,np,nlev)  :: dp3d
    real (kind=real_kind), dimension(np,np,nlev+1) :: dpnh_dp_i
    real (kind=real_kind), dimension(np,np,nlev)  :: dp, exner, vtheta_dp, Rstar
    real (kind=real_kind), dimension(np,np,nlev) :: dpscm
    real (kind=real_kind) :: cp_star1,cp_star2,qval_t1,qval_t2
    real (kind=real_kind) :: Qt,dt
    real (kind=real_kind), dimension(nlev,pcnst) :: stateQin1, stateQin2, stateQin_qfcst
    real (kind=real_kind), dimension(nlev,pcnst) :: forecast_q
    real (kind=real_kind), dimension(nlev) :: dummy1, dummy2, forecast_t, forecast_u, forecast_v
    real (kind=real_kind), dimension(nlev) :: tdiff_dyn, qdiff_dyn, temp_tend
    real (kind=real_kind), dimension(npsq,nlev) :: tdiff_out, qdiff_out
    real (kind=real_kind) :: forecast_ps
    real (kind=real_kind) :: temperature(np,np,nlev)
    logical :: wet

    integer:: icount

    nm_f = 1
    if (t_before_advance) then
       t1=tl%nm1
       t2=tl%n0
    else
       t1=tl%n0
       t2=tl%np1
    endif
    
    call TimeLevel_Qdp(tl,qsplit,tlQdp)

    ! Settings for traditional SCM run
    nelemd_todo = 1
    np_todo = 1
    
    if (scm_multcols) then
      nelemd_todo = nelemd
      np_todo = np
    endif

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif

    ! collect stats from dycore for analysis
#ifdef MODEL_THETA_L
    if (dp_crm) then
      call crm_resolved_turb(elem,hvcoord,hybrid,t1,nelemd_todo,np_todo)
    endif
#endif

    do ie=1,nelemd_todo

      do k=1,nlev
        p(:,:,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(:,:,t1)
        dpscm(:,:,k) = elem(ie)%state%dp3d(:,:,k,t1)
      end do

      dt=dtime

#ifdef MODEL_THETA_L
      ! If using the theta-l dycore then need to get the exner function
      !   and reference levels "dp", so we can convert the SCM forecasted
      !   temperature back potential temperature on reference levels.
      dp(:,:,:) = elem(ie)%state%dp3d(:,:,:,t1)
      call pnh_and_exner_from_eos(hvcoord,elem(ie)%state%vtheta_dp(:,:,:,t1),&
          dp,elem(ie)%state%phinh_i(:,:,:,t1),pnh,exner,dpnh_dp_i)
#endif

    ! Get temperature from dynamics state
    call get_temperature(elem(ie),temperature,hvcoord,t1)

      do j=1,np_todo
        do i=1,np_todo

          stateQin_qfcst(:nlev,:pcnst) = elem(ie)%state%Q(i,j,:nlev,:pcnst)
          stateQin1(:nlev,:pcnst) = stateQin_qfcst(:nlev,:pcnst)
          stateQin2(:nlev,:pcnst) = stateQin_qfcst(:nlev,:pcnst)
          forecast_q(:nlev,:pcnst) = stateQin_qfcst(:nlev,:pcnst)

          if (.not. use_3dfrc .or. dp_crm) then
            temp_tend(:) = 0.0_real_kind
          else
            temp_tend(:) = elem(ie)%derived%fT(i,j,:)
          endif
          dummy2(:) = 0.0_real_kind
          forecast_ps = elem(ie)%state%ps_v(i,j,t1)
          
#ifdef MODEL_THETA_L
          ! At first time step the tendency term is set to the
          !  initial condition temperature profile.  Do NOT use
          !  as it will cause temperature profile to blow up.
          if ( is_first_step() ) then
            temp_tend(:) = 0.0_real_kind
          endif
#endif          

          call forecast(begchunk,elem(ie)%state%ps_v(i,j,t1),&
            elem(ie)%state%ps_v(i,j,t1),forecast_ps,forecast_u,&
            elem(ie)%state%v(i,j,1,:,t1),elem(ie)%state%v(i,j,1,:,t1),&
            forecast_v,elem(ie)%state%v(i,j,2,:,t1),&
            elem(ie)%state%v(i,j,2,:,t1),forecast_t,&
            temperature(i,j,:),temperature(i,j,:),&
            forecast_q,stateQin2,stateQin1,dt,temp_tend,dummy2,dummy2,&
            stateQin_qfcst,p(i,j,:),stateQin1,1,&
            tdiff_dyn,qdiff_dyn)         

            ! Update the q related arrays.  NOTE that Qdp array must
            !  be updated first to ensure exact restarts
            do m=1,pcnst
              ! Update the Qdp array
              elem(ie)%state%Qdp(i,j,:nlev,m,tlQdp) = &
                forecast_q(:nlev,m) * dpscm(i,j,:nlev)
              ! Update the Q array
              elem(ie)%state%Q(i,j,:nlev,m) = &
                elem(ie)%state%Qdp(i,j,:nlev,m,tlQdp)/dpscm(i,j,:nlev)
            enddo

#ifdef MODEL_THETA_L
          ! If running theta-l model then the forecast temperature needs
          !   to be converted back to potential temperature on reference levels, 
          !   which is what dp_coupling expects
          call get_R_star(Rstar,elem(ie)%state%Q(:,:,:,1))
          elem(ie)%state%vtheta_dp(i,j,:,t1) = (forecast_t(:)*Rstar(i,j,:)*dp(i,j,:))/&
                (Rgas*exner(i,j,:))
#else
          elem(ie)%state%T(i,j,:,t1) = forecast_t(:)
#endif
          elem(ie)%state%v(i,j,1,:,t1) = forecast_u(:)
          elem(ie)%state%v(i,j,2,:,t1) = forecast_v(:)

          tdiff_out(i+(j-1)*np,:)=tdiff_dyn(:)
          qdiff_out(i+(j-1)*np,:)=qdiff_dyn(:)
  
        enddo
      enddo
      
      if (scm_multcols) then
        call outfld('TDIFF',tdiff_out,npsq,ie)
        call outfld('QDIFF',qdiff_out,npsq,ie)
      else
        call outfld('TDIFF',tdiff_dyn,plon,begchunk)
        call outfld('QDIFF',qdiff_dyn,plon,begchunk)
      endif
    
    enddo

    if ((iop_nudge_tq .or. iop_nudge_uv) .and. dp_crm) then
      ! If running in a doubly periodic CRM mode, then nudge the domain
      !   based on the domain mean and observed quantities of T, Q, u, and v
      call iop_domain_relaxation(elem,hvcoord,hybrid,t1,dp,nelemd_todo,np_todo,dt)
    endif

end subroutine apply_SC_forcing


subroutine iop_domain_relaxation(elem,hvcoord,hybrid,t1,dp,nelemd_todo,np_todo,dt)

  ! Subroutine intended to nudge a uniformly forced grid with heterogenous
  !   surface (i.e. doubly-periodic SCREAM).

  use scamMod
  use dimensions_mod, only : np, np, nlev, npsq, nelem
  use parallel_mod, only: global_shared_buf, global_shared_sum
  use global_norms_mod, only: wrap_repro_sum
  use hybvcoord_mod, only : hvcoord_t
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t
  use time_mod
  use physical_constants, only : Cp, Rgas

  ! Input/Output variables
  type (element_t)     , intent(inout), target :: elem(:)
  type (hvcoord_t)                  :: hvcoord
  type(hybrid_t),             intent(in) :: hybrid
  integer, intent(in) :: nelemd_todo, np_todo, t1
  real (kind=real_kind), intent(in):: dt

  ! Local variables
  real (kind=real_kind), dimension(np,np,nlev+1) :: dpnh_dp_i
  real (kind=real_kind), dimension(nlev) :: domain_q, domain_t, domain_u, domain_v, rtau
  real (kind=real_kind), dimension(nlev) :: relax_t, relax_q, relax_u, relax_v, iop_pres
  real (kind=real_kind), dimension(np,np,nlev) :: temperature, Rstar, pnh, exner, dp
  integer :: ie, i, j, k

  ! Compute pressure for IOP observations
  do k=1,nlev
    iop_pres(k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*psobs
  end do

  ! Compute the domain means for temperature, moisture, u, and v
  do k=1,nlev
    do ie=1,nelemd_todo
      ! Get absolute temperature
      call get_temperature(elem(ie),temperature,hvcoord,t1)

      ! Initialize the global buffer
      global_shared_buf(ie,1) = 0.0_real_kind
      global_shared_buf(ie,2) = 0.0_real_kind
      global_shared_buf(ie,3) = 0.0_real_kind
      global_shared_buf(ie,4) = 0.0_real_kind

      ! Sum each variable for each level
      global_shared_buf(ie,1) = global_shared_buf(ie,1) + SUM(elem(ie)%state%Q(:,:,k,1))
      global_shared_buf(ie,2) = global_shared_buf(ie,2) + SUM(temperature(:,:,k))
      global_shared_buf(ie,3) = global_shared_buf(ie,3) + SUM(elem(ie)%state%v(:,:,1,k,t1))
      global_shared_buf(ie,4) = global_shared_buf(ie,4) + SUM(elem(ie)%state%v(:,:,2,k,t1))
    enddo

    ! Compute the domain mean
    call wrap_repro_sum(nvars=4, comm=hybrid%par%comm)
    domain_q(k) = global_shared_sum(1)/(dble(nelem)*dble(np)*dble(np))
    domain_t(k) = global_shared_sum(2)/(dble(nelem)*dble(np)*dble(np))
    domain_u(k) = global_shared_sum(3)/(dble(nelem)*dble(np)*dble(np))
    domain_v(k) = global_shared_sum(4)/(dble(nelem)*dble(np)*dble(np))
  enddo

  ! Initialize relaxation arrays
  do k=1,nlev
    relax_q(k) = 0._real_kind
    relax_t(k) = 0._real_kind
    relax_u(k) = 0._real_kind
    relax_v(k) = 0._real_kind
  enddo

  ! Compute the relaxation for each level.  This is the relaxation based
  !   on the domain mean and observed quantity from IOP file.
  do k=1,nlev

    ! Define nudging timescale
    rtau(k) = iop_nudge_tscale
    rtau(k) = max(dt,rtau(k))

    ! Compute relaxation for winds
    relax_u(k) = -(domain_u(k) - uobs(k))/rtau(k)
    relax_v(k) = -(domain_v(k) - vobs(k))/rtau(k)

    ! Restrict nudging of T and Q to certain levels if requested by user
    ! pmidm1 variable is in unitis of [Pa], while iop_nudge_tq_low/high
    !   is in units of [hPa], thus convert iop_nudge_tq_low/high
    if (iop_pres(k) .le. iop_nudge_tq_low*100._real_kind .and. &
      iop_pres(k) .ge. iop_nudge_tq_high*100._real_kind) then

      ! compute relaxation for each variable
      !  units are [unit/s] (i.e. K/s for temperature)
      relax_q(k) = -(domain_q(k) - qobs(k))/rtau(k)
      relax_t(k) = -(domain_t(k) - tobs(k))/rtau(k)
    endif

  enddo

  ! Now apply the nudging tendency for each variable
  do ie=1,nelemd_todo
#ifdef MODEL_THETA_L
    dp(:,:,:) = elem(ie)%state%dp3d(:,:,:,t1)
    ! Get updated Rstar to convert to density weighted potential temperature
    call get_R_star(Rstar,elem(ie)%state%Q(:,:,:,1))
    ! Get updated exner
    call pnh_and_exner_from_eos(hvcoord,elem(ie)%state%vtheta_dp(:,:,:,t1),&
           dp,elem(ie)%state%phinh_i(:,:,:,t1),pnh,exner,dpnh_dp_i)
#endif
    call get_temperature(elem(ie),temperature,hvcoord,t1)
    do j=1,np_todo
      do i=1,np_todo
        do k=1,nlev

          if (iop_nudge_tq) then
            ! Apply vapor relaxation
            elem(ie)%state%Q(i,j,k,1) = elem(ie)%state%Q(i,j,k,1) + relax_q(k) * dt

            ! Apply temperature relaxation
            temperature(i,j,k) = temperature(i,j,k) + relax_t(k) * dt

            ! Update temperature appropriately based on dycore used
#ifdef MODEL_THETA_L
            elem(ie)%state%vtheta_dp(i,j,k,t1) = (temperature(i,j,k)*Rstar(i,j,k)*dp(i,j,k))/&
              (Rgas*exner(i,j,k))
#else
            elem(ie)%state%T(i,j,k,t1) = temperature(i,j,k)
#endif
          endif

          if (iop_nudge_uv) then
            ! Apply u and v relaxation
            elem(ie)%state%v(i,j,1,k,t1) = elem(ie)%state%v(i,j,1,k,t1) + relax_u(k) * dt
            elem(ie)%state%v(i,j,2,k,t1) = elem(ie)%state%v(i,j,2,k,t1) + relax_v(k) * dt
          endif

        enddo
      enddo
    enddo
  enddo

end subroutine iop_domain_relaxation

#ifdef MODEL_THETA_L
subroutine crm_resolved_turb(elem,hvcoord,hybrid,t1,&
                              nelemd_todo,np_todo)

  ! Subroutine intended to compute various resolved turbulence statistics
  ! (done so to be consistent with SHOC's definition of each for ease of
  ! comparison) for DP CRM runs.

  use scamMod
  use dimensions_mod, only : np, np, nlev, nlevp, npsq, nelem
  use parallel_mod, only: global_shared_buf, global_shared_sum
  use global_norms_mod, only: wrap_repro_sum
  use hybvcoord_mod, only : hvcoord_t
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t
  use time_mod
  use physical_constants, only : Cp, Rgas
  use cam_history, only: outfld
  use physconst, only: latvap

  ! Input/Output variables
  type (element_t)     , intent(inout), target :: elem(:)
  type (hvcoord_t)                  :: hvcoord
  type(hybrid_t),             intent(in) :: hybrid
  integer, intent(in) :: nelemd_todo, np_todo, t1

  ! Local variables
  real (kind=real_kind), dimension(nlev) :: domain_thetal, domain_qw, domain_u, domain_v
  real (kind=real_kind) :: thetal(np,np,nlev), qw(np,np,nlev)
  real (kind=real_kind) :: rho(np,np,nlev), temperature(np,np,nlev)
  real (kind=real_kind) :: wthl_res(np,np,nlev+1), wqw_res(np,np,nlev+1)
  real (kind=real_kind) :: w2_res(np,np,nlev+1), w3_res(np,np,nlev+1)
  real (kind=real_kind) :: u2_res(np,np,nlev), v2_res(np,np,nlev)
  real (kind=real_kind) :: thl2_res(np,np,nlev), qw2_res(np,np,nlev), qwthl_res(np,np,nlev)
  real (kind=real_kind) :: pres(nlev)
  integer :: ie, i, j, k

  ! Compute pressure for IOP observations
  do k=1,nlev
    pres(k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*psobs
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Find the domain mean of thetal, qw, u, and v
  do k=1,nlev
    do ie=1,nelemd_todo

      ! Compute potential temperature
      call get_pottemp(elem(ie),thetal,hvcoord,t1,1)

      ! Compute liquid water potential temperature
      thetal(:,:,k) = thetal(:,:,k) - (latvap/Cp)*(elem(ie)%state%Q(1:np,1:np,k,2))
      ! Compute total water
      qw(:,:,k) = elem(ie)%state%Q(1:np,1:np,k,1) + elem(ie)%state%Q(1:np,1:np,k,2)

      ! Initialize the global buffer
      global_shared_buf(ie,1) = 0.0_real_kind
      global_shared_buf(ie,2) = 0.0_real_kind
      global_shared_buf(ie,3) = 0.0_real_kind
      global_shared_buf(ie,4) = 0.0_real_kind

      ! Sum each variable for each level
      global_shared_buf(ie,1) = global_shared_buf(ie,1) + SUM(thetal(:,:,k))
      global_shared_buf(ie,2) = global_shared_buf(ie,2) + SUM(qw(:,:,k))
      global_shared_buf(ie,3) = global_shared_buf(ie,3) + SUM(elem(ie)%state%v(:,:,1,k,t1))
      global_shared_buf(ie,4) = global_shared_buf(ie,4) + SUM(elem(ie)%state%v(:,:,2,k,t1))
    enddo

    ! Compute the domain mean
    call wrap_repro_sum(nvars=4, comm=hybrid%par%comm)
    domain_thetal(k) = global_shared_sum(1)/(dble(nelem)*dble(np)*dble(np))
    domain_qw(k) = global_shared_sum(2)/(dble(nelem)*dble(np)*dble(np))
    domain_u(k) = global_shared_sum(3)/(dble(nelem)*dble(np)*dble(np))
    domain_v(k) = global_shared_sum(4)/(dble(nelem)*dble(np)*dble(np))
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Compute the flux
  do ie=1,nelemd_todo

    ! Compute potential temperature
    call get_pottemp(elem(ie),thetal,hvcoord,t1,1)
    call get_temperature(elem(ie),temperature,hvcoord,t1)

    ! Compute liquid water potential temperature
    do k=1,nlev
      thetal(:,:,k) = thetal(:,:,k) - (latvap/Cp)*(elem(ie)%state%Q(1:np,1:np,k,2))
      qw(:,:,k) = elem(ie)%state%Q(1:np,1:np,k,1) + elem(ie)%state%Q(1:np,1:np,k,2)
    enddo

    ! Compute air density
    do k=1,nlev
      rho(:,:,k) = pres(k)/(Rgas*temperature(:,:,k))
    enddo

    ! Compute fluxes
    wthl_res(:,:,:) = 0.0_real_kind ! initialize
    wqw_res(:,:,:) = 0.0_real_kind
    do k=2,nlev
      wthl_res(:,:,k) = elem(ie)%state%w_i(:,:,k,t1) * 0.5_real_kind * &
                        ((thetal(:,:,k)-domain_thetal(k))+&
                        (thetal(:,:,k-1)-domain_thetal(k-1)))
      wqw_res(:,:,k) = elem(ie)%state%w_i(:,:,k,t1) * 0.5_real_kind * &
                        ((qw(:,:,k)-domain_qw(k))+(qw(:,:,k-1)-domain_qw(k-1)))

      ! convert to W/m2
      wthl_res(:,:,k) = wthl_res(:,:,k) * 0.5_real_kind*(rho(:,:,k)+rho(:,:,k-1)) * Cp
      wqw_res(:,:,k) = wqw_res(:,:,k) * 0.5_real_kind*(rho(:,:,k)+rho(:,:,k-1)) * latvap
    enddo

    ! Compute variances and covariances on interface grid
    w2_res(:,:,:) = 0.0_real_kind
    w3_res(:,:,:) = 0.0_real_kind
    do k=1,nlev+1
      w2_res(:,:,k) = elem(ie)%state%w_i(:,:,k,t1)**2
      w3_res(:,:,k) = elem(ie)%state%w_i(:,:,k,t1)**3
    enddo

    ! Compute variances and covariances on midpoint grid
    u2_res(:,:,:) = 0.0_real_kind
    v2_res(:,:,:) = 0.0_real_kind
    thl2_res(:,:,:) = 0.0_real_kind
    qw2_res(:,:,:) = 0.0_real_kind
    qwthl_res(:,:,:) = 0.0_real_kind
    do k=1,nlev
      u2_res(:,:,k) = (elem(ie)%state%V(:,:,1,k,t1) - domain_u(k))**2
      v2_res(:,:,k) = (elem(ie)%state%V(:,:,2,k,t1) - domain_v(k))**2
      thl2_res(:,:,k) = (thetal(:,:,k)-domain_thetal(k))**2
      qw2_res(:,:,k) = (qw(:,:,k)-domain_qw(k))**2
      qwthl_res(:,:,k) = (thetal(:,:,k)-domain_thetal(k))* &
                         (qw(:,:,k)-domain_qw(k))
    enddo

    call outfld('WTHL_RES',wthl_res,npsq,ie)
    call outfld('WQW_RES',wqw_res,npsq,ie)
    call outfld('W2_RES',w2_res,npsq,ie)
    call outfld('W3_RES',w3_res,npsq,ie)
    call outfld('U2_RES',u2_res,npsq,ie)
    call outfld('V2_RES',v2_res,npsq,ie)
    call outfld('THL2_RES',thl2_res,npsq,ie)
    call outfld('QW2_RES',qw2_res,npsq,ie)
    call outfld('QWTHL_RES',qwthl_res,npsq,ie)

  enddo

end subroutine crm_resolved_turb
#endif  // MODEL_THETA_L

end module se_single_column_mod
