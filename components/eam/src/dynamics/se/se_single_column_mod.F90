module se_single_column_mod
!--------------------------------------------------------
! 
! Module for the SE single column model

use element_mod, only: element_t
use scamMod
use constituents, only: cnst_get_ind
use dimensions_mod, only: nelemd, np
use time_manager, only: get_nstep, dtime, is_first_step
use ppgrid, only: begchunk
#ifdef MODEL_THETA_L
use element_ops, only: get_R_star
use eos, only: pnh_and_exner_from_eos
#endif
use element_ops, only: get_temperature

implicit none

public scm_setinitial
public scm_setfield
public apply_SC_forcing

!=========================================================================
contains
!=========================================================================

subroutine scm_setinitial(elem)

  implicit none

  type(element_t), intent(inout) :: elem(:)

  integer i, j, k, ie, thelev
  integer inumliq, inumice, icldliq, icldice

  if (.not. use_replay .and. get_nstep() .eq. 0) then
    call cnst_get_ind('NUMLIQ', inumliq, abort=.false.)
    call cnst_get_ind('NUMICE', inumice, abort=.false.)
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
              if (have_omega) elem(ie)%derived%omega_p(i,j,k) = wfld(k)
            enddo

          endif

        enddo
      enddo
    enddo
  endif

end subroutine scm_setinitial

subroutine scm_setfield(elem,iop_update_phase1)

  implicit none

  logical, intent(in) :: iop_update_phase1
  type(element_t), intent(inout) :: elem(:)

#ifndef MODEL_THETA_L

  integer i, j, k, ie

  do ie=1,nelemd
    if (have_ps .and. use_replay .and. .not. iop_update_phase1) elem(ie)%state%ps_v(:,:,:) = psobs
    if (have_ps .and. .not. use_replay) elem(ie)%state%ps_v(:,:,:) = psobs
    do i=1, PLEV
      if (have_omega .and. iop_update_phase1) elem(ie)%derived%omega_p(:,:,i)=wfld(i)  !     set t to tobs at first
    end do
  end do

#endif

end subroutine scm_setfield

subroutine apply_SC_forcing(elem,hvcoord,tl,n,t_before_advance,nets,nete)
! 
    use scamMod, only: single_column, use_3dfrc
    use kinds, only : real_kind
    use dimensions_mod, only : np, np, nlev, npsq
    use control_mod, only : use_cpstar, qsplit
    use hybvcoord_mod, only : hvcoord_t
    use element_mod, only : element_t
    use physical_constants, only : Cp, Rgas, cpwater_vapor
    use time_mod
    use constituents, only: pcnst
    use time_manager, only: get_nstep
    use shr_const_mod, only: SHR_CONST_PI

    integer :: t1,t2,n,nets,nete,pp
    type (element_t)     , intent(inout), target :: elem(:)
    type (hvcoord_t)                  :: hvcoord
    type (TimeLevel_t), intent(in)       :: tl
    logical :: t_before_advance, do_column_scm
    real(kind=real_kind), parameter :: rad2deg = 180.0_real_kind / SHR_CONST_PI

    integer :: ie,k,i,j,t,nm_f
    real (kind=real_kind), dimension(np,np,nlev)  :: dpt1,dpt2   ! delta pressure
    real (kind=real_kind), dimension(np,np)  :: E
    real (kind=real_kind), dimension(np,np)  :: suml,suml2,v1,v2
    real (kind=real_kind), dimension(np,np,nlev)  :: sumlk, suml2k
    real (kind=real_kind), dimension(np,np,nlev)  :: p,T_v,phi, pnh
    real (kind=real_kind), dimension(np,np,nlev+1) :: dpnh_dp_i
    real (kind=real_kind), dimension(np,np,nlev)  :: dp, exner, vtheta_dp, Rstar
    real (kind=real_kind), dimension(np,np,nlev) :: dpscm
    real (kind=real_kind) :: cp_star1,cp_star2,qval_t1,qval_t2
    real (kind=real_kind) :: Qt,dt
    real (kind=real_kind), dimension(nlev,pcnst) :: stateQin1, stateQin2, stateQin_qfcst
    real (kind=real_kind), dimension(nlev,pcnst) :: forecast_q
    real (kind=real_kind), dimension(nlev) :: dummy1, dummy2, forecast_t, forecast_u, forecast_v
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

    ! For SCM only one column is considered
    ie = 1
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif

    do k=1,nlev
      p(:,:,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(:,:,t1)
    end do

#ifdef MODEL_THETA_L
      ! If using the theta-l dycore then need to get the exner function
      !   and floating levels "dp", so we can convert the SCM forecasted
      !   temperature back potential temperature on floating levels.
      dp(:,:,:) = elem(ie)%state%dp3d(:,:,:,t1)
      call pnh_and_exner_from_eos(hvcoord,elem(ie)%state%vtheta_dp(:,:,:,t1),&
          dp,elem(ie)%state%phinh_i(:,:,:,t1),pnh,exner,dpnh_dp_i) 
#endif

    dt=dtime

    ! Get temperature from dynamics state
    call get_temperature(elem(ie),temperature,hvcoord,t1)

    ! For SCM we only consider one point
    i=1
    j=1

    stateQin_qfcst(:,:) = elem(ie)%state%Q(i,j,:,:)
    stateQin1(:,:) = stateQin_qfcst(:,:)
    stateQin2(:,:) = stateQin_qfcst(:,:)        

    if (.not. use_3dfrc) then
      dummy1(:) = 0.0_real_kind
    else
      dummy1(:) = elem(ie)%derived%fT(i,j,:)
    endif
    dummy2(:) = 0.0_real_kind
    forecast_ps = elem(ie)%state%ps_v(i,j,t1)

#ifdef MODEL_THETA_L
    ! At first time step the forcing term is set to the
    !  initial condition temperature profile.  Do NOT use
    !  as it will cause temperature profile to blow up.
    if ( is_first_step() ) then
      dummy1(:) = 0.0
    endif
#endif

    call forecast(begchunk,elem(ie)%state%ps_v(i,j,t1),&
           elem(ie)%state%ps_v(i,j,t1),forecast_ps,forecast_u,&
           elem(ie)%state%v(i,j,1,:,t1),elem(ie)%state%v(i,j,1,:,t1),&
           forecast_v,elem(ie)%state%v(i,j,2,:,t1),&
           elem(ie)%state%v(i,j,2,:,t1),forecast_t,&
#ifdef MODEL_THETA_L  
           temperature(i,j,:),temperature(i,j,:),&
#else
           elem(ie)%state%T(i,j,:,t1),elem(ie)%state%T(i,j,:,t1),&
#endif
           forecast_q,stateQin2,stateQin1,dt,dummy1,dummy2,dummy2,&
           stateQin_qfcst,p(i,j,:),stateQin1,1)         

    elem(ie)%state%Q(i,j,:,:) = forecast_q(:,:)

#ifdef MODEL_THETA_L
    ! If running theta-l model then the forecaste temperature needs
    !   to be converted back to potential temperature on floating levels.
    call get_R_star(Rstar,elem(ie)%state%Q(:,:,:,1))
    elem(ie)%state%vtheta_dp(i,j,:,t1) = (forecast_t(:)*Rstar(i,j,:)*dp(i,j,:))/&
                (Rgas*exner(i,j,:))
#else
    elem(ie)%state%T(i,j,:,t1) = forecast_t(:)           
#endif
    elem(ie)%state%v(i,j,1,:,t1) = forecast_u(:)
    elem(ie)%state%v(i,j,2,:,t1) = forecast_v(:)

    end subroutine apply_SC_forcing

end module se_single_column_mod
