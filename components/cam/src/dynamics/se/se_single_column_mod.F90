module se_single_column_mod
!--------------------------------------------------------
! 
! Module for the SE single column model

use element_mod, only: element_t
use scamMod
use constituents, only: cnst_get_ind
use dimensions_mod, only: nelemd, np
use time_manager, only: get_nstep, dtime

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

  if (.not. use_camiop .and. get_nstep() .eq. 0) then
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
              tobs(k)=elem(ie)%state%T(i,j,k,1)
              qobs(k)=elem(ie)%state%Q(i,j,k,1)
            enddo
          else
            tobs(:)=elem(ie)%state%T(i,j,:,1)
            qobs(:)=elem(ie)%state%Q(i,j,:,1)
          endif

          if (get_nstep() .eq. 0) then
            do k=thelev, PLEV
              if (have_t) elem(ie)%state%T(i,j,k,1)=tobs(k)
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
              if (have_omega) elem(ie)%derived%omega_p(i,j,k) = wfldh(k)
            enddo

          endif

        enddo
      enddo
    enddo
  endif

end subroutine scm_setinitial

subroutine scm_setfield(elem)

  implicit none

  type(element_t), intent(inout) :: elem(:)

  integer i, j, k, ie

  do ie=1,nelemd
    if (have_ps) elem(ie)%state%ps_v(:,:,:) = psobs 
    do i=1, PLEV
      if (have_omega) elem(ie)%derived%omega_p(:,:,i)=wfldh(i)  !     set t to tobs at first
    end do
  end do

end subroutine scm_setfield

subroutine apply_SC_forcing(elem,hvcoord,tl,n,t_before_advance,nets,nete)
! 
    use scamMod, only: single_column, use_3dfrc
    use kinds, only : real_kind
    use dimensions_mod, only : np, np, nlev, npsq
    use control_mod, only : use_cpstar, qsplit
    use hybvcoord_mod, only : hvcoord_t
    use element_mod, only : element_t
    use physical_constants, only : Cp, cpwater_vapor
    use time_mod
    use constituents, only: pcnst
    use time_manager, only: get_nstep
    use shr_const_mod, only: SHR_CONST_PI

    integer :: t1,t2,n,nets,nete,pp
    type (element_t)     , intent(inout), target :: elem(:)
    type (hvcoord_t)                  :: hvcoord
    type (TimeLevel_t), intent(in)       :: tl
    logical :: t_before_advance, do_column_scm
    real(kind=real_kind), parameter :: rad2deg = 180.0 / SHR_CONST_PI

    integer :: ie,k,i,j,t,nm_f
    real (kind=real_kind), dimension(np,np,nlev)  :: dpt1,dpt2   ! delta pressure
    real (kind=real_kind), dimension(np,np)  :: E
    real (kind=real_kind), dimension(np,np)  :: suml,suml2,v1,v2
    real (kind=real_kind), dimension(np,np,nlev)  :: sumlk, suml2k
    real (kind=real_kind), dimension(np,np,nlev)  :: p,T_v,phi
    real (kind=real_kind) :: cp_star1,cp_star2,qval_t1,qval_t2
    real (kind=real_kind) :: Qt,dt
    real (kind=real_kind), dimension(nlev,pcnst) :: stateQin1, stateQin2, stateQin_qfcst
    real (kind=real_kind), dimension(nlev,pcnst) :: forecast_q
    real (kind=real_kind), dimension(nlev) :: dummy1, dummy2, forecast_t, forecast_u, forecast_v
    real (kind=real_kind) :: forecast_ps
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

    !   IE   Cp*dpdn*T  + (Cpv-Cp) Qdpdn*T
    !        Cp*dpdn(n)*T(n+1) + (Cpv-Cp) Qdpdn(n)*T(n+1)
    !        [Cp + (Cpv-Cp) Q(n)] *dpdn(n)*T(n+1) 

    ie=1

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif

    do k=1,nlev
      p(:,:,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(:,:,t1)
    end do

    dt=dtime

    i=1
    j=1

    stateQin_qfcst(:,:) = elem(ie)%state%Q(i,j,:,:)
    stateQin1(:,:) = stateQin_qfcst(:,:)
    stateQin2(:,:) = stateQin_qfcst(:,:)        

    if (.not. use_3dfrc) then
      dummy1(:) = 0.0
    else
      dummy1(:) = elem(ie)%derived%fT(i,j,:)
    endif
    dummy2(:) = 0.0
    forecast_ps = elem(ie)%state%ps_v(i,j,t1)

    call forecast(97,elem(ie)%state%ps_v(i,j,t1),&
           elem(ie)%state%ps_v(i,j,t1),forecast_ps,forecast_u,&
           elem(ie)%state%v(i,j,1,:,t1),elem(ie)%state%v(i,j,1,:,t1),&
           forecast_v,elem(ie)%state%v(i,j,2,:,t1),&
           elem(ie)%state%v(i,j,2,:,t1),forecast_t,&
           elem(ie)%state%T(i,j,:,t1),elem(ie)%state%T(i,j,:,t1),&
           forecast_q,stateQin2,stateQin1,dt,dummy1,dummy2,dummy2,&
           stateQin_qfcst,p(i,j,:),stateQin1,1)         

    elem(ie)%state%T(i,j,:,t1) = forecast_t(:)
    elem(ie)%state%v(i,j,1,:,t1) = forecast_u(:)
    elem(ie)%state%v(i,j,2,:,t1) = forecast_v(:)
    elem(ie)%state%Q(i,j,:,:) = forecast_q(:,:)

    end subroutine apply_SC_forcing

end module se_single_column_mod
