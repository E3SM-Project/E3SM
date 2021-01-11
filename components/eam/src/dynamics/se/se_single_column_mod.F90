module se_single_column_mod
!--------------------------------------------------------
! 
! Module for the SE single column model

use element_mod, only: element_t
use scamMod
use constituents, only: cnst_get_ind
use dimensions_mod, only: nelemd, np
use time_manager, only: get_nstep, dtime
use ppgrid, only: begchunk
use pmgrid, only: plev, plon
use parallel_mod,            only: par

implicit none

public scm_setinitial
public scm_setfield
public scm_broadcast
public apply_SC_forcing

!=========================================================================
contains
!=========================================================================

subroutine scm_setinitial(elem)

  implicit none

  type(element_t), intent(inout) :: elem(:)

  integer i, j, k, ie, thelev
  integer inumliq, inumice, icldliq, icldice

  if (.not. use_replay .and. get_nstep() .eq. 0 .and. par%dynproc) then
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
              !  If IOP mode we do NOT want to write over the dy-core vertical 
              !    velocity with the large-scale one.  wfld is used in forecast.F90
              !    for the compuation of the large-scale subsidence.
              if (have_omega .and. .not. iop_mode) elem(ie)%derived%omega_p(i,j,k) = wfld(k)
            enddo

          endif

        enddo
      enddo
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
    if (have_ps .and. .not. iop_update_phase1) elem(ie)%state%ps_v(:,:,:) = psobs 
    do i=1, PLEV
      ! If IOP mode do NOT write over dycore vertical velocity
      if ((have_omega .and. iop_update_phase1) .and. .not. iop_mode) elem(ie)%derived%omega_p(:,:,i)=wfld(i)  !     set t to tobs at first
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
    use cam_history, only: outfld
    use shr_const_mod, only: SHR_CONST_PI

    integer :: t1,t2,n,nets,nete,pp
    type (element_t)     , intent(inout), target :: elem(:)
    type (hvcoord_t)                  :: hvcoord
    type (TimeLevel_t), intent(in)       :: tl
    logical :: t_before_advance, do_column_scm
    real(kind=real_kind), parameter :: rad2deg = 180.0_real_kind / SHR_CONST_PI

    integer :: ie,k,i,j,t,nm_f
    integer :: nelemd_todo, np_todo
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
    real (kind=real_kind), dimension(nlev) :: tdiff_dyn, qdiff_dyn
    real (kind=real_kind), dimension(npsq,nlev) :: tdiff_out, qdiff_out
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
    
    ! Settings for traditional SCM run
    nelemd_todo = 1
    np_todo = 1
    
    if (iop_mode) then
      nelemd_todo = nelemd
      np_todo = np
    endif

    !   IE   Cp*dpdn*T  + (Cpv-Cp) Qdpdn*T
    !        Cp*dpdn(n)*T(n+1) + (Cpv-Cp) Qdpdn(n)*T(n+1)
    !        [Cp + (Cpv-Cp) Q(n)] *dpdn(n)*T(n+1) 

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif

    do ie=1,nelemd_todo

      do k=1,nlev
        p(:,:,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(:,:,t1)
      end do

      dt=dtime

      do j=1,np_todo
        do i=1,np_todo

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

          call forecast(begchunk,elem(ie)%state%ps_v(i,j,t1),&
            elem(ie)%state%ps_v(i,j,t1),forecast_ps,forecast_u,&
            elem(ie)%state%v(i,j,1,:,t1),elem(ie)%state%v(i,j,1,:,t1),&
            forecast_v,elem(ie)%state%v(i,j,2,:,t1),&
            elem(ie)%state%v(i,j,2,:,t1),forecast_t,&
#ifdef MODEL_THETA_L	    
            elem(ie)%derived%FT(i,j,:),elem(ie)%derived%FT(i,j,:),&
#else
            elem(ie)%state%T(i,j,:,t1),elem(ie)%state%T(i,j,:,t1),&
#endif
            forecast_q,stateQin2,stateQin1,dt,dummy1,dummy2,dummy2,&
            stateQin_qfcst,p(i,j,:),stateQin1,1,&
            tdiff_dyn,qdiff_dyn)         

#ifdef MODEL_THETA_L
          elem(ie)%derived%FT(i,j,:) = forecast_t(:)
#else
          elem(ie)%state%T(i,j,:,t1) = forecast_t(:)
#endif
          elem(ie)%state%v(i,j,1,:,t1) = forecast_u(:)
          elem(ie)%state%v(i,j,2,:,t1) = forecast_v(:)
          elem(ie)%state%Q(i,j,:,:) = forecast_q(:,:)
  
          tdiff_out(i+(j-1)*np,:)=tdiff_dyn(:)
          qdiff_out(i+(j-1)*np,:)=qdiff_dyn(:)
  
        enddo
      enddo
      
      if (iop_mode) then
        call outfld('TDIFF',tdiff_out,npsq,ie)
        call outfld('QDIFF',qdiff_out,npsq,ie)
      else
        call outfld('TDIFF',tdiff_dyn,plon,begchunk)
        call outfld('QDIFF',qdiff_dyn,plon,begchunk)
      endif
    
    enddo

    end subroutine apply_SC_forcing

end module se_single_column_mod
