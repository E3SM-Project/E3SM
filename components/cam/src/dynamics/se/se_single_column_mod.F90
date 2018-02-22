module se_single_column_mod
!--------------------------------------------------------
! 
! Module for the SE single column model

use element_mod, only: element_t
use scamMod
use constituents, only: cnst_get_ind
use dimensions_mod, only: nelemd, np
use time_manager, only: get_nstep

implicit none

public scm_setinitial
public scm_setfield

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
              if (have_omega) elem(ie)%derived%omega_p(i,j,k) = wfld(k)
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
    if (have_ps) elem(ie)%state%ps_v(:,:,1) = psobs 
    do i=1, PLEV
      if (have_omega) elem(ie)%derived%omega_p(:,:,i)=wfld(i)  !     set t to tobs at first
    end do
  end do

end subroutine scm_setfield

end module se_single_column_mod
