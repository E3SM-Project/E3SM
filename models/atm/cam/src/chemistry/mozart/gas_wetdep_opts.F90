!-----------------------------------------------------------------------
! Reads namelist options for gas-phase wet deposition
!
! Created by Francis Vitt -- 22 Apr 2011
!-----------------------------------------------------------------------
module gas_wetdep_opts

  use constituents, only : pcnst
  use cam_logfile,  only : iulog
  use constituents, only : pcnst
  use spmd_utils,   only : masterproc
  use abortutils,   only : endrun

  implicit none

  character(len=8) :: gas_wetdep_list(pcnst) = ' '
  character(len=3) :: gas_wetdep_method = 'MOZ'
  integer :: gas_wetdep_cnt = 0

contains

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  subroutine gas_wetdep_readnl(nlfile)

    use abortutils,      only: endrun
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
#ifdef SPMD
    use mpishorthand,    only: mpichar, mpicom
#endif

    implicit none

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    integer :: unitn, i, ierr

    namelist /wetdep_inparm/ gas_wetdep_list
    namelist /wetdep_inparm/ gas_wetdep_method

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'wetdep_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, wetdep_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun('mo_neu_wetdep->wetdep_readnl: ERROR reading wetdep_inparm namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    call mpibcast (gas_wetdep_list, len(gas_wetdep_list(1))*pcnst, mpichar, 0, mpicom)
    call mpibcast (gas_wetdep_method, len(gas_wetdep_method), mpichar, 0, mpicom)
#endif

    gas_wetdep_cnt = 0
    do i = 1,pcnst
       if ( len_trim(gas_wetdep_list(i)) > 0 ) then
          gas_wetdep_cnt = gas_wetdep_cnt + 1
       endif
    enddo

    if (( gas_wetdep_cnt>0 ).and.( .not.(gas_wetdep_method=='MOZ' .or. gas_wetdep_method=='NEU') )) then
       call endrun('gas_wetdep_readnl; gas_wetdep_method must be set to either MOZ or NEU')
    endif

  end subroutine gas_wetdep_readnl

end module gas_wetdep_opts
