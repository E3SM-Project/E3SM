module dshr_util_mod

  use ESMF         , only : ESMF_LOGERR_PASSTHRU, ESMF_LogFoundError
  use shr_kind_mod , only : r8=>shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl, cxx=>shr_kind_cxx

  implicit none
  public

  integer, parameter      :: memdebug_level=1
  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine memcheck(string, level, mastertask)

    ! input/output variables
    character(len=*) , intent(in) :: string
    integer          , intent(in) :: level
    logical          , intent(in) :: mastertask

    ! local variables
    integer :: ierr
    integer, external :: GPTLprint_memusage
    !-----------------------------------------------------------------------

    if ((mastertask .and. memdebug_level > level) .or. memdebug_level > level+1) then
       ierr = GPTLprint_memusage(string)
    endif

  end subroutine memcheck

  !===============================================================================
  logical function chkerr(rc, line, file)
    integer, intent(in) :: rc
    integer, intent(in) :: line
    character(len=*), intent(in) :: file
    integer :: lrc
    !-----------------------------------------------------------------------
    chkerr = .false.
    lrc = rc
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
       chkerr = .true.
    endif
  end function chkerr

end module dshr_util_mod
