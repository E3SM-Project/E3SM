module shr_nuopc_utils_mod

  implicit none
  private

  public :: shr_nuopc_memcheck
  public :: shr_nuopc_utils_ChkErr
  public :: shr_nuopc_log_clock_advance

  integer     , parameter :: memdebug_level=1
  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine shr_nuopc_memcheck(string, level, mastertask)
    character(len=*), intent(in) :: string
    integer, intent(in) :: level
    logical, intent(in) :: mastertask
    integer :: ierr
    integer, external :: GPTLprint_memusage
    if((mastertask .and. memdebug_level > level) .or. memdebug_level > level+1) then
       ierr = GPTLprint_memusage(string)
    endif
  end subroutine shr_nuopc_memcheck

!===============================================================================

  logical function shr_nuopc_utils_ChkErr(rc, line, file, mpierr)

    use mpi , only : MPI_ERROR_STRING, MPI_MAX_ERROR_STRING, MPI_SUCCESS
    use ESMF, only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_LOGMSG_INFO
    use ESMF, only : ESMF_FAILURE, ESMF_LogWrite

    integer, intent(in) :: rc
    integer, intent(in) :: line

    character(len=*), intent(in) :: file
    logical, optional, intent(in) :: mpierr

    character(MPI_MAX_ERROR_STRING) :: lstring
    integer :: dbrc, lrc, len, ierr

    shr_nuopc_utils_ChkErr = .false.
    lrc = rc
    if (present(mpierr) .and. mpierr) then
       if (rc == MPI_SUCCESS) return
       call MPI_ERROR_STRING(rc, lstring, len, ierr)
       call ESMF_LogWrite("ERROR: "//trim(lstring), ESMF_LOGMSG_INFO, line=line, file=file, rc=dbrc)
       lrc = ESMF_FAILURE
    endif

    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
      shr_nuopc_utils_ChkErr = .true.
    endif

  end function shr_nuopc_utils_ChkErr

!===============================================================================

  subroutine shr_nuopc_log_clock_advance(clock, component, logunit)
    use ESMF, only : ESMF_Clock, ESMF_ClockPrint
    use med_constants_mod, only : CL

    type(ESMF_Clock) :: clock
    character(len=*), intent(in) :: component
    integer, intent(in) :: logunit

    character(len=CL) :: cvalue, prestring
    integer :: rc

    write(prestring, *) "------>Advancing ",trim(component)," from: "
    call ESMF_ClockPrint(clock, options="currTime", unit=cvalue, &
         preString=trim(prestring), rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    write(logunit, *) trim(cvalue)

    call ESMF_ClockPrint(clock, options="stopTime", unit=cvalue, &
         preString="--------------------------------> to: ", rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    write(logunit, *) trim(cvalue)

  end subroutine shr_nuopc_log_clock_advance

end module shr_nuopc_utils_mod
