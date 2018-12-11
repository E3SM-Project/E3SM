module shr_nuopc_utils_mod

  use shr_sys_mod    , only : shr_nuopc_abort => shr_sys_abort
  use shr_string_mod , only : shr_nuopc_string_listGetName => shr_string_listGetName

  implicit none
  private
  public :: shr_nuopc_abort, shr_nuopc_string_listGetName
  public :: shr_nuopc_memcheck
  public :: shr_nuopc_get_component_instance
  public :: shr_nuopc_set_component_logging
  public :: shr_nuopc_utils_ChkErr
  public :: shr_nuopc_log_clock_advance

  integer, parameter :: memdebug_level=1
  character(*),parameter :: u_FILE_u = __FILE__

contains
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

  subroutine shr_nuopc_get_component_instance(gcomp, inst_suffix, inst_index)
    use ESMF, only : ESMF_SUCCESS, ESMF_GridComp
    use NUOPC, only : NUOPC_CompAttributeGet

    type(ESMF_GridComp) :: gcomp
    character(len=*), intent(out) :: inst_suffix
    integer, intent(out) :: inst_index
    integer :: rc
    logical :: isPresent
    character(len=4) :: cvalue

    call NUOPC_CompAttributeGet(gcomp, name="inst_suffix", isPresent=isPresent, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name="inst_suffix", value=inst_suffix, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
       cvalue = inst_suffix(2:)
       read(cvalue, *) inst_index
    else
       inst_suffix = ""
       inst_index=1
    endif

  end subroutine shr_nuopc_get_component_instance

  subroutine shr_nuopc_set_component_logging(gcomp, mastertask, logunit, shrlogunit, shrloglev)
    use ESMF, only : ESMF_GridComp, ESMF_VM, ESMF_VMGet, ESMF_GridCompGet
    use NUOPC, only : NUOPC_CompAttributeGet
    use med_constants_mod, only : shr_file_getunit, shr_file_getLogUnit, shr_file_getLogLevel
    use med_constants_mod, only : shr_file_setLogLevel, CL, shr_file_setlogunit

    type(ESMF_GridComp) :: gcomp
    logical, intent(in) :: mastertask
    integer, intent(out) :: logunit
    integer, intent(out) :: shrlogunit
    integer, intent(out) :: shrloglev

    character(len=CL) :: diro
    character(len=CL) :: logfile
    integer :: rc
    shrlogunit = 6
    if (mastertask) then
       call NUOPC_CompAttributeGet(gcomp, name="diro", value=diro, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeGet(gcomp, name="logfile", value=logfile, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
       open(newunit=logunit,file=trim(diro)//"/"//trim(logfile))
    else
       logUnit = 6
    endif
    call shr_file_setLogUnit (logunit)
  end subroutine shr_nuopc_set_component_logging

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

  !-----------------------------------------------------------------------------
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
