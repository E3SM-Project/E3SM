module util

  !-----------------------------------------------------------------------------
  ! CustomFieldDictionaryProto utility module
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  
  implicit none
  
  private
  
  public FieldDictionaryLog
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  subroutine FieldDictionaryLog(label, iofmt, rc)

    character(len=*),      intent(in)            :: label
    type(ESMF_IOFmt_Flag), intent(in),  optional :: iofmt
    integer,               intent(out), optional :: rc

    integer                    :: ibeg, iend, length
    character(len=ESMF_MAXSTR) :: sep, title
    type(NUOPC_FreeFormat)     :: freeFormat

    if (present(rc)) rc = ESMF_SUCCESS

    write(sep,'(64("="))')

    ! build section separator with title
    title = "> Begin Field Dictionary: " // trim(label) // " <"
    length = len_trim(title)

    ! center title within separator
    ibeg = max((64 - length)/2,1)
    iend = min(ibeg+length-1,ESMF_MAXSTR)
    sep(ibeg:iend) = title(1:iend-ibeg+1)

    call ESMF_LogWrite(sep, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_FieldDictionaryEgest(freeFormat, iofmt=iofmt, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_FreeFormatLog(freeFormat, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! build section separator with title
    title = "> End Field Dictionary: " // trim(label) // " <"
    length = len_trim(title)

    ! align with opening title
    sep = ""
    write(sep,'(64("="))')
    iend = min(ibeg+length-1,ESMF_MAXSTR)
    sep(ibeg:iend) = title(1:iend-ibeg+1)

    call ESMF_LogWrite(sep, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine FieldDictionaryLog

end module
