module t_drv_timers_mod

  use perf_mod

contains

  !===============================================================================

  subroutine t_drvstartf(string,cplrun,cplcom,budget,barrier)

    implicit none

    character(len=*),intent(in) :: string
    logical,intent(in),optional :: cplrun
    logical,intent(in),optional :: cplcom
    logical,intent(in),optional :: budget
    integer,intent(in),optional :: barrier

    character(len=128) :: strbar
    character(len=*),parameter :: strcpl = 'DRIVER_CPL_RUN'
    character(len=*),parameter :: strcom = 'DRIVER_CPL_COMM'
    character(len=*),parameter :: strbud = 'DRIVER_BUDGET'
    logical :: lcplrun,lcplcom,lbudget

    !-------------------------------------------------------------------------------

    lcplrun  = .false.
    lcplcom  = .false.
    lbudget  = .false.
    if (present(cplrun)) then
       lcplrun = cplrun
    endif
    if (present(cplcom)) then
       lcplcom = cplcom
    endif
    if (present(budget)) then
       lbudget = budget
    endif

    if (present(barrier)) then
       strbar = trim(string)//'_BARRIER'
       call t_barrierf (trim(strbar), barrier)
    endif

    if (lcplrun) then
       call t_startf   (trim(strcpl))
    endif

    if (lcplcom) then
       call t_startf   (trim(strcom))
    endif

    if (lbudget) then
       call t_startf   (trim(strbud))
    endif

    call t_startf   (trim(string))

  end subroutine t_drvstartf

  !===============================================================================

  subroutine t_drvstopf(string,cplrun,cplcom,budget)

    implicit none

    character(len=*),intent(in) :: string
    logical,intent(in),optional :: cplrun
    logical,intent(in),optional :: cplcom
    logical,intent(in),optional :: budget

    character(len=128) :: strbar
    character(len=*),parameter :: strcpl = 'DRIVER_CPL_RUN'
    character(len=*),parameter :: strcom = 'DRIVER_CPL_COMM'
    character(len=*),parameter :: strbud = 'DRIVER_BUDGET'
    logical :: lcplrun,lcplcom,lbudget

    !-------------------------------------------------------------------------------

    lcplrun = .false.
    lcplcom = .false.
    lbudget = .false.
    if (present(cplrun)) then
       lcplrun = cplrun
    endif
    if (present(cplcom)) then
       lcplcom = cplcom
    endif
    if (present(budget)) then
       lbudget = budget
    endif

    !  strbar = trim(string)//'_BARRIER'

    call t_stopf   (trim(string))

    if (lbudget) then
       call t_stopf   (trim(strbud))
    endif

    if (lcplrun) then
       call t_stopf   (trim(strcpl))
    endif

    if (lcplcom) then
       call t_stopf   (trim(strcom))
    endif

  end subroutine t_drvstopf

end module t_drv_timers_mod
