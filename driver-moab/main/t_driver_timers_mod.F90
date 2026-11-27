module t_drv_timers_mod

  use perf_mod
  integer, private :: cpl_run_hash=0, cpl_comm_hash=0, cpl_budget_hash=0
  character(len=*),parameter :: strcpl = 'CPL:RUN'
  character(len=*),parameter :: strcom = 'CPL:COMM'
  character(len=*),parameter :: strbud = 'CPL:BUDGET'

contains

  !===============================================================================

  subroutine t_drvstartf(string,cplrun,cplcom,budget,barrier, hashint)

    implicit none

    character(len=*),intent(in) :: string
    logical,intent(in),optional :: cplrun
    logical,intent(in),optional :: cplcom
    logical,intent(in),optional :: budget
    integer,intent(in),optional :: barrier
    integer,intent(inout), optional :: hashint

    character(len=128) :: strbar

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
       call t_startf   (trim(strcpl), cpl_run_hash)
       call t_adj_detailf(+1)
    endif

    if (lcplcom) then
       call t_startf   (trim(strcom), cpl_comm_hash)
       call t_adj_detailf(+1)
    endif

    if (lbudget) then
       call t_startf   (trim(strbud), cpl_budget_hash)
       call t_adj_detailf(+1)
    endif

    call t_startf   (trim(string),hashint)
    call t_adj_detailf(+1)

  end subroutine t_drvstartf

  !===============================================================================

  subroutine t_drvstopf(string,cplrun,cplcom,budget,hashint)

    implicit none

    character(len=*),intent(in) :: string
    logical,intent(in),optional :: cplrun
    logical,intent(in),optional :: cplcom
    logical,intent(in),optional :: budget
    integer, intent(in), optional :: hashint
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

    call t_adj_detailf(-1)
    call t_stopf   (trim(string), hashint)

    if (lbudget) then
       call t_adj_detailf(-1)
       call t_stopf   (trim(strbud), cpl_budget_hash)
    endif

    if (lcplrun) then
       call t_adj_detailf(-1)
       call t_stopf   (trim(strcpl), cpl_run_hash)
    endif

    if (lcplcom) then
       call t_adj_detailf(-1)
       call t_stopf   (trim(strcom),cpl_comm_hash)
    endif

  end subroutine t_drvstopf

  !===============================================================================
  ! MOAB sub-phase timers.
  !
  ! These bracket individual iMOAB operations (tag send / receive / free sender
  ! buffers / apply weights / tag get-set) so the cost of a coupler phase can be
  ! split into communication versus local work.  They are opened one timing
  ! detail level below their enclosing phase timer, so they are controlled by
  ! the existing prof_inparm namelist knob profile_detail_limit: lower it below
  ! the enclosing phase's detail level to switch them all off without a rebuild.
  !===============================================================================

  subroutine t_mbstartf(string)

    implicit none
    character(len=*),intent(in) :: string
    !-------------------------------------------------------------------------------

    call t_adj_detailf(+1)
    call t_startf (trim(string))

  end subroutine t_mbstartf

  !===============================================================================

  subroutine t_mbstopf(string)

    implicit none
    character(len=*),intent(in) :: string
    !-------------------------------------------------------------------------------

    call t_stopf (trim(string))
    call t_adj_detailf(-1)

  end subroutine t_mbstopf

end module t_drv_timers_mod
