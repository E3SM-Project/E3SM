module rdycore_varctl
  !
  use shr_kind_mod, only: r8 => shr_kind_r8, SHR_KIND_CL, SHR_KIND_CS
  use shr_sys_mod , only: shr_sys_abort
  !
  implicit none
  public :: rdycore_varctl_set

  private
  save

  integer , parameter, public ::  iundef = -9999999
  real(r8), parameter, public ::  rundef = -9999999._r8
  integer , parameter, public ::  fname_len = SHR_KIND_CL   ! max length of file names in this module

  character(len=256), public :: caseid = ' ' ! case id
  character(len=256), public :: ctitle = ' ' ! case title

  integer, public, parameter :: nsrStartup  = 0        ! Startup from initial conditions
  integer, public, parameter :: nsrContinue = 1        ! Continue from restart files
  integer, public, parameter :: nsrBranch   = 2        ! Branch from restart files
  integer, public :: nsrest = iundef ! type of run: 0 = initial run; 1 = restart; 3 = branch

  logical, private :: rdycore_varctl_isset = .false.

  character(len=256), public :: rpntdir = '.'            ! directory name for local restart pointer file
  character(len=256), public :: rpntfil = 'rpointer.rof' ! file name for local restart pointer file

  character(len=256), public :: hostname = ' ' ! Hostname of machine running on
  character(len=256), public :: username = ' ' ! username of user running program

  integer, public :: iulog = 6 ! "stdout" log file unit number, default is 6

  ! instance control
  integer, public :: inst_index
  character(len=16), public :: inst_name
  character(len=16), public :: inst_suffix

contains

  !---------------------------------------------------------------------------
  subroutine rdycore_varctl_Set(caseid_in, ctitle_in, nsrest_in, hostname_in, username_in)
    !
    ! !DESCRIPTION:
    ! Set input control variables.
    !
    ! !ARGUMENTS:
    character(len=256), optional, intent(IN) :: caseid_in   ! case id
    character(len=256), optional, intent(IN) :: ctitle_in   ! case title
    integer,            optional, intent(IN) :: nsrest_in   ! 0: initial run. 1: restart: 3: branch
    character(len=256), optional, intent(IN) :: hostname_in ! hostname running on
    character(len=256), optional, intent(IN) :: username_in ! username running job

    if ( rdycore_varctl_isset )then
       call shr_sys_abort(' ERROR:: control variables already set, cannot call this routine')
    end if

    if ( present(caseid_in   ) ) caseid   = caseid_in
    if ( present(ctitle_in   ) ) ctitle   = ctitle_in
    if ( present(nsrest_in   ) ) nsrest   = nsrest_in
    if ( present(username_in ) ) username = username_in
    if ( present(hostname_in ) ) hostname = hostname_in

    rdycore_varctl_isset = .true.

  end subroutine rdycore_varctl_Set

end module rdycore_varctl
