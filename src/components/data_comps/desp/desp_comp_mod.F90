module desp_comp_mod

! !USES:

  use shr_sys_mod,     only: shr_sys_abort, shr_sys_flush
  use shr_kind_mod,    only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8
  use shr_kind_mod,    only: CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_file_mod,    only: shr_file_getunit, shr_file_freeunit, shr_file_setio
  use shr_file_mod,    only: shr_file_getlogunit, shr_file_getloglevel
  use shr_file_mod,    only: shr_file_setlogunit, shr_file_setloglevel
  use shr_mpi_mod,     only: shr_mpi_bcast
  use esmf,            only: ESMF_Clock
  use perf_mod,        only: t_startf, t_stopf, t_barrierf

  use shr_strdata_mod, only: shr_strdata_type, shr_strdata_advance
  use shr_strdata_mod, only: shr_strdata_pioinit

  use seq_timemgr_mod, only: seq_timemgr_EClockGetData, seq_timemgr_RestartAlarmIsOn
  use seq_comm_mct,    only: seq_comm_inst, seq_comm_name, seq_comm_suffix

  implicit none
  private

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public                     :: desp_comp_init
  public                     :: desp_comp_run
  public                     :: desp_comp_final

!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

  integer, parameter, public :: atm_ind = 1
  integer, parameter, public :: lnd_ind = 2
  integer, parameter, public :: ice_ind = 3
  integer, parameter, public :: ocn_ind = 4
  integer, parameter, public :: glc_ind = 5
  integer, parameter, public :: rof_ind = 6
  integer, parameter, public :: wav_ind = 7
  integer, parameter, public :: cpl_ind = 8
  integer, parameter, public :: max_ind = cpl_ind

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  character(len=CS)           :: myModelName = 'esp'       ! user defined model name
  integer(IN)                 :: mpicom
  integer(IN)                 :: COMPID                    ! mct comp id
  integer(IN)                 :: my_task                   ! my task in mpi communicator mpicom
  integer(IN)                 :: npes                      ! total number of tasks
  integer(IN),      parameter :: master_task=0             ! task number of master task
  integer(IN)                 :: logunit                   ! logging unit number
  integer(IN)                 :: loglevel
  integer                     :: inst_index                ! number of current instance (ie. 1)
  character(len=16)           :: inst_name                 ! fullname of current instance (ie. "lnd_0001")
  character(len=16)           :: inst_suffix               ! char string associated with instance
                                                           ! (ie. "_0001" or "")
  character(len=CL)           :: desp_mode                 ! mode of operation

  character(len=3), parameter :: comp_names(8) = (/ 'atm', 'lnd', 'ice',      &
       'ocn', 'glc', 'rof', 'wav', 'drv' /)
  character(len=*), parameter :: rpprefix  = 'rpointer.'
  character(len=*), parameter :: rpfile    = rpprefix//'esp'
  character(len=*), parameter :: nullstr   = 'undefined'
  character(len=*), parameter :: null_mode = 'NULL'     ! Take no action
  character(len=*), parameter :: noop_mode = 'NOCHANGE' ! Do not change data
  character(len=*), parameter :: test_mode = 'DATATEST' ! Modify restart data

  integer,          parameter :: NOERR       =  0
  integer,          parameter :: BAD_ID      = -1
  integer,          parameter :: NO_RPOINTER = -2
  integer,          parameter :: NO_RFILE    = -3
  integer,          parameter :: NO_RPREAD   = -4
  integer,          parameter :: NO_READ     = -5
  integer,          parameter :: NO_WRITE    = -6

  type(shr_strdata_type)      :: SDESP

!--------------------------------------------------------------------------
! Private interface
!--------------------------------------------------------------------------
interface get_restart_filenames
  module procedure get_restart_filenames_a
  module procedure get_restart_filenames_s
end interface get_restart_filenames

  SAVE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: desp_comp_init
!
! !DESCRIPTION:
!     initialize data esp model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

  subroutine desp_comp_init(EClock, espid, mpicom_in, phase, read_restart,    &
       esp_present, esp_prognostic)
    use pio,         only: iosystem_desc_t
    use shr_pio_mod, only: shr_pio_getiosys, shr_pio_getiotype

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock), intent(in)  :: EClock
    integer,          intent(in)  :: espid
    integer,          intent(in)  :: mpicom_in
    integer,          intent(in)  :: phase
    logical,          intent(in)  :: read_restart
    logical,          intent(out) :: esp_present    ! flag
    logical,          intent(out) :: esp_prognostic ! flag
    
    !EOP
    
    !--- local variables ---
    integer(IN)                    :: ierr       ! error code
    integer(IN)                    :: shrlogunit ! original log unit
    integer(IN)                    :: shrloglev  ! original log level
    integer(IN)                    :: nunit      ! unit number

    type(iosystem_desc_t), pointer :: iosystem

    character(len=CL)              :: fileName   ! generic file name

    character(len=CL)              :: rest_file  ! restart filename
    character(len=CL)              :: restfilm   ! model restart file namelist
    logical                        :: exists     ! filename existance
    integer(IN)                    :: info_debug ! logging level
    integer(IN)                    :: nu         ! unit number
    integer(IN)                    :: CurrentYMD ! model date
    integer(IN)                    :: CurrentTOD ! model sec into model date
    integer(IN)                    :: stepno     ! step number
    character(len=CL)              :: calendar   ! calendar type

    !----- define namelist -----
    namelist / desp_nml /                                                     &
         desp_mode, info_debug, restfilm

    !--- formats ---
    character(*), parameter :: F00   = "('(desp_comp_init) ',8a)"
    character(*), parameter :: F0L   = "('(desp_comp_init) ',a, l2)"
    character(*), parameter :: F01   = "('(desp_comp_init) ',a,5i8)"
    character(*), parameter :: F02   = "('(desp_comp_init) ',a,4es13.6)"
    character(*), parameter :: F03   = "('(desp_comp_init) ',a,i8,a)"
    character(*), parameter :: F04   = "('(desp_comp_init) ',2a,2i8,'s')"
    character(*), parameter :: F05   = "('(desp_comp_init) ',a,2f10.4)"
    character(*), parameter :: F90   = "('(desp_comp_init) ',73('='))"
    character(*), parameter :: F91   = "('(desp_comp_init) ',73('-'))"
    character(*), parameter :: subName = "(desp_comp_init) "
!-------------------------------------------------------------------------------

    call t_startf('DESP_INIT')

    !------------------------------------------------------------------------
    ! Initialize module variables from inputs
    !------------------------------------------------------------------------
    COMPID       = espid
    mpicom       = mpicom_in

    !------------------------------------------------------------------------
    ! Initialize output variables
    !------------------------------------------------------------------------
    esp_present    = .false.
    esp_prognostic = .false.

    inst_name   = seq_comm_name(COMPID)
    inst_index  = seq_comm_inst(COMPID)
    inst_suffix = seq_comm_suffix(COMPID)

    if (phase == 1) then
      ! Determine communicator groups and sizes
      call mpi_comm_rank(mpicom, my_task, ierr)
      call mpi_comm_size(mpicom, npes, ierr)

      !--- open log file ---
      if (my_task == master_task) then
        logUnit = shr_file_getUnit()
        call shr_file_setIO('esp_modelio.nml'//trim(inst_suffix),logUnit)
      else
        logUnit = 6
      end if

      !------------------------------------------------------------------------
      ! Reset shr logging to my log file
      !------------------------------------------------------------------------
      call shr_file_getLogUnit (shrlogunit)
      call shr_file_getLogLevel(shrloglev)
      call shr_file_setLogUnit (logUnit)

      !------------------------------------------------------------------------
      ! Read desp_in
      !------------------------------------------------------------------------

      call t_startf('desp_readnml')

      filename   = "desp_in"//trim(inst_suffix)
      desp_mode  = trim(nullstr)
      info_debug = -1
      restfilm   = trim(nullstr)
      if (my_task == master_task) then
        nunit = shr_file_getUnit() ! get unused unit number
        open (nunit,file=trim(filename),status="old",action="read")
        read (nunit,nml=desp_nml,iostat=ierr)
        close(nunit)

        call shr_file_freeUnit(nunit)
        if (ierr > 0) then
          write(logunit,F01) 'ERROR: reading input namelist, '//trim(filename)//' iostat=',ierr
          call shr_sys_abort(subName//': namelist read error '//trim(filename))
        end if
        write(logunit,F01)' info_debug  = ',info_debug
        write(logunit,F00)' restfilm    = ',trim(restfilm)
        write(logunit,F01) 'inst_index  = ',inst_index
        write(logunit,F00) 'inst_name   = ',trim(inst_name)
        write(logunit,F00) 'inst_suffix = ',trim(inst_suffix)
        call shr_sys_flush(logunit)
      end if
      call shr_mpi_bcast(desp_mode,  mpicom, 'desp_mode')
      call shr_mpi_bcast(info_debug, mpicom, 'info_debug')
      call shr_mpi_bcast(restfilm,   mpicom, 'restfilm')

      rest_file = trim(restfilm)
      loglevel = info_debug
      call shr_file_setLogLevel(loglevel)

      !------------------------------------------------------------------------
      ! Initialize PIO
      !------------------------------------------------------------------------

      iosystem => shr_pio_getiosys(trim(inst_name))
      call shr_strdata_pioinit(SDESP, iosystem, shr_pio_getiotype(trim(inst_name)))

      !------------------------------------------------------------------------
      ! Validate mode
      !------------------------------------------------------------------------

      if ( (trim(desp_mode) == null_mode) .or.                                 &
           (trim(desp_mode) == noop_mode) .or.                                 &
           (trim(desp_mode) == test_mode)) then

        if (my_task == master_task) then
          write(logunit,F00) ' desp mode = ',trim(desp_mode)
          call shr_sys_flush(logunit)
        end if
        if (trim(desp_mode) /= null_mode) then
          esp_present = .true.
        end if
      else
        write(logunit,F00) ' ERROR illegal esp mode = ',trim(desp_mode)
        call shr_sys_abort(subName//' Illegal ESP mode = '//trim(desp_mode))
      end if

      call t_stopf('desp_readnml')

      !------------------------------------------------------------------------
      ! Read restart
      !------------------------------------------------------------------------

      if (read_restart) then
        if (trim(rest_file) == trim(nullstr)) then
          if (my_task == master_task) then
            write(logunit,F00) ' restart filename from rpointer'
            call shr_sys_flush(logunit)
            inquire(file=trim(rpfile)//trim(inst_suffix),exist=exists)
            if (.not. exists) then
              write(logunit,F00) ' ERROR: rpointer file does not exist'
              call shr_sys_abort(trim(subname)//' ERROR: rpointer file missing')
            end if
            nu = shr_file_getUnit()
            open(nu, file=trim(rpfile)//trim(inst_suffix), form='formatted')
            read(nu,'(a)') rest_file
            close(nu)
            call shr_file_freeUnit(nu)
            inquire(file=trim(rest_file),exist=exists)
          end if
          call shr_mpi_bcast(rest_file,mpicom,'rest_file')
        else
          ! use namelist already read
          if (my_task == master_task) then
            write(logunit,F00) ' restart filename from namelist '
            call shr_sys_flush(logunit)
            inquire(file=trim(rest_file),exist=exists)
          end if
        end if

        call shr_mpi_bcast(exists, mpicom, 'exists')

        if (exists) then
          if (my_task == master_task) then
            write(logunit,F00) ' reading ',trim(rest_file)
          end if
          ! Read any restart info here? --XXgoldyXX
        else
          if (my_task == master_task) then
            write(logunit,F00) ' file not found, skipping ',trim(rest_file)
          end if
        end if
        call shr_sys_flush(logunit)
      end if

      if (read_restart) then
        call seq_timemgr_EClockGetData( EClock, curr_ymd=CurrentYMD, curr_tod=CurrentTOD)
        call seq_timemgr_EClockGetData( EClock, calendar=calendar )
      else
        call seq_timemgr_EClockGetData( EClock, stepno=stepno )
      end if
    else
      call shr_sys_abort(trim(subname)//' DESP initialization only has phase 1')
    end if

    !--------------------------------------------------------------------------
    ! Reset shr logging to original values
    !--------------------------------------------------------------------------

    if (my_task == master_task) then
      write(logunit,F00) 'desp_comp_init done'
    end if
    call shr_sys_flush(logunit)

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

    call t_stopf('DESP_INIT')

  end subroutine desp_comp_init

  !=============================================================================
  !BOP==========================================================================
  !
  ! !IROUTINE: desp_comp_run
  !
  ! !DESCRIPTION:
  !     run method for data esp model
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ---------------------------------------------------------------

  subroutine desp_comp_run(EClock, case_name, pause_sig, atm_resume,          &
       lnd_resume, rof_resume, ocn_resume, ice_resume, glc_resume,            &
       wav_resume, cpl_resume)
    use seq_comm_mct, only: num_inst_atm, num_inst_lnd, num_inst_rof
    use seq_comm_mct, only: num_inst_ocn, num_inst_ice, num_inst_glc
    use seq_comm_mct, only: num_inst_wav
    use esp_utils,    only: esp_pio_modify_variable

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock),  intent(in)    :: EClock
    character(len=*),  intent(in)    :: case_name
    logical,           intent(in)    :: pause_sig(max_ind)
    character(len=CL), intent(inout) :: atm_resume(num_inst_atm)
    character(len=CL), intent(inout) :: lnd_resume(num_inst_lnd)
    character(len=CL), intent(inout) :: rof_resume(num_inst_rof)
    character(len=CL), intent(inout) :: ocn_resume(num_inst_ocn)
    character(len=CL), intent(inout) :: ice_resume(num_inst_ice)
    character(len=CL), intent(inout) :: glc_resume(num_inst_glc)
    character(len=CL), intent(inout) :: wav_resume(num_inst_wav)
    character(len=CL), intent(inout) :: cpl_resume

                                                              !--- local ---
    integer(IN)                      :: CurrentYMD            ! model date
    integer(IN)                      :: CurrentTOD            ! model sec into model date
    integer(IN)                      :: yy,mm,dd              ! year month day
    integer(IN)                      :: ind                   ! loop index
    integer(IN)                      :: inst                  ! loop index
    integer(IN)                      :: errcode               ! error return code
    character(len=CL), allocatable   :: rfilenames(:)         ! Restart filenames
    integer(IN)                      :: shrlogunit, shrloglev ! original log unit and level
    integer(IN)                      :: idt                   ! integer timestep
    real(R8)                         :: dt                    ! timestep
    logical                          :: write_restart         ! restart now
    logical                          :: var_found             ! var on file
    character(len=CL)                :: rest_file             ! restart_file
    integer(IN)                      :: nu                    ! unit number
    integer(IN)                      :: stepno                ! step number
    character(len=CL)                :: calendar              ! calendar type
    character(len=CS)                :: varname

    character(len=*), parameter      :: F00   = "('(desp_comp_run) ',8a)"
    character(len=*), parameter      :: F04   = "('(desp_comp_run) ',2a,2i8,'s')"
    character(len=*), parameter      :: subName = "(desp_comp_run) "
    !--------------------------------------------------------------------------

    call t_startf('DESP_RUN')

    call t_startf('desp_run1')

    !--------------------------------------------------------------------------
    ! Make sure output variables are set
    !--------------------------------------------------------------------------
    atm_resume(:) = ' '
    lnd_resume(:) = ' '
    rof_resume(:) = ' '
    ocn_resume(:) = ' '
    ice_resume(:) = ' '
    glc_resume(:) = ' '
    wav_resume(:) = ' '
    cpl_resume    = ' '

    !--------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !--------------------------------------------------------------------------
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logUnit)
    call shr_file_setLogLevel(loglevel)

    call seq_timemgr_EClockGetData( EClock, curr_ymd=CurrentYMD, curr_tod=CurrentTOD)
    call seq_timemgr_EClockGetData( EClock, curr_yr=yy, curr_mon=mm, curr_day=dd)
    call seq_timemgr_EClockGetData( EClock, stepno=stepno, dtime=idt)
    call seq_timemgr_EClockGetData( EClock, calendar=calendar)
    dt = idt * 1.0_r8
    write_restart = seq_timemgr_RestartAlarmIsOn(EClock)

    call t_stopf('desp_run1')

    !--------------------
    ! ADVANCE ESP
    !--------------------

    call t_barrierf('desp_BARRIER',mpicom)
    call t_startf('desp')

    if (trim(desp_mode) /= null_mode) then
      call t_startf('desp_strdata_advance')
      call shr_strdata_advance(SDESP, currentYMD, currentTOD, mpicom, 'desp')
      call t_stopf('desp_strdata_advance')
    end if

    call t_startf('desp_mode')

    ! Find the active components and their restart files
    ! Note hard-coded variable names are just for testing. This could be
    !   changed if this feature comes to regular use
    do ind = 1, max_ind
      if (pause_sig(ind)) then
        select case (ind)
        case(atm_ind)
          call get_restart_filenames(ind, atm_resume, errcode)
          allocate(rfilenames(size(atm_resume)))
          rfilenames = atm_resume
          varname = 'T'
        case(lnd_ind)
          call get_restart_filenames(ind, lnd_resume, errcode)
          allocate(rfilenames(size(lnd_resume)))
          rfilenames = lnd_resume
          varname = 'T'
        case(ice_ind)
          call get_restart_filenames(ind, ice_resume, errcode)
          allocate(rfilenames(size(ice_resume)))
          rfilenames = ice_resume
          varname = 'T'
        case(ocn_ind)
          call get_restart_filenames(ind, ocn_resume, errcode)
          allocate(rfilenames(size(ocn_resume)))
          rfilenames = ocn_resume
          varname = 'T'
        case(glc_ind)
          call get_restart_filenames(ind, glc_resume, errcode)
          allocate(rfilenames(size(glc_resume)))
          rfilenames = glc_resume
          varname = 'T'
        case(rof_ind)
          call get_restart_filenames(ind, rof_resume, errcode)
          allocate(rfilenames(size(rof_resume)))
          rfilenames = rof_resume
          varname = 'T'
        case(wav_ind)
          call get_restart_filenames(ind, wav_resume, errcode)
          allocate(rfilenames(size(wav_resume)))
          rfilenames = wav_resume
          varname = 'T'
        case(cpl_ind)
          call get_restart_filenames(ind, cpl_resume, errcode)
          allocate(rfilenames(1))
          rfilenames(1) = cpl_resume
          varname = 'a2x_ax_Sa_tbot'
        case default
          call shr_sys_abort(subname//'Unrecognized ind')
        end select
        ! Just die on errors for now
        select case (errcode)
        case(NO_RPOINTER)
          call shr_sys_abort(subname//'Missing rpointer file for '//comp_names(ind))
        case(NO_RPREAD)
          call shr_sys_abort(subname//'Cannot read rpointer file for '//comp_names(ind))
        case(NO_RFILE)
          call shr_sys_abort(subname//'Missing restart file for '//comp_names(ind))
        case(NO_READ)
          call shr_sys_abort(subname//'No restart read access for '//comp_names(ind))
        case(NO_WRITE)
          call shr_sys_abort(subname//'No restart write access for '//comp_names(ind))
          ! No default case needed, just fall through
        end select
        select case (trim(desp_mode))
        case(noop_mode)
          ! Find the correct restart files but do not change them.
          do inst = 1, size(rfilenames)
            if ((my_task == master_task) .and. (loglevel > 0)) then
              write(logunit, *) subname, 'Found restart file ',trim(rfilenames(inst))
            end if
          end do
        case(test_mode)
          ! Find the correct restart files and 'tweak' them
          do inst = 1, size(rfilenames)
            if ((my_task == master_task) .and. (loglevel > 0)) then
              write(logunit, *) subname, 'Found restart file ',trim(rfilenames(inst))
            end if
            call esp_pio_modify_variable(COMPID, mpicom, rfilenames(inst), varname, var_found)
          end do
        case (null_mode)
          ! Since DESP is not 'present' for this mode, we should not get here.
          call shr_sys_abort(subname//'DESP should not run in "NULL" mode')
        end select
        if (allocated(rfilenames)) then
          deallocate(rfilenames)
        end if
      end if
    end do

    call t_stopf('desp_mode')

    if (write_restart) then
      call t_startf('desp_restart')
      write(rest_file,"(2a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)")                     &
           trim(case_name), '.desp'//trim(inst_suffix)//'.r.',                &
           yy,'-',mm,'-',dd,'-',currentTOD,'.nc'
      if (my_task == master_task) then
        nu = shr_file_getUnit()
        open(nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
        write(nu,'(a)') rest_file
        close(nu)
        call shr_file_freeUnit(nu)
      end if

      if (my_task == master_task) then
        write(logunit,F04) ' writing ',trim(rest_file), currentYMD, currentTOD
      end if
      ! Write any restart info here? -- XXgoldyXX
      call shr_sys_flush(logunit)
      call t_stopf('desp_restart')

    end if

    call t_stopf('desp')

    !--------------------------------------------------------------------------
    ! Log output for model date
    ! Reset shr logging to original values
    !--------------------------------------------------------------------------

    call t_startf('desp_run2')
    if ((loglevel > 1) .and. (my_task == master_task)) then
      write(logunit,F04) trim(myModelName),': model date ', CurrentYMD,CurrentTOD
      call shr_sys_flush(logunit)
    end if

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)
    call t_stopf('desp_run2')

    call t_stopf('DESP_RUN')

  end subroutine desp_comp_run

  !============================================================================
  !BOP ========================================================================
  !
  ! !IROUTINE: desp_comp_final
  !
  ! !DESCRIPTION:
  !     finalize method for data esp model
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ---------------------------------------------------------------
  !
  subroutine desp_comp_final()

    !EOP

    !--- formats ---
    character(*), parameter :: F00   = "('(desp_comp_final) ',8a)"
    character(*), parameter :: F91   = "('(desp_comp_final) ',73('-'))"
    character(*), parameter :: subName = "(desp_comp_final) "

    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------

    call t_startf('DESP_FINAL')

    if (my_task == master_task) then
      write(logunit,F91)
      write(logunit,F00) trim(myModelName),': end of main integration loop'
      write(logunit,F91)
    end if

    call t_stopf('DESP_FINAL')

  end subroutine desp_comp_final
  !============================================================================
  !============================================================================

  subroutine get_restart_filenames_a(comp_ind, filenames, retcode)
    use seq_comm_mct, only: ATMID, LNDID, OCNID, ICEID, GLCID, ROFID
    use seq_comm_mct, only: WAVID, CPLID, seq_comm_suffix
    use shr_file_mod, only: shr_file_getUnit, shr_file_freeUnit

    ! Dummy arguments
    integer,           intent(in)    :: comp_ind ! Internal comp. type index
    character(len=CL), intent(inout) :: filenames(:)
    integer,           intent(out)   :: retcode

    ! Local variables
    integer                          :: num_inst
    integer                          :: ind
    integer                          :: ierr
    integer                          :: unitn
    integer,          allocatable    :: ids(:)
    logical                          :: file_exists
    character(len=8)                 :: file_read
    character(len=CS)                :: rpointer_name
    character(len=CL)                :: errmsg
    character(len=*), parameter      :: subname = "(desp_restart_filenames) "

    retcode = NOERR
    filenames = ' '
    num_inst = size(filenames)
    allocate(ids(num_inst))
    select case (comp_ind)
      case(atm_ind)
        ids = ATMID
      case(lnd_ind)
        ids = LNDID
      case(ice_ind)
        ids = ICEID
      case(ocn_ind)
        ids = OCNID
      case(glc_ind)
        ids = GLCID
      case(rof_ind)
        ids = ROFID
      case(wav_ind)
        ids = WAVID
      case(cpl_ind)
        ids = CPLID
      case default
        call shr_sys_abort(subname//'Unrecognized comp_ind')
    end select
    rpointer_name = rpprefix//comp_names(comp_ind)

    do ind = 1, num_inst
      rpointer_name = rpprefix//comp_names(comp_ind)//trim(seq_comm_suffix(ids(ind)))
      if (my_task == master_task) then
        inquire(file=rpointer_name, EXIST=file_exists)
        if (.not. file_exists) then
          retcode = NO_RPOINTER
        else
          unitn = shr_file_getUnit()
          if (loglevel > 0) then
            write(logunit,"(3A)") subname,"read rpointer file ", rpointer_name
          end if
          open(unitn, file=rpointer_name, form='FORMATTED', status='old',     &
               action='READ', iostat=ierr)
          if (ierr /= 0) then
            write(errmsg, '(a,i0)') 'rpointer file open returns error condition, ',ierr
            call shr_sys_abort(subname//trim(errmsg))
          end if
          read(unitn,'(a)', iostat=ierr) filenames(ind)
          if (ierr /= 0) then
            write(errmsg, '(a,i0)') 'rpointer file read returns error condition, ',ierr
            call shr_sys_abort(subname//trim(errmsg))
          end if
          close(unitn)
          call shr_file_freeUnit(unitn)
          inquire(file=filenames(ind), EXIST=file_exists)
          if (.not. file_exists) then
            retcode = NO_RFILE
            ! No else
          end if
        end if
      end if
      ! Broadcast what we learned about files
      call shr_mpi_bcast(retcode,  mpicom, 'rpointer status')
      call shr_mpi_bcast(filenames(ind),  mpicom, 'filenames(ind)')
    end do

    if (allocated(ids)) then
      deallocate(ids)
    end if

  end subroutine get_restart_filenames_a

  subroutine get_restart_filenames_s(comp_ind, filename, retcode)

    ! Dummy arguments
    integer,           intent(in)    :: comp_ind ! Internal comp. type index
    character(len=CL), intent(inout) :: filename
    integer,           intent(out)   :: retcode

    ! Local variable
    character(len=CL)                :: filenames(1)

    call get_restart_filenames(comp_ind, filenames, retcode)
    filename = filenames(1)
  end subroutine get_restart_filenames_s

end module desp_comp_mod
