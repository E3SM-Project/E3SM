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
  use seq_timemgr_mod, only: seq_timemgr_EClockGetData
  use seq_timemgr_mod, only: seq_timemgr_RestartAlarmIsOn
  use seq_comm_mct,    only: num_inst_cpl => num_inst_driver

  ! Used to link esp components across multiple drivers
  use seq_comm_mct,    only: global_comm

  implicit none
  private

#include <mpif.h>

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public                     :: desp_comp_init
  public                     :: desp_comp_run
  public                     :: desp_comp_final
  public                     :: desp_bcast_res_files

  !--------------------------------------------------------------------------
  ! Public module data
  !--------------------------------------------------------------------------
  integer,          public, parameter :: desp_num_comps = 8
  character(len=3), public, parameter :: comp_names(desp_num_comps) =         &
       (/ 'atm', 'lnd', 'ice', 'ocn', 'glc', 'rof', 'wav', 'drv' /)

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  character(len=CS)           :: myModelName = 'esp'    ! user defined model name
  integer(IN)                 :: mpicom
  integer(IN)                 :: COMPID                 ! mct comp id
  integer(IN)                 :: my_task                ! my task in mpicom
  integer(IN)                 :: npes                   ! total number of tasks
  integer(IN),      parameter :: master_task=0          ! task number of master task
  integer(IN)                 :: global_numpes          ! #PEs in global_commm
  integer(IN)                 :: global_mype            ! My rank in global_comm
  integer(IN)                 :: logunit                ! logging unit number
  integer(IN)                 :: loglevel               ! logging level
  integer(IN)                 :: inst_index             ! number of current instance (ie. 1)
  character(len=16)           :: inst_name              ! fullname of current instance (ie. "esp_0001")
  character(len=16)           :: inst_suffix            ! char string associated with instance (ie. "_0001" or "")
  character(len=CL)           :: desp_mode              ! mode of operation
  logical                     :: bcast_filenames = .false.
  character(len=*), parameter :: rpprefix  = 'rpointer.'
  character(len=*), parameter :: rpfile    = rpprefix//'esp'
  character(len=*), parameter :: nullstr   = 'undefined'
  character(len=*), parameter :: null_mode = 'NULL'     ! Take no action
  character(len=*), parameter :: noop_mode = 'NOCHANGE' ! Do not change data
  character(len=*), parameter :: test_mode = 'DATATEST' ! Modify restart data

  integer,          parameter :: NOERR       =  0
  integer,          parameter :: BAD_ID      = -1
  integer,          parameter :: NO_READ     = -5
  integer,          parameter :: NO_WRITE    = -6

  type(shr_strdata_type)      :: SDESP

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !============================================================================
  subroutine desp_comp_init(EClock, espid, mpicom_in, phase, read_restart,    &
       inst_name_in, inst_index_in, inst_suffix_in, esp_present, esp_prognostic)

    ! !DESCRIPTION: initialize data esp model

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock),  intent(in)  :: EClock
    integer,           intent(in)  :: espid
    integer,           intent(in)  :: mpicom_in
    integer,           intent(in)  :: phase
    logical,           intent(in)  :: read_restart
    integer,           intent(in)  :: inst_index_in  ! (e.g., 1)
    character(len=16), intent(in)  :: inst_name_in   ! (e.g. "exp_0001")
    character(len=16), intent(in)  :: inst_suffix_in ! (e.g. "_0001" or "")
    logical,           intent(out) :: esp_present    ! flag
    logical,           intent(out) :: esp_prognostic ! flag

                                                     !--- local variables ---
    integer(IN)                    :: ierr           ! error code
    integer(IN)                    :: shrlogunit     ! original log unit
    integer(IN)                    :: shrloglev      ! original log level
    integer(IN)                    :: nunit          ! unit number

    character(len=CL)              :: fileName       ! generic file name

    character(len=CL)              :: rest_file      ! restart filename
    character(len=CL)              :: restfilm       ! rest_file namelist entry
    logical                        :: exists         ! filename existance
    integer(IN)                    :: nu             ! unit number
    integer(IN)                    :: CurrentYMD     ! model date
    integer(IN)                    :: CurrentTOD     ! model sec into model date
    integer(IN)                    :: stepno         ! step number
    character(len=CL)              :: calendar       ! calendar type

    !----- define namelist -----
    namelist / desp_nml /                                                     &
         desp_mode, restfilm

    !--- formats ---
    character(len=*), parameter :: subName = "(desp_comp_init) "
    character(len=*), parameter :: F00     = "('"//subName//"',8a)"
    character(len=*), parameter :: F01     = "('"//subName//"',a,5i8)"
    character(len=*), parameter :: F04     = "('"//subName//"',2(a,i0))"
    !-------------------------------------------------------------------------------

    call t_startf('DESP_INIT')

    !------------------------------------------------------------------------
    ! Initialize module variables from inputs
    !------------------------------------------------------------------------
    COMPID      = espid
    mpicom      = mpicom_in
    inst_name   = inst_name_in
    inst_index  = inst_index_in
    inst_suffix = inst_suffix_in

    !------------------------------------------------------------------------
    ! Initialize output variables
    !------------------------------------------------------------------------
    esp_present    = .false.
    esp_prognostic = .false.

    if (phase == 1) then
       ! Determine communicator groups and sizes
       call mpi_comm_rank(mpicom, my_task, ierr)
       call mpi_comm_size(mpicom, npes, ierr)
       call mpi_comm_rank(global_comm, global_mype, ierr)
       call mpi_comm_size(global_comm, global_numpes, ierr)

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
          write(logunit,F00) 'restfilm    = ',trim(restfilm)
          write(logunit,F01) 'inst_index  = ',inst_index
          write(logunit,F00) 'inst_name   = ',trim(inst_name)
          write(logunit,F00) 'inst_suffix = ',trim(inst_suffix)
          write(logunit,F04) 'root global rank ',global_mype,' of ',global_numpes
          call shr_sys_flush(logunit)
       end if
       call shr_mpi_bcast(desp_mode,  mpicom, 'desp_mode')
       call shr_mpi_bcast(restfilm,   mpicom, 'restfilm')

       rest_file = trim(restfilm)
       loglevel = 1 ! could be shrloglev
       call shr_file_setLogLevel(loglevel)

       !------------------------------------------------------------------------
       ! Initialize PIO
       !------------------------------------------------------------------------

       call shr_strdata_pioinit(SDESP, COMPID)

       !------------------------------------------------------------------------
       ! Validate mode
       !------------------------------------------------------------------------

       if ( (trim(desp_mode) == null_mode) .or.                                 &
            (trim(desp_mode) == noop_mode) .or.                                 &
            (trim(desp_mode) == test_mode)) then

          if (my_task == master_task) then
             write(logunit,F00) 'desp mode = ',trim(desp_mode)
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
                write(logunit,F00) ' restart filename from rpointer = ',trim(rpfile)
                call shr_sys_flush(logunit)
                inquire(file=trim(rpfile)//trim(inst_suffix),exist=exists)
                if (exists) then
                   nu = shr_file_getUnit()
                   open(nu, file=trim(rpfile)//trim(inst_suffix), form='formatted')
                   read(nu,'(a)') rest_file
                   close(nu)
                   call shr_file_freeUnit(nu)
                   inquire(file=trim(rest_file),exist=exists)
                endif
             endif
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
  subroutine desp_comp_run(EClock, case_name, pause_sig, atm_resume,          &
       lnd_resume, rof_resume, ocn_resume, ice_resume, glc_resume,            &
       wav_resume, cpl_resume)

    ! !DESCRIPTION: run method for data esp model

    use seq_comm_mct, only: num_inst_atm, num_inst_lnd, num_inst_rof
    use seq_comm_mct, only: num_inst_ocn, num_inst_ice, num_inst_glc
    use seq_comm_mct, only: num_inst_wav
    use esp_utils,    only: esp_pio_modify_variable
    use shr_cal_mod, only: shr_cal_ymdtod2string

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock), intent(in) :: EClock
    character(len=*), intent(in) :: case_name
    logical,          intent(in) :: pause_sig(desp_num_comps)
    character(len=CL), pointer   :: atm_resume(:)
    character(len=CL), pointer   :: lnd_resume(:)
    character(len=CL), pointer   :: rof_resume(:)
    character(len=CL), pointer   :: ocn_resume(:)
    character(len=CL), pointer   :: ice_resume(:)
    character(len=CL), pointer   :: glc_resume(:)
    character(len=CL), pointer   :: wav_resume(:)
    character(len=CL), pointer   :: cpl_resume(:)

    !--- local ---
    integer(IN)                      :: CurrentYMD    ! model date
    integer(IN)                      :: CurrentTOD    ! model sec into model date
    integer(IN)                      :: yy,mm,dd      ! year month day
    integer(IN)                      :: ind           ! loop index
    integer(IN)                      :: inst          ! loop index
    integer(IN)                      :: errcode       ! error return code
    character(len=CL)                :: rfilename     ! Restart filenames
    integer(IN)                      :: shrlogunit    ! original log unit
    integer(IN)                      :: shrloglev     ! original log level
    integer(IN)                      :: idt           ! integer timestep
    real(R8)                         :: dt            ! timestep
    logical                          :: write_restart ! restart now
    logical                          :: var_found     ! var on file
    character(len=CL)                :: rest_file     ! restart_file
    integer(IN)                      :: nu            ! unit number
    integer(IN)                      :: stepno        ! step number
    character(len=CL)                :: calendar      ! calendar type
    character(len=CS)                :: varname
    character(len=18)                :: date_str
    character(len=*), parameter      :: subName = "(desp_comp_run) "
    character(len=*), parameter      :: F00   = "('"//subName//"',8a)"
    character(len=*), parameter      :: F04   = "('"//subName//"',2a,2i8,'s')"
    !--------------------------------------------------------------------------

    call t_startf('DESP_RUN')

    call t_startf('desp_run1')

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

    if (.not. ANY(pause_sig)) then
       if ( (my_task == master_task) .and.                                     &
            ((loglevel > 1) .or. (trim(desp_mode) == test_mode))) then
          write(logunit, '(2a,i6.4,"-",i2.2,"-",i2.2,"-",i5.5)') subname,       &
               'WARNING: No pause signals found at ',yy,mm,dd,CurrentTOD
       end if
    end if

    errcode = NOERR
    ! Find the active components and their restart files
    ! Note hard-coded variable names are just for testing. This could be
    !   changed if this feature comes to regular use
    do ind = 1, desp_num_comps
       if (pause_sig(ind)) then
          select case (comp_names(ind))
          case('atm')
             rfilename = atm_resume(inst_index)
             varname = 'T'
          case('lnd')
             rfilename = lnd_resume(inst_index)
             varname = 'T'
          case('ice')
             rfilename = ice_resume(inst_index)
             varname = 'T'
          case('ocn')
             rfilename = ocn_resume(inst_index)
             varname = 'PSURF_CUR'
          case('glc')
             rfilename = glc_resume(inst_index)
             varname = 'T'
          case('rof')
             rfilename = rof_resume(inst_index)
             varname = 'T'
          case('wav')
             rfilename = wav_resume(inst_index)
             varname = 'T'
          case('drv')
             ! The driver is special, there may only be one (not multi-driver)
             rfilename = cpl_resume(min(inst_index,size(cpl_resume,1)))
             varname = 'x2oacc_ox_Foxx_swnet'
          case default
             call shr_sys_abort(subname//'Unrecognized ind')
          end select
          ! Just die on errors for now
          select case (errcode)
          case(NO_READ)
             call shr_sys_abort(subname//'No restart read access for '//comp_names(ind))
          case(NO_WRITE)
             call shr_sys_abort(subname//'No restart write access for '//comp_names(ind))
             ! No default case needed, just fall through
          end select
          select case (trim(desp_mode))
          case(noop_mode)
             ! Find the correct restart files but do not change them.
             if ((my_task == master_task) .and. (loglevel > 0)) then
                write(logunit, *) subname, 'Found restart file ',trim(rfilename)
             end if
          case(test_mode)
             ! Find the correct restart files and 'tweak' them
             if ((my_task == master_task) .and. (loglevel > 0)) then
                write(logunit, *) subname, 'Found restart file ',trim(rfilename)
             end if
             call esp_pio_modify_variable(COMPID, mpicom, rfilename, varname, var_found)
             if (.not. var_found) then
                call shr_sys_abort(subname//'Variable, '//trim(varname)//', not found on '//rfilename)
             end if
          case (null_mode)
             ! Since DESP is not 'present' for this mode, we should not get here.
             call shr_sys_abort(subname//'DESP should not run in "NULL" mode')
          end select
       end if
    end do

    call t_stopf('desp_mode')

    if (write_restart) then
       call t_startf('desp_restart')
       call shr_cal_ymdtod2string(date_str, yy,mm,dd,currentTOD)
       write(rest_file,"(6a)")                     &
            trim(case_name), '.desp',trim(inst_suffix),'.r.',                &
            trim(date_str),'.nc'
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
  subroutine desp_comp_final()

    ! !DESCRIPTION: finalize method for data esp model
    !--- formats ---
    character(len=*), parameter :: subName = "(desp_comp_final) "
    character(len=*), parameter :: F00   = "('"//subName//"',8a)"
    character(len=*), parameter :: F91   = "('"//subName//"',73('-'))"
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
  logical function desp_bcast_res_files(oneletterid)
     character(len=1), intent(in) :: oneletterid

     desp_bcast_res_files = bcast_filenames
  end function desp_bcast_res_files

  !============================================================================

end module desp_comp_mod
