#ifdef AIX
@PROCESS ALIAS_SIZE(805306368)
#endif
module dwav_comp_mod

  ! !USES:

  use esmf
  use mct_mod
  use perf_mod
  use shr_pcdf_mod
  use shr_sys_mod
  use shr_kind_mod      , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_file_mod      , only: shr_file_getunit, shr_file_freeunit
  use shr_mpi_mod       , only: shr_mpi_bcast
  use shr_strdata_mod   , only: shr_strdata_type, shr_strdata_pioinit, shr_strdata_init
  use shr_strdata_mod   , only: shr_strdata_print, shr_strdata_restRead
  use shr_strdata_mod   , only: shr_strdata_advance, shr_strdata_restWrite
  use shr_dmodel_mod    , only: shr_dmodel_gsmapcreate, shr_dmodel_rearrGGrid
  use shr_dmodel_mod    , only: shr_dmodel_translate_list, shr_dmodel_translateAV_list, shr_dmodel_translateAV
  use seq_timemgr_mod   , only: seq_timemgr_EClockGetData, seq_timemgr_RestartAlarmIsOn

  use dwav_shr_mod   , only: datamode       ! namelist input
  use dwav_shr_mod   , only: decomp         ! namelist input
  use dwav_shr_mod   , only: rest_file      ! namelist input
  use dwav_shr_mod   , only: rest_file_strm ! namelist input
  use dwav_shr_mod   , only: nullstr

  ! !PUBLIC TYPES:
  implicit none
  save
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: dwav_comp_init
  public :: dwav_comp_run
  public :: dwav_comp_final

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  character(CS)              :: myModelName = 'wav'   ! user defined model name
  character(len=*),parameter :: rpfile = 'rpointer.wav'
  type(mct_rearr)            :: rearr

  !--------------------------------------------------------------------------
  integer(IN),parameter :: ktrans = 3
  character(12),parameter  :: avofld(1:ktrans) =  (/"Sw_lamult   ","Sw_ustokes  ","Sw_vstokes  "/)
  character(12),parameter  :: avifld(1:ktrans) =  (/"lamult      ","ustokes     ","vstokes     "/)
  !--------------------------------------------------------------------------

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !===============================================================================
  subroutine dwav_comp_init(Eclock, x2w, w2x, &
       SDWAV, gsmap, ggrid, mpicom, compid, my_task, master_task, &
       inst_suffix, inst_name, logunit, read_restart, &
       seq_flds_w2x_fields, seq_flds_x2w_fields)

    ! !DESCRIPTION: initialize dwav model
    use pio        , only : iosystem_desc_t
    use shr_pio_mod, only : shr_pio_getiosys, shr_pio_getiotype
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)       , intent(in)    :: EClock
    type(mct_aVect)        , intent(inout) :: x2w, w2x            ! input/output attribute vectors
    type(shr_strdata_type) , intent(inout) :: SDWAV               ! model
    type(mct_gsMap)        , pointer       :: gsMap               ! model global seg map (output)
    type(mct_gGrid)        , pointer       :: ggrid               ! model ggrid (output)
    integer(IN)            , intent(in)    :: mpicom              ! mpi communicator
    integer(IN)            , intent(in)    :: compid              ! mct comp id
    integer(IN)            , intent(in)    :: my_task             ! my task in mpi communicator mpicom
    integer(IN)            , intent(in)    :: master_task         ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix         ! char string associated with instance
    character(len=*)       , intent(in)    :: inst_name           ! fullname of current instance (ie. "wav_0001")
    integer(IN)            , intent(in)    :: logunit             ! logging unit number
    logical                , intent(in)    :: read_restart        ! start from restart
    character(len=*)       , intent(in)    :: seq_flds_x2w_fields ! fields to mediator
    character(len=*)       , intent(in)    :: seq_flds_w2x_fields ! fields from mediator

    !--- local variables ---
    integer(IN)   :: n,k       ! generic counters
    integer(IN)   :: ierr      ! error code
    integer(IN)   :: lsize     ! local size
    logical       :: exists    ! file existance
    integer(IN)   :: nu        ! unit number
    character(CL) :: calendar  ! model calendar
    type(iosystem_desc_t), pointer :: wav_pio_subsystem

    !--- formats ---
    character(*), parameter :: F00   = "('(dwav_comp_init) ',8a)"
    character(*), parameter :: F0L   = "('(dwav_comp_init) ',a, l2)"
    character(*), parameter :: F01   = "('(dwav_comp_init) ',a,5i8)"
    character(*), parameter :: F02   = "('(dwav_comp_init) ',a,4es13.6)"
    character(*), parameter :: F03   = "('(dwav_comp_init) ',a,i8,a)"
    character(*), parameter :: F04   = "('(dwav_comp_init) ',2a,2i8,'s')"
    character(*), parameter :: F05   = "('(dwav_comp_init) ',a,2f10.4)"
    character(*), parameter :: F90   = "('(dwav_comp_init) ',73('='))"
    character(*), parameter :: F91   = "('(dwav_comp_init) ',73('-'))"
    character(*), parameter :: subName = "(dwav_comp_init) "
    !-------------------------------------------------------------------------------

    call t_startf('DWAV_INIT')

    !----------------------------------------------------------------------------
    ! Initialize pio
    !----------------------------------------------------------------------------

    wav_pio_subsystem => shr_pio_getiosys(trim(inst_name))
    call shr_strdata_pioinit(SDWAV, wav_pio_subsystem, shr_pio_getiotype(trim(inst_name)))

    !----------------------------------------------------------------------------
    ! Initialize SDWAV
    !----------------------------------------------------------------------------

    call t_startf('dwav_strdata_init')

    call seq_timemgr_EClockGetData( EClock, calendar=calendar )

    ! NOTE: shr_strdata_init calls shr_dmodel_readgrid which reads the data model
    ! grid and from that computes SDWAV%gsmap and SDWAV%ggrid. DWAV%gsmap is created
    ! using the decomp '2d1d' (1d decomp of 2d grid)

    call shr_strdata_init(SDWAV,mpicom,compid,name='wav', calendar=calendar)

    if (my_task == master_task) then
       call shr_strdata_print(SDWAV,'SDWAV data')
    endif

    call t_stopf('dwav_strdata_init')

    !----------------------------------------------------------------------------
    ! Initialize MCT global seg map, 1d decomp
    !----------------------------------------------------------------------------

    call t_startf('dwav_initgsmaps')

    if (my_task == master_task) write(logunit,F00) ' initialize gsmaps'
    call shr_sys_flush(logunit)

    ! create a data model global seqmap (gsmap) given the data model global grid sizes
    ! NOTE: gsmap is initialized using the decomp read in from the docn_in namelist
    ! (which by default is "1d")
    call shr_dmodel_gsmapcreate(gsmap,SDWAV%nxg*SDWAV%nyg,compid,mpicom,decomp)
    lsize = mct_gsmap_lsize(gsmap,mpicom)

    ! create a rearranger from the data model SDOCN%gsmap to gsmap
    call mct_rearr_init(SDWAV%gsmap,gsmap,mpicom,rearr)

    write(logunit,*)'lsize= ',lsize
    call shr_sys_flush(logunit)

    call t_stopf('dwav_initgsmaps')

    !----------------------------------------------------------------------------
    ! Initialize MCT domain
    !----------------------------------------------------------------------------

    call t_startf('dwav_initmctdom')
    !write(logunit,F00)' dwav_initmctdom...'

    if (my_task == master_task) write(logunit,F00) 'copy domains'
    call shr_sys_flush(logunit)

    call shr_dmodel_rearrGGrid(SDWAV%grid, ggrid, gsmap, rearr, mpicom)

    call t_stopf('dwav_initmctdom')

    !----------------------------------------------------------------------------
    ! Initialize MCT attribute vectors
    !----------------------------------------------------------------------------

    call t_startf('dwav_initmctavs')
    !write(logunit,F00)' dwav_initmctavs...'

    if (my_task == master_task) write(logunit,F00) 'allocate AVs'
    call shr_sys_flush(logunit)

    call mct_avect_init(w2x, rlist=seq_flds_w2x_fields, lsize=lsize)
    call mct_avect_zero(w2x)

    call mct_avect_init(x2w, rlist=seq_flds_x2w_fields, lsize=lsize)
    call mct_avect_zero(x2w)

    call t_stopf('dwav_initmctavs')

    !----------------------------------------------------------------------------
    ! Read restart
    !----------------------------------------------------------------------------

    if (read_restart) then
       if (trim(rest_file) == trim(nullstr) .and. trim(rest_file_strm) == trim(nullstr)) then
          if (my_task == master_task) then
             write(logunit,F00) ' restart filenames from rpointer'
             call shr_sys_flush(logunit)
             inquire(file=trim(rpfile)//trim(inst_suffix),exist=exists)
             if (.not.exists) then
                write(logunit,F00) ' ERROR: rpointer file does not exist'
                call shr_sys_abort(trim(subname)//' ERROR: rpointer file missing')
             endif
             nu = shr_file_getUnit()
             open(nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
             read(nu,'(a)') rest_file
             read(nu,'(a)') rest_file_strm
             close(nu)
             call shr_file_freeUnit(nu)
             inquire(file=trim(rest_file_strm),exist=exists)
          endif
          call shr_mpi_bcast(rest_file,mpicom,'rest_file')
          call shr_mpi_bcast(rest_file_strm,mpicom,'rest_file_strm')
       else
          ! use namelist already read
          if (my_task == master_task) then
             write(logunit,F00) ' restart filenames from namelist '
             call shr_sys_flush(logunit)
             inquire(file=trim(rest_file_strm),exist=exists)
          endif
       endif
       call shr_mpi_bcast(exists,mpicom,'exists')
       if (exists) then
          if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file_strm)
          call shr_strdata_restRead(trim(rest_file_strm),SDWAV,mpicom)
       else
          if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file_strm)
       endif
       call shr_sys_flush(logunit)
    endif

    !----------------------------------------------------------------------------
    ! Set initial wav state, needed for CCSM atm initialization
    !----------------------------------------------------------------------------

    call t_adj_detailf(+2)
    call dwav_comp_run(EClock, x2w, w2x, &
         SDWAV, gsmap, ggrid, mpicom, compid, my_task, master_task, &
         inst_suffix, logunit)
    call t_adj_detailf(-2)

    if (my_task == master_task) write(logunit, F00) 'dwav_comp_init done'
    call shr_sys_flush(logunit)

    call t_stopf('DWAV_INIT')

  end subroutine dwav_comp_init

  !===============================================================================
  subroutine dwav_comp_run(EClock, x2w, w2x, &
       SDWAV, gsmap, ggrid, mpicom, compid, my_task, master_task, &
       inst_suffix, logunit, case_name)

    use shr_cal_mod, only : shr_cal_ymdtod2string
    ! !DESCRIPTION:  run method for dwav model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)       , intent(in)    :: EClock
    type(mct_aVect)        , intent(inout) :: x2w
    type(mct_aVect)        , intent(inout) :: w2x
    type(shr_strdata_type) , intent(inout) :: SDWAV
    type(mct_gsMap)        , pointer       :: gsMap
    type(mct_gGrid)        , pointer       :: ggrid
    integer(IN)            , intent(in)    :: mpicom           ! mpi communicator
    integer(IN)            , intent(in)    :: compid           ! mct comp id
    integer(IN)            , intent(in)    :: my_task          ! my task in mpi communicator mpicom
    integer(IN)            , intent(in)    :: master_task      ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix      ! char string associated with instance
    integer(IN)            , intent(in)    :: logunit          ! logging unit number
    character(CL)          , intent(in), optional :: case_name ! case name

    !--- local ---
    integer(IN)   :: CurrentYMD            ! model date
    integer(IN)   :: CurrentTOD            ! model sec into model date
    integer(IN)   :: yy,mm,dd              ! year month day
    integer(IN)   :: n                     ! indices
    integer(IN)   :: idt                   ! integer timestep
    real(R8)      :: dt                    ! timestep
    integer(IN)   :: nu                    ! unit number
    logical       :: write_restart         ! restart now
    character(len=18) :: date_str

    character(*), parameter :: F00   = "('(dwav_comp_run) ',8a)"
    character(*), parameter :: F04   = "('(dwav_comp_run) ',2a,2i8,'s')"
    character(*), parameter :: subName = "(dwav_comp_run) "
    !-------------------------------------------------------------------------------

    call t_startf('DWAV_RUN')

    call t_startf('dwav_run1')
    call seq_timemgr_EClockGetData( EClock, curr_ymd=CurrentYMD, curr_tod=CurrentTOD)
    call seq_timemgr_EClockGetData( EClock, curr_yr=yy, curr_mon=mm, curr_day=dd)
    call seq_timemgr_EClockGetData( EClock, dtime=idt)
    dt = idt * 1.0_r8
    write_restart = seq_timemgr_RestartAlarmIsOn(EClock)
    call t_stopf('dwav_run1')

    !--------------------
    ! UNPACK
    !--------------------

    call t_startf('dwav_unpack')
    ! Nothing to be done for now
    call t_stopf('dwav_unpack')

    !--------------------
    ! ADVANCE WAV
    !--------------------

    call t_barrierf('dwav_BARRIER',mpicom)
    call t_startf('dwav')

    call t_startf('dwav_strdata_advance')
    call shr_strdata_advance(SDWAV,currentYMD,currentTOD,mpicom,'dwav')
    call t_stopf('dwav_strdata_advance')

    !--- copy all fields from streams to w2x as default ---
    call t_barrierf('dwav_scatter_BARRIER',mpicom)
    call t_startf('dwav_scatter')
    do n = 1,SDWAV%nstreams
       call shr_dmodel_translateAV(SDWAV%avs(n),w2x,avifld,avofld,rearr)
    enddo
    call t_stopf('dwav_scatter')

    !-------------------------------------------------
    ! Determine data model behavior based on the mode
    !-------------------------------------------------

    call t_startf('datamode')
    select case (trim(datamode))

    case('COPYALL')
       ! do nothing extra

    end select

    call t_stopf('datamode')

    !--------------------
    ! Write restart
    !--------------------

    if (write_restart) then
       call t_startf('dwav_restart')
       call shr_cal_ymdtod2string(date_str, yy,mm,dd,currentTOD)
       write(rest_file,"(6a)") &
            trim(case_name), '.dwav',trim(inst_suffix),'.r.', &
            trim(date_str),'.nc'
       write(rest_file_strm,"(6a)") &
            trim(case_name), '.dwav',trim(inst_suffix),'.rs1.', &
            trim(date_str),'.bin'
       if (my_task == master_task) then
          nu = shr_file_getUnit()
          open(nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
          write(nu,'(a)') rest_file
          write(nu,'(a)') rest_file_strm
          close(nu)
          call shr_file_freeUnit(nu)
       endif
       if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file_strm),currentYMD,currentTOD
       call shr_strdata_restWrite(trim(rest_file_strm),SDWAV,mpicom,trim(case_name),'SDWAV strdata')
       call shr_sys_flush(logunit)
       call t_stopf('dwav_restart')
    endif

    call t_stopf('dwav')

    !----------------------------------------------------------------------------
    ! Log output for model date
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call t_startf('dwav_run2')
    if (my_task == master_task) then
       write(logunit,F04) trim(myModelName),': model date ', CurrentYMD,CurrentTOD
       call shr_sys_flush(logunit)
    end if
    call t_stopf('dwav_run2')

    call t_stopf('DWAV_RUN')

  end subroutine dwav_comp_run

  !===============================================================================
  subroutine dwav_comp_final(my_task, master_task, logunit)

    ! !DESCRIPTION:  finalize method for dwav model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    integer(IN) , intent(in) :: my_task     ! my task in mpi communicator mpicom
    integer(IN) , intent(in) :: master_task ! task number of master task
    integer(IN) , intent(in) :: logunit     ! logging unit number

    !--- formats ---
    character(*), parameter :: F00   = "('(dwav_comp_final) ',8a)"
    character(*), parameter :: F91   = "('(dwav_comp_final) ',73('-'))"
    character(*), parameter :: subName = "(dwav_comp_final) "
    !-------------------------------------------------------------------------------

    call t_startf('DWAV_FINAL')

    if (my_task == master_task) then
       write(logunit,F91)
       write(logunit,F00) trim(myModelName),': end of main integration loop'
       write(logunit,F91)
    end if

    call t_stopf('DWAV_FINAL')

  end subroutine dwav_comp_final
  !===============================================================================

end module dwav_comp_mod
