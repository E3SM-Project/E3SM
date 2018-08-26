#ifdef AIX
@PROCESS ALIAS_SIZE(805306368)
#endif
module dwav_comp_mod

  ! !USES:

  use NUOPC                 , only : NUOPC_Advertise
  use ESMF                  , only : ESMF_State, ESMF_SUCCESS
  use perf_mod              , only : t_startf, t_stopf
  use perf_mod              , only : t_adj_detailf, t_barrierf
  use mct_mod               , only : mct_rearr, mct_gsmap_lsize, mct_rearr_init, mct_gsmap, mct_ggrid
  use mct_mod               , only : mct_avect, mct_avect_indexRA, mct_avect_zero, mct_aVect_nRattr
  use mct_mod               , only : mct_avect_init, mct_avect_lsize, mct_avect_clean, mct_aVect
  use shr_sys_mod           , only : shr_sys_abort
  use shr_kind_mod          , only : IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_kind_mod          , only : CXX=>SHR_KIND_CXX 
  use shr_string_mod        , only : shr_string_listGetName
  use shr_sys_mod           , only : shr_sys_abort
  use shr_file_mod          , only : shr_file_getunit, shr_file_freeunit
  use shr_mpi_mod           , only : shr_mpi_bcast
  use shr_strdata_mod       , only : shr_strdata_type, shr_strdata_pioinit, shr_strdata_init
  use shr_strdata_mod       , only : shr_strdata_print, shr_strdata_restRead
  use shr_strdata_mod       , only : shr_strdata_advance, shr_strdata_restWrite
  use shr_dmodel_mod        , only : shr_dmodel_gsmapcreate, shr_dmodel_rearrGGrid, shr_dmodel_translateAV
  use shr_cal_mod           , only : shr_cal_datetod2string
  use shr_nuopc_scalars_mod , only : flds_scalar_name
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
  use dshr_nuopc_mod        , only : fld_list_type
  use dshr_nuopc_mod        , only : dshr_fld_add
  use dwav_shr_mod          , only : datamode       ! namelist input
  use dwav_shr_mod          , only : decomp         ! namelist input
  use dwav_shr_mod          , only : rest_file      ! namelist input
  use dwav_shr_mod          , only : rest_file_strm ! namelist input
  use dwav_shr_mod          , only : nullstr

  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: dwav_comp_advertise
  public :: dwav_comp_init
  public :: dwav_comp_run
  public :: dwav_comp_final

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  type(mct_rearr)             :: rearr
  character(len=CS), pointer  :: avifld(:) ! character array for field names coming from streams
  character(len=CS), pointer  :: avofld(:) ! character array for field names to be sent/received from mediator
  character(len=CXX)          :: flds_w2x_mod
  character(len=CXX)          :: flds_x2w_mod
  integer                     :: dbug = 2  ! debug level (higher is more)
  character(len=*), parameter :: rpfile = 'rpointer.wav'
  character(*)    , parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine dwav_comp_advertise(importState, exportState, &
       wav_present, wav_prognostic, &
       fldsFrWav_num, fldsFrWav, fldsToWav_num, fldsToWav, &
       flds_w2x, flds_x2w, rc)

    ! 1. determine export and import fields to advertise to mediator
    ! 2. determine translation of fields from streams to export/import fields

    ! input/output arguments
    type(ESMF_State)                   :: importState
    type(ESMF_State)                   :: exportState
    logical              , intent(in)  :: wav_present
    logical              , intent(in)  :: wav_prognostic
    integer              , intent(out) :: fldsFrWav_num
    type (fld_list_type) , intent(out) :: fldsFrWav(:)
    integer              , intent(out) :: fldsToWav_num
    type (fld_list_type) , intent(out) :: fldsToWav(:)
    character(len=*)     , intent(out) :: flds_w2x
    character(len=*)     , intent(out) :: flds_x2w
    integer              , intent(out) :: rc

    ! local variables 
    integer :: n
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (.not. wav_present) return

    !-------------------
    ! export fields
    !-------------------
    
    ! scalar fields that need to be advertised

    fldsFrWav_num=1
    fldsFrWav(1)%stdname = trim(flds_scalar_name)

    ! export fields that have a corresponding stream field

    call dshr_fld_add(data_fld="lamult", data_fld_array=avifld, model_fld="Sw_lamult", model_fld_array=avofld, &
         model_fld_concat=flds_w2x, fldlist_num=fldsFrWav_num, fldlist=fldsFrWav)

    call dshr_fld_add(data_fld="ustokes", data_fld_array=avifld, model_fld="Sw_ustokes", model_fld_array=avofld, &
         model_fld_concat=flds_w2x, fldlist_num=fldsFrWav_num, fldlist=fldsFrWav)

    call dshr_fld_add(data_fld="vstokes", data_fld_array=avifld, model_fld="Sw_vstokes", model_fld_array=avofld, &
         model_fld_concat=flds_w2x, fldlist_num=fldsFrWav_num, fldlist=fldsFrWav)

    !-------------------
    ! advertise export state
    !-------------------

    do n = 1,fldsFrWav_num
       call NUOPC_Advertise(exportState, standardName=fldsFrWav(n)%stdname, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    !-------------------
    ! Save flds_w2x and flds_x2w as module variables for use in debugging
    !-------------------

    flds_x2w_mod = trim(flds_x2w)
    flds_w2x_mod = trim(flds_w2x)

  end subroutine dwav_comp_advertise

  !===============================================================================

  subroutine dwav_comp_init(x2w, w2x, &
       w2x_fields, x2w_fields, &
       SDWAV, gsmap, ggrid, mpicom, compid, my_task, master_task, &
       inst_suffix, inst_name, logunit, read_restart, calendar)

    ! !DESCRIPTION: initialize dwav model

    ! !INPUT/OUTPUT PARAMETERS:
    type(mct_aVect)        , intent(inout) :: x2w, w2x            ! input/output attribute vectors
    character(len=*)       , intent(in)    :: x2w_fields          ! fields to mediator
    character(len=*)       , intent(in)    :: w2x_fields          ! fields from mediator
    type(shr_strdata_type) , intent(inout) :: SDWAV               ! model
    type(mct_gsMap)        , pointer       :: gsMap               ! model global seg map (output)
    type(mct_gGrid)        , pointer       :: ggrid               ! model ggrid (output)
    integer                , intent(in)    :: mpicom              ! mpi communicator
    integer                , intent(in)    :: compid              ! mct comp id
    integer                , intent(in)    :: my_task             ! my task in mpi communicator mpicom
    integer                , intent(in)    :: master_task         ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix         ! char string associated with instance
    character(len=*)       , intent(in)    :: inst_name           ! fullname of current instance (ie. "wav_0001")
    integer                , intent(in)    :: logunit             ! logging unit number
    logical                , intent(in)    :: read_restart        ! start from restart
    character(len=*)       , intent(in)    :: calendar         ! calendar type

    !--- local variables ---
    integer       :: n,k       ! generic counters
    integer       :: lsize     ! local size
    logical       :: exists    ! file existance
    integer       :: nu        ! unit number
    character(*), parameter :: F00   = "('(dwav_comp_init) ',8a)"
    character(*), parameter :: subName = "(dwav_comp_init) "
    !-------------------------------------------------------------------------------

    call t_startf('DWAV_INIT')

    !----------------------------------------------------------------------------
    ! Initialize pio
    !----------------------------------------------------------------------------

    call shr_strdata_pioinit(SDWAV, compid)

    !----------------------------------------------------------------------------
    ! Initialize SDWAV
    !----------------------------------------------------------------------------

    call t_startf('dwav_strdata_init')

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

    ! create a data model global seqmap (gsmap) given the data model global grid sizes
    ! NOTE: gsmap is initialized using the decomp read in from the dwav_in namelist
    ! (which by default is "1d")
    call shr_dmodel_gsmapcreate(gsmap,SDWAV%nxg*SDWAV%nyg,compid,mpicom,decomp)
    lsize = mct_gsmap_lsize(gsmap,mpicom)

    ! create a rearranger from the data model SDWAV%gsmap to gsmap
    call mct_rearr_init(SDWAV%gsmap,gsmap,mpicom,rearr)
    call t_stopf('dwav_initgsmaps')

    !----------------------------------------------------------------------------
    ! Initialize MCT domain
    !----------------------------------------------------------------------------

    call t_startf('dwav_initmctdom')
    if (my_task == master_task) write(logunit,F00) 'copy domains'
    call shr_dmodel_rearrGGrid(SDWAV%grid, ggrid, gsmap, rearr, mpicom)
    call t_stopf('dwav_initmctdom')

    !----------------------------------------------------------------------------
    ! Initialize MCT attribute vectors
    !----------------------------------------------------------------------------

    if (my_task == master_task) write(logunit,F00) 'allocate AVs'

    call mct_avect_init(w2x, rlist=w2x_fields, lsize=lsize)
    call mct_avect_zero(w2x)
    call mct_avect_init(x2w, rlist=x2w_fields, lsize=lsize)
    call mct_avect_zero(x2w)

    !----------------------------------------------------------------------------
    ! Read restart
    !----------------------------------------------------------------------------

    if (read_restart) then
       if (trim(rest_file) == trim(nullstr) .and. trim(rest_file_strm) == trim(nullstr)) then
          if (my_task == master_task) then
             write(logunit,F00) ' restart filenames from rpointer'
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
    endif

    call t_stopf('DWAV_INIT')

  end subroutine dwav_comp_init

  !===============================================================================

  subroutine dwav_comp_run(x2w, w2x, &
       SDWAV, gsmap, ggrid, mpicom, my_task, master_task, &
       inst_suffix, logunit, read_restart, write_restart, &
       target_ymd, target_tod, case_name)

    ! DESCRIPTION:  run method for dwav model

    ! input/output parameters:
    type(mct_aVect)        , intent(inout) :: x2w
    type(mct_aVect)        , intent(inout) :: w2x
    type(shr_strdata_type) , intent(inout) :: SDWAV
    type(mct_gsMap)        , pointer       :: gsMap
    type(mct_gGrid)        , pointer       :: ggrid
    integer                , intent(in)    :: mpicom           ! mpi communicator
    integer                , intent(in)    :: my_task          ! my task in mpi communicator mpicom
    integer                , intent(in)    :: master_task      ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix      ! char string associated with instance
    integer                , intent(in)    :: logunit          ! logging unit number
    logical                , intent(in)    :: read_restart     ! start from restart
    logical                , intent(in)    :: write_restart    ! write restart
    integer(IN)            , intent(in)    :: target_ymd
    integer(IN)            , intent(in)    :: target_tod
    character(CL)          , intent(in), optional :: case_name ! case name

    !--- local ---
    integer       :: n                     ! indices
    integer       :: idt                   ! integer timestep
    integer       :: nu                    ! unit number
    character(len=18) :: date_str

    character(*), parameter :: F00   = "('(dwav_comp_run) ',8a)"
    character(*), parameter :: F04   = "('(dwav_comp_run) ',2a,2i8,'s')"
    character(*), parameter :: subName = "(dwav_comp_run) "
    !-------------------------------------------------------------------------------

    call t_startf('DWAV_RUN')

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
    call shr_strdata_advance(SDWAV,target_ymd,target_tod,mpicom,'dwav')
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
       call shr_cal_datetod2string(date_str, target_ymd, target_tod)
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
       if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file),target_ymd,target_tod
       call shr_strdata_restWrite(trim(rest_file_strm),SDWAV,mpicom,trim(case_name),'SDWAV strdata')
       call t_stopf('dwav_restart')
    endif

    call t_stopf('dwav')

    !----------------------------------------------------------------------------
    ! Log output for model date
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    if (my_task == master_task) then
       write(logunit,*) 'dwav: model date ', target_ymd,target_tod
    end if

    call t_stopf('DWAV_RUN')

  end subroutine dwav_comp_run

  !===============================================================================

  subroutine dwav_comp_final(my_task, master_task, logunit)

    ! !DESCRIPTION:  finalize method for dwav model

    ! !INPUT/OUTPUT PARAMETERS:
    integer , intent(in) :: my_task     ! my task in mpi communicator mpicom
    integer , intent(in) :: master_task ! task number of master task
    integer , intent(in) :: logunit     ! logging unit number

    !--- formats ---
    character(*), parameter :: F00   = "('(dwav_comp_final) ',8a)"
    character(*), parameter :: F91   = "('(dwav_comp_final) ',73('-'))"
    character(*), parameter :: subName = "(dwav_comp_final) "
    !-------------------------------------------------------------------------------

    if (my_task == master_task) then
       write(logunit,F91)
       write(logunit,F00) 'dwav: end of main integration loop'
       write(logunit,F91)
    end if

  end subroutine dwav_comp_final

end module dwav_comp_mod
