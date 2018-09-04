#ifdef AIX
@PROCESS ALIAS_SIZE(805306368)
#endif
module drof_comp_mod

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
  use drof_shr_mod          , only : datamode       ! namelist input
  use drof_shr_mod          , only : decomp         ! namelist input
  use drof_shr_mod          , only : rest_file      ! namelist input
  use drof_shr_mod          , only : rest_file_strm ! namelist input
  use drof_shr_mod          , only : nullstr

  !
  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: drof_comp_advertise
  public :: drof_comp_init
  public :: drof_comp_run

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  type(mct_rearr)             :: rearr
  character(len=CS), pointer  :: avifld(:) ! character array for field names coming from streams
  character(len=CS), pointer  :: avofld(:) ! character array for field names to be sent/received from mediator
  character(len=CXX)          :: flds_r2x_mod
  character(len=CXX)          :: flds_x2r_mod
  character(len=*), parameter :: rpfile = 'rpointer.rof'
  character(*)    , parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine drof_comp_advertise(importState, exportState, &
       rof_present, rof_prognostic, &
       fldsFrRof_num, fldsFrRof, fldsToRof_num, fldsToRof, &
       flds_r2x, flds_x2r, rc)

    ! 1. determine export and import fields to advertise to mediator
    ! 2. determine translation of fields from streams to export/import fields

    ! input/output arguments
    type(ESMF_State)                   :: importState
    type(ESMF_State)                   :: exportState
    logical              , intent(in)  :: rof_present
    logical              , intent(in)  :: rof_prognostic
    integer              , intent(out) :: fldsFrRof_num
    type (fld_list_type) , intent(out) :: fldsFrRof(:)
    integer              , intent(out) :: fldsToRof_num
    type (fld_list_type) , intent(out) :: fldsToRof(:)
    character(len=*)     , intent(out) :: flds_r2x
    character(len=*)     , intent(out) :: flds_x2r
    integer              , intent(out) :: rc

    ! local variables 
    integer :: n
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (.not. rof_present) return

    !-------------------
    ! export fields
    !-------------------
    
    ! scalar fields that need to be advertised

    fldsFrRof_num=1
    fldsFrRof(1)%stdname = trim(flds_scalar_name)

    ! export fields that have a corresponding stream field

    call dshr_fld_add(data_fld="rofl", data_fld_array=avifld, model_fld="Forr_rofl", model_fld_array=avofld, &
         model_fld_concat=flds_r2x, fldlist_num=fldsFrRof_num, fldlist=fldsFrRof)

    call dshr_fld_add(data_fld="rofi", data_fld_array=avifld, model_fld="Forr_rofi", model_fld_array=avofld, &
         model_fld_concat=flds_r2x, fldlist_num=fldsFrRof_num, fldlist=fldsFrRof)

    !-------------------
    ! advertise export state
    !-------------------

    do n = 1,fldsFrRof_num
       call NUOPC_Advertise(exportState, standardName=fldsFrRof(n)%stdname, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    !-------------------
    ! Save flds_r2x and flds_x2r as module variables for use in debugging
    !-------------------

    flds_x2r_mod = trim(flds_x2r)
    flds_r2x_mod = trim(flds_r2x)

  end subroutine drof_comp_advertise

  !===============================================================================

  subroutine drof_comp_init(x2r, r2x, &
       SDROF, gsmap, ggrid, mpicom, compid, my_task, master_task, &
       inst_suffix, inst_name, logunit, read_restart, &
       target_ymd, target_tod, calendar)

    ! !DESCRIPTION: initialize drof model

    ! !INPUT/OUTPUT PARAMETERS:
    type(mct_aVect)        , intent(inout) :: x2r, r2x     ! input/output attribute vectors
    type(shr_strdata_type) , intent(inout) :: SDROF        ! model shr_strdata instance (output)
    type(mct_gsMap)        , pointer       :: gsMap        ! model global seg map (output)
    type(mct_gGrid)        , pointer       :: ggrid        ! model ggrid (output)
    integer(IN)            , intent(in)    :: mpicom       ! mpi communicator
    integer(IN)            , intent(in)    :: compid       ! mct comp id
    integer(IN)            , intent(in)    :: my_task      ! my task in mpi communicator mpicom
    integer(IN)            , intent(in)    :: master_task  ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix  ! char string associated with instance
    character(len=*)       , intent(in)    :: inst_name    ! fullname of current instance (ie. "lnd_0001")
    integer(IN)            , intent(in)    :: logunit      ! logging unit number
    logical                , intent(in)    :: read_restart ! start from restart
    integer                , intent(in)    :: target_ymd   ! model date
    integer                , intent(in)    :: target_tod   ! model sec into model date
    character(len=*)       , intent(in)    :: calendar     ! model calendar

    !--- local variables ---
    integer(IN)   :: n,k     ! generic counters
    integer(IN)   :: ierr    ! error code
    integer(IN)   :: lsize   ! local size
    logical       :: exists  ! file existance logical
    logical       :: exists1 ! file existance logical
    integer(IN)   :: nu      ! unit number

    !--- formats ---
    character(*), parameter :: F00   = "('(drof_comp_init) ',8a)"
    character(*), parameter :: subName = "(drof_comp_init) "
    !-------------------------------------------------------------------------------

    call t_startf('DROF_INIT')

    !----------------------------------------------------------------------------
    ! Initialize pio
    !----------------------------------------------------------------------------

    call shr_strdata_pioinit(SDROF, COMPID)

    !----------------------------------------------------------------------------
    ! Initialize SDROF
    !----------------------------------------------------------------------------

    call t_startf('drof_strdata_init')

    ! NOTE: shr_strdata_init calls shr_dmodel_readgrid which reads the data model
    ! grid and from that computes SDROF%gsmap and SDROF%ggrid. DROF%gsmap is created
    ! using the decomp '2d1d' (1d decomp of 2d grid)

    call shr_strdata_init(SDROF, mpicom, compid, name='rof', calendar=calendar)

    if (my_task == master_task) then
       call shr_strdata_print(SDROF,'SDROF data')
    endif

    call t_stopf('drof_strdata_init')

    !----------------------------------------------------------------------------
    ! Initialize MCT global seg map, 1d decomp
    !----------------------------------------------------------------------------

    call t_startf('drof_initgsmaps')
    if (my_task == master_task) write(logunit,F00) ' initialize gsmaps'

    ! create a data model global seqmap (gsmap) given the data model global grid sizes
    ! NOTE: gsmap is initialized using the decomp read in from the drof_in namelist
    ! (which by default is "1d")
    call shr_dmodel_gsmapcreate(gsmap, SDROF%nxg*SDROF%nyg, compid, mpicom, decomp)
    lsize = mct_gsmap_lsize(gsmap,mpicom)

    ! create a rearranger from the data model SDROF%gsmap to gsmap
    call mct_rearr_init(SDROF%gsmap, gsmap, mpicom, rearr)
    call t_stopf('drof_initgsmaps')

    !----------------------------------------------------------------------------
    ! Initialize MCT domain
    !----------------------------------------------------------------------------

    call t_startf('drof_initmctdom')
    if (my_task == master_task) write(logunit,F00) 'copy domains'
    call shr_dmodel_rearrGGrid(SDROF%grid, ggrid, gsmap, rearr, mpicom)
    call t_stopf('drof_initmctdom')

    !----------------------------------------------------------------------------
    ! Initialize MCT attribute vectors
    !----------------------------------------------------------------------------

    call t_startf('drof_initmctavs')
    if (my_task == master_task) write(logunit,F00) 'allocate AVs'

    call mct_aVect_init(x2r, rList=flds_x2r_mod, lsize=lsize)
    call mct_aVect_zero(x2r)

    call mct_aVect_init(r2x, rList=flds_r2x_mod, lsize=lsize)
    call mct_aVect_zero(r2x)
    call t_stopf('drof_initmctavs')

    !-------------------------------------------------
    ! Determine data model behavior based on the mode
    !-------------------------------------------------

    call t_startf('drof_datamode')
    select case (trim(datamode))

    case('COPYALL')
       ! do nothing extra

    end select
    call t_stopf('drof_datamode')

    !----------------------------------------------------------------------------
    ! Read restart
    !----------------------------------------------------------------------------

    if (read_restart) then
       exists  = .false.
       exists1 = .false.
       if (trim(rest_file)      == trim(nullstr) .and. &
           trim(rest_file_strm) == trim(nullstr)) then
          if (my_task == master_task) then
             write(logunit,F00) ' restart filenames from rpointer = ',trim(rpfile)
             inquire(file=trim(rpfile)//trim(inst_suffix),exist=exists)
             if (exists) then
                nu = shr_file_getUnit()
                open(nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
                read(nu,'(a)') rest_file
                read(nu,'(a)') rest_file_strm
                close(nu)
                call shr_file_freeUnit(nu)
                inquire(file=trim(rest_file_strm),exist=exists)
                inquire(file=trim(rest_file),exist=exists1)
             endif
          endif
          call shr_mpi_bcast(rest_file,mpicom,'rest_file')
          call shr_mpi_bcast(rest_file_strm,mpicom,'rest_file_strm')
       else
          ! use namelist already read
          if (my_task == master_task) then
             write(logunit,F00) ' restart filenames from namelist '
             inquire(file=trim(rest_file_strm),exist=exists)
          endif
       end if

       call shr_mpi_bcast(exists,mpicom,'exists')
       call shr_mpi_bcast(exists1,mpicom,'exists1')

       ! if (exists1) then
       !    if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file)
       !    call shr_pcdf_readwrite('read',SDROF%pio_subsystem, SDROF%io_type, &
       !         trim(rest_file),mpicom,gsmap=gsmap,rf1=water,rf1n='water',io_format=SDROF%io_format)
       ! else
       !    if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file)
       ! endif

       if (exists) then
          if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file_strm)
          call shr_strdata_restRead(trim(rest_file_strm),SDROF,mpicom)
       else
          if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file_strm)
       endif
    end if

    !----------------------------------------------------------------------------
    ! Set initial rof state
    !----------------------------------------------------------------------------

    call t_adj_detailf(+2)

    call drof_comp_run(&
         x2r=x2r, & 
         r2x=r2x, &
         SDROF=SDROF, &
         gsmap=gsmap, &
         ggrid=ggrid, &
         mpicom=mpicom, &
         compid=compid, &
         my_task=my_task, &
         master_task=master_task, &
         inst_suffix=inst_suffix, &
         logunit=logunit, &
         read_restart=read_restart, &
         write_restart=.false., &
         target_ymd=target_ymd, &
         target_tod=target_tod)

    if (my_task == master_task) write(logunit,F00) 'drof_comp_init done'

    call t_adj_detailf(-2)

    call t_stopf('DROF_INIT')

  end subroutine drof_comp_init

  !===============================================================================

  subroutine drof_comp_run(x2r, r2x, &
       SDROF, gsmap, ggrid, mpicom, compid, my_task, master_task, &
       inst_suffix, logunit, read_restart, write_restart, &
       target_ymd, target_tod, case_name)

    ! !DESCRIPTION:  run method for drof model

    ! !INPUT/OUTPUT PARAMETERS:
    type(mct_aVect)        , intent(inout) :: x2r
    type(mct_aVect)        , intent(inout) :: r2x
    type(shr_strdata_type) , intent(inout) :: SDROF
    type(mct_gsMap)        , pointer       :: gsMap
    type(mct_gGrid)        , pointer       :: ggrid
    integer(IN)            , intent(in)    :: mpicom           ! mpi communicator
    integer(IN)            , intent(in)    :: compid           ! mct comp id
    integer(IN)            , intent(in)    :: my_task          ! my task in mpi communicator mpicom
    integer(IN)            , intent(in)    :: master_task      ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix      ! char string associated with instance
    integer(IN)            , intent(in)    :: logunit          ! logging unit number
    logical                , intent(in)    :: read_restart     ! start from restart
    logical                , intent(in)    :: write_restart    ! write restart
    integer(IN)            , intent(in)    :: target_ymd       ! model date
    integer(IN)            , intent(in)    :: target_tod       ! model sec into model date
    character(CL)          , intent(in), optional :: case_name ! case name

    !--- local ---
    integer(IN)   :: n  ! indices
    integer(IN)   :: nf ! fields loop index
    integer(IN)   :: nu ! unit number
    character(len=18) :: date_str

    character(*), parameter :: F00   = "('(drof_comp_run) ',8a)"
    character(*), parameter :: F04   = "('(drof_comp_run) ',2a,2i8,'s')"
    character(*), parameter :: subName = "(drof_comp_run) "
    !-------------------------------------------------------------------------------

    call t_startf('DROF_RUN')

    !--------------------
    ! UNPACK
    !--------------------

    ! do nothing currently

    !--------------------
    ! ADVANCE ROF
    !--------------------

    call t_barrierf('drof_r_BARRIER',mpicom)
    call t_startf('drof_r')

    call t_startf('drof_r_strdata_advance')
    call shr_strdata_advance(SDROF, target_ymd, target_tod, mpicom, 'drof_r')
    call t_stopf('drof_r_strdata_advance')

    !--- copy streams to r2x ---
    call t_barrierf('drof_r_scatter_BARRIER', mpicom)
    call t_startf('drof_r_scatter')
    do n = 1,SDROF%nstreams
       call shr_dmodel_translateAV(SDROF%avs(n), r2x, avifld, avofld, rearr)
    enddo
    call t_stopf('drof_r_scatter')

    ! zero out "special values"
    do nf = 1, mct_avect_nRattr(r2x)
       do n = 1, mct_avect_lsize(r2x)
          if (abs(r2x%rAttr(nf,n)) > 1.0e28) then
             r2x%rAttr(nf,n) = 0.0_r8
          end if
       enddo
    enddo

    call t_stopf('drof_r')

    !-------------------------------------------------
    ! Determine data model behavior based on the mode
    !-------------------------------------------------

    call t_startf('drof_datamode')
    select case (trim(datamode))

    case('COPYALL')
       ! do nothing extra

    end select
    call t_stopf('drof_datamode')

    !--------------------
    ! Write restart
    !--------------------

    if (write_restart) then
       call t_startf('drof_restart')
       call shr_cal_datetod2string(date_str, target_ymd, target_tod)
       write(rest_file,"(6a)") &
            trim(case_name), '.drof',trim(inst_suffix),'.r.', &
            trim(date_str),'.nc'
       write(rest_file_strm,"(6a)") &
            trim(case_name), '.drof',trim(inst_suffix),'.rs1.', &
            trim(date_str),'.bin'
       if (my_task == master_task) then
          nu = shr_file_getUnit()
          open(nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
          write(nu,'(a)') rest_file
          write(nu,'(a)') rest_file_strm
          close(nu)
          call shr_file_freeUnit(nu)
       endif
       if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file_strm),target_ymd,target_tod
       call shr_strdata_restWrite(trim(rest_file_strm), SDROF, mpicom, trim(case_name), 'SDROF strdata')
       call t_stopf('drof_restart')
    end if

    !----------------------------------------------------------------------------
    ! Log output for model date
    !----------------------------------------------------------------------------

    if (my_task == master_task) then
       write(logunit,*) ' drof: model date ', target_ymd,target_tod
    end if

    call t_stopf('DROF_RUN')

  end subroutine drof_comp_run

end module drof_comp_mod
