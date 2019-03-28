#ifdef AIX
@PROCESS ALIAS_SIZE(805306368)
#endif
module drof_comp_mod

  ! !USES:
  use NUOPC                 , only : NUOPC_Advertise
  use ESMF                  , only : ESMF_State, ESMF_SUCCESS, ESMF_State
  use ESMF                  , only : ESMF_Mesh, ESMF_DistGrid, ESMF_MeshGet, ESMF_DistGridGet
  use perf_mod              , only : t_startf, t_stopf
  use perf_mod              , only : t_adj_detailf, t_barrierf
  use mct_mod               , only : mct_gsmap_init
  use mct_mod               , only : mct_avect, mct_avect_indexRA, mct_avect_zero, mct_aVect_nRattr
  use mct_mod               , only : mct_avect_init, mct_avect_lsize
  use shr_sys_mod           , only : shr_sys_abort
  use med_constants_mod     , only : R8, CS, CL, CXX
  use shr_string_mod        , only : shr_string_listGetName
  use shr_sys_mod           , only : shr_sys_abort
  use shr_file_mod          , only : shr_file_getunit, shr_file_freeunit
  use shr_mpi_mod           , only : shr_mpi_bcast
  use shr_strdata_mod       , only : shr_strdata_init_model_domain
  use shr_strdata_mod       , only : shr_strdata_init_streams
  use shr_strdata_mod       , only : shr_strdata_init_mapping
  use shr_strdata_mod       , only : shr_strdata_type, shr_strdata_pioinit
  use shr_strdata_mod       , only : shr_strdata_print, shr_strdata_restRead
  use shr_strdata_mod       , only : shr_strdata_advance, shr_strdata_restWrite
  use shr_dmodel_mod        , only : shr_dmodel_translateAV
  use shr_cal_mod           , only : shr_cal_calendarname
  use shr_cal_mod           , only : shr_cal_datetod2string
  use shr_nuopc_scalars_mod , only : flds_scalar_name
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
  use dshr_nuopc_mod        , only : fld_list_type
  use dshr_nuopc_mod        , only : dshr_fld_add
  use drof_shr_mod          , only : datamode       ! namelist input
  use drof_shr_mod          , only : rest_file      ! namelist input
  use drof_shr_mod          , only : rest_file_strm ! namelist input
  use drof_shr_mod          , only : nullstr

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
       SDROF, mpicom, compid, my_task, master_task, &
       inst_suffix, logunit, read_restart, &
       target_ymd, target_tod, calendar, mesh)

    ! !DESCRIPTION: initialize drof model

    ! input/output arguments
    type(mct_aVect)        , intent(inout) :: x2r, r2x     ! input/output attribute vectors
    type(shr_strdata_type) , intent(inout) :: SDROF        ! model shr_strdata instance (output)
    integer                , intent(in)    :: mpicom       ! mpi communicator
    integer                , intent(in)    :: compid       ! mct comp id
    integer                , intent(in)    :: my_task      ! my task in mpi communicator mpicom
    integer                , intent(in)    :: master_task  ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix  ! char string associated with instance
    integer                , intent(in)    :: logunit      ! logging unit number
    logical                , intent(in)    :: read_restart ! start from restart
    integer                , intent(in)    :: target_ymd   ! model date
    integer                , intent(in)    :: target_tod   ! model sec into model date
    character(len=*)       , intent(in)    :: calendar     ! model calendar
    type(ESMF_Mesh)        , intent(in)    :: mesh         ! ESMF docn mesh

    ! local variables
    integer                      :: n,k     ! generic counters
    integer                      :: ierr    ! error code
    integer                      :: lsize   ! local size
    logical                      :: exists  ! file existance logical
    logical                      :: exists1 ! file existance logical
    integer                      :: nu      ! unit number
    type(ESMF_DistGrid)          :: distGrid
    integer, allocatable, target :: gindex(:)
    integer                      :: rc
    logical                      :: write_restart
    integer                      :: dimCount
    integer                      :: tileCount
    integer                      :: deCount
    integer                      :: gsize
    integer, allocatable         :: elementCountPTile(:)
    integer, allocatable         :: indexCountPDE(:,:)
    integer                      :: spatialDim
    integer                      :: numOwnedElements
    real(R8), pointer            :: ownedElemCoords(:)
    integer                      :: index_lat, index_lon
    real(R8), pointer            :: xc(:), yc(:)       ! arrays of model latitudes and longitudes
    character(*), parameter      :: F00   = "('(drof_comp_init) ',8a)"
    character(*), parameter      :: subName = "(drof_comp_init) "
    !-------------------------------------------------------------------------------

    call t_startf('DROF_INIT')

    !----------------------------------------------------------------------------
    ! Initialize pio
    !----------------------------------------------------------------------------

    call shr_strdata_pioinit(SDROF, COMPID)

    !----------------------------------------------------------------------------
    ! Create a data model global segmap
    !----------------------------------------------------------------------------

    call t_startf('drof_strdata_init')

    if (my_task == master_task) write(logunit,F00) ' initialize SDROF gsmap'

    ! obtain the distgrid from the mesh that was read in
    call ESMF_MeshGet(Mesh, elementdistGrid=distGrid, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determin local size on my processor
    call ESMF_distGridGet(distGrid, localDe=0, elementCount=lsize, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine global index space for my processor
    allocate(gindex(lsize))
    call ESMF_distGridGet(distGrid, localDe=0, seqIndexList=gindex, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine global size of distgrid
    call ESMF_distGridGet(distGrid, dimCount=dimCount, deCount=deCount, tileCount=tileCount, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    allocate(elementCountPTile(tileCount))
    call ESMF_distGridGet(distGrid, elementCountPTile=elementCountPTile, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    gsize = 0
    do n = 1,size(elementCountPTile)
       gsize = gsize + elementCountPTile(n)
    end do
    deallocate(elementCountPTile)

    ! create the data model gsmap given the local size, global size and gindex
    call mct_gsMap_init( SDROF%gsmap, gindex, mpicom, compid, lsize, gsize)
    deallocate(gindex)

    !----------------------------------------------------------------------------
    ! Initialize SDROF
    !----------------------------------------------------------------------------

    ! The call to shr_strdata_init_model_domain creates the SDROF%gsmap which
    ! is a '2d1d' decommp (1d decomp of 2d grid) and also create SDROF%grid

    SDROF%calendar = trim(shr_cal_calendarName(trim(calendar)))

    call shr_strdata_init_model_domain(SDROF, mpicom, compid, my_task, gsmap=SDROF%gsmap)

    if (my_task == master_task) then
       call shr_strdata_print(SDROF,'SDROF data')
    endif

    ! obtain mesh lats and lons
    call ESMF_MeshGet(mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numOwnedElements))
    allocate(xc(numOwnedElements), yc(numOwnedElements))
    call ESMF_MeshGet(mesh, ownedElemCoords=ownedElemCoords)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (numOwnedElements /= lsize) then
       call shr_sys_abort('ERROR: numOwnedElements is not equal to lsize')
    end if
    do n = 1,lsize
       xc(n) = ownedElemCoords(2*n-1)
       yc(n) = ownedElemCoords(2*n)
    end do

    ! error check that mesh lats and lons correspond to those on the input domain file
    index_lon = mct_aVect_indexRA(SDROF%grid%data,'lon')
    do n = 1, lsize
       if (abs( SDROF%grid%data%rattr(index_lon,n) - xc(n)) > 1.e-10) then
          write(6,*)'ERROR: lon diff = ',abs(SDROF%grid%data%rattr(index_lon,n) -  xc(n)),' too large'
          call shr_sys_abort()
       end if
       !SDROF%grid%data%rattr(index_lon,n) = xc(n) ! overwrite ggrid with mesh data
    end do
    index_lat = mct_aVect_indexRA(SDROF%grid%data,'lat')
    do n = 1, lsize
       if (abs( SDROF%grid%data%rattr(index_lat,n) -  yc(n)) > 1.e-10) then
          write(6,*)'ERROR: lat diff = ',abs(SDROF%grid%data%rattr(index_lat,n) -  yc(n)),' too large'
          call shr_sys_abort()
       end if
       !SDROF%grid%data%rattr(index_lat,n) = yc(n) ! overwrite ggrid with mesh data
    end do

    !----------------------------------------------------------------------------
    ! Initialize SDROF attributes for streams and mapping of streams to model domain
    !----------------------------------------------------------------------------

    call shr_strdata_init_streams(SDROF, compid, mpicom, my_task)
    call shr_strdata_init_mapping(SDROF, compid, mpicom, my_task)

    call t_stopf('drof_strdata_init')

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
       !         trim(rest_file),mpicom,gsmap=SDROF%gsmap,rf1=water,rf1n='water',io_format=SDROF%io_format)
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

    write_restart=.false.
    call drof_comp_run(x2r, r2x, &
         SDROF, mpicom, compid, my_task, master_task, &
         inst_suffix, logunit, read_restart, write_restart, &
         target_ymd, target_tod)

    if (my_task == master_task) write(logunit,F00) 'drof_comp_init done'

    call t_adj_detailf(-2)

    call t_stopf('DROF_INIT')

  end subroutine drof_comp_init

  !===============================================================================

  subroutine drof_comp_run(x2r, r2x, &
       SDROF, mpicom, compid, my_task, master_task, &
       inst_suffix, logunit, read_restart, write_restart, &
       target_ymd, target_tod, case_name)

    ! !DESCRIPTION:  run method for drof model

    ! input/output arguments
    type(mct_aVect)        , intent(inout) :: x2r
    type(mct_aVect)        , intent(inout) :: r2x
    type(shr_strdata_type) , intent(inout) :: SDROF
    integer                , intent(in)    :: mpicom           ! mpi communicator
    integer                , intent(in)    :: compid           ! mct comp id
    integer                , intent(in)    :: my_task          ! my task in mpi communicator mpicom
    integer                , intent(in)    :: master_task      ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix      ! char string associated with instance
    integer                , intent(in)    :: logunit          ! logging unit number
    logical                , intent(in)    :: read_restart     ! start from restart
    logical                , intent(in)    :: write_restart    ! write restart
    integer                , intent(in)    :: target_ymd       ! model date
    integer                , intent(in)    :: target_tod       ! model sec into model date
    character(len=*)       , intent(in), optional :: case_name ! case name

    ! local
    integer                 :: n  ! indices
    integer                 :: nf ! fields loop index
    integer                 :: nu ! unit number
    character(len=18)       :: date_str
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
       call shr_dmodel_translateAV(SDROF%avs(n), r2x, avifld, avofld)
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

    call t_stopf('DROF_RUN')

  end subroutine drof_comp_run

end module drof_comp_mod
