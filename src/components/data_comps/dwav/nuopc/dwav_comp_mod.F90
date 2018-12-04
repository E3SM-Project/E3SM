#ifdef AIX
@PROCESS ALIAS_SIZE(805306368)
#endif
module dwav_comp_mod

  ! !USES:

  use NUOPC                 , only : NUOPC_Advertise
  use ESMF                  , only : ESMF_State, ESMF_SUCCESS, ESMF_STATE
  use ESMF                  , only : ESMF_Mesh, ESMF_DistGrid, ESMF_MeshGet, ESMF_DistGridGet
  use perf_mod              , only : t_startf, t_stopf, t_adj_detailf, t_barrierf
  use mct_mod               , only : mct_gsmap_init
  use mct_mod               , only : mct_avect, mct_avect_indexRA, mct_avect_zero, mct_aVect_nRattr
  use mct_mod               , only : mct_avect_init, mct_avect_lsize
  use shr_sys_mod           , only : shr_sys_abort
  use shr_kind_mod          , only : IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_kind_mod          , only : CXX=>SHR_KIND_CXX
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
  use dwav_shr_mod          , only : datamode       ! namelist input
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

  character(len=CS), pointer  :: avifld(:) ! character array for field names coming from streams
  character(len=CS), pointer  :: avofld(:) ! character array for field names to be sent/received from mediator
  character(len=CXX)          :: flds_w2x_mod
  character(len=CXX)          :: flds_x2w_mod
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
       SDWAV, mpicom, compid, my_task, master_task, &
       inst_suffix, logunit, read_restart, &
       target_ymd, target_tod, calendar, mesh)

    ! !DESCRIPTION: initialize dwav model

    ! !INPUT/OUTPUT PARAMETERS:
    type(mct_aVect)        , intent(inout) :: x2w, w2x     ! input/output attribute vectors
    type(shr_strdata_type) , intent(inout) :: SDWAV        ! model
    integer                , intent(in)    :: mpicom       ! mpi communicator
    integer                , intent(in)    :: compid       ! mct comp id
    integer                , intent(in)    :: my_task      ! my task in mpi communicator mpicom
    integer                , intent(in)    :: master_task  ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix  ! char string associated with instance
    integer                , intent(in)    :: logunit      ! logging unit number
    logical                , intent(in)    :: read_restart ! start from restart
    integer                , intent(in)    :: target_ymd   ! model date
    integer                , intent(in)    :: target_tod   ! model sec into model date
    character(len=*)       , intent(in)    :: calendar     ! calendar type
    type(ESMF_Mesh)        , intent(in)    :: mesh         ! ESMF docn mesh

    !--- local variables ---
    integer                      :: n,k       ! generic counters
    integer                      :: lsize     ! local size
    logical                      :: exists    ! file existance
    integer                      :: nu        ! unit number
    logical                      :: write_restart
    type(ESMF_DistGrid)          :: distGrid
    integer, allocatable, target :: gindex(:)
    integer                      :: dimCount
    integer                      :: tileCount
    integer                      :: deCount
    integer                      :: gsize
    integer, allocatable         :: elementCountPTile(:)
    integer, allocatable         :: indexCountPDE(:,:)
    integer                      :: spatialDim
    integer                      :: numOwnedElements
    real(R8), pointer            :: ownedElemCoords(:)
    integer                      :: klat, klon, kfrac   ! AV indices
    real(R8), pointer            :: domlon(:),domlat(:) ! ggrid domain lats and lots
    real(R8), pointer            :: xc(:), yc(:)        ! mesh lats and lons
    real(r8)                     :: maxerr, err
    integer                      :: maxn
    real(r8)                     :: tolerance = 1.e-4
    integer                      :: rc
    character(*), parameter      :: F00   = "('(dwav_comp_init) ',8a)"
    character(*), parameter      :: subName = "(dwav_comp_init) "
    !-------------------------------------------------------------------------------

    call t_startf('DWAV_INIT')

    !----------------------------------------------------------------------------
    ! Initialize pio
    !----------------------------------------------------------------------------

    call shr_strdata_pioinit(SDWAV, compid)

    !----------------------------------------------------------------------------
    ! Create a data model global segmap
    !----------------------------------------------------------------------------

    call t_startf('dwav_strdata_init')

    if (my_task == master_task) write(logunit,F00) ' initialize SDWAV gsmap'

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
    call mct_gsMap_init( SDWAV%gsmap, gindex, mpicom, compid, lsize, gsize)
    deallocate(gindex)

    !----------------------------------------------------------------------------
    ! Initialize SDWAV
    !----------------------------------------------------------------------------

    ! The call to shr_strdata_init_model_domain creates the SDWAV%gsmap which
    ! is a '2d1d' decommp (1d decomp of 2d grid) and also create SDWAV%grid

    SDWAV%calendar = trim(shr_cal_calendarName(trim(calendar)))

    call shr_strdata_init_model_domain(SDWAV, mpicom, compid, my_task, gsmap=SDWAV%gsmap)

    if (my_task == master_task) then
       call shr_strdata_print(SDWAV,'SDWAV data')
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
    allocate(domlon(lsize), domlat(lsize))

    maxerr = 0._r8
    maxn = 0
    klon = mct_aVect_indexRA(SDWAV%grid%data,'lon')
    do n = 1, lsize
       domlon(n) = SDWAV%grid%data%rattr(klon,n)
       err = abs(domlon(n) -  xc(n))
       if (err > maxerr) then
          maxerr = err
          maxn = n
       end if
       !SDWAV%grid%data%rattr(klon,n) = xc(n)
    end do
    if (maxerr > 0._r8) then
       write(6,100) maxn, domlon(maxn), xc(maxn), maxerr
100    format('WARNING: DWAV n, dom_lon, mesh_lon, max_diff_lon = ',i6,2(f21.13,3x),d21.5)
    end if
    if (maxerr > tolerance) then
       write(6,*)'ERROR: diff of dom_lon and mesh_lon is greater than tolerance of ',tolerance
       call shr_sys_abort()
    end if

    maxerr = 0._r8
    maxn = 0
    klat = mct_aVect_indexRA(SDWAV%grid%data,'lat')
    do n = 1, lsize
       domlat(n) = SDWAV%grid%data%rattr(klat,n)
       err = abs(domlat(n) -  yc(n))
       if ( err > maxerr ) then
          maxerr = err
          maxn = n
       end if
       !SDWAV%grid%data%rattr(klat,n) = yc(n)
    end do
    if (maxerr > 0._r8) then
       write(6,101) maxn, domlat(maxn), yc(maxn), maxerr
101    format('WARNING: DWAV n, dom_lat, mesh_lat, max_diff_lat = ',i6,2(f21.13,3x),d21.5)
    end if
    if (maxerr > tolerance) then
       write(6,*)'ERROR: diff of dom_lat and mesh_lat is greater than tolerance of ',tolerance
       call shr_sys_abort()
    end if

    deallocate(domlon, domlat)

    !----------------------------------------------------------------------------
    ! Initialize SDLND attributes for streams and mapping of streams to model domain
    !----------------------------------------------------------------------------

    call shr_strdata_init_streams(SDWAV, compid, mpicom, my_task)
    call shr_strdata_init_mapping(SDWAV, compid, mpicom, my_task)

    call t_stopf('dwav_strdata_init')

    !----------------------------------------------------------------------------
    ! Initialize MCT attribute vectors
    !----------------------------------------------------------------------------

    if (my_task == master_task) write(logunit,F00) 'allocate AVs'

    call mct_avect_init(w2x, rlist=flds_w2x_mod, lsize=lsize)
    call mct_avect_zero(w2x)
    call mct_avect_init(x2w, rlist=flds_x2w_mod, lsize=lsize)
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

    !----------------------------------------------------------------------------
    ! Set initial wav state
    !----------------------------------------------------------------------------

    write_restart = .false.
    call dwav_comp_run(x2w, w2x, &
         SDWAV, mpicom, my_task, master_task, &
         inst_suffix, logunit, read_restart, write_restart, &
         target_ymd, target_tod)

    if (my_task == master_task) then
       write(logunit,F00) 'dwav_comp_init done'
    end if

    call t_stopf('DWAV_INIT')

  end subroutine dwav_comp_init

  !===============================================================================

  subroutine dwav_comp_run(x2w, w2x, &
       SDWAV, mpicom, my_task, master_task, &
       inst_suffix, logunit, read_restart, write_restart, &
       target_ymd, target_tod, case_name)

    ! DESCRIPTION:  run method for dwav model

    ! input/output parameters:
    type(mct_aVect)        , intent(inout) :: x2w
    type(mct_aVect)        , intent(inout) :: w2x
    type(shr_strdata_type) , intent(inout) :: SDWAV
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
    integer                 :: n                     ! indices
    integer                 :: idt                   ! integer timestep
    integer                 :: nu                    ! unit number
    character(len=18)       :: date_str
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
       call shr_dmodel_translateAV(SDWAV%avs(n), w2x, avifld, avofld)
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
       call shr_cal_datetod2string(date_str,  target_ymd,  target_tod)
       write(rest_file, "(6a)") &
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
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------
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
