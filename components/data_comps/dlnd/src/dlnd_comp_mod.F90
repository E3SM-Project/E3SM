#ifdef AIX
@PROCESS ALIAS_SIZE(805306368)
#endif
module dlnd_comp_mod

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
  use shr_cal_mod       , only: shr_cal_ymdtod2string
  use seq_timemgr_mod   , only: seq_timemgr_EClockGetData, seq_timemgr_RestartAlarmIsOn
  use glc_elevclass_mod , only: glc_get_num_elevation_classes, glc_elevclass_as_string

  use dlnd_shr_mod   , only: datamode       ! namelist input
  use dlnd_shr_mod   , only: decomp         ! namelist input
  use dlnd_shr_mod   , only: rest_file      ! namelist input
  use dlnd_shr_mod   , only: rest_file_strm ! namelist input
  use dlnd_shr_mod   , only: domain_fracname ! namelist input
  use dlnd_shr_mod   , only: nullstr

#ifdef HAVE_MOAB
  use seq_comm_mct,     only : mlnid ! id of moab lnd app
  use iso_c_binding
#endif
  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: dlnd_comp_init
  public :: dlnd_comp_run
  public :: dlnd_comp_final

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  !--- other ---
  character(CS) :: myModelName = 'lnd'   ! user defined model name
  character(len=*),parameter :: rpfile = 'rpointer.lnd'
  type(mct_rearr) :: rearr

  !--------------------------------------------------------------------------
  !--- names of fields ---
  integer(IN),parameter :: fld_len = 12       ! max character length of fields in avofld & avifld
  integer(IN),parameter :: nflds_nosnow = 29

  ! fields other than snow fields:
  character(fld_len),parameter  :: avofld_nosnow(1:nflds_nosnow) = &
       (/ "Sl_t        ","Sl_tref     ","Sl_qref     ","Sl_avsdr    ","Sl_anidr    ", &
       "Sl_avsdf    ","Sl_anidf    ","Sl_snowh    ","Fall_taux   ","Fall_tauy   ", &
       "Fall_lat    ","Fall_sen    ","Fall_lwup   ","Fall_evap   ","Fall_swnet  ", &
       "Sl_landfrac ","Sl_fv       ","Sl_ram1     ","Flrl_demand ",                &
        "Flrl_rofsur ","Flrl_rofgwl ","Flrl_rofsub ","Flrl_rofdto ","Flrl_rofi   ", &
       "Fall_flxdst1","Fall_flxdst2","Fall_flxdst3","Fall_flxdst4", "Flrl_rofmud "/)

  character(fld_len),parameter  :: avifld_nosnow(1:nflds_nosnow) = &
       (/ "t           ","tref        ","qref        ","avsdr       ","anidr       ", &
       "avsdf       ","anidf       ","snowh       ","taux        ","tauy        ", &
       "lat         ","sen         ","lwup        ","evap        ","swnet       ", &
       "lfrac       ","fv          ","ram1        ","demand      ",                &
        "rofsur      ","rofgwl      ","rofsub      ","rofdto      ","rofi        ", &
       "flddst1     ","flxdst2     ","flxdst3     ","flxdst4     " ,"rofmud      "/)

  integer(IN), parameter :: nflds_snow  = 3   ! number of snow fields in each elevation class
  integer(IN), parameter :: nec_len    = 2   ! length of elevation class index in field names
  ! for these snow fields, the actual field names will have the elevation class index at
  ! the end (e.g., Sl_tsrf01, tsrf01)

  character(fld_len-nec_len),parameter :: avofld_snow(nflds_snow) = &
       (/"Sl_tsrf  ", "Sl_topo  ", "Flgl_qice"/)

  character(fld_len-nec_len),parameter :: avifld_snow(nflds_snow) = &
       (/"tsrf", "topo", "qice"/)

  ! all fields:
  character(fld_len),dimension(:),allocatable :: avofld
  character(fld_len),dimension(:),allocatable :: avifld
  !--------------------------------------------------------------------------

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !===============================================================================
  subroutine dlnd_comp_init(Eclock, x2l, l2x, &
       seq_flds_x2l_fields, seq_flds_l2x_fields, &
       SDLND, gsmap, ggrid, mpicom, compid, my_task, master_task, &
       inst_suffix, inst_name, logunit, read_restart, &
       scmMode, scmlat, scmlon)

    ! !DESCRIPTION: initialize dlnd model
#ifdef HAVE_MOAB
       use iMOAB, only: iMOAB_DefineTagStorage, iMOAB_GetDoubleTagStorage, &
                        iMOAB_SetIntTagStorage, iMOAB_SetDoubleTagStorage, &
                        iMOAB_ResolveSharedEntities, iMOAB_CreateVertices, &
                        iMOAB_GetMeshInfo, iMOAB_UpdateMeshInfo, iMOAB_WriteMesh
#endif
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)       , intent(in)    :: EClock
    type(mct_aVect)        , intent(inout) :: x2l, l2x            ! input/output attribute vectors
    character(len=*)       , intent(in)    :: seq_flds_x2l_fields ! fields from mediator
    character(len=*)       , intent(in)    :: seq_flds_l2x_fields ! fields to mediator
    type(shr_strdata_type) , intent(inout) :: SDLND               ! model shr_strdata instance (output)
    type(mct_gsMap)        , pointer       :: gsMap               ! model global seg map (output)
    type(mct_gGrid)        , pointer       :: ggrid               ! model ggrid (output)
    integer(IN)            , intent(in)    :: mpicom              ! mpi communicator
    integer(IN)            , intent(in)    :: compid              ! mct comp id
    integer(IN)            , intent(in)    :: my_task             ! my task in mpi communicator mpicom
    integer(IN)            , intent(in)    :: master_task         ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix         ! char string associated with instance
    character(len=*)       , intent(in)    :: inst_name           ! fullname of current instance (ie. "lnd_0001")
    integer(IN)            , intent(in)    :: logunit             ! logging unit number
    logical                , intent(in)    :: read_restart        ! start from restart
    logical                , intent(in)    :: scmMode             ! single column mode
    real(R8)               , intent(in)    :: scmLat              ! single column lat
    real(R8)               , intent(in)    :: scmLon              ! single column lon

    !--- local variables ---
    integer(IN)        :: n,k       ! generic counters
    integer(IN)        :: ierr      ! error code
    integer(IN)        :: lsize     ! local size
    logical            :: exists    ! file existance
    integer(IN)        :: nu        ! unit number
    character(CL)      :: calendar  ! model calendar
    integer(IN)        :: glc_nec   ! number of elevation classes
    integer(IN)        :: nflds_glc_nec  ! number of snow fields separated by elev class
    integer(IN)        :: field_num ! field number
    character(nec_len) :: nec_str   ! elevation class, as character string
    character(*), parameter :: domain_fracname_unset = 'null'

#ifdef HAVE_MOAB
    character*400  tagname
    real(R8) latv, lonv
    integer iv, tagindex, ilat, ilon 
    real(R8), allocatable, target :: data(:)
    integer(IN), pointer :: idata(:)   ! temporary
    real(R8), dimension(:), allocatable :: moab_vert_coords  ! temporary
#ifdef MOABDEBUG
    character*100 outfile, wopts
#endif
#endif

    !--- formats ---
    character(*), parameter :: F00   = "('(dlnd_comp_init) ',8a)"
    character(*), parameter :: F0L   = "('(dlnd_comp_init) ',a, l2)"
    character(*), parameter :: F01   = "('(dlnd_comp_init) ',a,5i8)"
    character(*), parameter :: F02   = "('(dlnd_comp_init) ',a,4es13.6)"
    character(*), parameter :: F03   = "('(dlnd_comp_init) ',a,i8,a)"
    character(*), parameter :: F05   = "('(dlnd_comp_init) ',a,2f10.4)"
    character(*), parameter :: F90   = "('(dlnd_comp_init) ',73('='))"
    character(*), parameter :: F91   = "('(dlnd_comp_init) ',73('-'))"
    character(*), parameter :: subName = "(dlnd_comp_init) "
    !-------------------------------------------------------------------------------

    call t_startf('DLND_INIT')

    !----------------------------------------------------------------------------
    ! Initialize pio
    !----------------------------------------------------------------------------

    call shr_strdata_pioinit(SDLND, COMPID)

    !----------------------------------------------------------------------------
    ! Initialize SDLND
    !----------------------------------------------------------------------------

    call t_startf('dlnd_strdata_init')

    call seq_timemgr_EClockGetData( EClock, calendar=calendar )

    ! NOTE: shr_strdata_init calls shr_dmodel_readgrid which reads the data model
    ! grid and from that computes SDLND%gsmap and SDLND%ggrid. DLND%gsmap is created
    ! using the decomp '2d1d' (1d decomp of 2d grid)

    if (scmmode) then
       if (my_task == master_task) &
            write(logunit,F05) ' scm lon lat = ',scmlon,scmlat
       if (domain_fracname == domain_fracname_unset) then
          call shr_strdata_init(SDLND,mpicom,compid,name='lnd', &
               scmmode=scmmode,scmlon=scmlon,scmlat=scmlat, calendar=calendar)
       else
          call shr_strdata_init(SDLND,mpicom,compid,name='lnd', &
               scmmode=scmmode,scmlon=scmlon,scmlat=scmlat, calendar=calendar, &
               dmodel_domain_fracname_from_stream=domain_fracname)
       end if
    else
       if (domain_fracname == domain_fracname_unset) then
          call shr_strdata_init(SDLND,mpicom,compid,name='lnd', calendar=calendar)
       else
          call shr_strdata_init(SDLND,mpicom,compid,name='lnd', calendar=calendar, &
               dmodel_domain_fracname_from_stream=domain_fracname)
       end if
    endif

    if (my_task == master_task) then
       call shr_strdata_print(SDLND,'SDLND data')
    endif

    call t_stopf('dlnd_strdata_init')

    !----------------------------------------------------------------------------
    ! Build avofld & avifld
    !----------------------------------------------------------------------------

    glc_nec = glc_get_num_elevation_classes()
    if (glc_nec == 0) then
       nflds_glc_nec = 0
    else
       nflds_glc_nec = (glc_nec+1)*nflds_snow
    end if

    ! Start with non-snow fields
    allocate(avofld(nflds_nosnow + nflds_glc_nec))
    allocate(avifld(nflds_nosnow + nflds_glc_nec))
    avofld(1:nflds_nosnow) = avofld_nosnow
    avifld(1:nflds_nosnow) = avifld_nosnow
    field_num = nflds_nosnow

    ! Append each snow field
    if (glc_nec > 0) then
       do k = 1, nflds_snow
          do n = 0, glc_nec
             ! nec_str will be something like '02' or '10'
             nec_str = glc_elevclass_as_string(n)

             field_num = field_num + 1
             avofld(field_num) = trim(avofld_snow(k))//nec_str
             avifld(field_num) = trim(avifld_snow(k))//nec_str
          end do
       end do
    end if

    !----------------------------------------------------------------------------
    ! Initialize MCT global seg map, 1d decomp
    !----------------------------------------------------------------------------

    call t_startf('dlnd_initgsmaps')
    if (my_task == master_task) write(logunit,F00) ' initialize gsmaps'
    call shr_sys_flush(logunit)

    ! create a data model global seqmap (gsmap) given the data model global grid sizes
    ! NOTE: gsmap is initialized using the decomp read in from the dlnd_in namelist
    ! (which by default is "1d")
    call shr_dmodel_gsmapcreate(gsmap,SDLND%nxg*SDLND%nyg,compid,mpicom,decomp)
    lsize = mct_gsmap_lsize(gsmap,mpicom)

    ! create a rearranger from the data model DLND%gsmap to gsmap
    call mct_rearr_init(SDLND%gsmap, gsmap, mpicom, rearr)

    call t_stopf('dlnd_initgsmaps')

    !----------------------------------------------------------------------------
    ! Initialize MCT domain
    !----------------------------------------------------------------------------

    call t_startf('dlnd_initmctdom')
    if (my_task == master_task) write(logunit,F00) 'copy domains'
    call shr_sys_flush(logunit)

    call shr_dmodel_rearrGGrid(SDLND%grid, ggrid, gsmap, rearr, mpicom)

    call t_stopf('dlnd_initmctdom')

#ifdef HAVE_MOAB
    ilat = mct_aVect_indexRA(ggrid%data,'lat')
    ilon = mct_aVect_indexRA(ggrid%data,'lon')
    allocate(moab_vert_coords(lsize*3))
    do iv = 1, lsize
       lonv = ggrid%data%rAttr(ilon, iv) * SHR_CONST_PI/180.
       latv = ggrid%data%rAttr(ilat, iv) * SHR_CONST_PI/180.
       moab_vert_coords(3*iv-2)=COS(latv)*COS(lonv)
       moab_vert_coords(3*iv-1)=COS(latv)*SIN(lonv)
       moab_vert_coords(3*iv  )=SIN(latv)
    enddo
 
    ! create the vertices with coordinates from MCT domain
    ierr = iMOAB_CreateVertices(mlnid, lsize*3, 3, moab_vert_coords)
    if (ierr .ne. 0)  &
       call shr_sys_abort('Error: fail to create MOAB vertices in data lnd model')
 
    tagname='GLOBAL_ID'//C_NULL_CHAR
    ierr = iMOAB_DefineTagStorage(mlnid, tagname, &
                                  0, & ! dense, integer
                                  1, & ! number of components
                                  tagindex )
    if (ierr .ne. 0)  &
       call shr_sys_abort('Error: fail to retrieve GLOBAL_ID tag ')
 
    ! get list of global IDs for Dofs
    call mct_gsMap_orderedPoints(gsMap, my_task, idata)
 
    ierr = iMOAB_SetIntTagStorage ( mlnid, tagname, lsize, &
                                     0, & ! vertex type
                                     idata)
    if (ierr .ne. 0)  &
       call shr_sys_abort('Error: fail to set GLOBAL_ID tag ')
 
    ierr = iMOAB_ResolveSharedEntities( mlnid, lsize, idata );
    if (ierr .ne. 0)  &
       call shr_sys_abort('Error: fail to resolve shared entities')
 
    deallocate(moab_vert_coords)
    deallocate(idata)
 
    ierr = iMOAB_UpdateMeshInfo( mlnid )
    if (ierr .ne. 0)  &
       call shr_sys_abort('Error: fail to update mesh info ')
 
    allocate(data(lsize))
    ierr = iMOAB_DefineTagStorage( mlnid, "area:aream:frac:mask"//C_NULL_CHAR, &
                                     1, & ! dense, double
                                     1, & ! number of components
                                     tagindex )
    if (ierr > 0 )  &
       call shr_sys_abort('Error: fail to create tag: area:aream:frac:mask' )
 
    data(:) = ggrid%data%rAttr(mct_aVect_indexRA(ggrid%data,'area'),:)
    tagname='area'//C_NULL_CHAR
    ierr = iMOAB_SetDoubleTagStorage ( mlnid, tagname, lsize, &
                                     0, & ! set data on vertices
                                     data)
    if (ierr > 0 )  &
       call shr_sys_abort('Error: fail to get area tag ')
 
    ! set the same data for aream (model area) as area
    ! data(:) = ggrid%data%rAttr(mct_aVect_indexRA(ggrid%data,'aream'),:)
    tagname='aream'//C_NULL_CHAR
    ierr = iMOAB_SetDoubleTagStorage ( mlnid, tagname, lsize, &
                                     0, & ! set data on vertices
                                     data)
    if (ierr > 0 )  &
       call shr_sys_abort('Error: fail to set aream tag ')
 
    data(:) = ggrid%data%rAttr(mct_aVect_indexRA(ggrid%data,'mask'),:)
    tagname='mask'//C_NULL_CHAR
    ierr = iMOAB_SetDoubleTagStorage ( mlnid, tagname, lsize, &
                                     0, & ! set data on vertices
                                     data)
    if (ierr > 0 )  &
       call shr_sys_abort('Error: fail to set mask tag ')
 
    data(:) = ggrid%data%rAttr(mct_aVect_indexRA(ggrid%data,'frac'),:)
    tagname='frac'//C_NULL_CHAR
    ierr = iMOAB_SetDoubleTagStorage ( mlnid, tagname, lsize, &
                                     0, & ! set data on vertices
                                     data)
    if (ierr > 0 )  &
       call shr_sys_abort('Error: fail to set frac tag ')
 
    deallocate(data)
 
    ! define tags
    ierr = iMOAB_DefineTagStorage( mlnid, trim(seq_flds_x2l_fields)//C_NULL_CHAR, &
                                     1, & ! dense, double
                                     1, & ! number of components
                                     tagindex )
    if (ierr > 0 )  &
       call shr_sys_abort('Error: fail to create seq_flds_x2l_fields tags ')
 
    ierr = iMOAB_DefineTagStorage( mlnid, trim(seq_flds_l2x_fields)//C_NULL_CHAR, &
                                     1, & ! dense, double
                                     1, & ! number of components
                                     tagindex )
    if (ierr > 0 )  &
       call shr_sys_abort('Error: fail to create seq_flds_l2x_fields tags ')
#ifdef MOABDEBUG
       !      debug test
    outfile = 'LndDataMesh.h5m'//C_NULL_CHAR
    wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR !
       !      write out the mesh file to disk
    ierr = iMOAB_WriteMesh(mlnid, trim(outfile), trim(wopts))
    if (ierr .ne. 0) then
       call shr_sys_abort(subname//' ERROR in writing data mesh lnd ')
    endif
#endif
#endif
    !----------------------------------------------------------------------------
    ! Initialize MCT attribute vectors
    !----------------------------------------------------------------------------

    call t_startf('dlnd_initmctavs')
    if (my_task == master_task) write(logunit,F00) 'allocate AVs'
    call shr_sys_flush(logunit)

    call mct_aVect_init(l2x, rList=seq_flds_l2x_fields, lsize=lsize)
    call mct_aVect_zero(l2x)

    call mct_aVect_init(x2l, rList=seq_flds_x2l_fields, lsize=lsize)
    call mct_aVect_zero(x2l)

    call t_stopf('dlnd_initmctavs')

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
       !if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file)
       !call shr_pcdf_readwrite('read',trim(rest_file),mpicom,gsmap,rf1=somtp,rf1n='somtp')
       if (exists) then
          if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file_strm)
          call shr_strdata_restRead(trim(rest_file_strm),SDLND,mpicom)
       else
          if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file_strm)
       endif
       call shr_sys_flush(logunit)
    endif

    !----------------------------------------------------------------------------
    ! Set initial lnd state
    !----------------------------------------------------------------------------

    call t_adj_detailf(+2)
    call dlnd_comp_run(EClock, x2l, l2x, &
         SDLND, gsmap, ggrid, mpicom, compid, my_task, master_task, &
         inst_suffix, logunit)
    call t_adj_detailf(-2)

    if (my_task == master_task) write(logunit,F00) 'dlnd_comp_init done'
    call shr_sys_flush(logunit)

    call t_stopf('DLND_INIT')

  end subroutine dlnd_comp_init

  !===============================================================================
  subroutine dlnd_comp_run(EClock, x2l, l2x, &
       SDLND, gsmap, ggrid, mpicom, compid, my_task, master_task, &
       inst_suffix, logunit, case_name)

    ! !DESCRIPTION:  run method for dlnd model
#ifdef HAVE_MOAB
#ifdef MOABDEBUG
    use iMOAB, only: iMOAB_WriteMesh
#endif
    use seq_flds_mod    , only: seq_flds_l2x_fields 
    use seq_flds_mod    , only: moab_set_tag_from_av
#endif

    implicit none
    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)       , intent(in)    :: EClock
    type(mct_aVect)        , intent(inout) :: x2l
    type(mct_aVect)        , intent(inout) :: l2x
    type(shr_strdata_type) , intent(inout) :: SDLND
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
#ifdef HAVE_MOAB
    real(R8), allocatable, target :: datam(:)
    type(mct_list) :: temp_list
    integer :: size_list, index_list, lsize
    type(mct_string)    :: mctOStr  !
    character*400  tagname, mct_field
#ifdef MOABDEBUG
    integer  :: cur_dlnd_stepno, ierr
    character*100 outfile, wopts, lnum
#endif
#endif

    character(*), parameter :: F00   = "('(dlnd_comp_run) ',8a)"
    character(*), parameter :: F04   = "('(dlnd_comp_run) ',2a,2i8,'s')"
    character(*), parameter :: subName = "(dlnd_comp_run) "
    !-------------------------------------------------------------------------------

    call t_startf('DLND_RUN')

    call t_startf('dlnd_run1')
    call seq_timemgr_EClockGetData( EClock, curr_ymd=CurrentYMD, curr_tod=CurrentTOD)
    call seq_timemgr_EClockGetData( EClock, curr_yr=yy, curr_mon=mm, curr_day=dd)
    write_restart = seq_timemgr_RestartAlarmIsOn(EClock)
    call t_stopf('dlnd_run1')

    !--------------------
    ! UNPACK
    !--------------------

    call t_startf('dlnd_unpack')
    ! Nothing to be done for now
    call t_stopf('dlnd_unpack')

    !--------------------
    ! ADVANCE LAND
    !--------------------

    call t_barrierf('dlnd_BARRIER',mpicom)
    call t_startf('dlnd')

    call t_startf('dlnd_strdata_advance')
    call shr_strdata_advance(SDLND,currentYMD,currentTOD,mpicom,'dlnd')
    call t_stopf('dlnd_strdata_advance')

    call t_barrierf('dlnd_scatter_BARRIER',mpicom)
    call t_startf('dlnd_scatter')
    do n = 1,SDLND%nstreams
       call shr_dmodel_translateAV(SDLND%avs(n), l2x, avifld, avofld, rearr)
    enddo
    call t_stopf('dlnd_scatter')

    call t_stopf('dlnd')

    !-------------------------------------------------
    ! Determine data model behavior based on the mode
    !-------------------------------------------------

    call t_startf('dlnd_datamode')
    select case (trim(datamode))

    case('COPYALL')
       ! do nothing extra

    end select
    call t_stopf('dlnd_datamode')

    !--------------------
    ! Write restart
    !--------------------

    if (write_restart) then
       call t_startf('dlnd_restart')
       call shr_cal_ymdtod2string(date_str, yy,mm,dd,currentTOD)
       write(rest_file,"(6a)") &
            trim(case_name), '.dlnd',trim(inst_suffix),'.r.', &
            trim(date_str),'.nc'
       write(rest_file_strm,"(6a)") &
            trim(case_name), '.dlnd',trim(inst_suffix),'.rs1.', &
            trim(date_str),'.bin'
       if (my_task == master_task) then
          nu = shr_file_getUnit()
          open(nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
          write(nu,'(a)') rest_file
          write(nu,'(a)') rest_file_strm
          close(nu)
          call shr_file_freeUnit(nu)
       endif
       !if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file),currentYMD,currentTOD
       !call shr_pcdf_readwrite('write',trim(rest_file),mpicom,gsmap,clobber=.true., &
       !   rf1=somtp,rf1n='somtp')
       if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file_strm),currentYMD,currentTOD
       call shr_strdata_restWrite(trim(rest_file_strm),SDLND,mpicom,trim(case_name),'SDLND strdata')
       call shr_sys_flush(logunit)
       call t_stopf('dlnd_restart')
    endif

    !----------------------------------------------------------------------------
    ! Log output for model date
    !----------------------------------------------------------------------------

    call t_startf('dlnd_run2')
    if (my_task == master_task) then
       write(logunit,F04) trim(myModelName),': model date ', CurrentYMD,CurrentTOD
       call shr_sys_flush(logunit)
    end if
    call t_stopf('dlnd_run2')

    call t_stopf('DLND_RUN')

#ifdef HAVE_MOAB
    lsize = mct_avect_lsize(l2x) ! is it the same as mct_avect_lsize(avstrm) ?
    allocate(datam(lsize)) ! 
    call mct_list_init(temp_list ,seq_flds_l2x_fields)
    size_list=mct_list_nitem (temp_list)
    do index_list = 1, size_list
      call mct_list_get(mctOStr,index_list,temp_list)
      mct_field = mct_string_toChar(mctOStr)
      tagname= trim(mct_field)//C_NULL_CHAR
      call moab_set_tag_from_av(tagname, l2x, index_list, mlnid, datam, lsize) ! loop over all a2x fields, not just a few
    enddo
    call mct_list_clean(temp_list)
    deallocate(datam) ! maybe we should keep it around, deallocate at the final only?

#ifdef MOABDEBUG
    call seq_timemgr_EClockGetData( EClock, stepno=cur_dlnd_stepno )
    write(lnum,"(I0.2)")cur_dlnd_stepno
    outfile = 'dlnd_comp_run_'//trim(lnum)//'.h5m'//C_NULL_CHAR
    wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
    ierr = iMOAB_WriteMesh(mlnid, outfile, wopts)
    if (ierr > 0 )  then
       write(logunit,*) 'Failed to write data lnd component state '
    endif
#endif
#endif

  end subroutine dlnd_comp_run

  !===============================================================================
  subroutine dlnd_comp_final(my_task, master_task, logunit)

    ! !DESCRIPTION:  finalize method for dlnd model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    integer(IN) , intent(in) :: my_task     ! my task in mpi communicator mpicom
    integer(IN) , intent(in) :: master_task ! task number of master task
    integer(IN) , intent(in) :: logunit     ! logging unit number

    !--- formats ---
    character(*), parameter :: F00   = "('(dlnd_comp_final) ',8a)"
    character(*), parameter :: F91   = "('(dlnd_comp_final) ',73('-'))"
    character(*), parameter :: subName = "(dlnd_comp_final) "
    !-------------------------------------------------------------------------------

    call t_startf('DLND_FINAL')

    if (my_task == master_task) then
       write(logunit,F91)
       write(logunit,F00) trim(myModelName),': end of main integration loop'
       write(logunit,F91)
    end if

    call t_stopf('DLND_FINAL')

  end subroutine dlnd_comp_final
  !===============================================================================

end module dlnd_comp_mod
