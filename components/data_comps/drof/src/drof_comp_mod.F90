#ifdef AIX
@PROCESS ALIAS_SIZE(805306368)
#endif
module drof_comp_mod

  ! !USES:

  use esmf
  use mct_mod
  use perf_mod
  use shr_sys_mod     , only: shr_sys_flush, shr_sys_abort
  use shr_kind_mod    , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_file_mod    , only: shr_file_getunit, shr_file_freeunit
  use shr_mpi_mod     , only: shr_mpi_bcast
  use shr_strdata_mod , only: shr_strdata_type, shr_strdata_pioinit, shr_strdata_init
  use shr_strdata_mod , only: shr_strdata_print, shr_strdata_restRead
  use shr_strdata_mod , only: shr_strdata_advance, shr_strdata_restWrite
  use shr_dmodel_mod  , only: shr_dmodel_gsmapcreate, shr_dmodel_rearrGGrid
  use shr_dmodel_mod  , only: shr_dmodel_translate_list, shr_dmodel_translateAV_list, shr_dmodel_translateAV
  use shr_cal_mod     , only: shr_cal_ymdtod2string
  use seq_timemgr_mod , only: seq_timemgr_EClockGetData, seq_timemgr_RestartAlarmIsOn

  use drof_shr_mod   , only: datamode       ! namelist input
  use drof_shr_mod   , only: decomp         ! namelist input
  use drof_shr_mod   , only: rest_file      ! namelist input
  use drof_shr_mod   , only: rest_file_strm ! namelist input
  use drof_shr_mod   , only: nullstr
#ifdef HAVE_MOAB
  use seq_comm_mct,     only : mrofid ! id of moab rof app
  use seq_comm_mct,     only : mbrof_data ! turn on if the data rof
  use iso_c_binding
#endif
  !
  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: drof_comp_init
  public :: drof_comp_run
  public :: drof_comp_final

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  character(CS) :: myModelName = 'rof'   ! user defined model name
  character(len=*), parameter :: rpfile = 'rpointer.rof'
  type(mct_rearr) :: rearr

  !--------------------------------------------------------------------------
  integer(IN),parameter :: ktrans = 2

  character(12),parameter  :: avofld(1:ktrans) = &
       (/ "Forr_rofl   ","Forr_rofi   "/)

  character(12),parameter  :: avifld(1:ktrans) = &
       (/ "rofl        ","rofi        "/)
  !--------------------------------------------------------------------------

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !===============================================================================
  subroutine drof_comp_init(Eclock, x2r, r2x, &
       seq_flds_x2r_fields, seq_flds_r2x_fields, &
       SDROF, gsmap, ggrid, mpicom, compid, my_task, master_task, &
       inst_suffix, inst_name, logunit, read_restart)

#ifdef HAVE_MOAB
    use iMOAB, only: iMOAB_DefineTagStorage, iMOAB_GetDoubleTagStorage, &
                     iMOAB_SetIntTagStorage, iMOAB_SetDoubleTagStorage, &
                     iMOAB_ResolveSharedEntities, iMOAB_CreateVertices, &
                     iMOAB_GetMeshInfo, iMOAB_UpdateMeshInfo, iMOAB_WriteMesh
#endif
    ! !DESCRIPTION: initialize drof model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)       , intent(in)    :: EClock
    type(mct_aVect)        , intent(inout) :: x2r, r2x            ! input/output attribute vectors
    character(len=*)       , intent(in)    :: seq_flds_x2r_fields ! fields from mediator
    character(len=*)       , intent(in)    :: seq_flds_r2x_fields ! fields to mediator
    type(shr_strdata_type) , intent(inout) :: SDROF               ! model shr_strdata instance (output)
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

    !--- local variables ---
    integer(IN)   :: lsize       ! local size
    logical       :: exists      ! file existance logical
    integer(IN)   :: nu          ! unit number
    character(CL) :: calendar ! model calendar
#ifdef HAVE_MOAB
    character*400  tagname
    real(R8) latv, lonv
    integer iv, tagindex, ilat, ilon, ierr  !, arrsize, nfields
    real(R8), allocatable, target :: data(:)
    integer(IN), pointer :: idata(:)   ! temporary
    real(r8), dimension(:), allocatable :: moab_vert_coords  ! temporary
    !real(R8), allocatable, target :: vtags_zero(:, :)
#ifdef MOABDEBUG
    character*100 outfile, wopts
#endif
#endif
    !--- formats ---
    character(*), parameter :: F00   = "('(drof_comp_init) ',8a)"
    character(*), parameter :: F0L   = "('(drof_comp_init) ',a, l2)"
    character(*), parameter :: F01   = "('(drof_comp_init) ',a,5i8)"
    character(*), parameter :: F02   = "('(drof_comp_init) ',a,4es13.6)"
    character(*), parameter :: F03   = "('(drof_comp_init) ',a,i8,a)"
    character(*), parameter :: F05   = "('(drof_comp_init) ',a,2f10.4)"
    character(*), parameter :: F90   = "('(drof_comp_init) ',73('='))"
    character(*), parameter :: F91   = "('(drof_comp_init) ',73('-'))"
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

    call seq_timemgr_EClockGetData( EClock, calendar=calendar )

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
    call shr_sys_flush(logunit)

    ! create a data model global seqmap (gsmap) given the data model global grid sizes
    ! NOTE: gsmap is initialized using the decomp read in from the docn_in namelist
    ! (which by default is "1d")
    call shr_dmodel_gsmapcreate(gsmap, SDROF%nxg*SDROF%nyg, compid, mpicom, decomp)
    lsize = mct_gsmap_lsize(gsmap,mpicom)

    ! create a rearranger from the data model SDOCN%gsmap to gsmap
    call mct_rearr_init(SDROF%gsmap, gsmap, mpicom, rearr)

    call t_stopf('drof_initgsmaps')

    !----------------------------------------------------------------------------
    ! Initialize MCT domain
    !----------------------------------------------------------------------------

    call t_startf('drof_initmctdom')
    if (my_task == master_task) write(logunit,F00) 'copy domains'
    call shr_sys_flush(logunit)

    call shr_dmodel_rearrGGrid(SDROF%grid, ggrid, gsmap, rearr, mpicom)

    call t_stopf('drof_initmctdom')

! copy from data atm ; just need mrofid

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
   ierr = iMOAB_CreateVertices(mrofid, lsize*3, 3, moab_vert_coords)
   if (ierr .ne. 0)  &
      call shr_sys_abort('Error: fail to create MOAB vertices in land model')

   tagname='GLOBAL_ID'//C_NULL_CHAR
   ierr = iMOAB_DefineTagStorage(mrofid, tagname, &
                                 0, & ! dense, integer
                                 1, & ! number of components
                                 tagindex )
   if (ierr .ne. 0)  &
      call shr_sys_abort('Error: fail to retrieve GLOBAL_ID tag ')

   ! get list of global IDs for Dofs
   call mct_gsMap_orderedPoints(gsMap, my_task, idata)

   ierr = iMOAB_SetIntTagStorage ( mrofid, tagname, lsize, &
                                    0, & ! vertex type
                                    idata)
   if (ierr .ne. 0)  &
      call shr_sys_abort('Error: fail to set GLOBAL_ID tag ')

   ierr = iMOAB_ResolveSharedEntities( mrofid, lsize, idata );
   if (ierr .ne. 0)  &
      call shr_sys_abort('Error: fail to resolve shared entities')

   deallocate(moab_vert_coords)
   deallocate(idata)

   ierr = iMOAB_UpdateMeshInfo( mrofid )
   if (ierr .ne. 0)  &
      call shr_sys_abort('Error: fail to update mesh info ')

   allocate(data(lsize))
   ierr = iMOAB_DefineTagStorage( mrofid, "area:aream:frac:mask"//C_NULL_CHAR, &
                                    1, & ! dense, double
                                    1, & ! number of components
                                    tagindex )
   if (ierr > 0 )  &
      call shr_sys_abort('Error: fail to create tag: area:aream:frac:mask' )

   data(:) = ggrid%data%rAttr(mct_aVect_indexRA(ggrid%data,'area'),:)
   tagname='area'//C_NULL_CHAR
   ierr = iMOAB_SetDoubleTagStorage ( mrofid, tagname, lsize, &
                                    0, & ! set data on vertices
                                    data)
   if (ierr > 0 )  &
      call shr_sys_abort('Error: fail to get area tag ')

   ! set the same data for aream (model area) as area
   ! data(:) = ggrid%data%rAttr(mct_aVect_indexRA(ggrid%data,'aream'),:)
   tagname='aream'//C_NULL_CHAR
   ierr = iMOAB_SetDoubleTagStorage ( mrofid, tagname, lsize, &
                                    0, & ! set data on vertices
                                    data)
   if (ierr > 0 )  &
      call shr_sys_abort('Error: fail to set aream tag ')

   data(:) = ggrid%data%rAttr(mct_aVect_indexRA(ggrid%data,'mask'),:)
   tagname='mask'//C_NULL_CHAR
   ierr = iMOAB_SetDoubleTagStorage ( mrofid, tagname, lsize, &
                                    0, & ! set data on vertices
                                    data)
   if (ierr > 0 )  &
      call shr_sys_abort('Error: fail to set mask tag ')

   data(:) = ggrid%data%rAttr(mct_aVect_indexRA(ggrid%data,'frac'),:)
   tagname='frac'//C_NULL_CHAR
   ierr = iMOAB_SetDoubleTagStorage ( mrofid, tagname, lsize, &
                                    0, & ! set data on vertices
                                    data)
   if (ierr > 0 )  &
      call shr_sys_abort('Error: fail to set frac tag ')

   deallocate(data)

   ! define tags
   ierr = iMOAB_DefineTagStorage( mrofid, trim(seq_flds_x2r_fields)//C_NULL_CHAR, &
                                    1, & ! dense, double
                                    1, & ! number of components
                                    tagindex )
   if (ierr > 0 )  &
      call shr_sys_abort('Error: fail to create seq_flds_x2r_fields tags ')

   ierr = iMOAB_DefineTagStorage( mrofid, trim(seq_flds_r2x_fields)//C_NULL_CHAR, &
                                    1, & ! dense, double
                                    1, & ! number of components
                                    tagindex )
   if (ierr > 0 )  &
      call shr_sys_abort('Error: fail to create seq_flds_r2x_fields tags ')
   mbrof_data = .true. ! will have effects 
#ifdef MOABDEBUG
      !      debug test
   outfile = 'RofDataMesh.h5m'//C_NULL_CHAR
   wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR !
      !      write out the mesh file to disk
   ierr = iMOAB_WriteMesh(mrofid, trim(outfile), trim(wopts))
   if (ierr .ne. 0) then
      call shr_sys_abort(subname//' ERROR in writing data mesh rof ')
   endif
#endif
#endif
    !----------------------------------------------------------------------------
    ! Initialize MCT attribute vectors
    !----------------------------------------------------------------------------

    call t_startf('drof_initmctavs')
    if (my_task == master_task) write(logunit,F00) 'allocate AVs'
    call shr_sys_flush(logunit)

    call mct_aVect_init(x2r, rList=seq_flds_x2r_fields, lsize=lsize)
    call mct_aVect_zero(x2r)

    call mct_aVect_init(r2x, rList=seq_flds_r2x_fields, lsize=lsize)
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
       end if
       call shr_mpi_bcast(exists,mpicom,'exists')
       if (exists) then
          if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file_strm)
          call shr_strdata_restRead(trim(rest_file_strm),SDROF,mpicom)
       else
          if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file_strm)
       endif
       call shr_sys_flush(logunit)
    end if

    !----------------------------------------------------------------------------
    ! Set initial rof state
    !----------------------------------------------------------------------------

    call t_adj_detailf(+2)
    call drof_comp_run(EClock, x2r, r2x, &
         SDROF, gsmap, ggrid, mpicom, compid, my_task, master_task, &
         inst_suffix, logunit)
    call t_adj_detailf(-2)

    if (my_task == master_task) write(logunit,F00) 'drof_comp_init done'
    call shr_sys_flush(logunit)

    call t_stopf('DROF_INIT')

  end subroutine drof_comp_init

  !===============================================================================
  subroutine drof_comp_run(EClock, x2r, r2x, &
       SDROF, gsmap, ggrid, mpicom, compid, my_task, master_task, &
       inst_suffix, logunit, case_name)

#ifdef MOABDEBUG
    use iMOAB, only: iMOAB_WriteMesh
#endif
#ifdef HAVE_MOAB
    use seq_flds_mod    , only: seq_flds_r2x_fields 
    use seq_flds_mod    , only: moab_set_tag_from_av
#endif
    ! !DESCRIPTION:  run method for drof model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)       , intent(in)    :: EClock
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
    character(CL)          , intent(in), optional :: case_name ! case name

    !--- local ---
    integer(IN)   :: CurrentYMD        ! model date
    integer(IN)   :: CurrentTOD        ! model sec into model date
    integer(IN)   :: yy,mm,dd          ! year month day
    integer(IN)   :: n                 ! indices
    integer(IN)   :: nf                ! fields loop index
    integer(IN)   :: lsize             ! size of attr vect
    logical       :: write_restart     ! restart now
    integer(IN)   :: nu                ! unit number
    integer(IN)   :: nflds_r2x
    character(len=18) :: date_str
#ifdef HAVE_MOAB
    real(R8), allocatable, target :: datam(:)
    type(mct_list) :: temp_list
    integer :: size_list, index_list
    type(mct_string)    :: mctOStr  !
    character*400  tagname, mct_field
#ifdef MOABDEBUG
    integer  :: cur_drof_stepno, ierr
    character*100 outfile, wopts, lnum
#endif

#endif
    character(*), parameter :: F00   = "('(drof_comp_run) ',8a)"
    character(*), parameter :: F04   = "('(drof_comp_run) ',2a,2i8,'s')"
    character(*), parameter :: subName = "(drof_comp_run) "
    !-------------------------------------------------------------------------------

    call t_startf('DROF_RUN')

    call t_startf('drof_run1')

    call seq_timemgr_EClockGetData( EClock, curr_ymd=CurrentYMD, curr_tod=CurrentTOD)
    call seq_timemgr_EClockGetData( EClock, curr_yr=yy, curr_mon=mm, curr_day=dd)
    write_restart = seq_timemgr_RestartAlarmIsOn(EClock)

    call t_stopf('drof_run1')

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

    call shr_strdata_advance(SDROF, currentYMD, currentTOD, mpicom, 'drof_r')
    call t_stopf('drof_r_strdata_advance')

    !--- copy streams to r2x ---
    call t_barrierf('drof_r_scatter_BARRIER', mpicom)
    call t_startf('drof_r_scatter')
    do n = 1,SDROF%nstreams
       call shr_dmodel_translateAV(SDROF%avs(n), r2x, avifld, avofld, rearr)
    enddo
    call t_stopf('drof_r_scatter')

    ! zero out "special values"
    lsize = mct_avect_lsize(r2x)
    nflds_r2x = mct_avect_nRattr(r2x)
    do nf=1,nflds_r2x
       do n=1,lsize
          if (abs(r2x%rAttr(nf,n)) > 1.0e28) then
             r2x%rAttr(nf,n) = 0.0_r8
          end if
          ! write(6,*)'crrentymd, currenttod, nf,n,r2x= ',currentymd, currenttod, nf,n,r2x%rattr(nf,n)
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
       call shr_cal_ymdtod2string(date_str, yy, mm, dd, currentTOD)
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
       if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file_strm),currentYMD,currentTOD
       call shr_strdata_restWrite(trim(rest_file_strm), SDROF, mpicom, trim(case_name), 'SDROF strdata')
       call shr_sys_flush(logunit)
       call t_stopf('drof_restart')
    end if

    !----------------------------------------------------------------------------
    ! Log output for model date
    !----------------------------------------------------------------------------
#ifdef HAVE_MOAB
    lsize = mct_avect_lsize(r2x) ! is it the same as mct_avect_lsize(avstrm) ?
    allocate(datam(lsize)) ! 
    call mct_list_init(temp_list ,seq_flds_r2x_fields)
    size_list=mct_list_nitem (temp_list)
    do index_list = 1, size_list
      call mct_list_get(mctOStr,index_list,temp_list)
      mct_field = mct_string_toChar(mctOStr)
      tagname= trim(mct_field)//C_NULL_CHAR
      call moab_set_tag_from_av(tagname, r2x, index_list, mrofid, datam, lsize) ! loop over all a2x fields, not just a few
    enddo
    call mct_list_clean(temp_list)
    deallocate(datam) ! maybe we should keep it around, deallocate at the final only?

#ifdef MOABDEBUG
    call seq_timemgr_EClockGetData( EClock, stepno=cur_drof_stepno )
    write(lnum,"(I0.2)")cur_drof_stepno
    outfile = 'drof_comp_run_'//trim(lnum)//'.h5m'//C_NULL_CHAR
    wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
    ierr = iMOAB_WriteMesh(mrofid, outfile, wopts)
    if (ierr > 0 )  then
       write(logunit,*) 'Failed to write data rof component state '
    endif
#endif

#endif

    call t_startf('drof_run2')
    if (my_task == master_task) then
       write(logunit,F04) trim(myModelName),': model date ', CurrentYMD,CurrentTOD
       call shr_sys_flush(logunit)
    end if
    call t_stopf('drof_run2')

    call t_stopf('DROF_RUN')

  end subroutine drof_comp_run

  !===============================================================================
  subroutine drof_comp_final(my_task, master_task, logunit)

    ! !DESCRIPTION:  finalize method for drof model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    integer(IN) , intent(in) :: my_task     ! my task in mpi communicator mpicom
    integer(IN) , intent(in) :: master_task ! task number of master task
    integer(IN) , intent(in) :: logunit     ! logging unit number

    !--- formats ---
    character(*), parameter :: F00   = "('(drof_comp_final) ',8a)"
    character(*), parameter :: F91   = "('(drof_comp_final) ',73('-'))"
    character(*), parameter :: subName = "(drof_comp_final) "
    !-------------------------------------------------------------------------------

    call t_startf('DROF_FINAL')

    if (my_task == master_task) then
       write(logunit,F91)
       write(logunit,F00) trim(myModelName),': end of main integration loop'
       write(logunit,F91)
    end if

    call t_stopf('DROF_FINAL')

  end subroutine drof_comp_final
  !===============================================================================

end module drof_comp_mod
