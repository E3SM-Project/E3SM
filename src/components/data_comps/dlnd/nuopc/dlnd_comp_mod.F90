#ifdef AIX
@PROCESS ALIAS_SIZE(805306368)
#endif

module dlnd_comp_mod

  ! !USES:
  use mct_mod
  use perf_mod          , only: t_startf, t_stopf, t_barrierf, t_adj_detailf
  use shr_pcdf_mod      , only: shr_pcdf_readwrite 
  use shr_kind_mod      , only: R8=>SHR_KIND_R8
  use shr_file_mod      , only: shr_file_getunit, shr_file_freeunit
  use shr_mpi_mod       , only: shr_mpi_bcast
  use shr_strdata_mod   , only: shr_strdata_type, shr_strdata_pioinit, shr_strdata_init
  use shr_strdata_mod   , only: shr_strdata_print, shr_strdata_restRead
  use shr_strdata_mod   , only: shr_strdata_advance, shr_strdata_restWrite
  use shr_dmodel_mod    , only: shr_dmodel_gsmapcreate, shr_dmodel_rearrGGrid
  use shr_dmodel_mod    , only: shr_dmodel_translate_list, shr_dmodel_translateAV_list, shr_dmodel_translateAV
  use shr_cal_mod       , only: shr_cal_datetod2string
  use glc_elevclass_mod , only: glc_get_num_elevation_classes, glc_elevclass_as_string

  use dlnd_shr_mod      , only: datamode       ! namelist input
  use dlnd_shr_mod      , only: decomp         ! namelist input
  use dlnd_shr_mod      , only: rest_file      ! namelist input
  use dlnd_shr_mod      , only: rest_file_strm ! namelist input
  use dlnd_shr_mod      , only: domain_fracname ! namelist input
  use dlnd_shr_mod      , only: nullstr

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
  character(len=3)           :: myModelName = 'lnd' ! user defined model name
  type(mct_rearr)            :: rearr               ! MCT rearranger
  integer                    :: kf                  ! index for frac in AV
  real(R8),pointer           :: lfrac(:)            ! land frac
  character(len=*),parameter :: rpfile = 'rpointer.lnd'

  !--- names of fields ---
  integer    , parameter :: fld_len      = 12  ! max character length of fields in avofld & avifld
  integer    , parameter :: nflds_nosnow = 0   ! for now no fields without snow
  integer    , parameter :: nflds_snow   = 3   ! number of snow fields in each elevation class
  integer    , parameter :: nec_len      = 2   ! length of elevation class index in field names

  ! for these snow fields, the actual field names will have the elevation class index at
  ! the end (e.g., Sl_tsrf01, tsrf01)
  character(fld_len-nec_len),parameter :: avofld_snow(nflds_snow) = (/"Sl_tsrf", "Sl_topo", "Flgl_qice"/)
  character(fld_len-nec_len),parameter :: avifld_snow(nflds_snow) = (/"tsrf"   , "topo"   , "qice"     /)

  ! all fields:
  character(fld_len),dimension(:),allocatable :: avofld
  character(fld_len),dimension(:),allocatable :: avifld

!===============================================================================
contains
!===============================================================================

  subroutine dlnd_comp_init(x2l, l2x, &
       x2l_fields, l2x_fields, &
       SDLND, gsmap, ggrid, mpicom, compid, my_task, master_task, &
       inst_suffix, inst_name, logunit, read_restart, &
       scmMode, scmlat, scmlon, &
       calendar, current_ymd, current_tod)

    ! !DESCRIPTION: initialize dlnd model

    ! !INPUT/OUTPUT PARAMETERS:
    type(mct_aVect)        , intent(inout) :: x2l, l2x     ! input/output attribute vectors
    character(len=*)       , intent(in)    :: x2l_fields   ! fields from mediator
    character(len=*)       , intent(in)    :: l2x_fields   ! fields to mediator
    type(shr_strdata_type) , intent(inout) :: SDLND        ! model shr_strdata instance (output)
    type(mct_gsMap)        , pointer       :: gsMap        ! model global seg map (output)
    type(mct_gGrid)        , pointer       :: ggrid        ! model ggrid (output)
    integer                , intent(in)    :: mpicom       ! mpi communicator
    integer                , intent(in)    :: compid       ! mct comp id
    integer                , intent(in)    :: my_task      ! my task in mpi communicator mpicom
    integer                , intent(in)    :: master_task  ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix  ! char string associated with instance
    character(len=*)       , intent(in)    :: inst_name    ! fullname of current instance (ie. "lnd_0001")
    integer                , intent(in)    :: logunit      ! logging unit number
    logical                , intent(in)    :: read_restart ! start from restart
    logical                , intent(in)    :: scmMode      ! single column mode
    real(R8)               , intent(in)    :: scmLat       ! single column lat
    real(R8)               , intent(in)    :: scmLon       ! single column lon
    character(len=*)       , intent(in)    :: calendar     ! calendar name
    integer                , intent(in)    :: current_ymd  ! model date
    integer                , intent(in)    :: current_tod  ! model sec into model date

    !--- local variables ---
    integer            :: n,k             ! generic counters
    integer            :: lsize           ! local size
    logical            :: exists          ! file existance
    integer            :: nu              ! unit number
    integer            :: glc_nec         ! number of elevation classes
    integer            :: nflds_glc_nec   ! number of snow fields separated by elev class
    integer            :: field_num       ! field number
    character(nec_len) :: nec_str         ! elevation class, as character string
    integer            :: kfrac           ! AV index
    logical            :: write_restart

    !--- formats ---
    character(*), parameter :: F00   = "('(dlnd_comp_init) ',8a)"
    character(*), parameter :: F05   = "('(dlnd_comp_init) ',a,2f10.4)"
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

    ! NOTE: shr_strdata_init calls shr_dmodel_readgrid which reads the data model
    ! grid and from that computes SDLND%gsmap and SDLND%ggrid. DLND%gsmap is created
    ! using the decomp '2d1d' (1d decomp of 2d grid)

    if (scmmode) then
       if (my_task == master_task) &
            write(logunit,F05) ' scm lon lat = ',scmlon,scmlat
       call shr_strdata_init(SDLND,mpicom,compid,name='lnd', &
            scmmode=scmmode,scmlon=scmlon,scmlat=scmlat, calendar=calendar, &
            dmodel_domain_fracname_from_stream=domain_fracname)
    else
       call shr_strdata_init(SDLND,mpicom,compid,name='lnd', calendar=calendar, &
            dmodel_domain_fracname_from_stream=domain_fracname)
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

    ! For now only use snow fields
    allocate(avofld(nflds_nosnow + nflds_glc_nec))
    allocate(avifld(nflds_nosnow + nflds_glc_nec))

    ! Append each snow field
    field_num = nflds_nosnow
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

    ! create a data model global seqmap (gsmap) given the data model global grid sizes
    ! NOTE: gsmap is initialized using the decomp read in from the dlnd_in namelist
    ! (which by default is "1d")
    call shr_dmodel_gsmapcreate(gsmap,SDLND%nxg*SDLND%nyg,compid,mpicom,decomp)
    lsize = mct_gsmap_lsize(gsmap,mpicom)

    ! create a rearranger from the data model SDLND%gsmap to gsmap
    call mct_rearr_init(SDLND%gsmap, gsmap, mpicom, rearr)

    call t_stopf('dlnd_initgsmaps')

    !----------------------------------------------------------------------------
    ! Initialize MCT domain
    !----------------------------------------------------------------------------

    call t_startf('dlnd_initmctdom')
    if (my_task == master_task) write(logunit,F00) 'copy domains'

    call shr_dmodel_rearrGGrid(SDLND%grid, ggrid, gsmap, rearr, mpicom)
    kfrac = mct_aVect_indexRA(ggrid%data,'frac')
    allocate(lfrac(lsize))
    lfrac(:) = ggrid%data%rAttr(kfrac,:)

    call t_stopf('dlnd_initmctdom')

    !----------------------------------------------------------------------------
    ! Initialize MCT attribute vectors
    !----------------------------------------------------------------------------

    if (my_task == master_task) write(logunit,F00) 'allocate AVs'

    call mct_aVect_init(l2x, rList=l2x_fields, lsize=lsize)
    call mct_aVect_zero(l2x)
    call mct_aVect_init(x2l, rList=x2l_fields, lsize=lsize)
    call mct_aVect_zero(x2l)

    kf = mct_aVect_indexRA(l2x, 'Sl_lfrin', perrwith='quiet')

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
       !if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file)
       !call shr_pcdf_readwrite('read',trim(rest_file),mpicom,gsmap,rf1=somtp,rf1n='somtp')
       if (exists) then
          if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file_strm)
          call shr_strdata_restRead(trim(rest_file_strm),SDLND,mpicom)
       else
          if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file_strm)
       endif
    endif

    !----------------------------------------------------------------------------
    ! Set initial lnd state
    !----------------------------------------------------------------------------

    call t_adj_detailf(+2)

    write_restart = .false.
    call dlnd_comp_run(x2l, l2x, &
         SDLND, gsmap, ggrid, mpicom, compid, my_task, master_task, &
         inst_suffix, logunit, read_restart, write_restart, &
         current_ymd, current_tod)

    call t_adj_detailf(-2)

    if (my_task == master_task) then
       write(logunit,F00) 'dlnd_comp_init done'
    end if

    call t_stopf('DLND_INIT')

  end subroutine dlnd_comp_init

  !===============================================================================

  subroutine dlnd_comp_run(x2l, l2x, &
       SDLND, gsmap, ggrid, mpicom, compid, my_task, master_task, &
       inst_suffix, logunit, read_restart, write_restart, &
       target_ymd, target_tod, case_name)

    ! !DESCRIPTION:  run method for dlnd model

    ! !INPUT/OUTPUT PARAMETERS:
    type(mct_aVect)        , intent(inout) :: x2l
    type(mct_aVect)        , intent(inout) :: l2x
    type(shr_strdata_type) , intent(inout) :: SDLND
    type(mct_gsMap)        , pointer       :: gsMap
    type(mct_gGrid)        , pointer       :: ggrid
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

    !--- local ---
    integer       :: n                     ! indices
    integer       :: nu                    ! unit number
    integer       :: lsize                 ! local size
    character(len=18) :: date_str

    character(*), parameter :: F00   = "('(dlnd_comp_run) ',8a)"
    character(*), parameter :: F04   = "('(dlnd_comp_run) ',2a,2i8,'s')"
    character(*), parameter :: subName = "(dlnd_comp_run) "
    !-------------------------------------------------------------------------------

    call t_startf('DLND_RUN')

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
    if (kf /= 0) then
       lsize = mct_avect_lsize(l2x)
       do n = 1,lsize
          l2x%rAttr(kf,n) = lfrac(n)
       enddo
    end if
    call shr_strdata_advance(SDLND,target_ymd,target_tod,mpicom,'dlnd')
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
       call shr_cal_datetod2string(date_str, target_ymd, target_tod)
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
       if (my_task == master_task) then
          write(logunit,F04) ' writing ',trim(rest_file_strm),target_ymd,target_tod
       end if
       call shr_strdata_restWrite(trim(rest_file_strm),SDLND,mpicom,trim(case_name),'SDLND strdata')
       call shr_sys_flush(logunit)
       call t_stopf('dlnd_restart')
    endif

    !----------------------------------------------------------------------------
    ! Log output for model date
    !----------------------------------------------------------------------------

    if (my_task == master_task) then
       write(logunit,F04) trim(myModelName),': model date ', target_ymd,target_tod
    end if

    call t_stopf('DLND_RUN')

  end subroutine dlnd_comp_run

  !===============================================================================

  subroutine dlnd_comp_final(my_task, master_task, logunit)

    ! !DESCRIPTION:  finalize method for dlnd model

    ! !INPUT/OUTPUT PARAMETERS:
    integer , intent(in) :: my_task     ! my task in mpi communicator mpicom
    integer , intent(in) :: master_task ! task number of master task
    integer , intent(in) :: logunit     ! logging unit number

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
