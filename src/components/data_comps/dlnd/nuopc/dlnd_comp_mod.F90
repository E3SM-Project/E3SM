#ifdef AIX
@PROCESS ALIAS_SIZE(805306368)
#endif

module dlnd_comp_mod

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
  use glc_elevclass_mod     , only : glc_get_num_elevation_classes, glc_elevclass_as_string, glc_elevclass_init 
  use dlnd_shr_mod          , only : datamode        ! namelist input
  use dlnd_shr_mod          , only : decomp          ! namelist input
  use dlnd_shr_mod          , only : rest_file       ! namelist input
  use dlnd_shr_mod          , only : rest_file_strm  ! namelist input
  use dlnd_shr_mod          , only : domain_fracname ! namelist input
  use dlnd_shr_mod          , only : nullstr

  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: dlnd_comp_advertise
  public :: dlnd_comp_init
  public :: dlnd_comp_run

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  type(mct_rearr)             :: rearr
  character(len=CS), pointer  :: avifld(:)           ! char array field names coming from streams
  character(len=CS), pointer  :: avofld(:)           ! char array field names to be sent/recd from med
  character(len=CXX)          :: flds_l2x_mod
  character(len=CXX)          :: flds_x2l_mod
  integer                     :: kf                  ! index for frac in AV
  real(R8), pointer           :: lfrac(:)            ! land frac
  character(len=*), parameter :: rpfile = 'rpointer.lnd'
  integer         , parameter :: nec_len = 2         ! length of elevation class index in field names
  character(*)    , parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine dlnd_comp_advertise(importState, exportState, &
       lnd_present, lnd_prognostic, glc_nec, &
       fldsFrLnd_num, fldsFrLnd, fldsToLnd_num, fldsToLnd, &
       flds_l2x, flds_x2l, rc)

    ! 1. determine export and import fields to advertise to mediator
    ! 2. determine translation of fields from streams to export/import fields

    ! input/output arguments
    type(ESMF_State)                   :: importState
    type(ESMF_State)                   :: exportState
    integer              , intent(in)  :: glc_nec
    logical              , intent(in)  :: lnd_present
    logical              , intent(in)  :: lnd_prognostic
    integer              , intent(out) :: fldsFrLnd_num
    type (fld_list_type) , intent(out) :: fldsFrLnd(:)
    integer              , intent(out) :: fldsToLnd_num
    type (fld_list_type) , intent(out) :: fldsToLnd(:)
    character(len=*)     , intent(out) :: flds_l2x
    character(len=*)     , intent(out) :: flds_x2l
    integer              , intent(out) :: rc

    ! local variables 
    integer :: n
    character(nec_len) :: nec_str         ! elevation class, as character string
    character(len=CS)  :: data_fld_name
    character(len=CS)  :: model_fld_name
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (.not. lnd_present) return

    !-------------------
    ! export fields
    !-------------------
    
    ! scalar fields that need to be advertised

    fldsFrLnd_num=1
    fldsFrLnd(1)%stdname = trim(flds_scalar_name)

    call dshr_fld_add(model_fld="Sl_lfrin", model_fld_concat=flds_l2x, model_fld_index=kf, &
         fldlist_num=fldsFrLnd_num, fldlist=fldsFrLnd)

    ! The actual snow field names will have the elevation class index at the end (e.g., Sl_tsrf01, tsrf01)
    call glc_elevclass_init(glc_nec)
    if (glc_nec > 0) then
       do n = 0, glc_nec
          nec_str = glc_elevclass_as_string(n)

          data_fld_name  = "tsrf" // nec_str
          model_fld_name = "Sl_tsrf" // nec_str
          call dshr_fld_add(data_fld=trim(data_fld_name), data_fld_array=avifld, &
               model_fld=trim(model_fld_name), model_fld_array=avofld, &
               model_fld_concat=flds_l2x, fldlist_num=fldsFrLnd_num, fldlist=fldsFrLnd)

          data_fld_name  = "topo" // nec_str
          model_fld_name = "Sl_topo" // nec_str
          call dshr_fld_add(data_fld=trim(data_fld_name), data_fld_array=avifld, &
               model_fld=trim(model_fld_name), model_fld_array=avofld, &
               model_fld_concat=flds_l2x, fldlist_num=fldsFrLnd_num, fldlist=fldsFrLnd)

          data_fld_name  = "qice" // nec_str
          model_fld_name = "Flgl_qice" // nec_str
          call dshr_fld_add(data_fld=trim(data_fld_name), data_fld_array=avifld, &
               model_fld=trim(model_fld_name), model_fld_array=avofld, &
               model_fld_concat=flds_l2x, fldlist_num=fldsFrLnd_num, fldlist=fldsFrLnd)
       end do
    end if

    ! Non snow fields that nead to be added if dlnd is in cplhist mode
    ! "Sl_t        "
    ! "Sl_tref     "
    ! "Sl_qref     "
    ! "Sl_avsdr    "
    ! "Sl_anidr    "
    ! "Sl_avsdf    "
    ! "Sl_anidf    "
    ! "Sl_snowh    "
    ! "Fall_taux   "
    ! "Fall_tauy   "
    ! "Fall_lat    "
    ! "Fall_sen    "
    ! "Fall_lwup   "
    ! "Fall_evap   "
    ! "Fall_swnet  "
    ! "Sl_landfrac "
    ! "Sl_fv       "
    ! "Sl_ram1     "
    ! "Fall_flxdst1"
    ! "Fall_flxdst2"
    ! "Fall_flxdst3"
    ! "Fall_flxdst4"
    
    do n = 1,fldsFrLnd_num
       call NUOPC_Advertise(exportState, standardName=fldsFrLnd(n)%stdname, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    !-------------------
    ! Save flds_l2x and flds_x2l as module variables for use in debugging
    !-------------------

    flds_x2l_mod = trim(flds_x2l)
    flds_l2x_mod = trim(flds_l2x)

  end subroutine dlnd_comp_advertise

  !===============================================================================

  subroutine dlnd_comp_init(x2l, l2x, &
       SDLND, gsmap, ggrid, mpicom, compid, my_task, master_task, &
       inst_suffix, inst_name, logunit, read_restart, &
       scmMode, scmlat, scmlon, &
       calendar, current_ymd, current_tod)

    ! !DESCRIPTION: initialize dlnd model

    ! !INPUT/OUTPUT PARAMETERS:
    type(mct_aVect)        , intent(inout) :: x2l, l2x     ! input/output attribute vectors
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

    call shr_strdata_pioinit(SDLND, compid)

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

    allocate(lfrac(lsize))
    kfrac = mct_aVect_indexRA(ggrid%data,'frac')
    lfrac(:) = ggrid%data%rAttr(kfrac,:)

    call t_stopf('dlnd_initmctdom')

    !----------------------------------------------------------------------------
    ! Initialize MCT attribute vectors
    !----------------------------------------------------------------------------

    if (my_task == master_task) write(logunit,F00) 'allocate AVs'

    call mct_aVect_init(l2x, rList=flds_l2x_mod, lsize=lsize)
    call mct_aVect_zero(l2x)
    call mct_aVect_init(x2l, rList=flds_x2l_mod, lsize=lsize)
    call mct_aVect_zero(x2l)

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
         SDLND, gsmap, ggrid, mpicom, my_task, master_task, &
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
       SDLND, gsmap, ggrid, mpicom, my_task, master_task, &
       inst_suffix, logunit, read_restart, write_restart, &
       target_ymd, target_tod, case_name)

    ! !DESCRIPTION:  run method for dlnd model

    ! input/output variables:
    type(mct_aVect)        , intent(inout) :: x2l
    type(mct_aVect)        , intent(inout) :: l2x
    type(shr_strdata_type) , intent(inout) :: SDLND
    type(mct_gsMap)        , pointer       :: gsMap
    type(mct_gGrid)        , pointer       :: ggrid
    integer                , intent(in)    :: mpicom           ! mpi communicator
    integer                , intent(in)    :: my_task          ! my task in mpi communicator mpicom
    integer                , intent(in)    :: master_task      ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix      ! char string associated with instance
    integer                , intent(in)    :: logunit          ! logging unit number
    logical                , intent(in)    :: read_restart     ! start from restart
    logical                , intent(in)    :: write_restart    ! write restart
    integer                , intent(in)    :: target_ymd       ! model date
    integer                , intent(in)    :: target_tod       ! model sec into model date
    character(len=*)       , intent(in), optional :: case_name ! case name

    ! local variables
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
    lsize = mct_avect_lsize(l2x)
    do n = 1,lsize
       l2x%rAttr(kf,n) = lfrac(n)
    enddo
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
       call t_stopf('dlnd_restart')
    endif

    !----------------------------------------------------------------------------
    ! Log output for model date
    !----------------------------------------------------------------------------

    if (my_task == master_task) then
       write(logunit,*) ' dlnd: model date ', target_ymd,target_tod
    end if

    call t_stopf('DLND_RUN')

  end subroutine dlnd_comp_run

end module dlnd_comp_mod
