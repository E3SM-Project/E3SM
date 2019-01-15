#ifdef AIX
@PROCESS ALIAS_SIZE(805306368)
#endif
module docn_comp_mod

  ! !USES:
  use shr_pcdf_mod          , only : shr_pcdf_readwrite
  use NUOPC                 , only : NUOPC_Advertise
  use ESMF                  , only : ESMF_State, ESMF_SUCCESS, ESMF_State
  use ESMF                  , only : ESMF_Mesh, ESMF_DistGrid, ESMF_MeshGet, ESMF_DistGridGet
  use ESMF                  , only : ESMF_State, ESMF_LOGMSG_INFO, ESMF_LogWrite
  use perf_mod              , only : t_startf, t_stopf
  use perf_mod              , only : t_adj_detailf, t_barrierf
  use mct_mod               , only : mct_gsmap, mct_gsmap_init, mct_gsmap_lsize
  use mct_mod               , only : mct_avect, mct_avect_indexRA, mct_avect_zero, mct_aVect_nRattr
  use mct_mod               , only : mct_avect_init, mct_avect_lsize, mct_avect_clean
  use med_constants_mod     , only : R8, CS, CXX
  use shr_const_mod         , only : shr_const_cpsw, shr_const_rhosw, shr_const_TkFrz
  use shr_const_mod         , only : shr_const_TkFrzSw, shr_const_latice, shr_const_ocn_ref_sal
  use shr_const_mod         , only : shr_const_zsrflyr, shr_const_pi
  use shr_string_mod        , only : shr_string_listGetName
  use shr_sys_mod           , only : shr_sys_abort
  use shr_file_mod          , only : shr_file_getunit, shr_file_freeunit
  use shr_mpi_mod           , only : shr_mpi_bcast
  use shr_frz_mod           , only : shr_frz_freezetemp
  use shr_cal_mod           , only : shr_cal_calendarname
  use shr_cal_mod           , only : shr_cal_datetod2string
  use shr_nuopc_scalars_mod , only : flds_scalar_name
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
  use shr_strdata_mod       , only : shr_strdata_init_model_domain
  use shr_strdata_mod       , only : shr_strdata_init_streams
  use shr_strdata_mod       , only : shr_strdata_init_mapping
  use shr_strdata_mod       , only : shr_strdata_type, shr_strdata_pioinit
  use shr_strdata_mod       , only : shr_strdata_print, shr_strdata_restRead
  use shr_strdata_mod       , only : shr_strdata_advance, shr_strdata_restWrite
  use shr_dmodel_mod        , only : shr_dmodel_translateAV
  use dshr_nuopc_mod        , only : fld_list_type, dshr_fld_add
  use docn_shr_mod          , only : datamode       ! namelist input
  use docn_shr_mod          , only : aquap_option   ! derived from datamode namelist input
  use docn_shr_mod          , only : rest_file      ! namelist input
  use docn_shr_mod          , only : rest_file_strm ! namelist input
  use docn_shr_mod          , only : nullstr

  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public/Private interfaces
  !--------------------------------------------------------------------------

  public :: docn_comp_advertise
  public :: docn_comp_init
  public :: docn_comp_run

  private :: prescribed_sst

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  integer                    :: debug_import = 0      ! debug level (if > 0 will print all import fields)
  integer                    :: debug_export = 0      ! debug level (if > 0 will print all export fields)

  real(R8),parameter         :: cpsw    = shr_const_cpsw        ! specific heat of sea h2o ~ J/kg/K
  real(R8),parameter         :: rhosw   = shr_const_rhosw       ! density of sea water ~ kg/m^3
  real(R8),parameter         :: TkFrz   = shr_const_TkFrz       ! freezing point, fresh water (Kelvin)
  real(R8),parameter         :: TkFrzSw = shr_const_TkFrzSw     ! freezing point, sea   water (Kelvin)
  real(R8),parameter         :: latice  = shr_const_latice      ! latent heat of fusion
  real(R8),parameter         :: ocnsalt = shr_const_ocn_ref_sal ! ocean reference salinity

  integer                    :: kt,ks,ku,kv,kdhdx,kdhdy,kq,kswp ! field indices
  integer                    :: kswnet,klwup,klwdn
  integer                    :: ksen,klat,kmelth,ksnow,krofi
  integer                    :: kh,kqbot
  integer                    :: kmask, kfrac                    ! frac and mask field indices of docn domain
  integer                    :: ksomask                         ! So_omask field index

  type(mct_avect)            :: avstrm                          ! av of data created from all stream input
  character(len=CS), pointer :: avifld(:)                       ! names of fields in input streams
  character(len=CS), pointer :: avofld(:)                       ! local names of fields in input streams for import/export
  character(len=CS), pointer :: stifld(:)                       ! names of fields in input streams
  character(len=CS), pointer :: stofld(:)                       ! local names of fields in input streams for calculations
  character(CXX)             :: flds_strm = ''                  ! set in docn_comp_init
  character(len=CXX)         :: flds_o2x_mod                    ! set in docn_comp_advertise
  character(len=CXX)         :: flds_x2o_mod                    ! set in docn_comp_advertise
  logical                    :: ocn_prognostic_mod              ! set in docn_comp_advertise

  integer , pointer          :: imask(:)                        ! integer ocean mask
  real(R8), pointer          :: xc(:), yc(:)                    ! arrays of model latitudes and longitudes
  real(R8), pointer          :: somtp(:)                        ! SOM ocean temperature
  real(R8), pointer          :: tfreeze(:)                      ! SOM ocean freezing temperature

  logical                    :: firstcall = .true.              ! first call logical
  character(len=*),parameter :: rpfile = 'rpointer.ocn'         ! name of ocean ropinter file
  character(*),parameter     :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine docn_comp_advertise(importState, exportState, &
       ocn_present, ocn_prognostic, ocnrof_prognostic, &
       fldsFrOcn_num, fldsFrOcn, fldsToOcn_num, fldsToOcn, &
       flds_o2x, flds_x2o, rc)

    ! input/output arguments
    type(ESMF_State)     , intent(inout) :: importState
    type(ESMF_State)     , intent(inout) :: exportState
    logical              , intent(in)    :: ocn_present
    logical              , intent(in)    :: ocn_prognostic
    logical              , intent(in)    :: ocnrof_prognostic
    integer              , intent(out)   :: fldsToOcn_num
    integer              , intent(out)   :: fldsFrOcn_num
    type (fld_list_type) , intent(out)   :: fldsToOcn(:)
    type (fld_list_type) , intent(out)   :: fldsFrOcn(:)
    character(len=*)     , intent(out)   :: flds_o2x
    character(len=*)     , intent(out)   :: flds_x2o
    integer              , intent(out)   :: rc

    ! local variables
    integer         :: n
    !-------------------------------------------------------------------------------

    if (.not. ocn_present) return

    !--------------------------------
    ! export fields
    !--------------------------------

    fldsFrOcn_num=1
    fldsFrOcn(1)%stdname = trim(flds_scalar_name)

    ! export fields that have no corresponding stream field (computed internally)

    call dshr_fld_add(model_fld='So_omask', model_fld_concat=flds_o2x, model_fld_index=ksomask, &
         fldlist_num=fldsFrOcn_num, fldlist=fldsFrOcn)

    call dshr_fld_add(model_fld='Fioo_q', model_fld_concat=flds_o2x, model_fld_index=kq, &
         fldlist_num=fldsFrOcn_num, fldlist=fldsFrOcn)

    ! export fields that have a corresponding stream field

    call dshr_fld_add(data_fld='t', data_fld_array=avifld, model_fld='So_t', model_fld_array=avofld, &
         model_fld_concat=flds_o2x, model_fld_index=kt, fldlist_num=fldsFrOcn_num, fldlist=fldsFrOcn)

    call dshr_fld_add(data_fld='s', data_fld_array=avifld, model_fld='So_s', model_fld_array=avofld, &
         model_fld_concat=flds_o2x, model_fld_index=ks, fldlist_num=fldsFrOcn_num, fldlist=fldsFrOcn)

    call dshr_fld_add(data_fld='u', data_fld_array=avifld, model_fld='So_u', model_fld_array=avofld, &
         model_fld_concat=flds_o2x, model_fld_index=ku, fldlist_num=fldsFrOcn_num, fldlist=fldsFrOcn)

    call dshr_fld_add(data_fld='v', data_fld_array=avifld, model_fld='So_v', model_fld_array=avofld, &
         model_fld_concat=flds_o2x, model_fld_index=kv, fldlist_num=fldsFrOcn_num, fldlist=fldsFrOcn)

    call dshr_fld_add(data_fld='dhdx', data_fld_array=avifld, model_fld='So_dhdx', model_fld_array=avofld, &
         model_fld_concat=flds_o2x, model_fld_index=kdhdx, fldlist_num=fldsFrOcn_num, fldlist=fldsFrOcn)

    call dshr_fld_add(data_fld='dhdy', data_fld_array=avifld, model_fld='So_dhdy', model_fld_array=avofld, &
         model_fld_concat=flds_o2x, model_fld_index=kdhdy, fldlist_num=fldsFrOcn_num, fldlist=fldsFrOcn)

    !-------------------
    ! import fields (have no corresponding stream fields)
    !-------------------

    if (ocn_prognostic) then

       fldsToOcn_num=1
       fldsToOcn(1)%stdname = trim(flds_scalar_name)

       call dshr_fld_add(model_fld='Foxx_swnet', model_fld_concat=flds_x2o, model_fld_index=kswnet, &
            fldlist_num=fldsToOcn_num, fldlist=fldsToOcn)
       call dshr_fld_add(model_fld='Foxx_lwup',  model_fld_concat=flds_x2o, model_fld_index=klwup, &
            fldlist_num=fldsToOcn_num, fldlist=fldsToOcn)
       call dshr_fld_add(model_fld='Foxx_sen',   model_fld_concat=flds_x2o, model_fld_index=ksen, &
            fldlist_num=fldsToOcn_num, fldlist=fldsToOcn)
       call dshr_fld_add(model_fld='Foxx_lat',   model_fld_concat=flds_x2o, model_fld_index=klat, &
            fldlist_num=fldsToOcn_num, fldlist=fldsToOcn)
       call dshr_fld_add(model_fld='Faxa_lwdn',  model_fld_concat=flds_x2o, model_fld_index=klwdn, &
            fldlist_num=fldsToOcn_num, fldlist=fldsToOcn)
       call dshr_fld_add(model_fld='Faxa_snow',  model_fld_concat=flds_x2o, model_fld_index=ksnow, &
            fldlist_num=fldsToOcn_num, fldlist=fldsToOcn)
       call dshr_fld_add(model_fld='Fioi_melth', model_fld_concat=flds_x2o, model_fld_index=kmelth, &
            fldlist_num=fldsToOcn_num, fldlist=fldsToOcn)
       call dshr_fld_add(model_fld='Foxx_rofi', model_fld_concat=flds_x2o, model_fld_index=krofi, &
            fldlist_num=fldsToOcn_num, fldlist=fldsToOcn)
    end if

    !-------------------
    ! Advertise fields for import and export states
    !-------------------

    do n = 1,fldsFrOcn_num
       call NUOPC_Advertise(exportState, standardName=fldsFrOcn(n)%stdname, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite('(ocn_comp_nuopc):(InitializeAdvertise):Fr_ocn'//trim(fldsFrOcn(n)%stdname), &
            ESMF_LOGMSG_INFO)
    enddo

    if (ocn_prognostic) then
       do n = 1,fldsToOcn_num
          call NUOPC_Advertise(importState, standardName=fldsToOcn(n)%stdname, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_LogWrite('(ocn_comp_nuopc):(InitializeAdvertise):To_ocn'//trim(fldsToOcn(n)%stdname), &
               ESMF_LOGMSG_INFO)
       end do
    end if

    !-------------------
    ! Save flds_x2o and flds_o2x as module variables for use in debugging
    !-------------------

    flds_x2o_mod = trim(flds_x2o)
    flds_o2x_mod = trim(flds_o2x)
    ocn_prognostic_mod = ocn_prognostic

    !-------------------
    ! module character arrays stifld and stofld
    !-------------------

    ! - stifld is a character array of stream field names
    ! - stofld is a character array of data model field names that have a one-to-one correspondence with names in stifld
    ! - flds_strm is a colon delimited string of field names that is created from the field names in stofld for ONLY
    !   those field names that are available in the data streams present in SDOCN%sdatm
    ! - avstrm is an attribute vector created from flds_strm

    if (ocn_prognostic_mod) then
       call dshr_fld_add(data_fld="h"   , data_fld_array=stifld, model_fld="strm_h"   , model_fld_array=stofld)
       call dshr_fld_add(data_fld="qbot", data_fld_array=stifld, model_fld="strm_qbot", model_fld_array=stofld)
    end if

  end subroutine docn_comp_advertise

  !===============================================================================

  subroutine docn_comp_init(x2o, o2x, &
       SDOCN, mpicom, compid, my_task, master_task, &
       inst_suffix, logunit, read_restart, &
       scmMode, scmlat, scmlon, calendar, current_ymd, current_tod, modeldt, mesh)


    ! !DESCRIPTION: initialize docn model
    use pio        , only : iosystem_desc_t
    use shr_pio_mod, only : shr_pio_getiosys, shr_pio_getiotype

    ! --- input/output arguments ---
    type(mct_aVect)        , intent(inout) :: x2o, o2x       ! input/output attribute vectors
    type(shr_strdata_type) , intent(inout) :: SDOCN          ! model shr_strdata instance (output)
    integer                , intent(in)    :: mpicom         ! mpi communicator
    integer                , intent(in)    :: compid         ! mct comp id
    integer                , intent(in)    :: my_task        ! my task in mpi communicator mpicom
    integer                , intent(in)    :: master_task    ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix    ! char string associated with instance
    integer                , intent(in)    :: logunit        ! logging unit number
    logical                , intent(in)    :: read_restart   ! start from restart
    logical                , intent(in)    :: scmMode        ! single column mode
    real(R8)               , intent(in)    :: scmLat         ! single column lat
    real(R8)               , intent(in)    :: scmLon         ! single column lon
    character(len=*)       , intent(in)    :: calendar       ! model calendar type
    integer                , intent(in)    :: current_ymd    ! model date
    integer                , intent(in)    :: current_tod    ! model sec into model date
    integer                , intent(in)    :: modeldt        ! model time step
    type(ESMF_Mesh)        , intent(in)    :: mesh           ! ESMF docn mesh

    !--- local variables ---
    integer                        :: n,k      ! generic counters
    integer                        :: lsize    ! local size
    integer                        :: kfld     ! fld index
    integer                        :: cnt      ! counter
    logical                        :: exists   ! file existance
    logical                        :: exists1  ! file existance
    integer                        :: nu       ! unit number
    type(ESMF_DistGrid)            :: distGrid
    integer, allocatable, target   :: gindex(:)
    integer                        :: rc
    type(iosystem_desc_t), pointer :: ocn_pio_subsystem
    integer                        :: dimCount
    integer                        :: tileCount
    integer                        :: deCount
    integer                        :: gsize
    integer, allocatable           :: elementCountPTile(:)
    integer, allocatable           :: indexCountPDE(:,:)
    integer                        :: spatialDim
    integer                        :: numOwnedElements
    real(R8), pointer              :: ownedElemCoords(:)
    integer                        :: klat, klon
    character(*), parameter        :: F00   = "('(docn_comp_init) ',8a)"
    character(*), parameter        :: F05   = "('(docn_comp_init) ',a,2f10.4)"
    character(*), parameter        :: F06   = "('(docn_comp_init) ',a,f10.4)"
    character(*), parameter        :: subName = "(docn_comp_init) "
    !-------------------------------------------------------------------------------

    call t_startf('DOCN_INIT')

    !----------------------------------------------------------------------------
    ! Initialize pio
    !----------------------------------------------------------------------------

    call shr_strdata_pioinit(SDOCN, COMPID)

    !----------------------------------------------------------------------------
    ! Create a data model global seqmap
    !----------------------------------------------------------------------------

    call t_startf('docn_strdata_init')

    if (my_task == master_task) write(logunit,F00) ' initialize DOCN gsmap'

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
    call mct_gsMap_init( SDOCN%gsmap, gindex, mpicom, compid, lsize, gsize)
    deallocate(gindex)

    !----------------------------------------------------------------------------
    ! Initialize SDOCN model domain attributes
    !----------------------------------------------------------------------------

    ! The call to shr_strdata_init_model_domain creates the SDOCN%gsmap which
    ! is a '2d1d' decommp (1d decomp of 2d grid) and also create SDOCN%grid

    SDOCN%calendar = trim(shr_cal_calendarName(trim(calendar)))

    if (scmmode) then
       if (my_task == master_task) write(logunit,F05) ' scm lon lat = ',scmlon,scmlat
       call shr_strdata_init_model_domain(SDOCN, mpicom, compid, my_task, &
            scmmode=scmmode, scmlon=scmlon, scmlat=scmlat, gsmap=SDOCN%gsmap)
    else if (datamode == 'SST_AQUAPANAL' .or. datamode == 'SST_AQUAPFILE' .or. datamode == 'SOM_AQUAP') then
       call shr_strdata_init_model_domain(SDOCN, mpicom, compid, my_task, &
            reset_domain_mask=.true., gsmap=SDOCN%gsmap)
    else
       call shr_strdata_init_model_domain(SDOCN, mpicom, compid, my_task, gsmap=SDOCN%gsmap)
    end if

    if (my_task == master_task) then
       call shr_strdata_print(SDOCN,'SDOCN data')
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
    klon = mct_aVect_indexRA(SDOCN%grid%data,'lon')
    do n = 1, lsize
       if (abs( SDOCN%grid%data%rattr(klon,n) - xc(n)) > 1.e-4) then
          write(6,*)'ERROR: DOCN lon diff = ',abs(SDOCN%grid%data%rattr(klon,n) -  xc(n)),' too large'
          call shr_sys_abort()
       end if
       !SDOCN%grid%data%rattr(klon,n) = xc(n) ! overwrite ggrid with mesh data
       xc(n) = SDOCN%grid%data%rattr(klon,n)  ! overwrite mesh data with ggrid data
    end do
    klat = mct_aVect_indexRA(SDOCN%grid%data,'lat')
    do n = 1, lsize
       if (abs( SDOCN%grid%data%rattr(klat,n) -  yc(n)) > 1.e-4) then
          write(6,*)'ERROR: DOCN lat diff = ',abs(SDOCN%grid%data%rattr(klat,n) -  yc(n)),' too large'
          call shr_sys_abort()
       end if
       !SDOCN%grid%data%rattr(klat,n) = yc(n)
       yc(n) = SDOCN%grid%data%rattr(klat,n)
    end do

    ! determine module mask array (imask)
    allocate(imask(lsize))
    kmask = mct_aVect_indexRA(SDOCN%grid%data,'mask')
    imask(:) = nint(SDOCN%grid%data%rAttr(kmask,:))

    !----------------------------------------------------------------------------
    ! Initialize the SDOCN streams and mapping of streams to model domain
    !----------------------------------------------------------------------------

    call shr_strdata_init_streams(SDOCN, compid, mpicom, my_task)
    call shr_strdata_init_mapping(SDOCN, compid, mpicom, my_task)

    !----------------------------------------------------------------------------
    ! Allocate module arrays
    !----------------------------------------------------------------------------

    allocate(somtp(lsize))
    allocate(tfreeze(lsize))

    call t_stopf('docn_strdata_init')

    !----------------------------------------------------------------------------
    ! Initialize attribute vectors
    !----------------------------------------------------------------------------

    call t_startf('docn_initavs')
    if (my_task == master_task) write(logunit,F00) 'allocate AVs'

    call mct_aVect_init(o2x, rList=flds_o2x_mod, lsize=lsize)
    call mct_aVect_zero(o2x)

    kfrac = mct_aVect_indexRA(SDOCN%grid%data,'frac')
    o2x%rAttr(ksomask,:) = SDOCN%grid%data%rAttr(kfrac,:)

    if (ocn_prognostic_mod) then
       call mct_aVect_init(x2o, rList=flds_x2o_mod, lsize=lsize)
       call mct_aVect_zero(x2o)

       ! Initialize internal attribute vectors for optional streams
       ! Create the colon deliminted list flds_strm based on mapping the
       ! input stream fields from SDOCN%avs(n) with names in stifld to flds_strm with the names in stofld

       cnt = 0
       flds_strm = ''
       do n = 1,SDOCN%nstreams
          ! Loop over the field names in stifld
          do k = 1,size(stifld)
             ! Search input stream n for the field name stifld(k)
             kfld = mct_aVect_indexRA(SDOCN%avs(n), trim(stifld(k)), perrWith='quiet')
             if (kfld > 0) then
                cnt = cnt + 1
                ! Append the colon deliminted flds_strm with the mapped field name stofld(k)
                if (cnt == 1) then
                   flds_strm = trim(stofld(k))
                else
                   flds_strm = trim(flds_strm)//':'//trim(stofld(k))
                endif
             endif
          enddo
       enddo

       ! Initialize avstrm based on the active streams determined above
       if (my_task == master_task) write(logunit,F00) ' flds_strm = ',trim(flds_strm)
       call mct_aVect_init(avstrm, rList=flds_strm, lsize=lsize)
       call mct_aVect_zero(avstrm)

       ! Note: because the following needs to occur AFTER we determine the fields in
       ! flds_strm - the indices below CANNOT be set in the docn_comp_advertise phase

       ! Now set indices into these active streams
       kh    = mct_aVect_indexRA(avstrm,'strm_h'   , perrWith='quiet')
       kqbot = mct_aVect_indexRA(avstrm,'strm_qbot', perrWith='quiet')
    end if

    call t_stopf('docn_initavs')

    !----------------------------------------------------------------------------
    ! Read restart
    !----------------------------------------------------------------------------

    if (read_restart) then
       exists = .false.
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
       endif

       call shr_mpi_bcast(exists,mpicom,'exists')
       call shr_mpi_bcast(exists1,mpicom,'exists1')

       if (trim(datamode) == 'SOM' .or. trim(datamode) == 'SOM_AQUAP') then
          if (exists1) then
             if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file)
             call shr_pcdf_readwrite('read',SDOCN%pio_subsystem, SDOCN%io_type, &
                  trim(rest_file), mpicom, gsmap=SDOCN%gsmap, rf1=somtp, rf1n='somtp', &
                  io_format=SDOCN%io_format)
          else
             if (my_task == master_task) then
                write(logunit,F00) ' file not found, skipping ',trim(rest_file)
             end if
          endif
       endif

       if (exists) then
          if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file_strm)
          call shr_strdata_restRead(trim(rest_file_strm),SDOCN,mpicom)
       else
          if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file_strm)
       endif
    endif

    !----------------------------------------------------------------------------
    ! Set initial ocn state
    !----------------------------------------------------------------------------

    call t_adj_detailf(+2)

    call docn_comp_run(&
         x2o=x2o, &
         o2x=o2x, &
         SDOCN=SDOCN, &
         mpicom=mpicom, &
         compid=compid, &
         my_task=my_task, &
         master_task=master_task, &
         inst_suffix=inst_suffix, &
         logunit=logunit, &
         read_restart=read_restart, &
         write_restart=.false., &
         target_ymd=current_ymd, &
         target_tod=current_tod, &
         modeldt=modeldt)

    if (my_task == master_task) then
       write(logunit,F00) 'docn_comp_init done'
    end if

    call t_adj_detailf(-2)

    call t_stopf('DOCN_INIT')

  end subroutine docn_comp_init

  !===============================================================================

  subroutine docn_comp_run(x2o, o2x, &
       SDOCN, mpicom, compid, my_task, master_task, &
       inst_suffix, logunit, read_restart, write_restart, &
       target_ymd, target_tod, modeldt, case_name)

    ! !DESCRIPTION:  run method for docn model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(mct_aVect)        , intent(inout) :: x2o
    type(mct_aVect)        , intent(inout) :: o2x
    type(shr_strdata_type) , intent(inout) :: SDOCN
    integer                , intent(in)    :: mpicom        ! mpi communicator
    integer                , intent(in)    :: compid        ! mct comp id
    integer                , intent(in)    :: my_task       ! my task in mpi communicator mpicom
    integer                , intent(in)    :: master_task   ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix   ! char string associated with instance
    integer                , intent(in)    :: logunit       ! logging unit number
    logical                , intent(in)    :: read_restart  ! start from restart
    logical                , intent(in)    :: write_restart ! restart alarm is on
    integer                , intent(in)    :: target_ymd    ! model date
    integer                , intent(in)    :: target_tod    ! model sec into model date
    integer                , intent(in)    :: modeldt
    character(len=*)       , intent(in), optional :: case_name ! case name

    !--- local ---
    integer           :: n,nfld ! indices
    integer           :: lsize  ! size of attr vect
    real(R8)          :: dt     ! timestep
    integer           :: nu     ! unit number
    character(len=18) :: date_str
    character(len=CS) :: fldname

    real(R8), parameter :: &
         swp = 0.67_R8*(exp((-1._R8*shr_const_zsrflyr) /1.0_R8)) + 0.33_R8*exp((-1._R8*shr_const_zsrflyr)/17.0_R8)

    character(*), parameter :: F00   = "('(docn_comp_run) ',8a)"
    character(*), parameter :: F01   = "('(docn_comp_run) ',a, i7,2x,i5,2x,i5,2x,d21.14)"
    character(*), parameter :: F04   = "('(docn_comp_run) ',2a,2i8,'s')"
    character(*), parameter :: F0D   = "('(docn_comp_run) ',a, i7,2x,i5,2x,i5,2x,d21.14)"
    character(*), parameter :: subName = "(docn_comp_run) "
    !-------------------------------------------------------------------------------

    !--------------------
    ! Debug input
    !--------------------

    if (debug_import > 0 .and. my_task == master_task .and. ocn_prognostic_mod) then
       do nfld = 1, mct_aVect_nRAttr(x2o)
          call shr_string_listGetName(trim(flds_x2o_mod), nfld, fldname)
          do n = 1, mct_aVect_lsize(x2o)
             write(logunit,F0D)'import: ymd,tod,n  = '// trim(fldname),target_ymd, target_tod, &
                  n, x2o%rattr(nfld,n)
          end do
       end do
    end if

    call t_startf('DOCN_RUN')
    call t_barrierf('docn_BARRIER',mpicom)

    !--------------------
    ! ADVANCE OCN
    !--------------------

    call t_startf('docn')

    !--- defaults, copy all fields from streams to o2x ---

    lsize = mct_avect_lsize(o2x)
    do n = 1,lsize
       if (ksomask /= 0) then
          o2x%rAttr(ksomask, n) = SDOCN%grid%data%rAttr(kfrac,n)
       end if
       o2x%rAttr(kt   ,n) = TkFrz
       o2x%rAttr(ks   ,n) = ocnsalt
       o2x%rAttr(ku   ,n) = 0.0_R8
       o2x%rAttr(kv   ,n) = 0.0_R8
       o2x%rAttr(kdhdx,n) = 0.0_R8
       o2x%rAttr(kdhdy,n) = 0.0_R8
       o2x%rAttr(kq   ,n) = 0.0_R8
       if (kswp /= 0) then
          o2x%rAttr(kswp ,n) = swp
       end if
    enddo

    ! NOTE: for SST_AQUAPANAL, the docn buildnml sets the stream to "null"
    ! and thereby shr_strdata_advance does nothing

    call t_startf('docn_strdata_advance')
    call shr_strdata_advance(SDOCN, target_ymd, target_tod, mpicom, 'docn')
    call t_stopf('docn_strdata_advance')

    !--- copy streams to o2x ---
    call t_barrierf('docn_scatter_BARRIER', mpicom)
    call t_startf('docn_scatter')
    do n = 1, SDOCN%nstreams
       call shr_dmodel_translateAV(SDOCN%avs(n), o2x, avifld, avofld)
    enddo
    call t_stopf('docn_scatter')

    !-------------------------------------------------
    ! Determine data model behavior based on the mode
    !-------------------------------------------------

    call t_startf('docn_datamode')
    select case (trim(datamode))

    case('COPYALL')
       ! do nothing extra

    case('SSTDATA')
       lsize = mct_avect_lsize(o2x)
       do n = 1,lsize
          o2x%rAttr(kt   ,n) = o2x%rAttr(kt,n) + TkFrz
          o2x%rAttr(ks   ,n) = ocnsalt
          o2x%rAttr(ku   ,n) = 0.0_R8
          o2x%rAttr(kv   ,n) = 0.0_R8
          o2x%rAttr(kdhdx,n) = 0.0_R8
          o2x%rAttr(kdhdy,n) = 0.0_R8
          o2x%rAttr(kq   ,n) = 0.0_R8
          if (kswp /= 0) then
             o2x%rAttr(kswp ,n) = swp
          end if
       enddo

    case('SST_AQUAPANAL')
       lsize = mct_avect_lsize(o2x)
       ! Zero out the attribute vector before calling the prescribed_sst
       ! function - so this also zeroes out the So_omask if it is needed
       ! so need to re-introduce it
       do n = 1,lsize
          o2x%rAttr(:,n) = 0.0_r8
       end do
       call prescribed_sst(xc, yc, lsize, aquap_option, o2x%rAttr(kt,:))
       do n = 1,lsize
          o2x%rAttr(kt,n) = o2x%rAttr(kt,n) + TkFrz
          if (ksomask /= 0) then
             o2x%rAttr(ksomask, n) = SDOCN%grid%data%rAttr(kfrac,n)
          end if
       enddo

    case('SST_AQUAPFILE')
       lsize = mct_avect_lsize(o2x)
       do n = 1,lsize
          o2x%rAttr(kt   ,n) = o2x%rAttr(kt,n) + TkFrz
          o2x%rAttr(ks   ,n) = ocnsalt
          o2x%rAttr(ku   ,n) = 0.0_R8
          o2x%rAttr(kv   ,n) = 0.0_R8
          o2x%rAttr(kdhdx,n) = 0.0_R8
          o2x%rAttr(kdhdy,n) = 0.0_R8
          o2x%rAttr(kq   ,n) = 0.0_R8
          if (kswp /= 0) then
             o2x%rAttr(kswp ,n) = swp
          end if
       enddo

    case('IAF')
       lsize = mct_avect_lsize(o2x)
       do n = 1,lsize
          o2x%rAttr(kt   ,n) = o2x%rAttr(kt,n) + TkFrz
          o2x%rAttr(ks   ,n) = ocnsalt
          o2x%rAttr(ku   ,n) = 0.0_R8
          o2x%rAttr(kv   ,n) = 0.0_R8
          o2x%rAttr(kdhdx,n) = 0.0_R8
          o2x%rAttr(kdhdy,n) = 0.0_R8
          o2x%rAttr(kq   ,n) = 0.0_R8
          if (kswp /= 0) then
             o2x%rAttr(kswp ,n) = swp
          end if
       enddo

    case('SOM')
       lsize = mct_avect_lsize(o2x)
       do n = 1,SDOCN%nstreams
          call shr_dmodel_translateAV(SDOCN%avs(n),avstrm,stifld,stofld)
       enddo
       if (firstcall) then
          do n = 1,lsize
             if (.not. read_restart) then
                somtp(n) = o2x%rAttr(kt,n) + TkFrz
             endif
             o2x%rAttr(kt,n) = somtp(n)
             o2x%rAttr(kq,n) = 0.0_R8
          enddo
       else   ! firstcall
          tfreeze = shr_frz_freezetemp(o2x%rAttr(ks,:)) + TkFrz
          dt = modeldt * 1.0_R8
          do n = 1,lsize
             if (imask(n) /= 0) then
                !--- compute new temp ---
                o2x%rAttr(kt,n) = somtp(n) + &
                     (x2o%rAttr(kswnet,n) + &  ! shortwave
                      x2o%rAttr(klwup ,n) + &  ! longwave
                      x2o%rAttr(klwdn ,n) + &  ! longwave
                      x2o%rAttr(ksen  ,n) + &  ! sensible
                      x2o%rAttr(klat  ,n) + &  ! latent
                      x2o%rAttr(kmelth,n) - &  ! ice melt
                      avstrm%rAttr(kqbot ,n) - &  ! flux at bottom
                     (x2o%rAttr(ksnow,n)+x2o%rAttr(krofi,n))*latice) * &  ! latent by prec and roff
                     dt/(cpsw*rhosw* avstrm%rAttr(kh,n))
                !--- compute ice formed or melt potential ---
                o2x%rAttr(kq,n) = (tfreeze(n) - o2x%rAttr(kt,n))*(cpsw*rhosw*avstrm%rAttr(kh,n))/dt  ! ice formed q>0
                o2x%rAttr(kt,n) = max(tfreeze(n),o2x%rAttr(kt,n))                    ! reset temp
                somtp(n) = o2x%rAttr(kt,n)                                           ! save temp
             endif
          end do
       endif   ! firstcall

    case('SOM_AQUAP')
       lsize = mct_avect_lsize(o2x)
       do n = 1,SDOCN%nstreams
          call shr_dmodel_translateAV(SDOCN%avs(n),avstrm,stifld,stofld)
       enddo
       if (firstcall) then
          do n = 1,lsize
             if (.not. read_restart) then
                somtp(n) = o2x%rAttr(kt,n) + TkFrz
             endif
             o2x%rAttr(kt,n) = somtp(n)
             o2x%rAttr(kq,n) = 0.0_R8
          enddo
       else   ! firstcall
          tfreeze = shr_frz_freezetemp(o2x%rAttr(ks,:)) + TkFrz
          do n = 1,lsize
             !--- compute new temp ---
             o2x%rAttr(kt,n) = somtp(n) + &
                  (x2o%rAttr(kswnet,n) + &  ! shortwave
                  x2o%rAttr(klwup ,n) + &  ! longwave
                  x2o%rAttr(klwdn ,n) + &  ! longwave
                  x2o%rAttr(ksen  ,n) + &  ! sensible
                  x2o%rAttr(klat  ,n) + &  ! latent
                  x2o%rAttr(kmelth,n) - &  ! ice melt
                  avstrm%rAttr(kqbot ,n) - &  ! flux at bottom
                  (x2o%rAttr(ksnow,n)+x2o%rAttr(krofi,n))*latice) * &  ! latent by prec and roff
                  dt/(cpsw*rhosw*avstrm%rAttr(kh,n))
             !--- compute ice formed or melt potential ---
             o2x%rAttr(kq,n) = (tfreeze(n) - o2x%rAttr(kt,n))*(cpsw*rhosw*avstrm%rAttr(kh,n))/dt  ! ice formed q>0
             somtp(n) = o2x%rAttr(kt,n)                                        ! save temp
          enddo
       endif   ! firstcall

    end select

    call t_stopf('docn_datamode')

    !--------------------
    ! Debug output
    !--------------------

    if (debug_export > 1 .and. my_task == master_task) then
       do nfld = 1, mct_aVect_nRAttr(o2x)
          call shr_string_listGetName(trim(flds_o2x_mod), nfld, fldname)
          do n = 1, mct_aVect_lsize(o2x)
             write(logunit,F0D)'export: ymd,tod,n  = '// trim(fldname),target_ymd, target_tod, &
                  n, o2x%rattr(nfld,n)
          end do
       end do
    end if

    !--------------------
    ! Write restart
    !--------------------

    if (write_restart) then
       call t_startf('docn_restart')
       call shr_cal_datetod2string(date_str, target_ymd, target_tod)
       write(rest_file,"(6a)") &
            trim(case_name), '.docn',trim(inst_suffix),'.r.', &
            trim(date_str),'.nc'
       write(rest_file_strm,"(6a)") &
            trim(case_name), '.docn',trim(inst_suffix),'.rs1.', &
            trim(date_str),'.bin'
       if (my_task == master_task) then
          nu = shr_file_getUnit()
          open(nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
          write(nu,'(a)') rest_file
          write(nu,'(a)') rest_file_strm
          close(nu)
          call shr_file_freeUnit(nu)
       endif
       if (trim(datamode) == 'SOM' .or. trim(datamode) == 'SOM_AQUAP') then
          if (my_task == master_task) then
             write(logunit,F04) ' writing ',trim(rest_file),target_ymd,target_tod
          end if
          call shr_pcdf_readwrite('write', SDOCN%pio_subsystem, SDOCN%io_type,&
               trim(rest_file), mpicom, SDOCN%gsmap, clobber=.true., rf1=somtp,rf1n='somtp')
       endif
       if (my_task == master_task) then
          write(logunit,F04) ' writing ',trim(rest_file_strm),target_ymd,target_tod
       end if
       call shr_strdata_restWrite(trim(rest_file_strm), SDOCN, mpicom, trim(case_name), 'SDOCN strdata')
       call t_stopf('docn_restart')
    endif

    firstcall = .false.

    call t_stopf('docn')
    call t_stopf('DOCN_RUN')

  end subroutine docn_comp_run

  !===============================================================================

  subroutine prescribed_sst(xc, yc, lsize, sst_option, sst)

    real(R8)     , intent(in)    :: xc(:)  !degrees
    real(R8)     , intent(in)    :: yc(:)  !degrees
    integer      , intent(in)    :: lsize
    integer      , intent(in)    :: sst_option
    real(R8)     , intent(inout) :: sst(:)

    ! local
    integer  :: i
    real(r8) :: tmp, tmp1, pi
    real(r8) :: rlon(lsize), rlat(lsize)

    real(r8), parameter :: pio180 = SHR_CONST_PI/180._r8

    ! Parameters for zonally symmetric experiments
    real(r8), parameter ::   t0_max     = 27._r8
    real(r8), parameter ::   t0_min     = 0._r8
    real(r8), parameter ::   maxlat     = 60._r8*pio180
    real(r8), parameter ::   shift      = 5._r8*pio180
    real(r8), parameter ::   shift9     = 10._r8*pio180
    real(r8), parameter ::   shift10    = 15._r8*pio180

    ! Parameters for zonally asymmetric experiments
    real(r8), parameter ::   t0_max6    = 1._r8
    real(r8), parameter ::   t0_max7    = 3._r8
    real(r8), parameter ::   latcen     = 0._r8*pio180
    real(r8), parameter ::   loncen     = 0._r8*pio180
    real(r8), parameter ::   latrad6    = 15._r8*pio180
    real(r8), parameter ::   latrad8    = 30._r8*pio180
    real(r8), parameter ::   lonrad     = 30._r8*pio180
    !-------------------------------------------------------------------------------

    pi = SHR_CONST_PI

    ! convert xc and yc from degrees to radians

    rlon(:) = xc(:) * pio180
    rlat(:) = yc(:) * pio180

    ! Control

    if (sst_option < 1 .or. sst_option > 10) then
       call shr_sys_abort ('prescribed_sst: ERROR: sst_option must be between 1 and 10')
    end if

    if (sst_option == 1 .or. sst_option == 6 .or. sst_option == 7 .or. sst_option == 8) then
       do i = 1,lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else
             tmp = sin(rlat(i)*pi*0.5_r8/maxlat)
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if

    ! Flat

    if (sst_option == 2) then
       do i = 1,lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else
             tmp = sin(rlat(i)*pi*0.5_r8/maxlat)
             tmp = 1._r8 - tmp*tmp*tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if

    ! Qobs

    if (sst_option == 3) then
       do i = 1,lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else
             tmp = sin(rlat(i)*pi*0.5_r8/maxlat)
             tmp = (2._r8 - tmp*tmp*tmp*tmp - tmp*tmp)*0.5_r8
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if

    ! Peaked

    if (sst_option == 4) then
       do i = 1,lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else
             tmp = (maxlat - abs(rlat(i)))/maxlat
             tmp1 = 1._r8 - tmp
             sst(i) = t0_max*tmp + t0_min*tmp1
          end if
       end do
    end if

    ! Control-5N

    if (sst_option == 5) then
       do i = 1,lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else if (rlat(i) > shift) then
             tmp = sin((rlat(i)-shift)*pi*0.5_r8/(maxlat-shift))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          else
             tmp = sin((rlat(i)-shift)*pi*0.5_r8/(maxlat+shift))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if

    ! 1KEQ

    if (sst_option == 6) then
       do i = 1,lsize
          if (abs(rlat(i)-latcen) <= latrad6) then
             tmp1 = cos((rlat(i)-latcen)*pi*0.5_r8/latrad6)
             tmp1 = tmp1*tmp1
             tmp = abs(rlon(i)-loncen)
             tmp = min(tmp , 2._r8*pi-tmp)
             if(tmp <= lonrad) then
                tmp = cos(tmp*pi*0.5_r8/lonrad)
                tmp = tmp*tmp
                sst(i) = sst(i) + t0_max6*tmp*tmp1
             end if
          end if
       end do
    end if

    ! 3KEQ

    if (sst_option == 7) then
       do i = 1, lsize
          if (abs(rlat(i)-latcen) <= latrad6) then
             tmp1 = cos((rlat(i)-latcen)*pi*0.5_r8/latrad6)
             tmp1 = tmp1*tmp1
             tmp = abs(rlon(i)-loncen)
             tmp = min(tmp , 2._r8*pi-tmp)
             if (tmp <= lonrad) then
                tmp = cos(tmp*pi*0.5_r8/lonrad)
                tmp = tmp*tmp
                sst(i) = sst(i) + t0_max7*tmp*tmp1
             end if
          end if
       end do
    end if

    ! 3KW1

    if (sst_option == 8) then
       do i = 1, lsize
          if (abs(rlat(i)-latcen) <= latrad8) then
             tmp1 = cos((rlat(i)-latcen)*pi*0.5_r8/latrad8)
             tmp1 = tmp1*tmp1
             tmp = cos(rlon(i)-loncen)
             sst(i) = sst(i) + t0_max7*tmp*tmp1
          end if
       end do
    end if

    ! Control-10N

    if (sst_option == 9) then
       do i = 1, lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else if (rlat(i) > shift9) then
             tmp = sin((rlat(i)-shift9)*pi*0.5_r8/(maxlat-shift9))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          else
             tmp = sin((rlat(i)-shift9)*pi*0.5_r8/(maxlat+shift9))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if

    ! Control-15N

    if (sst_option == 10) then
       do i = 1, lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else if(rlat(i) > shift10) then
             tmp = sin((rlat(i)-shift10)*pi*0.5_r8/(maxlat-shift10))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          else
             tmp = sin((rlat(i)-shift10)*pi*0.5_r8/(maxlat+shift10))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if

  end subroutine prescribed_sst

end module docn_comp_mod
