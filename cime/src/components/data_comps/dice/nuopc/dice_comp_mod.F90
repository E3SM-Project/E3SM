#ifdef AIX
@PROCESS ALIAS_SIZE(805306368)
#endif

module dice_comp_mod

  ! !USES:
  use NUOPC                 , only : NUOPC_Advertise
  use ESMF                  , only : ESMF_State, ESMF_SUCCESS, ESMF_State
  use ESMF                  , only : ESMF_Mesh, ESMF_DistGrid, ESMF_MeshGet, ESMF_DistGridGet
  use perf_mod              , only : t_startf, t_stopf, t_adj_detailf, t_barrierf
  use mct_mod               , only : mct_gsmap_init
  use mct_mod               , only : mct_avect, mct_avect_indexRA, mct_avect_zero, mct_aVect_nRattr
  use mct_mod               , only : mct_avect_init, mct_avect_lsize
  use med_constants_mod     , only : R8, CS, CXX
  use shr_const_mod         , only : shr_const_pi, shr_const_spval, shr_const_tkfrz, shr_const_latice
  use shr_file_mod          , only : shr_file_getunit, shr_file_freeunit
  use shr_mpi_mod           , only : shr_mpi_bcast
  use shr_frz_mod           , only : shr_frz_freezetemp
  use shr_cal_mod           , only : shr_cal_calendarname
  use shr_cal_mod           , only : shr_cal_datetod2string
  use shr_string_mod        , only : shr_string_listGetName
  use shr_sys_mod           , only : shr_sys_abort
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
  use dice_shr_mod          , only : datamode       ! namelist input
  use dice_shr_mod          , only : rest_file      ! namelist input
  use dice_shr_mod          , only : rest_file_strm ! namelist input
  use dice_shr_mod          , only : flux_swpf      ! namelist input -short-wave penatration factor
  use dice_shr_mod          , only : flux_Qmin      ! namelist input -bound on melt rate
  use dice_shr_mod          , only : flux_Qacc      ! namelist input -activates water accumulation/melt wrt Q
  use dice_shr_mod          , only : flux_Qacc0     ! namelist input -initial water accumulation value
  use dice_shr_mod          , only : nullstr
  use dice_flux_atmice_mod  , only : dice_flux_atmice
  use shr_pcdf_mod

  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: dice_comp_advertise
  public :: dice_comp_init
  public :: dice_comp_run

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  integer                    :: debug_import = 0      ! debug level (if > 0 will print all import fields)
  integer                    :: debug_export = 0      ! debug level (if > 0 will print all export fields)

  real(R8),parameter         :: pi     = shr_const_pi      ! pi
  real(R8),parameter         :: spval  = shr_const_spval   ! flags invalid data
  real(R8),parameter         :: tFrz   = shr_const_tkfrz   ! temp of freezing
  real(R8),parameter         :: latice = shr_const_latice  ! latent heat of fusion
  real(R8),parameter         :: waterMax = 1000.0_R8       ! wrt iFrac comp & frazil ice (kg/m^2)

  !----- surface albedo constants ------
  real(R8),parameter         :: snwfrac = 0.286_R8 ! snow cover fraction ~ [0,1]
  real(R8),parameter         :: as_nidf = 0.950_R8 ! albedo: snow,near-infr,diffuse
  real(R8),parameter         :: as_vsdf = 0.700_R8 ! albedo: snow,visible  ,diffuse
  real(R8),parameter         :: as_nidr = 0.960_R8 ! albedo: snow,near-infr,direct
  real(R8),parameter         :: as_vsdr = 0.800_R8 ! albedo: snow,visible  ,direct
  real(R8),parameter         :: ai_nidf = 0.700_R8 ! albedo: ice, near-infr,diffuse
  real(R8),parameter         :: ai_vsdf = 0.500_R8 ! albedo: ice, visible  ,diffuse
  real(R8),parameter         :: ai_nidr = 0.700_R8 ! albedo: ice, near-infr,direct
  real(R8),parameter         :: ai_vsdr = 0.500_R8 ! albedo: ice, visible  ,direct
  real(R8),parameter         :: ax_nidf = ai_nidf*(1.0_R8-snwfrac) + as_nidf*snwfrac
  real(R8),parameter         :: ax_vsdf = ai_vsdf*(1.0_R8-snwfrac) + as_vsdf*snwfrac
  real(R8),parameter         :: ax_nidr = ai_nidr*(1.0_R8-snwfrac) + as_nidr*snwfrac
  real(R8),parameter         :: ax_vsdr = ai_vsdr*(1.0_R8-snwfrac) + as_vsdr*snwfrac

  integer                    :: km
  integer                    :: kswvdr,kswndr,kswvdf,kswndf,kq,kz,kua,kva,kptem,kshum,kdens,ktbot
  integer                    :: kiFrac,kt,kavsdr,kanidr,kavsdf,kanidf,kswnet,kmelth,kmeltw
  integer                    :: ksen,klat,klwup,kevap,ktauxa,ktauya,ktref,kqref,kswpen,ktauxo,ktauyo,ksalt
  integer                    :: ksalinity
  integer                    :: kbcpho, kbcphi, kflxdst
  integer                    :: kbcphidry, kbcphodry, kbcphiwet
  integer                    :: kocphidry, kocphodry, kocphiwet
  integer                    :: kdstdry1, kdstdry2, kdstdry3, kdstdry4
  integer                    :: kdstwet1, kdstwet2, kdstwet3, kdstwet4
  integer                    :: kiFrac_01,kswpen_iFrac_01 ! optional per thickness category fields
  integer                    :: index_lat, index_lon

  integer     , pointer      :: imask(:)
  real(R8)    , pointer      :: xc(:), yc(:)       ! arrays of model latitudes and longitudes
  real(R8)    , pointer      :: water(:)
  real(R8)    , pointer      :: tfreeze(:)
  !real(R8)   , pointer      :: ifrac0(:)

  character(len=CS), pointer :: avifld(:)
  character(len=CS), pointer :: avofld(:)
  character(len=CS), pointer :: strmifld(:)
  character(len=CS), pointer :: strmofld(:)
  character(len=CXX)         :: flds_strm = ''   ! colon deliminated string of field names
  character(len=CXX)         :: flds_i2x_mod
  character(len=CXX)         :: flds_x2i_mod

  logical                    :: firstcall = .true. ! first call logical
  character(len=*),parameter :: rpfile = 'rpointer.ice'
  character(*),parameter     :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine dice_comp_advertise(importState, exportState, &
       ice_present, ice_prognostic,  &
       fldsFrIce_num, fldsFrIce, fldsToIce_num, fldsToIce, &
       flds_i2x, flds_x2i, rc)

    ! input/output arguments
    type(ESMF_State)     , intent(inout) :: importState
    type(ESMF_State)     , intent(inout) :: exportState
    logical              , intent(in)    :: ice_present
    logical              , intent(in)    :: ice_prognostic
    integer              , intent(out)   :: fldsToIce_num
    integer              , intent(out)   :: fldsFrIce_num
    type (fld_list_type) , intent(out)   :: fldsToIce(:)
    type (fld_list_type) , intent(out)   :: fldsFrIce(:)
    character(len=*)     , intent(out)   :: flds_i2x
    character(len=*)     , intent(out)   :: flds_x2i
    integer              , intent(out)   :: rc

    ! local variables
    integer         :: n
    !-------------------------------------------------------------------------------

    if (.not. ice_present) return

    !--------------------------------
    ! export fields
    !--------------------------------

    fldsFrIce_num=1
    fldsFrIce(1)%stdname = trim(flds_scalar_name)

    ! export fields that have a corresponding stream field

    call dshr_fld_add(data_fld='ifrac', data_fld_array=avifld, model_fld='Si_ifrac', model_fld_array=avofld, &
         model_fld_concat=flds_i2x, model_fld_index=kiFrac, fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    ! export fields that have no corresponding stream field (computed internally)

    call dshr_fld_add(model_fld='Si_imask', model_fld_concat=flds_i2x, model_fld_index=km, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Si_t', model_fld_concat=flds_i2x, model_fld_index=kt, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Si_tref', model_fld_concat=flds_i2x, model_fld_index=ktref, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Si_qref', model_fld_concat=flds_i2x, model_fld_index=kqref, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Si_avsdr', model_fld_concat=flds_i2x, model_fld_index=kavsdr, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Si_anidr', model_fld_concat=flds_i2x, model_fld_index=kanidr, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Si_avsdf', model_fld_concat=flds_i2x, model_fld_index=kavsdf, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Si_anidf', model_fld_concat=flds_i2x, model_fld_index=kanidf, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Faii_swnet', model_fld_concat=flds_i2x, model_fld_index=kswnet, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Faii_sen', model_fld_concat=flds_i2x, model_fld_index=ksen, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Faii_lat', model_fld_concat=flds_i2x, model_fld_index=klat, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Faii_lwup', model_fld_concat=flds_i2x, model_fld_index=klwup, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Faii_evap', model_fld_concat=flds_i2x, model_fld_index=kevap, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Faii_taux', model_fld_concat=flds_i2x, model_fld_index=ktauxa, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Faii_tauy', model_fld_concat=flds_i2x, model_fld_index=ktauya, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Fioi_melth', model_fld_concat=flds_i2x, model_fld_index=kmelth, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Fioi_meltw', model_fld_concat=flds_i2x, model_fld_index=kmeltw, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Fioi_swpen', model_fld_concat=flds_i2x, model_fld_index=kswpen, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Fioi_taux', model_fld_concat=flds_i2x, model_fld_index=ktauxo, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Fioi_tauy', model_fld_concat=flds_i2x, model_fld_index=ktauyo, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Fioi_salt', model_fld_concat=flds_i2x, model_fld_index=ksalt, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Fioi_bcpho', model_fld_concat=flds_i2x, model_fld_index=kbcpho, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Fioi_bcphi', model_fld_concat=flds_i2x, model_fld_index=kbcphi, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    call dshr_fld_add(model_fld='Fioi_flxdst', model_fld_concat=flds_i2x, model_fld_index=kflxdst, &
         fldlist_num=fldsFrIce_num, fldlist=fldsFrIce)

    !-------------------
    ! import fields (have no corresponding stream fields)
    !-------------------

    if (ice_prognostic) then

       fldsToIce_num=1
       fldsToIce(1)%stdname = trim(flds_scalar_name)

       call dshr_fld_add(model_fld='Faxa_swvdr', model_fld_concat=flds_x2i, model_fld_index=kswvdr, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Faxa_swvdf', model_fld_concat=flds_x2i, model_fld_index=kswvdf, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Faxa_swndr', model_fld_concat=flds_x2i, model_fld_index=kswndr, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Faxa_swndf', model_fld_concat=flds_x2i, model_fld_index=kswndf, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Fioo_q', model_fld_concat=flds_x2i, model_fld_index=kq, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Sa_z', model_fld_concat=flds_x2i, model_fld_index=kz, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Sa_u', model_fld_concat=flds_x2i, model_fld_index=kua, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Sa_v', model_fld_concat=flds_x2i, model_fld_index=kva, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Sa_ptem', model_fld_concat=flds_x2i, model_fld_index=kptem, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Sa_shum', model_fld_concat=flds_x2i, model_fld_index=kshum, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Sa_dens', model_fld_concat=flds_x2i, model_fld_index=kdens, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Sa_tbot', model_fld_concat=flds_x2i, model_fld_index=ktbot, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='So_s', model_fld_concat=flds_x2i, model_fld_index=ksalinity, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Faxa_bcphidry', model_fld_concat=flds_x2i, model_fld_index=kbcphidry, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Faxa_bcphodry', model_fld_concat=flds_x2i, model_fld_index=kbcphodry, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Faxa_bcphiwet', model_fld_concat=flds_x2i, model_fld_index=kbcphiwet, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Faxa_ocphidry', model_fld_concat=flds_x2i, model_fld_index=kocphidry, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Faxa_ocphodry', model_fld_concat=flds_x2i, model_fld_index=kocphodry, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Faxa_ocphiwet', model_fld_concat=flds_x2i, model_fld_index=kocphiwet, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Faxa_dstdry1', model_fld_concat=flds_x2i, model_fld_index=kdstdry1, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Faxa_dstdry2', model_fld_concat=flds_x2i, model_fld_index=kdstdry2, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Faxa_dstdry3', model_fld_concat=flds_x2i, model_fld_index=kdstdry3, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Faxa_dstdry4', model_fld_concat=flds_x2i, model_fld_index=kdstdry4, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Faxa_dstwet1', model_fld_concat=flds_x2i, model_fld_index=kdstwet1, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Faxa_dstwet2', model_fld_concat=flds_x2i, model_fld_index=kdstwet2, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Faxa_dstwet3', model_fld_concat=flds_x2i, model_fld_index=kdstwet3, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

       call dshr_fld_add(model_fld='Faxa_dstwet4', model_fld_concat=flds_x2i, model_fld_index=kdstwet4, &
            fldlist_num=fldsToIce_num, fldlist=fldsToIce)

    end if

    do n = 1,fldsFrIce_num
       call NUOPC_Advertise(exportState, standardName=fldsFrIce(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    if (ice_prognostic) then
       do n = 1,fldsToIce_num
          call NUOPC_Advertise(importState, standardName=fldsToIce(n)%stdname, &
               TransferOfferGeomObject='will provide', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       enddo
    end if

    ! Save flds_x2i and flds_i2x as module variables for use in debugging

    flds_x2i_mod = trim(flds_x2i)
    flds_i2x_mod = trim(flds_i2x)

  end subroutine dice_comp_advertise

  !===============================================================================

  subroutine dice_comp_init(x2i, i2x, &
       flds_x2i_fields, flds_i2x_fields, flds_i2o_per_cat, &
       SDICE, mpicom, compid, my_task, master_task, &
       inst_suffix, inst_name, logunit, read_restart, &
       scmMode, scmlat, scmlon, calendar, mesh)

    ! !DESCRIPTION: initialize dice model

    ! input/output parameters:
    type(mct_aVect)        , intent(inout) :: x2i, i2x         ! input/output attribute vectors
    character(len=*)       , intent(in)    :: flds_x2i_fields  ! fields from mediator
    character(len=*)       , intent(in)    :: flds_i2x_fields  ! fields to mediator
    logical                , intent(in)    :: flds_i2o_per_cat ! .true. if select per ice thickness fields from ice
    type(shr_strdata_type) , intent(inout) :: SDICE            ! dice shr_strdata instance (output)
    integer                , intent(in)    :: mpicom           ! mpi communicator
    integer                , intent(in)    :: compid           ! mct comp id
    integer                , intent(in)    :: my_task          ! my task in mpi communicator mpicom
    integer                , intent(in)    :: master_task      ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix      ! char string associated with instance
    character(len=*)       , intent(in)    :: inst_name        ! fullname of current instance (ie. "lnd_0001")
    integer                , intent(in)    :: logunit          ! logging unit number
    logical                , intent(in)    :: read_restart     ! start from restart
    logical                , intent(in)    :: scmMode          ! single column mode
    real(R8)               , intent(in)    :: scmLat           ! single column lat
    real(R8)               , intent(in)    :: scmLon           ! single column lon
    character(len=*)       , intent(in)    :: calendar         ! calendar type
    type(ESMF_Mesh)        , intent(in)    :: mesh             ! ESMF dice mesh

    !--- local variables ---
    integer                      :: n,k            ! generic counters
    integer                      :: ierr           ! error code
    integer                      :: lsize          ! local size
    integer                      :: kfld           ! field reference
    logical                      :: exists,exists1 ! file existance logical
    integer                      :: nu             ! unit number
    type(ESMF_DistGrid)          :: distGrid
    integer, allocatable, target :: gindex(:)
    integer                      :: rc
    integer                      :: dimCount
    integer                      :: tileCount
    integer                      :: deCount
    integer                      :: gsize
    integer, allocatable         :: elementCountPTile(:)
    integer, allocatable         :: indexCountPDE(:,:)
    integer                      :: spatialDim
    integer                      :: numOwnedElements
    real(R8), pointer            :: ownedElemCoords(:)
    character(*), parameter      :: F00   = "('(dice_comp_init) ',8a)"
    character(*), parameter      :: F01   = "('(dice_comp_init) ',a,2f10.4)"
    character(*), parameter      :: subName = "(dice_comp_init) "
    !-------------------------------------------------------------------------------

    call t_startf('DICE_INIT')

    !----------------------------------------------------------------------------
    ! Initialize PIO
    !----------------------------------------------------------------------------

    call shr_strdata_pioinit(SDICE, compid)

    !----------------------------------------------------------------------------
    ! Create a data model global segmap
    !----------------------------------------------------------------------------

    call t_startf('dice_strdata_init')

    if (my_task == master_task) write(logunit,F00) ' initialize SDICE gsmap'

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
    call mct_gsMap_init( SDICE%gsmap, gindex, mpicom, compid, lsize, gsize)
    deallocate(gindex)

    !----------------------------------------------------------------------------
    ! Initialize SDICE
    !----------------------------------------------------------------------------

    ! The call to shr_strdata_init_model_domain creates the SDICE%gsmap which
    ! is a '2d1d' decommp (1d decomp of 2d grid) and also create SDICE%grid

    SDICE%calendar = trim(shr_cal_calendarName(trim(calendar)))

    if (scmmode) then
       if (my_task == master_task) write(logunit,F01) ' scm lon lat = ',scmlon,scmlat
       call shr_strdata_init_model_domain(SDICE, mpicom, compid, my_task, &
            scmmode=scmmode, scmlon=scmlon, scmlat=scmlat, gsmap=SDICE%gsmap)
    else
       call shr_strdata_init_model_domain(SDICE, mpicom, compid, my_task, gsmap=SDICE%gsmap)
    end if

    if (my_task == master_task) then
       call shr_strdata_print(SDICE,'SDICE data')
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
    index_lon = mct_aVect_indexRA(SDICE%grid%data,'lon')
    do n = 1, lsize
       if (abs( SDICE%grid%data%rattr(index_lon,n) - xc(n)) > 1.e-4) then
          write(6,*)'ERROR: lon diff = ',abs(SDICE%grid%data%rattr(index_lon,n) -  xc(n)),' too large'
          call shr_sys_abort()
       end if
       !SDICE%grid%data%rattr(index_lon,n) = xc(n) ! overwrite ggrid with mesh data
    end do
    index_lat = mct_aVect_indexRA(SDICE%grid%data,'lat')
    do n = 1, lsize
       if (abs( SDICE%grid%data%rattr(index_lat,n) -  yc(n)) > 1.e-4) then
          write(6,*)'ERROR: lat diff = ',abs(SDICE%grid%data%rattr(index_lat,n) -  yc(n)),' too large'
          call shr_sys_abort()
       end if
       !SDICE%grid%data%rattr(index_lat,n) = yc(n) ! overwrite ggrid with mesh data
    end do

    ! Note that the module array, imask, does not change after initialization
    allocate(imask(lsize))
    kfld = mct_aVect_indexRA(SDICE%grid%data,'mask')
    imask(:) = nint(SDICE%grid%data%rAttr(kfld,:))

    if (my_task == master_task) then
       call shr_strdata_print(SDICE,'SDICE data')
    endif

    call t_stopf('dice_strdata_init')

    !----------------------------------------------------------------------------
    ! Initialize SDICE attributes for streams and mapping of streams to model domain
    !----------------------------------------------------------------------------

    call shr_strdata_init_streams(SDICE, compid, mpicom, my_task)
    call shr_strdata_init_mapping(SDICE, compid, mpicom, my_task)

    !----------------------------------------------------------------------------
    ! Initialize MCT attribute vectors
    !----------------------------------------------------------------------------

    call t_startf('dice_initmctavs')
    if (my_task == master_task) write(logunit,F00) 'allocate AVs'

    call mct_aVect_init(i2x, rList=flds_i2x_fields, lsize=lsize)
    call mct_aVect_zero(i2x)

    ! optional per thickness category fields
    if (flds_i2o_per_cat) then
       kiFrac_01       = mct_aVect_indexRA(i2x,'Si_ifrac_01')
       kswpen_iFrac_01 = mct_aVect_indexRA(i2x,'PFioi_swpen_ifrac_01')
    end if

    call mct_aVect_init(x2i, rList=flds_x2i_fields, lsize=lsize)
    call mct_aVect_zero(x2i)

    allocate(water(lsize))
    allocate(tfreeze(lsize))
    ! allocate(iFrac0(lsize))

    if (km /= 0) then
       i2x%rAttr(km, :) = imask(:)
    end if

    call t_stopf('dice_initmctavs')

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

       if (exists1) then
          if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file)
          call shr_pcdf_readwrite('read',SDICE%pio_subsystem, SDICE%io_type, &
               trim(rest_file), mpicom, gsmap=SDICE%gsmap, rf1=water, rf1n='water', io_format=SDICE%io_format)
       else
          if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file)
       endif

       if (exists) then
          if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file_strm)
          call shr_strdata_restRead(trim(rest_file_strm),SDICE,mpicom)
       else
          if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file_strm)
       endif
    endif

    !----------------------------------------------------------------------------
    ! On initial call, x2i is unset, so set for use in run method
    !  These values should have no impact on the solution!!
    !----------------------------------------------------------------------------

    x2i%rAttr(kz,:)    = 10.0_R8
    x2i%rAttr(kua,:)   = 5.0_R8
    x2i%rAttr(kva,:)   = 5.0_R8
    x2i%rAttr(kptem,:) = 260.0_R8
    x2i%rAttr(ktbot,:) = 260.0_R8
    x2i%rAttr(kshum,:) = 0.0014_R8
    x2i%rAttr(kdens,:) = 1.3_R8

    call t_stopf('DICE_INIT')

  end subroutine dice_comp_init

  !===============================================================================

  subroutine dice_comp_run(x2i, i2x, flds_i2o_per_cat, &
       SDICE, mpicom, my_task, master_task, &
       inst_suffix, logunit, read_restart, write_restart, &
       calendar, modeldt, target_ymd, target_tod, cosArg, case_name )

    ! !DESCRIPTION: run method for dice model

    ! input/output parameters:
    type(mct_aVect)        , intent(inout) :: x2i
    type(mct_aVect)        , intent(inout) :: i2x
    logical                , intent(in)    :: flds_i2o_per_cat     ! .true. if select per ice thickness fields from ice
    type(shr_strdata_type) , intent(inout) :: SDICE
    integer                , intent(in)    :: mpicom               ! mpi communicator
    integer                , intent(in)    :: my_task              ! my task in mpi communicator mpicom
    integer                , intent(in)    :: master_task          ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix          ! char string associated with instance
    integer                , intent(in)    :: logunit              ! logging unit number
    logical                , intent(in)    :: read_restart         ! start from restart
    logical                , intent(in)    :: write_restart        ! restart now
    character(len=*)       , intent(in)    :: calendar
    integer                , intent(in)    :: modeldt
    integer                , intent(in)    :: target_ymd
    integer                , intent(in)    :: target_tod
    real(R8)               , intent(in)    :: cosarg               ! for setting ice temp pattern
    character(len=*)       , intent(in), optional :: case_name     ! case name

    !--- local ---
    integer           :: n,nfld            ! indices
    integer           :: lsize             ! size of attr vect
    real(R8)          :: dt                ! timestep
    integer           :: nu                ! unit number
    real(R8)          :: qmeltall          ! q that would melt all accumulated water
    character(len=CS) :: fldname
    character(len=18) :: date_str

    character(*), parameter :: F00   = "('(dice_comp_run) ',8a)"
    character(*), parameter :: F04   = "('(dice_comp_run) ',2a,2i8,'s')"
    character(*), parameter :: F0D   = "('(dice_comp_run) ',a, i7,2x,i5,2x,i5,2x,d21.14)"
    character(*), parameter :: subName = "(dice_comp_run) "
    !-------------------------------------------------------------------------------

    !--------------------
    ! Debug output
    !--------------------

    if (debug_import > 1 .and. my_task == master_task) then
       do nfld = 1, mct_aVect_nRAttr(x2i)
          call shr_string_listGetName(trim(flds_x2i_mod), nfld, fldname)
          do n = 1, mct_aVect_lsize(x2i)
             write(logunit,F0D)'import: ymd,tod,n  = '// trim(fldname),target_ymd, target_tod, &
                  n, x2i%rattr(nfld,n)
          end do
       end do
    end if

    !--------------------
    ! ADVANCE ICE
    !--------------------

    call t_startf('DICE_RUN')
    call t_barrierf('dice_BARRIER',mpicom)
    call t_startf('dice')

    dt = modeldt * 1.0_r8

    !--- copy all stream fields to i2x as default (avifld in streams -> avofld in i2x)

    if (trim(datamode) /= 'NULL') then
       call t_startf('dice_strdata_advance')
       call shr_strdata_advance(SDICE,target_ymd,target_tod,mpicom,'dice')
       call t_stopf('dice_strdata_advance')
       call t_barrierf('dice_scatter_BARRIER',mpicom)
       call t_startf('dice_scatter')
       do n = 1,SDICE%nstreams
          call shr_dmodel_translateAV(SDICE%avs(n),i2x,avifld,avofld)
       enddo
       call t_stopf('dice_scatter')
    else
       call mct_aVect_zero(i2x)
    endif

    !-------------------------------------------------
    ! Determine data model behavior based on the mode
    !-------------------------------------------------

    call t_startf('dice_datamode')
    select case (trim(datamode))

    case('COPYALL')
       ! do nothing extra

    case('SSTDATA')
       if (firstcall .and. .not. read_restart) then
          ! iFrac0 = iFrac  ! previous step's ice fraction
          water  = 0.0_R8 ! previous step's water accumulation
          where (i2x%rAttr(kiFrac,:) > 0.0_R8) water(:) = flux_Qacc0
       endif

       lsize = mct_avect_lsize(i2x)

       tfreeze = shr_frz_freezetemp(x2i%rAttr(ksalinity,:)) + tFrz ! convert to Kelvin

       do n = 1,lsize

          !--- fix erroneous iFrac ---
          i2x%rAttr(kiFrac,n) = min(1.0_R8,max(0.0_R8,i2x%rAttr(kiFrac,n)))

          !--- fabricate ice surface T, fix erroneous iFrac ---
          if ( yc(n) > 0.0_R8) then
             i2x%rAttr(kt,n) = 260.0_R8 + 10.0_R8*cos(cosArg)
          else
             i2x%rAttr(kt,n) = 260.0_R8 - 10.0_R8*cos(cosArg)
          end if

          !--- set albedos (constant) ---
          i2x%rAttr(kavsdr,n) = ax_vsdr
          i2x%rAttr(kanidr,n) = ax_nidr
          i2x%rAttr(kavsdf,n) = ax_vsdf
          i2x%rAttr(kanidf,n) = ax_nidf

          !--- swnet is sent to cpl as a diagnostic quantity only ---
          !--- newly recv'd swdn goes with previously sent albedo ---
          !--- but albedos are (currently) time invariant         ---
          i2x%rAttr(kswnet,n) = (1.0_R8 - i2x%rAttr(kavsdr,n))*x2i%rAttr(kswvdr,n) &
                              + (1.0_R8 - i2x%rAttr(kanidr,n))*x2i%rAttr(kswndr,n) &
                              + (1.0_R8 - i2x%rAttr(kavsdf,n))*x2i%rAttr(kswvdf,n) &
                              + (1.0_R8 - i2x%rAttr(kanidf,n))*x2i%rAttr(kswndf,n)

          !--- compute melt/freeze water balance, adjust iFrac  -------------
          if ( .not. flux_Qacc ) then ! Q accumulation option is OFF
             i2x%rAttr(kmelth,n) = min(x2i%rAttr(kq,n),0.0_R8 ) ! q<0 => melt potential
             i2x%rAttr(kmelth,n) = max(i2x%rAttr(kmelth,n),Flux_Qmin   ) ! limit the melt rate
             i2x%rAttr(kmeltw,n) =    -i2x%rAttr(kmelth,n)/latice   ! corresponding water flux

          else                                 ! Q accumulation option is ON
             !--------------------------------------------------------------
             ! 1a) Q<0 & iFrac > 0  =>  infinite supply of water to melt
             ! 1b) Q<0 & iFrac = 0  =>  melt accumulated water only
             ! 2a) Q>0 & iFrac > 0  =>  zero-out accumulated water
             ! 2b) Q>0 & iFrac = 0  =>  accumulated water
             !--------------------------------------------------------------
             if ( x2i%rAttr(kq,n) <  0.0_R8 ) then ! Q<0 => melt
                if (i2x%rAttr(kiFrac,n) > 0.0_R8 ) then
                   i2x%rAttr(kmelth,n) = i2x%rAttr(kiFrac,n)*max(x2i%rAttr(kq,n),Flux_Qmin)
                   i2x%rAttr(kmeltw,n) =    -i2x%rAttr(kmelth,n)/latice
                   !  water(n) = < don't change this value >
                else
                   Qmeltall   = -water(n)*latice/dt
                   i2x%rAttr(kmelth,n) = max(x2i%rAttr(kq,n), Qmeltall, Flux_Qmin )
                   i2x%rAttr(kmeltw,n) = -i2x%rAttr(kmelth,n)/latice
                   water(n) =  water(n) - i2x%rAttr(kmeltw,n)*dt
                end if
             else                       ! Q>0 => freeze
                if (i2x%rAttr(kiFrac,n) > 0.0_R8 ) then
                   i2x%rAttr(kmelth,n) = 0.0_R8
                   i2x%rAttr(kmeltw,n) = 0.0_R8
                   water(n) = 0.0_R8
                else
                   i2x%rAttr(kmelth,n) = 0.0_R8
                   i2x%rAttr(kmeltw,n) = 0.0_R8
                   water(n) = water(n) + dt*x2i%rAttr(kq,n)/latice
                end if
             end if

             if (water(n) < 1.0e-16_R8 ) water(n) = 0.0_R8

             !--- non-zero water => non-zero iFrac ---
             if (i2x%rAttr(kiFrac,n) <= 0.0_R8  .and.  water(n) > 0.0_R8) then
                i2x%rAttr(kiFrac,n) = min(1.0_R8,water(n)/waterMax)
                ! i2x%rAttr(kT,n) = tfreeze(n)     ! T can be above freezing?!?
             end if

             !--- cpl multiplies melth & meltw by iFrac ---
             !--- divide by iFrac here => fixed quantity flux (not per area) ---
             if (i2x%rAttr(kiFrac,n) > 0.0_R8) then
                i2x%rAttr(kiFrac,n) = max( 0.01_R8, i2x%rAttr(kiFrac,n)) ! min iFrac
                i2x%rAttr(kmelth,n) = i2x%rAttr(kmelth,n)/i2x%rAttr(kiFrac,n)
                i2x%rAttr(kmeltw,n) = i2x%rAttr(kmeltw,n)/i2x%rAttr(kiFrac,n)
             else
                i2x%rAttr(kmelth,n) = 0.0_R8
                i2x%rAttr(kmeltw,n) = 0.0_R8
             end if
          end if

          !--- modify T wrt iFrac: (iFrac -> 0) => (T -> tfreeze) ---
          i2x%rAttr(kt,n) = tfreeze(n) + i2x%rAttr(kiFrac,n)*(i2x%rAttr(kt,n)-tfreeze(n))

       end do

       ! compute atm/ice surface fluxes
       call dice_flux_atmice( &
            iMask              ,x2i%rAttr(kz,:)     ,x2i%rAttr(kua,:)    ,x2i%rAttr(kva,:)  , &
            x2i%rAttr(kptem,:) ,x2i%rAttr(kshum,:)  ,x2i%rAttr(kdens,:)  ,x2i%rAttr(ktbot,:), &
            i2x%rAttr(kt,:)    ,i2x%rAttr(ksen,:)   ,i2x%rAttr(klat,:)   ,i2x%rAttr(klwup,:), &
            i2x%rAttr(kevap,:) ,i2x%rAttr(ktauxa,:) ,i2x%rAttr(ktauya,:) ,i2x%rAttr(ktref,:), &
            i2x%rAttr(kqref,:) ,logunit )

       ! compute ice/oce surface fluxes (except melth & meltw, see above)
       do n=1,lsize
          if (iMask(n) == 0) then
             i2x%rAttr(kswpen,n) = spval
             i2x%rAttr(kmelth,n) = spval
             i2x%rAttr(kmeltw,n) = spval
             i2x%rAttr(ksalt ,n) = spval
             i2x%rAttr(ktauxo,n) = spval
             i2x%rAttr(ktauyo,n) = spval
             i2x%rAttr(kiFrac,n) = 0.0_R8
          else
             !--- penetrating short wave ---
             i2x%rAttr(kswpen,n) = max(0.0_R8, flux_swpf*i2x%rAttr(kswnet,n) ) ! must be non-negative

             !--- i/o surface stress ( = atm/ice stress) ---
             i2x%rAttr(ktauxo,n) = i2x%rAttr(ktauxa,n)
             i2x%rAttr(ktauyo,n) = i2x%rAttr(ktauya,n)

             !--- salt flux ---
             i2x%rAttr(ksalt ,n) = 0.0_R8
          end if
          if (km /= 0) then
             i2x%rAttr(km, n) = imask(n)
          end if
          ! !--- save ifrac for next timestep
          ! iFrac0(n) = i2x%rAttr(kiFrac,n)
       end do

       ! Compute outgoing aerosol fluxes
       do n = 1,lsize
          i2x%rAttr(kbcpho ,n) = x2i%rAttr(kbcphodry,n)
          i2x%rAttr(kbcphi ,n) = x2i%rAttr(kbcphidry,n) + x2i%rAttr(kbcphiwet,n)
          i2x%rAttr(kflxdst,n) = x2i%rAttr(kdstdry1,n) + x2i%rAttr(kdstwet1,n) &
                               + x2i%rAttr(kdstdry2,n) + x2i%rAttr(kdstwet2,n) &
                               + x2i%rAttr(kdstdry3,n) + x2i%rAttr(kdstwet3,n) &
                               + x2i%rAttr(kdstdry4,n) + x2i%rAttr(kdstwet4,n)
       end do

    end select

    !-------------------------------------------------
    ! optional per thickness category fields
    !-------------------------------------------------

    if (flds_i2o_per_cat) then
       do n=1,lsize
          i2x%rAttr(kiFrac_01,n)       = i2x%rAttr(kiFrac,n)
          i2x%rAttr(kswpen_iFrac_01,n) = i2x%rAttr(kswpen,n) * i2x%rAttr(kiFrac,n)
       end do
    end if

    call t_stopf('dice_datamode')

    !--------------------
    ! Debug output
    !--------------------

    if (debug_export > 1 .and. my_task == master_task) then
       do nfld = 1, mct_aVect_nRAttr(i2x)
          call shr_string_listGetName(trim(flds_i2x_mod), nfld, fldname)
          do n = 1, mct_aVect_lsize(i2x)
             write(logunit,F0D)'export: ymd,tod,n  = '// trim(fldname),target_ymd, target_tod, &
                  n, i2x%rattr(nfld,n)
          end do
       end do
    end if

    !--------------------
    ! Write restart
    !--------------------

    if (write_restart) then
       call t_startf('dice_restart')
       call shr_cal_datetod2string(date_str, target_ymd, target_tod)
       write(rest_file,"(6a)") &
            trim(case_name), '.dice',trim(inst_suffix),'.r.', &
            trim(date_str),'.nc'
       write(rest_file_strm,"(6a)") &
            trim(case_name), '.dice',trim(inst_suffix),'.rs1.', &
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
       call shr_pcdf_readwrite('write',SDICE%pio_subsystem, SDICE%io_type, &
            trim(rest_file), mpicom, SDICE%gsmap, clobber=.true., rf1=water, rf1n='water')
       if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file_strm),target_ymd,target_tod
       call shr_strdata_restWrite(trim(rest_file_strm),SDICE,mpicom,trim(case_name),'SDICE strdata')
       call t_stopf('dice_restart')
    endif

    call t_stopf('dice')

    firstcall = .false.

    call t_stopf('DICE_RUN')

  end subroutine dice_comp_run

end module dice_comp_mod
