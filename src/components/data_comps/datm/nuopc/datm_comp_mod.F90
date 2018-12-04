#ifdef AIX
  @PROCESS ALIAS_SIZE(805306368)
#endif
module datm_comp_mod

  ! !USES:
  use NUOPC                 , only : NUOPC_Advertise
  use ESMF                  , only : ESMF_State, ESMF_SUCCESS, ESMF_State
  use ESMF                  , only : ESMF_Mesh, ESMF_DistGrid, ESMF_MeshGet, ESMF_DistGridGet
  use perf_mod              , only : t_startf, t_stopf, t_adj_detailf, t_barrierf
  use mct_mod               , only : mct_gsmap_init
  use mct_mod               , only : mct_avect, mct_avect_indexRA, mct_avect_zero, mct_aVect_nRattr
  use mct_mod               , only : mct_avect_init, mct_avect_lsize
  use shr_const_mod         , only : SHR_CONST_SPVAL
  use shr_const_mod         , only : SHR_CONST_TKFRZ
  use shr_const_mod         , only : SHR_CONST_PI
  use shr_const_mod         , only : SHR_CONST_PSTD
  use shr_const_mod         , only : SHR_CONST_STEBOL
  use shr_const_mod         , only : SHR_CONST_RDAIR
  use med_constants_mod     , only : R8, CS, CL, CXX
  use shr_string_mod        , only : shr_string_listGetName
  use shr_sys_mod           , only : shr_sys_abort
  use shr_file_mod          , only : shr_file_getunit, shr_file_freeunit
  use shr_cal_mod           , only : shr_cal_calendarname
  use shr_cal_mod           , only : shr_cal_date2julian, shr_cal_datetod2string
  use shr_mpi_mod           , only : shr_mpi_bcast, shr_mpi_max
  use shr_precip_mod        , only : shr_precip_partition_rain_snow_ramp
  use shr_strdata_mod       , only : shr_strdata_init_model_domain
  use shr_strdata_mod       , only : shr_strdata_init_streams
  use shr_strdata_mod       , only : shr_strdata_init_mapping
  use shr_strdata_mod       , only : shr_strdata_type, shr_strdata_pioinit
  use shr_strdata_mod       , only : shr_strdata_print, shr_strdata_restRead
  use shr_strdata_mod       , only : shr_strdata_advance, shr_strdata_restWrite
  use shr_strdata_mod       , only : shr_strdata_setorbs
  use shr_dmodel_mod        , only : shr_dmodel_translate_list, shr_dmodel_translateAV_list
  use shr_nuopc_scalars_mod , only : flds_scalar_name
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
  use dshr_nuopc_mod        , only : fld_list_type
  use dshr_nuopc_mod        , only : dshr_fld_add
  use datm_shr_mod          , only : datm_shr_esat, datm_shr_CORE2getFactors
  use datm_shr_mod          , only : datamode       ! namelist input
  use datm_shr_mod          , only : wiso_datm      ! namelist input
  use datm_shr_mod          , only : rest_file      ! namelist input
  use datm_shr_mod          , only : rest_file_strm ! namelist input
  use datm_shr_mod          , only : factorfn       ! namelist input
  use datm_shr_mod          , only : iradsw         ! namelist input
  use datm_shr_mod          , only : nullstr
  use datm_shr_mod          , only : presaero

  ! !PUBLIC TYPES:

  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: datm_comp_advertise
  public :: datm_comp_init
  public :: datm_comp_run

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  integer                    :: debug_import = 0      ! debug level (if > 0 will print all import fields)
  integer                    :: debug_export = 0      ! debug level (if > 0 will print all export fields)

  real(R8)                   :: tbotmax               ! units detector
  real(R8)                   :: tdewmax               ! units detector
  real(R8)                   :: anidrmax              ! existance detector

  ! Attribute vectors field indices
  integer                    :: kz,ktopo,ku,kv,ktbot,kptem,kshum,kdens,kpbot,kpslv,klwdn
  integer                    :: krc,krl,ksc,ksl,kswndr,kswndf,kswvdr,kswvdf,kswnet
  integer                    :: kanidr,kanidf,kavsdr,kavsdf
  integer                    :: kshum_16O, kshum_18O, kshum_HDO
  integer                    :: krc_18O, krc_HDO
  integer                    :: krl_18O, krl_HDO
  integer                    :: ksc_18O, ksc_HDO
  integer                    :: ksl_18O, ksl_HDO
  integer                    :: stbot,swind,sz,spbot,sshum,stdew,srh,slwdn,sswdn,sswdndf,sswdndr
  integer                    :: sprecc,sprecl,sprecn,sco2p,sco2d,sswup,sprec,starcf
  integer                    :: srh_16O, srh_18O, srh_HDO, sprecn_16O, sprecn_18O, sprecn_HDO
  integer                    :: sprecsf
  integer                    :: sprec_af,su_af,sv_af,stbot_af,sshum_af,spbot_af,slwdn_af,sswdn_af

  type(mct_avect)            :: avstrm         ! av of data from stream
  character(len=CS), pointer :: avifld(:)      ! character array for field names coming from streams
  character(len=CS), pointer :: avofld(:)      ! character array for field names to be sent/received from mediator
  character(len=CS), pointer :: stifld(:)      ! character array for field names coming from streams
  character(len=CS), pointer :: stofld(:)      ! character array for field intermediate avs for calculations
  character(len=CL), pointer :: ilist_av(:)    ! input  character array for translation (avifld->avofld)
  character(len=CL), pointer :: olist_av(:)    ! output character array for translation (avifld->avofld)
  integer    ,       pointer :: count_av(:)    ! number of fields in translation (avifld->avofld)
  character(len=CL), pointer :: ilist_st(:)    ! input  character array for translation (stifld->strmofld)
  character(len=CL), pointer :: olist_st(:)    ! output character array for translation (stifld->strmofld)
  integer      ,     pointer :: count_st(:)    ! number of fields in translation (stifld->strmofld)
  character(len=CXX)         :: flds_strm = '' ! colon deliminated string of field names
  character(len=CXX)         :: flds_a2x_mod
  character(len=CXX)         :: flds_x2a_mod

  real(R8), pointer          :: xc(:), yc(:)   ! arrays of model latitudes and longitudes
  real(R8), pointer          :: windFactor(:)
  real(R8), pointer          :: winddFactor(:)
  real(R8), pointer          :: qsatFactor(:)

  real(R8),parameter         :: tKFrz  = SHR_CONST_TKFRZ
  real(R8),parameter         :: degtorad = SHR_CONST_PI/180.0_R8
  real(R8),parameter         :: pstd   = SHR_CONST_PSTD     ! standard pressure ~ Pa
  real(R8),parameter         :: stebol = SHR_CONST_STEBOL   ! Stefan-Boltzmann constant ~ W/m^2/K^4
  real(R8),parameter         :: rdair  = SHR_CONST_RDAIR    ! dry air gas constant   ~ J/K/kg
  real(R8),parameter         :: avg_c0 =  61.846_R8
  real(R8),parameter         :: avg_c1 =   1.107_R8
  real(R8),parameter         :: amp_c0 = -21.841_R8
  real(R8),parameter         :: amp_c1 =  -0.447_R8
  real(R8),parameter         :: phs_c0 =   0.298_R8
  real(R8),parameter         :: dLWarc =  -5.000_R8

  real(R8)           :: dTarc(12)
  data   dTarc      / 0.49_R8, 0.06_R8,-0.73_R8,  -0.89_R8,-0.77_R8,-1.02_R8, &
                     -1.99_R8,-0.91_R8, 1.72_R8,   2.30_R8, 1.81_R8, 1.06_R8/

  character(len=*),parameter :: rpfile = 'rpointer.atm'
  character(*),parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine datm_comp_advertise(importState, exportState, &
       atm_prognostic, &
       flds_wiso, flds_co2a, flds_co2b, flds_co2c, &
       fldsFrAtm_num, fldsFrAtm, fldsToAtm_num, fldsToAtm, &
       flds_a2x, flds_x2a, rc)

    ! 1. determine export and import fields to advertise to mediator
    ! 2. determine translation of fields from streams to export/import fields
    ! 3. determine module indices for attribute vectors

    ! input/output arguments
    type(ESMF_State)                   :: importState
    type(ESMF_State)                   :: exportState
    logical              , intent(in)  :: atm_prognostic
    logical              , intent(in)  :: flds_wiso      ! use case
    logical              , intent(in)  :: flds_co2a      ! use case
    logical              , intent(in)  :: flds_co2b      ! use case
    logical              , intent(in)  :: flds_co2c      ! use case
    integer              , intent(out) :: fldsFrAtm_num
    type (fld_list_type) , intent(out) :: fldsFrAtm(:)
    integer              , intent(out) :: fldsToAtm_num
    type (fld_list_type) , intent(out) :: fldsToAtm(:)
    character(len=*)     , intent(out) :: flds_a2x
    character(len=*)     , intent(out) :: flds_x2a
    integer              , intent(out) :: rc

    ! local variables
    integer :: n
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !-------------------
    ! export fields
    !-------------------

    ! scalar fields that need to be advertised

    fldsFrAtm_num=1
    fldsFrAtm(1)%stdname = trim(flds_scalar_name)

    ! export fields that have a corresponding stream field

    call dshr_fld_add(data_fld="topo", data_fld_array=avifld, model_fld="Sa_topo", model_fld_array=avofld, &
         model_fld_concat=flds_a2x, model_fld_index=ktopo, fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

    call dshr_fld_add(data_fld="z", data_fld_array=avifld, model_fld="Sa_z", model_fld_array=avofld, &
         model_fld_concat=flds_a2x, model_fld_index=kz, fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

    call dshr_fld_add(data_fld="u", data_fld_array=avifld, model_fld="Sa_u", model_fld_array=avofld, &
         model_fld_concat=flds_a2x, model_fld_index=ku, fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

    call dshr_fld_add(data_fld="v", data_fld_array=avifld, model_fld="Sa_v", model_fld_array=avofld, &
         model_fld_concat=flds_a2x, model_fld_index=kv, fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

    call dshr_fld_add(data_fld="ptem", data_fld_array=avifld, model_fld="Sa_ptem", model_fld_array=avofld, &
         model_fld_concat=flds_a2x, model_fld_index=kptem, fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

    call dshr_fld_add(data_fld="dens", data_fld_array=avifld, model_fld="Sa_dens", model_fld_array=avofld, &
         model_fld_concat=flds_a2x, model_fld_index=kdens, fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

    call dshr_fld_add(data_fld="pslv", data_fld_array=avifld, model_fld="Sa_pslv", model_fld_array=avofld, &
         model_fld_concat=flds_a2x, model_fld_index=kpslv, fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

    call dshr_fld_add(data_fld="rainc", data_fld_array=avifld, model_fld="Faxa_rainc", model_fld_array=avofld, &
         model_fld_concat=flds_a2x, model_fld_index=krc, fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

    call dshr_fld_add(data_fld="rainl", data_fld_array=avifld, model_fld="Faxa_rainl", model_fld_array=avofld, &
         model_fld_concat=flds_a2x, model_fld_index=krl, fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

    call dshr_fld_add(data_fld="snowc", data_fld_array=avifld, model_fld="Faxa_snowc", model_fld_array=avofld, &
         model_fld_concat=flds_a2x, model_fld_index=ksc, fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

    call dshr_fld_add(data_fld="snowl", data_fld_array=avifld, model_fld="Faxa_snowl", model_fld_array=avofld, &
         model_fld_concat=flds_a2x, model_fld_index=ksl, fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

    call dshr_fld_add(data_fld="swndr", data_fld_array=avifld, model_fld="Faxa_swndr", model_fld_array=avofld, &
         model_fld_concat=flds_a2x, model_fld_index=kswndr, fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

    call dshr_fld_add(data_fld="swvdr", data_fld_array=avifld, model_fld="Faxa_swvdr", model_fld_array=avofld, &
         model_fld_concat=flds_a2x, model_fld_index=kswvdr, fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

    call dshr_fld_add(data_fld="swndf", data_fld_array=avifld, model_fld="Faxa_swndf", model_fld_array=avofld, &
         model_fld_concat=flds_a2x, model_fld_index=kswndf, fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

    call dshr_fld_add(data_fld="swvdf", data_fld_array=avifld, model_fld="Faxa_swvdf", model_fld_array=avofld, &
         model_fld_concat=flds_a2x, model_fld_index=kswvdf, fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

    call dshr_fld_add(data_fld="swnet", data_fld_array=avifld, model_fld="Faxa_swnet", model_fld_array=avofld, &
         model_fld_concat=flds_a2x, model_fld_index=kswnet, fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

    ! export fields that have a corresponding stream field AND that have a corresponding internal field

    call dshr_fld_add(data_fld="tbot", data_fld_array=avifld, model_fld="Sa_tbot", model_fld_array=avofld, &
         model_fld_concat=flds_a2x, model_fld_index=ktbot , fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

    call dshr_fld_add(data_fld="pbot", data_fld_array=avifld, model_fld="Sa_pbot", model_fld_array=avofld, &
         model_fld_concat=flds_a2x, model_fld_index=kpbot , fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

    call dshr_fld_add(data_fld="shum", data_fld_array=avifld, model_fld="Sa_shum", model_fld_array=avofld, &
         model_fld_concat=flds_a2x, model_fld_index=kshum , fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

    call dshr_fld_add(data_fld="lwdn", data_fld_array=avifld, model_fld="Faxa_lwdn", model_fld_array=avofld, &
         model_fld_concat=flds_a2x, model_fld_index=klwdn , fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

    if (flds_co2a .or. flds_co2b .or. flds_co2c) then
       call dshr_fld_add(data_fld="co2prog", data_fld_array=avifld, model_fld="Sa_co2prog", model_fld_array=avofld, &
            model_fld_concat=flds_x2a,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="co2diag", data_fld_array=avifld, model_fld="Sa_co2diag", model_fld_array=avofld, &
            model_fld_concat=flds_x2a,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)
    end if

    if (presaero) then

       call dshr_fld_add(data_fld="bcphidry", data_fld_array=avifld, model_fld="Faxa_bcphidry", model_fld_array=avofld, &
            model_fld_concat=flds_a2x,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="bcphodry", data_fld_array=avifld, model_fld="Faxa_bcphodry", model_fld_array=avofld, &
            model_fld_concat=flds_a2x,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="bcphiwet", data_fld_array=avifld, model_fld="Faxa_bcphiwet", model_fld_array=avofld, &
            model_fld_concat=flds_a2x,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="ocphidry", data_fld_array=avifld, model_fld="Faxa_ocphidry", model_fld_array=avofld, &
            model_fld_concat=flds_a2x,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="ocphodry", data_fld_array=avifld, model_fld="Faxa_ocphodry", model_fld_array=avofld, &
            model_fld_concat=flds_a2x,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="ocphiwet", data_fld_array=avifld, model_fld="Faxa_ocphiwet", model_fld_array=avofld, &
            model_fld_concat=flds_a2x,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="dstwet1", data_fld_array=avifld, model_fld="Faxa_dstwet1", model_fld_array=avofld, &
            model_fld_concat=flds_a2x,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="dstwet2", data_fld_array=avifld, model_fld="Faxa_dstwet2", model_fld_array=avofld, &
            model_fld_concat=flds_a2x,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="dstwet3", data_fld_array=avifld, model_fld="Faxa_dstwet3", model_fld_array=avofld, &
            model_fld_concat=flds_a2x,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="dstwet4", data_fld_array=avifld, model_fld="Faxa_dstwet4", model_fld_array=avofld, &
            model_fld_concat=flds_a2x,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="dstdry1", data_fld_array=avifld, model_fld="Faxa_dstdry1", model_fld_array=avofld, &
            model_fld_concat=flds_a2x,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="dstdry2", data_fld_array=avifld, model_fld="Faxa_dstdry2", model_fld_array=avofld, &
            model_fld_concat=flds_a2x,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="dstdry3", data_fld_array=avifld, model_fld="Faxa_dstdry3", model_fld_array=avofld, &
            model_fld_concat=flds_a2x,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="dstdry4", data_fld_array=avifld, model_fld="Faxa_dstdry4", model_fld_array=avofld, &
            model_fld_concat=flds_a2x,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)
    end if

    ! isotopic forcing

    if (flds_wiso) then

       call dshr_fld_add(data_fld="rainc_18O", data_fld_array=avifld, model_fld="Faxa_rainc_18O", model_fld_array=avofld, &
            model_fld_concat=flds_a2x, model_fld_index=krc_18O,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="rainc_HDO", data_fld_array=avifld, model_fld="Faxa_rainc_HDO", model_fld_array=avofld, &
            model_fld_concat=flds_a2x, model_fld_index=krc_HDO,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="rainl_18O", data_fld_array=avifld, model_fld="Faxa_rainl_18O", model_fld_array=avofld, &
            model_fld_concat=flds_a2x, model_fld_index=krl_18O,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="rainl_HDO", data_fld_array=avifld, model_fld="Faxa_rainl_HDO", model_fld_array=avofld, &
            model_fld_concat=flds_a2x, model_fld_index=krl_HDO,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="snowc_18O", data_fld_array=avifld, model_fld="Faxa_snowc_18O", model_fld_array=avofld, &
            model_fld_concat=flds_a2x, model_fld_index=ksc_18O,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="snowc_HDO", data_fld_array=avifld, model_fld="Faxa_snowc_HDO", model_fld_array=avofld, &
            model_fld_concat=flds_a2x, model_fld_index=ksc_HDO,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="snowl_18O", data_fld_array=avifld, model_fld="Faxa_snowl_18O", model_fld_array=avofld, &
            model_fld_concat=flds_a2x, model_fld_index=ksl_18O,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="snowl_HDO", data_fld_array=avifld, model_fld="Faxa_snowl_HDO", model_fld_array=avofld, &
            model_fld_concat=flds_a2x, model_fld_index=ksl_HDO,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="shum_16O", data_fld_array=avifld, model_fld="Sa_shum_16O", model_fld_array=avofld, &
            model_fld_concat=flds_a2x, model_fld_index=kshum_16O,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="shum_18O", data_fld_array=avifld, model_fld="Sa_shum_18O", model_fld_array=avofld, &
            model_fld_concat=flds_a2x, model_fld_index=kshum_18O,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)

       call dshr_fld_add(data_fld="shum_HDO", data_fld_array=avifld, model_fld="Sa_shum_HDO", model_fld_array=avofld, &
            model_fld_concat=flds_a2x, model_fld_index=kshum_HDO,    fldlist_num=fldsFrAtm_num, fldlist=fldsFrAtm)
    end if

    !-------------------
    ! import fields (have no corresponding stream fields)
    !-------------------

    if (atm_prognostic) then

       fldsToAtm_num=1
       fldsToAtm(1)%stdname = trim(flds_scalar_name)

       ! The module indices set by the model_fld_index argument are used in the run phase
       call dshr_fld_add(model_fld="Sx_avsdr", model_fld_concat=flds_x2a, model_fld_index=kavsdr, &
            fldlist_num=fldsToAtm_num, fldlist=fldsToAtm)
       call dshr_fld_add(model_fld="Sx_anidr", model_fld_concat=flds_x2a, model_fld_index=kanidr, &
            fldlist_num=fldsToAtm_num, fldlist=fldsToAtm)
       call dshr_fld_add(model_fld="Sx_avsdf", model_fld_concat=flds_x2a, model_fld_index=kavsdf, &
            fldlist_num=fldsToAtm_num, fldlist=fldsToAtm)
       call dshr_fld_add(model_fld="Sx_anidf", model_fld_concat=flds_x2a, model_fld_index=kanidf, &
            fldlist_num=fldsToAtm_num, fldlist=fldsToAtm)

       call dshr_fld_add(model_fld="Sx_tref"       , model_fld_concat=flds_x2a, fldlist_num=fldsToAtm_num, fldlist=fldsToAtm)
       call dshr_fld_add(model_fld="Sx_qref"       , model_fld_concat=flds_x2a, fldlist_num=fldsToAtm_num, fldlist=fldsToAtm)
       call dshr_fld_add(model_fld="Sx_t"          , model_fld_concat=flds_x2a, fldlist_num=fldsToAtm_num, fldlist=fldsToAtm)
       call dshr_fld_add(model_fld="So_t"          , model_fld_concat=flds_x2a, fldlist_num=fldsToAtm_num, fldlist=fldsToAtm)
       call dshr_fld_add(model_fld="Sl_snowh"      , model_fld_concat=flds_x2a, fldlist_num=fldsToAtm_num, fldlist=fldsToAtm)
       call dshr_fld_add(model_fld="Sl_lfrac"      , model_fld_concat=flds_x2a, fldlist_num=fldsToAtm_num, fldlist=fldsToAtm)
       call dshr_fld_add(model_fld="Si_ifrac"      , model_fld_concat=flds_x2a, fldlist_num=fldsToAtm_num, fldlist=fldsToAtm)
       call dshr_fld_add(model_fld="So_ofrac"      , model_fld_concat=flds_x2a, fldlist_num=fldsToAtm_num, fldlist=fldsToAtm)
       call dshr_fld_add(model_fld="Faxx_taux"     , model_fld_concat=flds_x2a, fldlist_num=fldsToAtm_num, fldlist=fldsToAtm)
       call dshr_fld_add(model_fld="Faxx_tauy"     , model_fld_concat=flds_x2a, fldlist_num=fldsToAtm_num, fldlist=fldsToAtm)
       call dshr_fld_add(model_fld="Faxx_lat"      , model_fld_concat=flds_x2a, fldlist_num=fldsToAtm_num, fldlist=fldsToAtm)
       call dshr_fld_add(model_fld="Faxx_sen"      , model_fld_concat=flds_x2a, fldlist_num=fldsToAtm_num, fldlist=fldsToAtm)
       call dshr_fld_add(model_fld="Faxx_lwup"     , model_fld_concat=flds_x2a, fldlist_num=fldsToAtm_num, fldlist=fldsToAtm)
       call dshr_fld_add(model_fld="Faxx_evap"     , model_fld_concat=flds_x2a, fldlist_num=fldsToAtm_num, fldlist=fldsToAtm)
     ! call dshr_fld_add(model_fld="Fall_fco2_lnd" , model_fld_concat=flds_x2a ,fldlist_num=fldsToAtm_num, fldlist=fldsToAtm)
     ! call dshr_fld_add(model_fld="Faoo_fco2_ocn" , model_fld_concat=flds_x2a ,fldlist_num=fldsToAtm_num, fldlist=fldsToAtm)

    end if

    !-------------------
    ! advertise fields for import and export states
    !-------------------

    do n = 1,fldsFrAtm_num
       call NUOPC_Advertise(exportState, standardName=fldsFrAtm(n)%stdname, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    if (atm_prognostic) then
       do n = 1,fldsToAtm_num
          call NUOPC_Advertise(importState, standardName=fldsToAtm(n)%stdname, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end do
    end if

    !-------------------
    ! Save flds_x2a and flds_a2x as module variables for use in debugging
    !-------------------

    flds_x2a_mod = trim(flds_x2a)
    flds_a2x_mod = trim(flds_a2x)

    !-------------------
    ! module character arrays stifld and stofld
    !-------------------

    ! - stifld is a character array of stream field names
    ! - stofld is a character array of data model field names that have a one-to-one correspondence with names in stifld
    ! - flds_strm is a colon delimited string of field names that is created from the field names in stofld for ONLY
    !   those field names that are available in the data streams present in SDATM%sdatm
    ! - avstrm is an attribute vector created from flds_strm

    call dshr_fld_add(data_fld="wind"      , data_fld_array=stifld, model_fld="strm_wind"      , model_fld_array=stofld)
    call dshr_fld_add(data_fld="tdew"      , data_fld_array=stifld, model_fld="strm_tdew"      , model_fld_array=stofld)
    call dshr_fld_add(data_fld="tbot"      , data_fld_array=stifld, model_fld="strm_tbot"      , model_fld_array=stofld)
    call dshr_fld_add(data_fld="pbot"      , data_fld_array=stifld, model_fld="strm_pbot"      , model_fld_array=stofld)
    call dshr_fld_add(data_fld="shum"      , data_fld_array=stifld, model_fld="strm_shum"      , model_fld_array=stofld)
    call dshr_fld_add(data_fld="lwdn"      , data_fld_array=stifld, model_fld="strm_lwdn"      , model_fld_array=stofld)
    call dshr_fld_add(data_fld="wind"      , data_fld_array=stifld, model_fld="strm_wind"      , model_fld_array=stofld)
    call dshr_fld_add(data_fld="rh"        , data_fld_array=stifld, model_fld="strm_rh"        , model_fld_array=stofld)
    call dshr_fld_add(data_fld="swdn"      , data_fld_array=stifld, model_fld="strm_swdn"      , model_fld_array=stofld)
    call dshr_fld_add(data_fld="swdndf"    , data_fld_array=stifld, model_fld="strm_swdndf"    , model_fld_array=stofld)
    call dshr_fld_add(data_fld="swdndr"    , data_fld_array=stifld, model_fld="strm_swdndr"    , model_fld_array=stofld)
    call dshr_fld_add(data_fld="prec"      , data_fld_array=stifld, model_fld="strm_prec"      , model_fld_array=stofld)
    call dshr_fld_add(data_fld="precc"     , data_fld_array=stifld, model_fld="strm_precc"     , model_fld_array=stofld)
    call dshr_fld_add(data_fld="precl"     , data_fld_array=stifld, model_fld="strm_precl"     , model_fld_array=stofld)
    call dshr_fld_add(data_fld="precn"     , data_fld_array=stifld, model_fld="strm_precn"     , model_fld_array=stofld)
    call dshr_fld_add(data_fld="swup"      , data_fld_array=stifld, model_fld="strm_swup"      , model_fld_array=stofld)
    call dshr_fld_add(data_fld="tarcf"     , data_fld_array=stifld, model_fld="strm_tarcf"     , model_fld_array=stofld)

    ! water isotopes
    call dshr_fld_add(data_fld="rh_16O"    , data_fld_array=stifld, model_fld="strm_rh_16O"    , model_fld_array=stofld)
    call dshr_fld_add(data_fld="rh_18O"    , data_fld_array=stifld, model_fld="strm_rh_18O"    , model_fld_array=stofld)
    call dshr_fld_add(data_fld="rh_HDO"    , data_fld_array=stifld, model_fld="strm_rh_HDO"    , model_fld_array=stofld)
    call dshr_fld_add(data_fld="precn_16O" , data_fld_array=stifld, model_fld="strm_precn_16O" , model_fld_array=stofld)
    call dshr_fld_add(data_fld="precn_18O" , data_fld_array=stifld, model_fld="strm_precn_18O" , model_fld_array=stofld)
    call dshr_fld_add(data_fld="precn_HDO" , data_fld_array=stifld, model_fld="strm_precn_HDO" , model_fld_array=stofld)

    ! values for optional bias correction / anomaly forcing (add Sa_precsf for precip scale factor)
    call dshr_fld_add(data_fld="precsf"    , data_fld_array=stifld, model_fld="strm_precsf"    , model_fld_array=stofld)
    call dshr_fld_add(data_fld="prec_af"   , data_fld_array=stifld, model_fld="strm_prec_af"   , model_fld_array=stofld)
    call dshr_fld_add(data_fld="u_af"      , data_fld_array=stifld, model_fld="strm_u_af"      , model_fld_array=stofld)
    call dshr_fld_add(data_fld="v_af"      , data_fld_array=stifld, model_fld="strm_v_af"      , model_fld_array=stofld)
    call dshr_fld_add(data_fld="tbot_af"   , data_fld_array=stifld, model_fld="strm_tbot_af"   , model_fld_array=stofld)
    call dshr_fld_add(data_fld="pbot_af"   , data_fld_array=stifld, model_fld="strm_pbot_af"   , model_fld_array=stofld)
    call dshr_fld_add(data_fld="shum_af"   , data_fld_array=stifld, model_fld="strm_shum_af"   , model_fld_array=stofld)
    call dshr_fld_add(data_fld="swdn_af"   , data_fld_array=stifld, model_fld="strm_swdn_af"   , model_fld_array=stofld)
    call dshr_fld_add(data_fld="lwdn_af"   , data_fld_array=stifld, model_fld="strm_lwdn_af"   , model_fld_array=stofld)

    if (flds_co2a .or. flds_co2b .or. flds_co2c) then
       call dshr_fld_add(data_fld="co2prog", data_fld_array=stifld, model_fld="strm_co2prog"   , model_fld_array=stofld)
       call dshr_fld_add(data_fld="co2diag", data_fld_array=stifld, model_fld="strm_co2diag"   , model_fld_array=stofld)
    end if

  end subroutine datm_comp_advertise

  !===============================================================================

  subroutine datm_comp_init(x2a, a2x, &
       SDATM, mpicom, compid, my_task, master_task, &
       inst_suffix, inst_name, logunit, read_restart, &
       scmMode, scmlat, scmlon, &
       orbEccen, orbMvelpp, orbLambm0, orbObliqr, &
       calendar, modeldt, current_ymd, current_tod, current_mon, &
       atm_prognostic, mesh)

    use dshr_nuopc_mod, only : dshr_fld_add

    ! !DESCRIPTION: initialize data atm model

    ! !INPUT/OUTPUT PARAMETERS:
    type(mct_aVect)        , intent(inout) :: x2a
    type(mct_aVect)        , intent(inout) :: a2x
    type(shr_strdata_type) , intent(inout) :: SDATM          ! model shr_strdata instance (output)
    integer                , intent(in)    :: mpicom         ! mpi communicator
    integer                , intent(in)    :: compid         ! mct comp id
    integer                , intent(in)    :: my_task        ! my task in mpi communicator mpicom
    integer                , intent(in)    :: master_task    ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix    ! char string associated with instance
    character(len=*)       , intent(in)    :: inst_name      ! fullname of current instance (ie."lnd_0001")
    integer                , intent(in)    :: logunit        ! logging unit number
    logical                , intent(in)    :: read_restart   ! start from restart
    logical                , intent(in)    :: scmMode        ! single column mode
    real(R8)               , intent(in)    :: scmLat         ! single column lat
    real(R8)               , intent(in)    :: scmLon         ! single column lon
    real(R8)               , intent(in)    :: orbEccen       ! orb eccentricity (unit-less)
    real(R8)               , intent(in)    :: orbMvelpp      ! orb moving vernal eq (radians)
    real(R8)               , intent(in)    :: orbLambm0      ! orb mean long of perhelion (radians)
    real(R8)               , intent(in)    :: orbObliqr      ! orb obliquity (radians)
    character(len=*)       , intent(in)    :: calendar       ! calendar type
    integer                , intent(in)    :: modeldt        ! model time step
    integer                , intent(in)    :: current_ymd    ! model date
    integer                , intent(in)    :: current_tod    ! model sec into model date
    integer                , intent(in)    :: current_mon    ! model month
    logical                , intent(in)    :: atm_prognostic ! if true, need x2a data
    type(ESMF_Mesh)        , intent(inout) :: mesh

    !--- local variables ---
    integer                      :: n,k            ! generic counters
    integer                      :: lsize          ! local size
    integer                      :: kmask          ! field reference
    integer                      :: klon,klat      ! field reference
    integer                      :: kfld           ! fld index
    integer                      :: cnt            ! counter
    logical                      :: exists,exists1 ! filename existance
    integer                      :: nu             ! unit number
    integer                      :: stepno         ! step number
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
    character(*), parameter      :: F00   ="('(datm_comp_init) ',8a)"
    character(*), parameter      :: F01   ="('(datm_comp_init) ',a,2f10.4)"
    character(*), parameter      :: subName ="(datm_comp_init)"
    !-------------------------------------------------------------------------------

    call t_startf('DATM_INIT')

    !----------------------------------------------------------------------------
    ! Initialize PIO
    !----------------------------------------------------------------------------

    call shr_strdata_pioinit(SDATM, COMPID)

    !----------------------------------------------------------------------------
    ! Create a data model global seqmap
    !----------------------------------------------------------------------------

    call t_startf('datm_strdata_init')

    if (my_task == master_task) write(logunit,F00) ' initialize SDATM gsmap'

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
    call mct_gsMap_init( SDATM%gsmap, gindex, mpicom, compid, lsize, gsize)
    deallocate(gindex)

    !----------------------------------------------------------------------------
    ! Initialize SDATM
    !----------------------------------------------------------------------------

    ! The call to shr_strdata_init_model_domain creates the SDATM%gsmap which
    ! is a '2d1d' decommp (1d decomp of 2d grid) and also create SDATM%grid

    SDATM%calendar = trim(shr_cal_calendarName(trim(calendar)))

    if (scmmode) then
       if (my_task == master_task) write(logunit,F01) ' scm lon lat = ',scmlon,scmlat
       call shr_strdata_init_model_domain(SDATM, mpicom, compid, my_task, &
            scmmode=scmmode, scmlon=scmlon, scmlat=scmlat, gsmap=SDATM%gsmap)
    else
       call shr_strdata_init_model_domain(SDATM, mpicom, compid, my_task, gsmap=SDATM%gsmap)
    end if

    if (my_task == master_task) then
       call shr_strdata_print(SDATM,'SDATM data')
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
    klon = mct_aVect_indexRA(SDATM%grid%data,'lon')
    do n = 1, lsize
       if (abs( SDATM%grid%data%rattr(klon,n) - xc(n)) > 1.e-5) then
          write(6,*)'ERROR: DATM n, lon(domain), lon(mesh) = ',n, SDATM%grid%data%rattr(klon,n),xc(n)
          write(6,*)'ERROR: DATM lon diff = ',abs(SDATM%grid%data%rattr(klon,n) -  xc(n)),' too large'
          call shr_sys_abort()
       end if
       !SDATM%grid%data%rattr(klon,n) = xc(n) ! overwrite ggrid with mesh data
       xc(n) = SDATM%grid%data%rattr(klon,n)
    end do
    klat = mct_aVect_indexRA(SDATM%grid%data,'lat')
    do n = 1, lsize
       if (abs( SDATM%grid%data%rattr(klat,n) -  yc(n)) > 1.e-5) then
          write(6,*)'ERROR: DATM n, lat(domain), lat(mesh) = ',n,SDATM%grid%data%rattr(klat,n),yc(n)
          write(6,*)'ERROR: DATM lat diff = ',abs(SDATM%grid%data%rattr(klat,n) -  yc(n)),' too large'
          call shr_sys_abort()
       end if
       !SDATM%grid%data%rattr(klat,n) = yc(n) ! overwrite ggrid with mesh data
       yc(n) = SDATM%grid%data%rattr(klat,n)
    end do

    ! overwrite mask and frac
    k = mct_aVect_indexRA(SDATM%grid%data,'mask')
    SDATM%grid%data%rAttr(k,:) = 1.0_R8

    k = mct_aVect_indexRA(SDATM%grid%data,'frac')
    SDATM%grid%data%rAttr(k,:) = 1.0_R8

    if (my_task == master_task) then
       call shr_strdata_print(SDATM,'ATM data')
    endif

    !----------------------------------------------------------------------------
    ! Initialize SDATM attributes for streams and mapping of streams to model domain
    !----------------------------------------------------------------------------

    call shr_strdata_init_streams(SDATM, compid, mpicom, my_task)
    call shr_strdata_init_mapping(SDATM, compid, mpicom, my_task)

    !----------------------------------------------------------------------------
    ! allocate module arrays
    !----------------------------------------------------------------------------

    allocate(windFactor(lsize))
    allocate(winddFactor(lsize))
    allocate(qsatFactor(lsize))

    call t_stopf('datm_strdata_init')

    !----------------------------------------------------------------------------
    ! Initialize attribute vectors
    !----------------------------------------------------------------------------

    call t_startf('datm_initmctavs')
    if (my_task == master_task) write(logunit,F00) 'allocate AVs'

    call mct_aVect_init(a2x, rList=flds_a2x_mod, lsize=lsize)
    call mct_aVect_zero(a2x)
    call mct_aVect_init(x2a, rList=flds_x2a_mod, lsize=lsize)
    call mct_aVect_zero(x2a)

    ! Initialize internal attribute vectors for optional streams
    ! Create the colon deliminted list flds_strm based on mapping the
    ! input stream fields from SDATM%avs(n) to with the names in stifld to stofld

    cnt = 0
    flds_strm = ''
    do n = 1,SDATM%nstreams
       ! Loop over the field names in stifld
       do k = 1,size(stifld)
          ! Search the streams for the field name stifld(k)
          kfld = mct_aVect_indexRA(SDATM%avs(n), trim(stifld(k)), perrWith='quiet')
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
    ! flds_strm - the indices below CANNOT be set in the datm_comp_advertise phase

    ! Now set indices into these active streams
    stbot  = mct_aVect_indexRA(avstrm,'strm_tbot'   ,perrWith='quiet')
    swind  = mct_aVect_indexRA(avstrm,'strm_wind'   ,perrWith='quiet')
    sz     = mct_aVect_indexRA(avstrm,'strm_z'      ,perrWith='quiet')
    spbot  = mct_aVect_indexRA(avstrm,'strm_pbot'   ,perrWith='quiet')
    sshum  = mct_aVect_indexRA(avstrm,'strm_shum'   ,perrWith='quiet')
    stdew  = mct_aVect_indexRA(avstrm,'strm_tdew'   ,perrWith='quiet')
    srh    = mct_aVect_indexRA(avstrm,'strm_rh'     ,perrWith='quiet')
    slwdn  = mct_aVect_indexRA(avstrm,'strm_lwdn'   ,perrWith='quiet')
    sswdn  = mct_aVect_indexRA(avstrm,'strm_swdn'   ,perrWith='quiet')
    sswdndf= mct_aVect_indexRA(avstrm,'strm_swdndf' ,perrWith='quiet')
    sswdndr= mct_aVect_indexRA(avstrm,'strm_swdndr' ,perrWith='quiet')
    sprecc = mct_aVect_indexRA(avstrm,'strm_precc'  ,perrWith='quiet')
    sprecl = mct_aVect_indexRA(avstrm,'strm_precl'  ,perrWith='quiet')
    sprecn = mct_aVect_indexRA(avstrm,'strm_precn'  ,perrWith='quiet')
    sco2p  = mct_aVect_indexRA(avstrm,'strm_co2p'   ,perrWith='quiet')
    sco2d  = mct_aVect_indexRA(avstrm,'strm_co2d'   ,perrWith='quiet')
    sswup  = mct_aVect_indexRA(avstrm,'strm_swup'   ,perrWith='quiet')
    sprec  = mct_aVect_indexRA(avstrm,'strm_prec'   ,perrWith='quiet')
    starcf = mct_aVect_indexRA(avstrm,'strm_tarcf'  ,perrWith='quiet')

    ! anomaly forcing
    sprecsf  = mct_aVect_indexRA(avstrm,'strm_precsf'  ,perrWith='quiet')
    sprec_af = mct_aVect_indexRA(avstrm,'strm_prec_af' ,perrWith='quiet')
    su_af    = mct_aVect_indexRA(avstrm,'strm_u_af'    ,perrWith='quiet')
    sv_af    = mct_aVect_indexRA(avstrm,'strm_v_af'    ,perrWith='quiet')
    stbot_af = mct_aVect_indexRA(avstrm,'strm_tbot_af' ,perrWith='quiet')
    spbot_af = mct_aVect_indexRA(avstrm,'strm_pbot_af' ,perrWith='quiet')
    sshum_af = mct_aVect_indexRA(avstrm,'strm_shum_af' ,perrWith='quiet')
    sswdn_af = mct_aVect_indexRA(avstrm,'strm_swdn_af' ,perrWith='quiet')
    slwdn_af = mct_aVect_indexRA(avstrm,'strm_lwdn_af' ,perrWith='quiet')

    ! isotopic forcing
    if (wiso_datm) then
       sprecn_16O = mct_aVect_indexRA(avstrm,'strm_precn_16O',perrWith='quiet')
       sprecn_18O = mct_aVect_indexRA(avstrm,'strm_precn_18O',perrWith='quiet')
       sprecn_HDO = mct_aVect_indexRA(avstrm,'strm_precn_HDO',perrWith='quiet')
       ! Okay here to just use srh_18O and srh_HDO, because the forcing is (should)
       ! just be deltas, applied in CTSM to the base tracer
       srh_16O    = mct_aVect_indexRA(avstrm,'strm_rh_16O',perrWith='quiet')
       srh_18O    = mct_aVect_indexRA(avstrm,'strm_rh_18O',perrWith='quiet')
       srh_HDO    = mct_aVect_indexRA(avstrm,'strm_rh_HDO',perrWith='quiet')
    end if

    call t_stopf('datm_initmctavs')

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

       ! if (exists1) then
       !    if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file)
       !    call shr_pcdf_readwrite('read',SDATM%pio_subsystem, SDATM%io_type, &
       !         trim(rest_file),mpicom,gsmap=SDATM%gsmap,rf1=water,rf1n='water',io_format=SDATM%io_format)
       ! else
       !    if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file)
       ! endif

       if (exists) then
          if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file_strm)
          call shr_strdata_restRead(trim(rest_file_strm),SDATM,mpicom)
       else
          if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file_strm)
       endif
    endif

    !----------------------------------------------------------------------------
    ! Set initial atm state
    !----------------------------------------------------------------------------

    call t_adj_detailf(+2)
    call datm_comp_run(&
         x2a=x2a, &
         a2x=a2x, &
         SDATM=SDATM, &
         mpicom=mpicom, &
         compid=compid, &
         my_task=my_task, &
         master_task=master_task, &
         inst_suffix=inst_suffix, &
         logunit=logunit, &
         orbEccen=orbEccen, &
         orbMvelpp=orbMvelpp, &
         orbLambm0=orbLambm0, &
         orbObliqr=orbObliqr, &
         write_restart=.false., &
         target_ymd=current_ymd, &
         target_tod=current_tod, &
         target_mon=current_mon, &
         calendar=calendar, &
         modeldt=modeldt, &
         atm_prognostic=atm_prognostic)
    call t_adj_detailf(-2)

    call t_stopf('DATM_INIT')

  end subroutine datm_comp_init

  !===============================================================================

  subroutine datm_comp_run(x2a, a2x, &
       SDATM, mpicom, compid, my_task, master_task, &
       inst_suffix, logunit, &
       orbEccen, orbMvelpp, orbLambm0, orbObliqr, &
       write_restart, target_ymd, target_tod, target_mon, modeldt, calendar, &
       atm_prognostic, case_name)

    ! !DESCRIPTION: run method for datm model

    ! !INPUT/OUTPUT PARAMETERS:
    type(mct_aVect)        , intent(inout) :: x2a
    type(mct_aVect)        , intent(inout) :: a2x
    type(shr_strdata_type) , intent(inout) :: SDATM
    integer                , intent(in)    :: mpicom           ! mpi communicator
    integer                , intent(in)    :: compid           ! mct comp id
    integer                , intent(in)    :: my_task          ! my task in mpi communicator mpicom
    integer                , intent(in)    :: master_task      ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix      ! char string associated with instance
    integer                , intent(in)    :: logunit          ! logging unit number
    real(R8)               , intent(in)    :: orbEccen         ! orb eccentricity (unit-less)
    real(R8)               , intent(in)    :: orbMvelpp        ! orb moving vernal eq (radians)
    real(R8)               , intent(in)    :: orbLambm0        ! orb mean long of perhelion (radians)
    real(R8)               , intent(in)    :: orbObliqr        ! orb obliquity (radians)
    logical                , intent(in)    :: write_restart    ! restart alarm is on
    integer                , intent(in)    :: target_ymd       ! model date
    integer                , intent(in)    :: target_tod       ! model sec into model date
    integer                , intent(in)    :: target_mon       ! model month
    character(len=*)       , intent(in)    :: calendar         ! calendar type
    Integer                , intent(in)    :: modeldt          ! model time step
    logical                , intent(in)    :: atm_prognostic
    character(len=*)       , intent(in), optional :: case_name ! case name

    !--- local ---
    integer                 :: n,nfld            ! indices
    integer                 :: lsize             ! size of attr vect
    character(CL)           :: rest_file         ! restart_file
    character(CL)           :: rest_file_strm    ! restart_file
    integer                 :: nu                ! unit number
    integer                 :: eday              ! elapsed day
    real(R8)                :: rday              ! elapsed day
    real(R8)                :: cosFactor         ! cosine factor
    real(R8)                :: factor            ! generic/temporary correction factor
    real(R8)                :: avg_alb           ! average albedo
    real(R8)                :: tMin              ! minimum temperature
    character(len=18)       :: date_str
    character(len=CS)       :: fldname
    real(R8)                :: uprime,vprime,swndr,swndf,swvdr,swvdf,ratio_rvrf
    real(R8)                :: tbot,pbot,rtmp,vp,ea,e,qsat,frac,qsatT
    logical                 :: firstcall = .true.    ! first call logical
    character(*), parameter :: F00   = "('(datm_comp_run) ',8a)"
    character(*), parameter :: F04   = "('(datm_comp_run) ',2a,2i8,'s')"
    character(*), parameter :: F0D   = "('(datm_comp_run) ',a, i7,2x,i5,2x,i5,2x,d21.14)"
    character(*), parameter :: subName = "(datm_comp_run) "
    !-------------------------------------------------------------------------------

    !--------------------
    ! Debug output
    !--------------------

    if (debug_import > 0 .and. my_task == master_task .and. atm_prognostic) then
       do nfld = 1, mct_aVect_nRAttr(x2a)
          call shr_string_listGetName(trim(flds_x2a_mod), nfld, fldname)
          do n = 1, mct_aVect_lsize(x2a)
             write(logunit,F0D)'import: ymd,tod,n  = '// trim(fldname),target_ymd, target_tod, &
                  n, x2a%rattr(nfld,n)
          end do
       end do
    end if

    !--------------------
    ! ADVANCE ATM
    !--------------------

    call t_startf('DATM_RUN')
    call t_barrierf('datm_BARRIER',mpicom)
    call t_startf('datm')

    !--- set data needed for cosz t-interp method ---
    call shr_strdata_setOrbs(SDATM,orbEccen,orbMvelpp,orbLambm0,orbObliqr,modeldt)

    !--- copy all fields from streams to a2x as default ---
    call t_startf('datm_strdata_advance')
    call shr_strdata_advance(SDATM,target_ymd,target_tod,mpicom,'datm')
    call t_stopf('datm_strdata_advance')

    call t_barrierf('datm_scatter_BARRIER',mpicom)

    call t_startf('datm_scatter')
    if (firstcall) then
       allocate(ilist_av(SDATM%nstreams))
       allocate(olist_av(SDATM%nstreams))
       allocate(ilist_st(SDATM%nstreams))
       allocate(olist_st(SDATM%nstreams))
       allocate(count_av(SDATM%nstreams))
       allocate(count_st(SDATM%nstreams))
       do n = 1,SDATM%nstreams
          ! Obtain a smaller list for translate given the actual fields present in the streams
          ! This can only be done once the SDATM has been initialized
          call shr_dmodel_translate_list(   &
               avi=SDATM%avs(n),            & ! input  av
               avo=a2x,                     & ! output av
               avifld=avifld,               & ! input  field names for translation
               avofld=avofld,               & ! output field names for translation
               ilist=ilist_av(n),           & ! input  list for translation
               olist=olist_av(n),           & ! output list for translation
               cnt=count_av(n))               ! indices
       end do
       do n = 1,SDATM%nstreams
          call shr_dmodel_translate_list(   &
               avi=SDATM%avs(n),            & ! input av
               avo=avstrm,                  & ! output av
               avifld=stifld,               & ! input  field names for translation
               avofld=stofld,               & ! output field names for translation
               ilist=ilist_st(n),           & ! input  list for translation
               olist=olist_st(n),           & ! output list for translation
               cnt=count_st(n))               ! indices
       end do
    end if

    ! At this point DATM%avs(n) has been interpolated to the model
    ! grid and the model time

    ! Fill in a2x from ALL the streams in SDATM%avs(:)
    do n = 1,SDATM%nstreams
       if (count_av(n) > 0) then
          call shr_dmodel_translateAV_list( avi=SDATM%avs(n), avo=a2x, &
               ilist=ilist_av(n), olist=olist_av(n))
       end if
    enddo

    ! Fill in avstrm from ALL the streams in SDATM%avs(:)
    do n = 1,SDATM%nstreams
       if (count_st(n) > 0) then
          call shr_dmodel_translateAV_list( avi=SDATM%avs(n), avo=avstrm, &
               ilist=ilist_st(n), olist=olist_st(n))
       end if
    enddo
    call t_stopf('datm_scatter')

    !-------------------------------------------------
    ! Determine data model behavior based on the mode
    !-------------------------------------------------

    call t_startf('datm_datamode')
    select case (trim(datamode))

    case('COPYALL')
       ! do nothing extra

    case('CORE2_NYF','CORE2_IAF')
       if (firstcall) then
          if (sprec < 1 .or. sswdn < 1) then
             write(logunit,F00) 'ERROR: prec and swdn must be in streams for CORE2'
             call shr_sys_abort(trim(subname)//'ERROR: prec and swdn must be in streams for CORE2')
          endif
          if (trim(datamode) == 'CORE2_IAF' ) then
             if (starcf < 1 ) then
                write(logunit,F00) 'ERROR: tarcf must be in an input stream for CORE2_IAF'
                call shr_sys_abort(trim(subname)//'tarcf must be in an input stream for CORE2_IAF')
             endif
          endif
          call datm_shr_CORE2getFactors(factorFn,windFactor,winddFactor,qsatFactor, &
               mpicom,compid, SDATM%gsmap, SDATM%grid, SDATM%nxg, SDATM%nyg)
       endif
       call shr_cal_date2julian(target_ymd,target_tod,rday,calendar)
       rday = mod((rday - 1.0_R8),365.0_R8)
       cosfactor = cos((2.0_R8*SHR_CONST_PI*rday)/365 - phs_c0)

       lsize = mct_avect_lsize(a2x)
       do n = 1,lsize
          a2x%rAttr(kz,n) = 10.0_R8

          !--- correction to NCEP winds based on QSCAT ---
          uprime    = a2x%rAttr(ku,n)*windFactor(n)
          vprime    = a2x%rAttr(kv,n)*windFactor(n)
          a2x%rAttr(ku,n) = uprime*cos(winddFactor(n)*degtorad)- &
                            vprime*sin(winddFactor(n)*degtorad)
          a2x%rAttr(kv,n) = uprime*sin(winddFactor(n)*degtorad)+ &
                            vprime*cos(winddFactor(n)*degtorad)

          !--- density, tbot, & pslv taken directly from input stream, set pbot ---
          a2x%rAttr(kpbot,n) = a2x%rAttr(kpslv,n)

          !--- correction to NCEP Arctic & Antarctic air T & potential T ---
          if      ( yc(n) < -60.0_R8 ) then
             tMin = (avg_c0 + avg_c1*yc(n)) + (amp_c0 + amp_c1*yc(n))*cosFactor + tKFrz
             a2x%rAttr(ktbot,n) = max(a2x%rAttr(ktbot,n), tMin)
          else if ( yc(n) > 60.0_R8 ) then
             factor = MIN(1.0_R8, 0.1_R8*(yc(n)-60.0_R8) )
             a2x%rAttr(ktbot,n) = a2x%rAttr(ktbot,n) + factor * dTarc(target_mon)
          endif
          a2x%rAttr(kptem,n) = a2x%rAttr(ktbot,n)

          !---  correction to NCEP relative humidity for heat budget balance ---
          a2x%rAttr(kshum,n) = a2x%rAttr(kshum,n) + qsatFactor(n)

          !--- Dupont correction to NCEP Arctic air T  ---
          !--- don't correct during summer months (July-September)
          !--- ONLY correct when forcing year is 1997->2004
          if (trim(datamode) == 'CORE2_IAF' ) then
             a2x%rAttr(ktbot,n) = a2x%rAttr(ktbot,n) +  avstrm%rAttr(starcf,n)
             a2x%rAttr(kptem,n) = a2x%rAttr(ktbot,n)
          end if

          !-------------------------------------------------------------------------
          ! PRECIPITATION DATA
          !-------------------------------------------------------------------------

          avstrm%rAttr(sprec,n) = avstrm%rAttr(sprec,n)/86400.0_R8        ! convert mm/day to kg/m^2/s

          !  only correct satellite products, do not correct Serreze Arctic data
          if ( yc(n) < 58. ) then
             avstrm%rAttr(sprec,n) = avstrm%rAttr(sprec,n)*1.14168_R8
          endif
          if ( yc(n) >= 58. .and. yc(n) < 68. ) then
             factor = MAX(0.0_R8, 1.0_R8 - 0.1_R8*(yc(n)-58.0_R8) )
             avstrm%rAttr(sprec,n) = avstrm%rAttr(sprec,n)*(factor*(1.14168_R8 - 1.0_R8) + 1.0_R8)
          endif

          a2x%rAttr(krc,n) = 0.0_R8                    ! default zero
          a2x%rAttr(ksc,n) = 0.0_R8
          if (a2x%rAttr(ktbot,n) < tKFrz ) then        ! assign precip to rain/snow components
             a2x%rAttr(krl,n) = 0.0_R8
             a2x%rAttr(ksl,n) = avstrm%rAttr(sprec,n)
          else
             a2x%rAttr(krl,n) = avstrm%rAttr(sprec,n)
             a2x%rAttr(ksl,n) = 0.0_R8
          endif

          !-------------------------------------------------------------------------
          ! RADIATION DATA
          !-------------------------------------------------------------------------

          !--- fabricate required swdn components from net swdn ---
          a2x%rAttr(kswvdr,n) = avstrm%rAttr(sswdn,n)*(0.28_R8)
          a2x%rAttr(kswndr,n) = avstrm%rAttr(sswdn,n)*(0.31_R8)
          a2x%rAttr(kswvdf,n) = avstrm%rAttr(sswdn,n)*(0.24_R8)
          a2x%rAttr(kswndf,n) = avstrm%rAttr(sswdn,n)*(0.17_R8)

          !--- compute net short-wave based on LY08 latitudinally-varying albedo ---
          avg_alb = ( 0.069 - 0.011*cos(2.0_R8*yc(n)*degtorad ) )
          a2x%rAttr(kswnet,n) = avstrm%rAttr(sswdn,n)*(1.0_R8 - avg_alb)

          !--- corrections to GISS sswdn for heat budget balancing ---
          factor = 1.0_R8
          if      ( -60.0_R8 < yc(n) .and. yc(n) < -50.0_R8 ) then
             factor = 1.0_R8 - (yc(n) + 60.0_R8)*(0.05_R8/10.0_R8)
          else if ( -50.0_R8 < yc(n) .and. yc(n) <  30.0_R8 ) then
             factor = 0.95_R8
          else if (  30.0_R8 < yc(n) .and. yc(n) <  40._R8 ) then
             factor = 1.0_R8 - (40.0_R8 - yc(n))*(0.05_R8/10.0_R8)
          endif
          a2x%rAttr(kswnet,n) = a2x%rAttr(kswnet,n)*factor
          a2x%rAttr(kswvdr,n) = a2x%rAttr(kswvdr,n)*factor
          a2x%rAttr(kswndr,n) = a2x%rAttr(kswndr,n)*factor
          a2x%rAttr(kswvdf,n) = a2x%rAttr(kswvdf,n)*factor
          a2x%rAttr(kswndf,n) = a2x%rAttr(kswndf,n)*factor

          !--- correction to GISS lwdn in Arctic ---
          if ( yc(n) > 60._R8 ) then
             factor = MIN(1.0_R8, 0.1_R8*(yc(n)-60.0_R8) )
             a2x%rAttr(klwdn,n) = a2x%rAttr(klwdn,n) + factor * dLWarc
          endif

       enddo   ! lsize

    case('CORE_IAF_JRA')
       if (firstcall) then
          if (sprec < 1 .or. sswdn < 1) then
             write(logunit,F00) 'ERROR: prec and swdn must be in streams for CORE_IAF_JRA'
             call shr_sys_abort(trim(subname)//'ERROR: prec and swdn must be in streams for CORE_IAF_JRA')
          endif
          if (trim(datamode) == 'CORE_IAF_JRA' ) then
             if (starcf < 1 ) then
                write(logunit,F00) 'ERROR: tarcf must be in an input stream for CORE_IAF_JRA'
                call shr_sys_abort(trim(subname)//'tarcf must be in an input stream for CORE_IAF_JRA')
             endif
          endif
          if (trim(factorFn) == 'null') then
            windFactor = 1.0_R8
            winddFactor = 1.0_R8
            qsatFactor = 1.0_R8
          else
            call datm_shr_CORE2getFactors(factorFn,windFactor,winddFactor,qsatFactor, &
                 mpicom, compid, SDATM%gsmap, SDATM%grid, SDATM%nxg, SDATM%nyg)
          endif
       endif
       call shr_cal_date2julian(target_ymd,target_tod,rday,calendar)
       rday = mod((rday - 1.0_R8),365.0_R8)
       cosfactor = cos((2.0_R8*SHR_CONST_PI*rday)/365 - phs_c0)

       lsize = mct_avect_lsize(a2x)
       do n = 1,lsize
          a2x%rAttr(kz,n) = 10.0_R8

          !--- density, tbot, & pslv taken directly from input stream, set pbot ---
          a2x%rAttr(kpbot,n) = a2x%rAttr(kpslv,n)

          a2x%rAttr(kptem,n) = a2x%rAttr(ktbot,n)

          !--- density computation for JRA55 forcing ---
          a2x%rAttr(kdens,n) = a2x%rAttr(kpbot,n)/(rdair*a2x%rAttr(ktbot,n) &
                               *(1+0.608* a2x%rAttr(kshum,n)))

          !-------------------------------------------------------------------------
          ! PRECIPITATION DATA
          !-------------------------------------------------------------------------

          a2x%rAttr(krc,n) = 0.0_R8                    ! default zero
          a2x%rAttr(ksc,n) = 0.0_R8
          if (a2x%rAttr(ktbot,n) < tKFrz ) then        ! assign precip to rain/snow components
             a2x%rAttr(krl,n) = 0.0_R8
             a2x%rAttr(ksl,n) = avstrm%rAttr(sprec,n)
          else
             a2x%rAttr(krl,n) = avstrm%rAttr(sprec,n)
             a2x%rAttr(ksl,n) = 0.0_R8
          endif

          !-------------------------------------------------------------------------
          ! RADIATION DATA
          !-------------------------------------------------------------------------

          !--- fabricate required swdn components from net swdn ---
          a2x%rAttr(kswvdr,n) = avstrm%rAttr(sswdn,n)*(0.28_R8)
          a2x%rAttr(kswndr,n) = avstrm%rAttr(sswdn,n)*(0.31_R8)
          a2x%rAttr(kswvdf,n) = avstrm%rAttr(sswdn,n)*(0.24_R8)
          a2x%rAttr(kswndf,n) = avstrm%rAttr(sswdn,n)*(0.17_R8)

          !--- compute net short-wave based on LY08 latitudinally-varying albedo ---
          avg_alb = ( 0.069 - 0.011*cos(2.0_R8*yc(n)*degtorad ) )
          a2x%rAttr(kswnet,n) = avstrm%rAttr(sswdn,n)*(1.0_R8 - avg_alb)

       enddo   ! lsize

    case('CLMNCEP')
       if (firstcall) then
          if (swind < 1 .or. stbot < 1) then
             write(logunit,F00) ' ERROR: wind and tbot must be in streams for CLMNCEP'
             call shr_sys_abort(trim(subname)//' ERROR: wind and tbot must be in streams for CLMNCEP')
          endif
          rtmp = maxval(a2x%rAttr(ktbot,:))
          call shr_mpi_max(rtmp,tbotmax,mpicom,'datm_tbot',all=.true.)
          if (atm_prognostic) then
             rtmp = maxval(x2a%rAttr(kanidr,:))
             call shr_mpi_max(rtmp,anidrmax,mpicom,'datm_ani',all=.true.)
          else
             anidrmax = SHR_CONST_SPVAL ! see below for use
          end if
          if (stdew > 0) then
             rtmp = maxval(avstrm%rAttr(stdew,:))
             call shr_mpi_max(rtmp,tdewmax,mpicom,'datm_tdew',all=.true.)
          endif
          if (my_task == master_task) &
               write(logunit,*) trim(subname),' max values = ',tbotmax,tdewmax,anidrmax
       endif
       lsize = mct_avect_lsize(a2x)
       do n = 1,lsize
          !--- bottom layer height ---
          if (sz < 1) a2x%rAttr(kz,n) = 30.0_R8

          !--- temperature ---
          if (tbotmax < 50.0_R8) a2x%rAttr(ktbot,n) = a2x%rAttr(ktbot,n) + tkFrz
          ! Limit very cold forcing to 180K
          a2x%rAttr(ktbot,n) = max(180._r8, a2x%rAttr(ktbot,n))
          a2x%rAttr(kptem,n) = a2x%rAttr(ktbot,n)

          !--- pressure ---
          if (spbot < 1) a2x%rAttr(kpbot,n) = pstd
          a2x%rAttr(kpslv,n) = a2x%rAttr(kpbot,n)

          !--- u, v wind velocity ---
          a2x%rAttr(ku,n) = avstrm%rAttr(swind,n)/sqrt(2.0_R8)
          a2x%rAttr(kv,n) = a2x%rAttr(ku,n)

          !--- specific humidity ---
          tbot = a2x%rAttr(ktbot,n)
          pbot = a2x%rAttr(kpbot,n)
          if (sshum > 0) then
             e = datm_shr_esat(tbot,tbot)
             qsat = (0.622_R8 * e)/(pbot - 0.378_R8 * e)
             if (qsat < a2x%rAttr(kshum,n)) then
                a2x%rAttr(kshum,n) = qsat
             endif
          else if (srh > 0) then
             e = avstrm%rAttr(srh,n) * 0.01_R8 * datm_shr_esat(tbot,tbot)
             qsat = (0.622_R8 * e)/(pbot - 0.378_R8 * e)
             a2x%rAttr(kshum,n) = qsat
             if(wiso_datm) then
                ! isotopic forcing
                ! For tracer specific humidity, lnd_import_mct expects a delta, so
                ! just keep the delta from the input file - TW
                a2x%rAttr(kshum_16O,n) = avstrm%rAttr(srh_16O,n)
                a2x%rAttr(kshum_18O,n) = avstrm%rAttr(srh_18O,n)
                a2x%rAttr(kshum_HDO,n) = avstrm%rAttr(srh_HDO,n)
             end if
          else if (stdew > 0) then
             if (tdewmax < 50.0_R8) avstrm%rAttr(stdew,n) = avstrm%rAttr(stdew,n) + tkFrz
             e = datm_shr_esat(avstrm%rAttr(stdew,n),tbot)
             qsat = (0.622_R8 * e)/(pbot - 0.378_R8 * e)
             a2x%rAttr(kshum,n) = qsat
          else
             call shr_sys_abort(subname//'ERROR: cannot compute shum')
          endif

          !--- density ---
          vp = (a2x%rAttr(kshum,n)*pbot) / (0.622_R8 + 0.378_R8 * a2x%rAttr(kshum,n))
          a2x%rAttr(kdens,n) = (pbot - 0.378_R8 * vp) / (tbot*rdair)

          !--- downward longwave ---
          if (slwdn < 1) then
             e  = a2x%rAttr(kpslv,n) * a2x%rAttr(kshum,n) / (0.622_R8 + 0.378_R8 * a2x%rAttr(kshum,n))
             ea = 0.70_R8 + 5.95e-05_R8 * 0.01_R8 * e * exp(1500.0_R8/tbot)
             a2x%rAttr(klwdn,n) = ea * stebol * tbot**4
          endif

          !--- shortwave radiation ---
          if (sswdndf > 0 .and. sswdndr > 0) then
             a2x%rAttr(kswndr,n) = avstrm%rAttr(sswdndr,n) * 0.50_R8
             a2x%rAttr(kswvdr,n) = avstrm%rAttr(sswdndr,n) * 0.50_R8
             a2x%rAttr(kswndf,n) = avstrm%rAttr(sswdndf,n) * 0.50_R8
             a2x%rAttr(kswvdf,n) = avstrm%rAttr(sswdndf,n) * 0.50_R8
          elseif (sswdn > 0) then
             ! relationship between incoming NIR or VIS radiation and ratio of
             ! direct to diffuse radiation calculated based on one year's worth of
             ! hourly CAM output from CAM version cam3_5_55
             swndr = avstrm%rAttr(sswdn,n) * 0.50_R8
             ratio_rvrf =   min(0.99_R8,max(0.29548_R8 + 0.00504_R8*swndr  &
                  -1.4957e-05_R8*swndr**2 + 1.4881e-08_R8*swndr**3,0.01_R8))
             a2x%rAttr(kswndr,n) = ratio_rvrf*swndr
             swndf = avstrm%rAttr(sswdn,n) * 0.50_R8
             a2x%rAttr(kswndf,n) = (1._R8 - ratio_rvrf)*swndf

             swvdr = avstrm%rAttr(sswdn,n) * 0.50_R8
             ratio_rvrf =   min(0.99_R8,max(0.17639_R8 + 0.00380_R8*swvdr  &
                  -9.0039e-06_R8*swvdr**2 + 8.1351e-09_R8*swvdr**3,0.01_R8))
             a2x%rAttr(kswvdr,n) = ratio_rvrf*swvdr
             swvdf = avstrm%rAttr(sswdn,n) * 0.50_R8
             a2x%rAttr(kswvdf,n) = (1._R8 - ratio_rvrf)*swvdf
          else
             call shr_sys_abort(subName//'ERROR: cannot compute short-wave down')
          endif

          !--- swnet: a diagnostic quantity ---
          if (anidrmax < 1.0e-8 .or. anidrmax > SHR_CONST_SPVAL * 0.9_R8) then
             a2x%rAttr(kswnet,n) = 0.0_R8
          else
             a2x%rAttr(kswnet,n) = (1.0_R8-x2a%rAttr(kanidr,n))*a2x%rAttr(kswndr,n) + &
                                   (1.0_R8-x2a%rAttr(kavsdr,n))*a2x%rAttr(kswvdr,n) + &
                                   (1.0_R8-x2a%rAttr(kanidf,n))*a2x%rAttr(kswndf,n) + &
                                   (1.0_R8-x2a%rAttr(kavsdf,n))*a2x%rAttr(kswvdf,n)
          endif

          !--- rain and snow ---
          if (sprecc > 0 .and. sprecl > 0) then
             a2x%rAttr(krc,n) = avstrm%rAttr(sprecc,n)
             a2x%rAttr(krl,n) = avstrm%rAttr(sprecl,n)
          elseif (sprecn > 0) then
             a2x%rAttr(krc,n) = avstrm%rAttr(sprecn,n)*0.1_R8
             a2x%rAttr(krl,n) = avstrm%rAttr(sprecn,n)*0.9_R8
          else
             call shr_sys_abort(subName//'ERROR: cannot compute rain and snow')
          endif

          !--- split precip between rain & snow ---
          call shr_precip_partition_rain_snow_ramp(tbot, frac)
          a2x%rAttr(ksc,n) = max(0.0_R8, a2x%rAttr(krc,n)*(1.0_R8 - frac) )
          a2x%rAttr(ksl,n) = max(0.0_R8, a2x%rAttr(krl,n)*(1.0_R8 - frac) )
          a2x%rAttr(krc,n) = max(0.0_R8, a2x%rAttr(krc,n)*(         frac) )
          a2x%rAttr(krl,n) = max(0.0_R8, a2x%rAttr(krl,n)*(         frac) )

       enddo

    end select

    !----------------------------------------------------------
    ! bias correction / anomaly forcing ( start block )
    !----------------------------------------------------------

    ! modify atmospheric input fields if streams exist
    lsize = mct_avect_lsize(avstrm)

    ! bias correct precipitation relative to observed
    ! (via bias_correct nameslist option)
    if (sprecsf > 0) then
       do n = 1,lsize
          a2x%rAttr(ksc,n) = a2x%rAttr(ksc,n) * min(1.e2_r8,avstrm%rAttr(sprecsf,n))
          a2x%rAttr(ksl,n) = a2x%rAttr(ksl,n) * min(1.e2_r8,avstrm%rAttr(sprecsf,n))
          a2x%rAttr(krc,n) = a2x%rAttr(krc,n) * min(1.e2_r8,avstrm%rAttr(sprecsf,n))
          a2x%rAttr(krl,n) = a2x%rAttr(krl,n) * min(1.e2_r8,avstrm%rAttr(sprecsf,n))

       end do
    endif

    ! adjust atmospheric input fields if anomaly forcing streams exist
    ! (via anomaly_forcing namelist option)

    ! wind
    if (su_af > 0 .and. sv_af > 0) then
       do n = 1,lsize
          a2x%rAttr(ku,n) = a2x%rAttr(ku,n) + avstrm%rAttr(su_af,n)
          a2x%rAttr(kv,n) = a2x%rAttr(kv,n) + avstrm%rAttr(sv_af,n)
       end do
    endif

    ! specific humidity
    if (sshum_af > 0) then
       do n = 1,lsize
          a2x%rAttr(kshum,n) = a2x%rAttr(kshum,n) + avstrm%rAttr(sshum_af,n)

          ! avoid possible negative q values
          if(a2x%rAttr(kshum,n) < 0._r8) then
             a2x%rAttr(kshum,n) = 1.e-6_r8
          endif

       end do
    endif

    ! pressure
    if (spbot_af > 0) then
       do n = 1,lsize
          a2x%rAttr(kpbot,n) = a2x%rAttr(kpbot,n) + avstrm%rAttr(spbot_af,n)
       end do
    endif

    ! temperature
    if (stbot_af > 0) then
       do n = 1,lsize
          a2x%rAttr(ktbot,n) = a2x%rAttr(ktbot,n) + avstrm%rAttr(stbot_af,n)
       end do
    endif

    ! longwave
    if (slwdn_af > 0) then
       do n = 1,lsize
          a2x%rAttr(klwdn,n) = a2x%rAttr(klwdn,n) * avstrm%rAttr(slwdn_af,n)
       end do
    endif

    ! precipitation
    if (sprec_af > 0) then
       do n = 1,lsize
          a2x%rAttr(ksc,n) = a2x%rAttr(ksc,n) * avstrm%rAttr(sprec_af,n)
          a2x%rAttr(ksl,n) = a2x%rAttr(ksl,n) * avstrm%rAttr(sprec_af,n)
          a2x%rAttr(krc,n) = a2x%rAttr(krc,n) * avstrm%rAttr(sprec_af,n)
          a2x%rAttr(krl,n) = a2x%rAttr(krl,n) * avstrm%rAttr(sprec_af,n)
       enddo
    endif

    ! shortwave
    if (sswdn_af > 0) then
       do n = 1,lsize
          a2x%rAttr(kswndr,n) = a2x%rAttr(kswndr,n) * avstrm%rAttr(sswdn_af,n)
          a2x%rAttr(kswvdr,n) = a2x%rAttr(kswvdr,n) * avstrm%rAttr(sswdn_af,n)
          a2x%rAttr(kswndf,n) = a2x%rAttr(kswndf,n) * avstrm%rAttr(sswdn_af,n)
          a2x%rAttr(kswvdf,n) = a2x%rAttr(kswvdf,n) * avstrm%rAttr(sswdn_af,n)
       enddo
    endif
    !--------------------
    ! bias correction / anomaly forcing ( end block )
    !--------------------

    call t_stopf('datm_datamode')

    !--------------------
    ! Debug output
    !--------------------

    if (debug_export > 0 .and. my_task == master_task) then
       do nfld = 1, mct_aVect_nRAttr(a2x)
          call shr_string_listGetName(trim(flds_a2x_mod), nfld, fldname)
          do n = 1, mct_aVect_lsize(a2x)
             write(logunit,F0D)'export: ymd,tod,n  = '// trim(fldname),target_ymd, target_tod, &
                  n, a2x%rattr(nfld,n)
          end do
       end do
    end if

    !--------------------
    ! Write restart
    !--------------------

    if (write_restart) then
       call t_startf('datm_restart')
       call shr_cal_datetod2string(date_str, target_ymd, target_tod)

       write(rest_file,"(6a)") &
            trim(case_name), '.datm',trim(inst_suffix),'.r.', trim(date_str), '.nc'
       write(rest_file_strm,"(6a)") &
            trim(case_name), '.datm',trim(inst_suffix),'.rs1.', trim(date_str), '.bin'
       if (my_task == master_task) then
          nu = shr_file_getUnit()
          open(nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
          write(nu,'(a)') rest_file
          write(nu,'(a)') rest_file_strm
          close(nu)
          call shr_file_freeUnit(nu)
       endif

       if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file_strm),target_ymd,target_tod
       call shr_strdata_restWrite(trim(rest_file_strm),SDATM,mpicom,trim(case_name),'SDATM strdata')
       call t_stopf('datm_restart')
    endif

    firstcall = .false.

    call t_stopf('datm')
    call t_stopf('DATM_RUN')

  end subroutine datm_comp_run

end module datm_comp_mod
